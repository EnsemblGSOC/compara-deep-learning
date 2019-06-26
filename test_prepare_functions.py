import pandas as pd 
import progressbar
import json
import numpy as np 
import multiprocessing
from threading import Thread
from ete3 import Tree
from read_get_gene_seq import read_gene_seq,create_dict
from process_data import get_nearest_neighbors
from create_synteny_matrix import create_synteny_matrix_mul

def create_data_homology_ls(df,n,a,d,ld,ldg,cmap,cimap):
    buf=open("not_found.txt","w")
    lsy={} #dictionary which stores +/- n genes of the given gene by id. Each key is a gene id which corresponds to the one in center.
    for _,row in progressbar.progressbar(df.iterrows()):
        x=row[1]
        y=row[3]
        xs=row[2]
        ys=row[4]
        try:
            z=lsy[x]
        except:
            try:
                t2=d[xs.capitalize()]#see if the species exist in genomic maps
                xl,xr=get_nearest_neighbors(x,xs,n,a,d,ld,ldg,cmap,cimap)
                if len(xl)!=0:#check if neighboring genes were successfully found
                    lsy[x]=dict(b=xl,f=xr)                 
            except Exception as e:
                buf.write(x+"\t"+xs)
                buf.write("\n")
                continue
        try:
            z=lsy[y]
        except:
            try:
                t2=d[ys.capitalize()]
                yl,yr=get_nearest_neighbors(y,ys,n,a,d,ld,ldg,cmap,cimap)
                if len(yl)!=0:
                    lsy[y]=dict(b=yl,f=yr)
            except Exception as e:               
                buf.write(y+"\t"+ys)
                buf.write("\n")
                continue
    buf.close()
    return lsy

def group_seq_by_species(df,g_to_sp):
    sph=list(df[4])
    ghsp=list(df[3])
    sp=list(df[2])
    gsp=list(df[1])
    create_dict(gsp,sp,g_to_sp)
    create_dict(ghsp,sph,g_to_sp)
    return g_to_sp

def read_gene_sequences(df,lsy,data_dir,fname):

    """The basic idea here is to create a list/dictionary of all the genes by their species.
       Once the mapping is done, all the respective fasta sequence files are read by Species
       and the CDNA sequences for each gene in the species record are read and stored.
       Thus we don't have to read the same file multiple times."""
    grouped_genes={}
    gene_by_species_dict={}
    grouped_genes=group_seq_by_species(df,grouped_genes)        
    for i in df[2].unique():
        gene_by_species_dict[i]=[]
    for i in df[4].unique():
        gene_by_species_dict[i]=[]
            
    for x in progressbar.progressbar(lsy):
        try:
            species=grouped_genes[x]#get the species
        except:
            continue
        if x not in gene_by_species_dict[species]:#check if the gene already exists in the species dict or not.
            gene_by_species_dict[species].append(x)
        xl=lsy[x]['b']
        xr=lsy[x]['f']
        for gxl in xl:
            if gxl=="NULL_GENE":
                break
            if gxl not in gene_by_species_dict[species]:
                gene_by_species_dict[species].append(gxl)
        for gxr in xr:
            if gxr=="NULL_GENE":
                break
            if gxr not in gene_by_species_dict[species]:
                gene_by_species_dict[species].append(gxr)
        

    s=[x for x in gene_by_species_dict if len(gene_by_species_dict[x])!=0]#select those species only whose gene sequences we have to read.
    s=[x.capitalize() for x in s]
    
    data=read_gene_seq(data_dir,s,gene_by_species_dict)
    not_found={}
    for species in gene_by_species_dict:
        for gene in gene_by_species_dict[species]:
            try:
                _=data[gene]
            except:
                not_found[gene]=1

    with open("processed/not_found.json","w") as file:
        json.dump(not_found,file)
    with open("processed/"+fname+".json","w") as file:#save the data
        json.dump(data,file)
    return data

def intermediate_process(gene_seq,x,y,n,index,sl,sg,ind):
    smgtemp,smltemp=create_synteny_matrix_mul(gene_seq,x,y,n)
    if np.all(smgtemp==0):
        return
    sg.append(smgtemp)
    sl.append(smltemp)
    ind.append(index)


def synteny_matrix(gene_seq,hdf,lsy,n,enable_break,sg,sl,ind):
    #sg=[]
    #sl=[]
    t=0
    #ind=[]
    
    for index,row in progressbar.progressbar(hdf.iterrows()):
        g1=str(row[1])
        g2=str(row[3])
        x=[]
        y=[]
        t+=1
        try:
            temp=lsy[g1]
        except:
            continue
        try:
            temp=lsy[g2]
        except:
            continue
        for i in range(len(lsy[g1]['b'])-1,-1,-1):
            x.append(lsy[g1]['b'][i])
        x.append(g1)
        for k in lsy[g1]['f']:
            x.append(k)

        for i in range(len(lsy[g2]['b'])-1,-1,-1):
            y.append(lsy[g2]['b'][i])
        y.append(g2)
        for k in lsy[g2]['f']:
            y.append(k)
        
        assert(len(x)==len(y))
        assert(len(x)==(2*n+1))
        """smgtemp,smltemp=create_synteny_matrix_mul(gene_seq,x,y,2*n+1)
        if np.all(smgtemp==0):
            continue
        sg.append(smgtemp)
        sl.append(smltemp)
        ind.append(index)"""
        try:
            th=Thread(target=intermediate_process,name="TimeOutDetector",args=(gene_seq,x,y,2*n+1,index,sl,sg,ind,))
            th.start()
            th.join(30)
            if th.is_alive():
                th.join()
                print(row)
        except Exception as e:
            print(e)
        if t==5 and enable_break==1:
            break
    #print("Time Taken:",end-start)
    #print("Average Time:",(end-start)/len(sg))
    print(t)
    return np.array(sg),np.array(sl),np.array(ind)


def create_branch_length_padding(bl):
    maxlen=29
    for x in bl:
        for i in range(len(x),maxlen):
            x.append(0)

def create_tree_data(treename,df):
    t=Tree(treename)
    branch_lengths_s=[]
    branch_lengths_hs=[]
    dist=[]
    ns=[]
    nhs=[]
    for index,row in progressbar.progressbar(df.iterrows()):
        d=0
        x=row[2]
        y=row[4]
        bl=[]
        c=0
        mca=t.get_common_ancestor(x,y)
        node=t&x
        while node.up!=mca:
            d+=node.dist
            bl.append(node.dist)
            node=node.up
            c+=1
        ns.append(c)
        c=0
        branch_lengths_s.append(bl)
        bl=[]
        node=t&y
        while node.up!=mca:
            d+=node.dist
            bl.append(node.dist)
            node=node.up
            c+=1
        nhs.append(c)
        branch_lengths_hs.append(bl)
        dist.append(d)
    create_branch_length_padding(branch_lengths_s)
    create_branch_length_padding(branch_lengths_hs)
    return np.array(branch_lengths_s),np.array(branch_lengths_hs),np.array(dist),np.array(ns),np.array(nhs)
