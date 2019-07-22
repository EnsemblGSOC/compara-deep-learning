import pandas as pd 
import progressbar
import json
import numpy as np 
import time
import pickle
from ete3 import Tree
from read_get_gene_seq import read_gene_seq,create_dict
from process_data import get_nearest_neighbors
from create_synteny_matrix import create_synteny_matrix_mul,update
import requests
import json

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

def read_gene_sequences(df,lsy,data_dir,fname,name):

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

    with open("processed/not_found_"+name+"_.json","w") as file:
        json.dump(not_found,file)
    return data

def intermediate_process(gene_seq,x,y,n,index,sl,sg,ind):
    smgtemp,smltemp=create_synteny_matrix_mul(gene_seq,x,y,n)
    if np.all(smgtemp==0):
        return
    sg.append(smgtemp)
    sl.append(smltemp)
    ind.append(index)


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

def read_data(nop,name):
    smg=[]
    sml=[]
    indexes=[]
    for i in range(nop):
        try:
            with open("temp_"+name+"/thread_"+str(i+1)+"_smg.temp","rb") as file:
                smg=smg+pickle.load(file)
            with open("temp_"+name+"/thread_"+str(i+1)+"_sml.temp","rb") as file:
                sml=sml+pickle.load(file)
            with open("temp_"+name+"/thread_"+str(i+1)+"_indexes.temp","rb") as file:
                indexes=indexes+pickle.load(file)
        except:
            continue
    print(len(indexes))
    return smg,sml,indexes

def update_rest(data,name):
    gids={}
    with open("processed/not_found_"+name+"_.json","r") as file:
        gids=dict(json.load(file))

    gids=list(gids.keys())

    geneseq={}

    server = "https://rest.ensembl.org"
    ext = "/sequence/id?type=cds"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    for i in progressbar.progressbar(range(0,len(gids)-50,50)):
        ids=dict(ids=list(gids[i:i+50]))
        while(1):
            try:
                r = requests.post(server+ext, headers=headers, data=str(json.dumps(ids)))
                if not r.ok:
                    r.raise_for_status()
                gs=r.json()
                tgs={}
                for g in gs:
                    tgs[g["query"]]=g["seq"]
                geneseq.update(tgs)
                break
            except Exception as e:
                print("Error:",e)
                continue

    data.update(geneseq)
    for genes in gids:
        try:
            _=data[genes]
        except:
            print(genes)
            update(data,genes)

    print("Gene Sequences Updated Successfully")
    return data
