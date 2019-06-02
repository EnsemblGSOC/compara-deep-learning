import numpy as np
import requests
import edlib as ed
import pandas as pd
import time
import sys
from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA,TabularMSA,RNA

def update(gene_seq,gene):
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/"+str(gene)+"?"

    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    gene_seq[gene]=str(r.text)

def create_synteny_matrix_mul(gene_seq,g1,g2,n):
    for gene in g1:
        if gene=="NULL_GENE":
            continue
        try:
            temp=gene_seq[gene]
        except:
            print("Updating gene sequences for gene:",gene)
            update(gene_seq,gene)
    for gene in g2:
        if gene=="NULL_GENE":
            continue
        try:
            temp=gene_seq[gene]
        except:
            print("Updating gene sequences for gene:",gene)
            update(gene_seq,gene)
    #print(n)
    sm=np.zeros((n,n,2))
    sml=np.zeros((n,n,2))
    for i in range(n):
        if g1[i]=="NULL_GENE":
            continue
        for j in range(n):
            if g2[j]=="NULL_GENE":
                continue
            norm_len=(len(gene_seq[g1[i]])+len(gene_seq[g2[j]]))
            result = ed.align(gene_seq[g1[i]],gene_seq[g2[j]], mode="NW", task="distance")
            sm[i][j][0]=result["editDistance"]/(norm_len)
            result = ed.align(gene_seq[g1[i]],gene_seq[g2[j]][::-1], mode="NW", task="distance")
            sm[i][j][1]=result["editDistance"]/(norm_len)
            _,result,_=local_pairwise_align_ssw(DNA(gene_seq[g1[i]]),DNA(gene_seq[g2[j]]))
            sml[i][j][0]=result/(norm_len)
            _,result,_=local_pairwise_align_ssw(DNA(gene_seq[g1[i]]),DNA(gene_seq[g2[j]][::-1]))
            sml[i][j][1]=result/(norm_len)
    return sm,sml

def synteny_matrix(gene_seq,hdf,lsy,n):
    sg=[]
    sl=[]
    t=0
    ind=[]
    start=time.time()
    for index,row in hdf.iterrows():
        g1=str(row["gene_stable_id"])
        g2=str(row["homology_gene_stable_id"])
        x=[]
        y=[]
        try:
            temp=lsy[g1]
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
        smgtemp,smltemp=create_synteny_matrix_mul(gene_seq,x,y,2*n+1)
        sg.append(smgtemp)
        sl.append(smltemp)
        ind.append(index)
        t+=1
        if t==5:
            break
    end=time.time()
    print("Time Taken:",end-start)
    print("Average Time:",(end-start)/len(sg))
    return np.array(sg),np.array(sl),np.array(ind)
