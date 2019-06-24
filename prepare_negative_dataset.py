import sys
import numpy as np
import pandas as pd
import os
import gc
import progressbar
import random
import pickle
from read_data import read_data_homology,read_data_genome

arg=sys.argv
arg=arg[1:]

nos=int(arg[0])
seed=int(arg[1])
random.seed(seed)

a_h,d_h=read_data_homology("data_homology")
gmap={}
indexes_hd=[]
for df in progressbar.progressbar(a_h):
    indexes_hd.append(list(df.index.values))
    hgids=df.homology_gene_stable_id.unique()
    for h in hgids:
        gmap[h]=1

data={}
with open("genome_maps","rb") as file:
    data=pickle.load(file)
a=data["a"]
d=data["d"]
data={}
assert(len(a)==len(d))
indexes_gd=[list(df.index.values) for df in a]
d=list(d.keys())
ld=[]
for i in progressbar.progressbar(range(len(a))):
    df=a[i]
    ldg=[]
    for _,row in df.iterrows():
        gid=row.gene_id
        try:
            temp=gmap[gid]
        except:
            ldg.append(gid)
    ld.append(ldg)
assert(len(ld)==len(a))

negativesamp={}
nohd=len(a_h)
nogd=len(a)
rows=[]
for i in progressbar.progressbar(range(nos)):
    while(1):
        try:
            slh=random.randrange(nohd)#sample a homology database
            slg=random.randrange(nogd)#sample a gene annotation file
            slhd=a_h[slh]#select the homology database
            indexes=indexes_hd[slh]#select the respective indexes
            ldg=ld[slg]#select the given gene_id annotations
            ind1=random.randrange(len(ldg))#sample a gene
            g1=ldg[ind1]#get the gene id
            ind2=random.randrange(len(indexes))#sample a row
            row=slhd.loc[indexes[ind2]]#get the row from the database
            r={}
            try:
                _=negativesamp[row.gene_stable_id+g1]#check if they exist in the database
                _=negativesamp[g1+row.gene_stable_id]#check if they exist in the database
                continue
            except:
                r["gene_stable_id"]=row.gene_stable_id#add it to the row
                r["species"]=row.species#add speccies to the row
                r["homology_gene_stable_id"]=g1#add the gene_id to the row
                r["homology_species"]=d[slg]#add the gene species to the row
                rows.append(r)#add it to the rows dict
                negativesamp[row.gene_stable_id+g1]=1#add to the map so duplicate samples are avoided
                break 
        except:
            continue  

a_h=[]
gc.collect()

ndf=pd.DataFrame(rows)
print(ndf[0:10])
ndf.to_hdf("negative_dataset.h5",key="ndf",mode="w")
