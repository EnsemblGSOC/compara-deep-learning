import sys
import numpy as np
import pandas as pd
import os
import gc
import progressbar
import random
from read_data import read_data_homology,read_data_genome

arg=sys.argv
arg=arg[1:]

nos=int(arg[0])
seed=int(arg[1])
random.seed(seed)

a_h,d_h=read_data_homology("data_homology")
a=[]
d={}
a,d=read_data_genome("data",a,d)
assert(len(a)==len(d))
indexes_gd=[list(df.index.values) for df in a]
d=list(d.keys())

gmap={}
indexes_hd=[]
for df in progressbar.progressbar(a_h):
    indexes_hd.append(list(df.index.values))
    hgids=df.homology_gene_stable_id.unique()
    for h in hgids:
        gmap[h]=1



col_names=["gid","species","hgid","h_species"]
nohd=len(a_h)
nogd=len(a)
rows=[]
for i in progressbar.progressbar(range(nos)):
    while(1):
        try:
            slh=random.randrange(nohd)
            slg=random.randrange(nogd)
            slgd=a[slg]
            n1=d[slg]
            indexes=indexes_gd[slg]
            ind=random.randrange(len(indexes))
            try:
                _=gmap[slgd.loc[indexes[ind]].gene_id]
                continue
            except:
                slhd=a_h[slh]
                lid=indexes_hd[slh]
                ind_1=random.randrange(len(slhd))
                row=slhd.loc[lid[ind_1]]
                r={}
                r["gene_stable_id"]=row.gene_stable_id 
                r["species"]=row.species
                r["homology_gene_stable_id"]=slgd.loc[indexes[ind]].gene_id
                r["homology_species"]=n1  
                rows.append(r) 
                break 
        except:
            continue  

a_h=[]
gc.collect()

ndf=pd.DataFrame(rows)
print(ndf[0:10])
ndf.to_hdf("negative_dataset.h5",key="ndf",mode="w")
