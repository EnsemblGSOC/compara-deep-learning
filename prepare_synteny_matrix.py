import pandas as pd 
import numpy as np 
import os 
import sys
import progressbar
import json
from read_data import read_data_homology
from read_get_gene_seq import read_gene_sequences
from create_synteny_matrix import synteny_matrix

if not os.path.isdir("processed/synteny_matrices"):
    os.mkdir("processed/synteny_matrices")

arg=sys.argv
arg=arg[1:]

nos=int(arg[0])

lsy={}
a_h,d_h=read_data_homology("data_homology")
print(d_h)
d_h=list(d_h.keys())
print("Homology Data Read")
for i in range(len(a_h)):
    df=a_h[i]
    random_indexes=np.random.permutation(len(df))
    random_indexes=random_indexes[:nos]
    df=df.loc[random_indexes]
    assert(len(df)==nos)
    a_h[i]=df

with open("processed/neighbor_genes.json","r") as file:
    lsy=dict(json.load(file))
print(len(lsy))
print("Neighbor Genes Loaded")


gene_sequences=read_gene_sequences(a_h,lsy,"geneseq","gene_sequences")
print("Gene Sequences Loaded")


n=3
ndir="processed/synteny_matrices/"
nf1="synteny_matrices_global"
nf2="synteny_matrices_local"
nf3="indexes"
for i in range(len(a_h)):
    df=a_h[i]
    synteny_matrices_global,synteny_matrices_local,indexes=synteny_matrix(gene_sequences,df,lsy,n,0)
    np.save(ndir+str(d_h[i])+"_"+nf1,synteny_matrices_global)
    np.save(ndir+str(d_h[i])+"_"+nf2,synteny_matrices_local)
    np.save(ndir+str(d_h[i])+"_"+nf3,indexes)
        
print("Synteny Matrices Created Successfully :)")
