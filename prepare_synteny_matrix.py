import pandas as pd 
import numpy as np 
import os 
import sys
import progressbar
import json
import gc
from read_data import read_data_homology
from read_get_gene_seq import read_gene_sequences
from create_synteny_matrix import synteny_matrix
from access_data_rest import update_rest
from select import select,create_map_reverse

if not os.path.isdir("processed/synteny_matrices"):
    os.mkdir("processed/synteny_matrices")

arg=sys.argv
arg=arg[1:]

nos=int(arg[0])

with open("dist_matrix","r") as file:
    matrix=file.readlines()
matrix=[x.split("\t") for x in matrix]
matrix=[[float(y) for y in x] for x in matrix]
matrix=np.array(matrix)
with open("sp_names","r") as file:
    dname=file.readlines()
dname=[x.split("\n")[0] for x in dname]
spnmap,nspmap=create_map_reverse(dname)

lsy={}
a_h,d_h=read_data_homology("data_homology")
print(d_h)
d_h=list(d_h.keys())
d_h=[x.split()[0] for x in d_h]
print("Homology Data Read")
for i in range(len(a_h)):
    a_h[i]=select(a_h[i],nos,matrix,spnmap,nspmap,d_h[i])
    assert(len(a_h[i])==nos)
gc.collect()
print("Data Selected")

with open("processed/neighbor_genes.json","r") as file:
    lsy=dict(json.load(file))
print(len(lsy))
print("Neighbor Genes Loaded")

gene_sequences={}
gene_sequences=read_gene_sequences(a_h,lsy,"geneseq","gene_sequences")
print("Gene Sequences Loaded")

print("Going to update not found sequences:")
gene_sequences=update_rest(gene_sequences)

n=3
ndir="processed/synteny_matrices/"
nf1="synteny_matrices_global"
nf2="synteny_matrices_local"
nf3="indexes"
for i in range(len(a_h)):
    df=a_h[i]
    print(len(df))
    synteny_matrices_global,synteny_matrices_local,indexes=synteny_matrix(gene_sequences,df,lsy,n,0)
    np.save(ndir+str(d_h[i])+"_"+nf1,synteny_matrices_global)
    np.save(ndir+str(d_h[i])+"_"+nf2,synteny_matrices_local)
    np.save(ndir+str(d_h[i])+"_"+nf3,indexes)
    print(len(indexes))    
print("Synteny Matrices Created Successfully :)")
