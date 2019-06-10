import pandas as pd 
import numpy as np 
import os 
import sys
import progressbar
import json
from read_data import read_data_homology
from read_get_gene_seq import read_gene_sequences
from create_synteny_matrix import synteny_matrix

if not os.path.exists("processed/synteny_matrcies"):
    os.mkdir("processed/synteny_matrices")

arg=sys.argv
arg=arg[1:]

enable_break=0

if arg[-1]=="-test":
    enable_break=1

lsy={}
a_h,d_h=read_data_homology("data_homology")
print("Homology Data Read")
with open("processed/neighbor_genes.json","r") as file:
    lsy=dict(json.load(file))
print(len(lsy))
print("Neighbor Genes Loaded")
if enable_break==1:
    gene_sequences=read_gene_sequences(a_h,lsy,"geneseq","gene_sequences")
else:
    gene_sequences=read_gene_sequences(a_h[0],lsy,"geneseq","gene_sequences")
print("Gene Sequences Loaded")


if enable_break==1:
    save_after=10
else:
    save_after=250000
n=3
c=0
j=0
ndir="processed/synteny_matrices/"
nf1="synteny_matrices_global"
nf2="synteny_matrices_local"
nf3="indexes"
for df in a_h:
    j+=1
    for i in progressbar.progressbar(range(0,len(df)//save_after)):
        synteny_matrices_global,synteny_matrices_local,indexes=synteny_matrix(gene_sequences,df[i*save_after:(i+1)*save_after],lsy,n,enable_break)
        np.save(ndir+nf1+"_"+str(j)+str(c+1),synteny_matrices_global)
        np.save(ndir+nf2+"_"+str(j)+str(c+1),synteny_matrices_local)
        np.save(ndir+nf3+"_"+str(j)+str(c+1),indexes)
        c+=1
        if enable_break==1:
            break
print("Synteny Matrices Created Successfully :)")
