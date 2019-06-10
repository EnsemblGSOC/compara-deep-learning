import numpy as np 
import os
import json
from read_data import read_data_homology
import matplotlib.pyplot as plt
import seaborn as sns

a_h,_=read_data_homology("data_homology")

ndir="processed/synteny_matrices/"
nf1="synteny_matrices_global_"
nf2="synteny_matrices_local_"
nf3="indexes_"

lsy={}
with open("processed/neighbor_genes.json","r") as file:
    lsy=dict(json.load(file))

synteny_matrices_global=np.load(ndir+nf1+"11"+".npy")
synteny_matrices_local=np.load(ndir+nf2+"11"+".npy")
indexes=np.load(ndir+nf3+"11"+".npy")

df=a_h[0].loc[indexes]

inddict={}
for i in range(len(indexes)):
    inddict[indexes[i]]=i

while(1):
    i=int(input("Enter the index"))
    if i in inddict:
        print("Species",df.loc[i].species)
        print("Homology Species",df.loc[i].homology_species)
        print("Gene Stable Id:",df.loc[i].gene_stable_id)
        print("Homology Gene Stable Id:",df.loc[i].homology_gene_stable_id)
        print("Global aligned matrix:")

        g1=df.loc[i].gene_stable_id
        g2=df.loc[i].homology_gene_stable_id
        x=[]
        y=[]
        for n in range(len(lsy[g1]['b'])-1,-1,-1):
            x.append(lsy[g1]['b'][n])
        x.append(g1)
        for k in lsy[g1]['f']:
            x.append(k)

        for n in range(len(lsy[g2]['b'])-1,-1,-1):
            y.append(lsy[g2]['b'][n])
        y.append(g2)
        for k in lsy[g2]['f']:
            y.append(k)

        loc=inddict[i]
        sg=synteny_matrices_global[loc]
        sl=synteny_matrices_local[loc]
        for m in range(sg.shape[-1]):
            matrix=sg[:,:,m]
            print(matrix)
            hmap=sns.heatmap(matrix,xticklabels=y, yticklabels=x,annot=True)
            plt.show()

        print("Local Alignment Matrix")
        for m in range(sg.shape[-1]):
            matrix=sl[:,:,m]
            print(matrix)
            hmap=sns.heatmap(matrix,xticklabels=y, yticklabels=x,annot=True)
            plt.show()

            

    else:
        print("Index not found")


