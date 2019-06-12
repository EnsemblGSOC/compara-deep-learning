import numpy as np 
import os
import json
import gc
from read_data import read_data_homology
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

a_h,d_h=read_data_homology("data_homology")
d_h=list(d_h.keys())

ndir="processed/synteny_matrices/"
nf1="_synteny_matrices_global"
nf2="_synteny_matrices_local"
nf3="_indexes"

lsy={}
with open("processed/neighbor_genes.json","r") as file:
    lsy=dict(json.load(file))

for i in range(len(d_h)):
    print("{}.{}".format(i+1,d_h[i]))

while(1):
    try:
        ch=int(input("Enter your choice:"))
        synteny_matrices_global=np.load(ndir+str(d_h[ch-1])+nf1+".npy")
        synteny_matrices_local=np.load(ndir+str(d_h[ch-1])+nf2+".npy")
        indexes=np.load(ndir+str(d_h[ch-1])+nf3+".npy")
        break
    except:
        print("Choice invalid or incomplete files!!!!!. Try Another Index.")

df=a_h[ch-1].loc[indexes]
a_h=[]
gc.collect()

inddict={}
for i in range(len(indexes)):
    inddict[indexes[i]]=i

ng=["Levenshtein Distance","Levenshtein Distance Reverse"]
nl=["Local Alignment Score","Local Alignment Score Reverse"]

print(indexes)
while(1):
    try:
        i=int(input("Enter the index:"))
    except:
        break
    if i in inddict:
        font = {'family': 'sans-serif',
        'color':  'darkturquoise',
        'weight': 'heavy',
        'size': 20,
        }
        pdf=PdfPages(str(i)+".pdf")
        print("Species",df.loc[i].species)
        print("Homology Species",df.loc[i].homology_species)
        print("Gene Stable Id:",df.loc[i].gene_stable_id)
        print("Homology Gene Stable Id:",df.loc[i].homology_gene_stable_id)
        text="Species:"+df.loc[i].species
        text+="\n"+"Gene Stable Id:"+df.loc[i].gene_stable_id
        text+="\n"+"Homology Species:"+df.loc[i].homology_species
        text+="\n"+"Homology Gene Stable Id:"+df.loc[i].homology_gene_stable_id
        text+="\n"+"Homology Type:"+df.loc[i].homology_type
        fp=plt.figure(figsize=(10,10))
        fp.text(0.5,0.5,text,ha="center",fontdict=font)
        pdf.savefig()
        plt.close()
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
            fig,ax=plt.subplots(figsize=(8,8))
            ax.set_xlabel(str(df.loc[i].homology_species),fontsize=10)
            ax.set_ylabel(str(df.loc[i].species),fontsize=10)
            hmap=sns.heatmap(matrix,xticklabels=y, yticklabels=x,annot=True,ax=ax,linewidths=.5,cmap="YlGnBu",annot_kws={"size": 10})
            hmap.figure.subplots_adjust(left=0.33,bottom=0.33,right=0.79,top=0.79)
            ax.set_title(ng[m])
            #plt.text(1,0.5,text,size=10)
            pdf.savefig()
            plt.show()

        print("Local Alignment Matrix")
        for m in range(sg.shape[-1]):
            matrix=sl[:,:,m]
            print(matrix)
            fig,ax=plt.subplots(figsize=(8,8))
            ax.set_xlabel(df.loc[i].homology_species)
            ax.set_ylabel(df.loc[i].species)
            hmap=sns.heatmap(matrix,xticklabels=y, yticklabels=x,annot=True,ax=ax,linewidths=.5,cmap="YlGnBu",annot_kws={"size": 10})
            hmap.figure.subplots_adjust(left=0.27,bottom=0.29,right=0.92)
            ax.set_title(nl[m])
            pdf.savefig()
            plt.show()
        pdf.close()
    else:
        print("Index not found")


