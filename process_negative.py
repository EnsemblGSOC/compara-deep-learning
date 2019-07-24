import pandas as pd 
import numpy as np 
import os 
import sys
import progressbar
import json
import sys
from neighbor_genes import read_genome_maps
from process_data import create_data_homology_ls
from threads import Procerssrunner
from read_get_gene_seq import read_gene_sequences
from access_data_rest import update_rest
from prepare_synteny_matrix import read_data_synteny
from save_data import write_dict_json

def read_database_txt(filename):
    df=pd.read_csv(filename,sep="\t",header=None)
    df=df.drop(0,axis=1)
    df.columns=["gene_stable_id","species","homology_gene_stable_id","homology_species","wga","goc","homology_type"]
    return df

def main():
    arg=sys.argv
    a,d,ld,ldg,cmap,cimap=read_genome_maps()
    print("Genome Maps Loaded.")
    df=read_database_txt(arg[-2])
    nop=int(arg[-1])
    print("Data Read.")
    a_h=[]
    d_h=[]
    a_h.append(df)
    d_h.append(arg[-2].split(".")[0])
    n=3
    lsy=create_data_homology_ls(a_h,d_h,n,a,d,ld,ldg,cmap,cimap,0)
    write_dict_json("neighbor_genes_negative","processed",lsy)
    print("Neighbor Genes Found and Saved Successfully:)")
    gene_sequences=read_gene_sequences(a_h,lsy,"geneseq","gene_seq_negative")
    gene_sequences=update_rest(gene_sequences,"gene_seq_negative")
    ndir="processed/synteny_matrices/"
    nf1="synteny_matrices_global"
    nf2="synteny_matrices_local"
    nf3="indexes"
    for i in range(len(a_h)):
        df=a_h[i]
        part=len(df)//nop
        pr=Procerssrunner()
        pr.start_processes(nop,df,gene_sequences,lsy,part,n,d_h[i])
        smg,sml,indexes=read_data_synteny(nop,d_h[i])
        print(len(indexes))
        np.save(ndir+str(d_h[i])+"_"+nf1,smg)
        np.save(ndir+str(d_h[i])+"_"+nf2,sml)
        np.save(ndir+str(d_h[i])+"_"+nf3,indexes)
    print("Synteny Matrices Created Successfully :)")

if __name__=="__main__":
    main()

