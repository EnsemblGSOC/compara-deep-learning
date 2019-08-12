import json
import gc
import pandas as pd 
import numpy as np 
import pickle
import sys
import progressbar
import os
from neighbor_genes import read_genome_maps
from process_data import create_data_homology_ls
from read_get_gene_seq import read_gene_sequences
from access_data_rest import update_rest,update_rest_protein
from prepare_synteny_matrix import write_fasta
from process_data import create_map_list

def read_database(fname,dirname):
    df=pd.read_csv(dirname+"/"+fname,sep="\t",header=None)
    label_dict=dict(ortholog_one2one=1,
                    other_paralog=0,
                    non_homolog=2,
                    ortholog_one2many=1,
                    ortholog_many2many=1,
                    within_species_paralog=0,
                    gene_split=4)
    label=[]
    for _,row in df.iterrows():
        label.append(label_dict[row[7]])
    df=df.assign(label=label)
    df=df.drop(7,axis=1)
    df=df.drop(0,axis=1)
    df.columns=["gene_stable_id","species","homology_gene_stable_id","homology_species","goc","wga","label"]
    return df

def read_prediction_file_folder(dir_name):
    lf=os.listdir(dir_name)
    a_h=[]
    d_h=[]
    for x in progressbar.progressbar(lf):
        df=read_database(x,dir_name)
        a_h.append(df)
        d_h.append(x.split(".")[0])
    return a_h,d_h

def create_synteny_features(a_h,d_h,n,a,d,ld,ldg,cmap,cimap,name):    
    lsy=create_data_homology_ls(a_h,d_h,n,a,d,ld,ldg,cmap,cimap,0)
    protein_sequences=read_gene_sequences(a_h,lsy,"pro_seq","prediction_"+name)
    protein_sequences=update_rest_protein(protein_sequences,"prediction_"+name)
    write_fasta(protein_sequences,"prediction_"+name)
    print("Protein Sequences Loaded")

def main():
    arg=sys.argv
    dirname=arg[-1]
    a_h,d_h=read_prediction_file_folder(dirname)
    n=3

    a,d,ld,ldg,cmap,cimap=read_genome_maps()#read the genome maps
    print("Genome Maps Loaded.")

    create_synteny_features(a_h,d_h,n,a,d,ld,ldg,cmap,cimap,dirname)
if __name__=="__main__":
    main()