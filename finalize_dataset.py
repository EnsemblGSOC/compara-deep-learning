import pandas as pd
import numpy as np 
import json
import gc
import pickle
import os
import sys
from tree_data import create_tree_data
from process_negative import read_database_txt
from select_data import read_db_homology

def read_data_homology(dirname,nfname):
    lf=os.listdir(dirname)
    if len(lf)==0:
        print("No Files in the Directory!!!!!!!")
        sys.exit(1)
    a_h=[]
    d_h=[]
    for x in lf:
        df,n=read_db_homology(dirname,x)
        n=n.split()[0]
        try:
            indexes=np.load("processed/synteny_matrices/"+n+"_indexes.npy")
        except:
            print("Incomplete data for:",n)
        df=df.loc[indexes]
        a_h.append(df)
        d_h.append(n)
    #read the negative dataset
    df=read_database_txt(nfname)
    indexes=np.load("processed/synteny_matrices/"+nfname.split(".")[0]+"_indexes.npy")
    df=df.loc[indexes]
    a_h.append(df)
    d_h.append(nfname.split(".")[0])
    return a_h,d_h

def prepare_features(a_h,d_h,sptree,label):
    rows=[]
    smg_name="_synteny_matrices_global.npy"
    sml_name="_synteny_matrices_local.npy"
    smi_name="_indexes.npy"
    dir_name="processed/synteny_matrices/"
    for i in range(len(a_h)):
        df=a_h[i]
        n=d_h[i]
        try:
            smg=np.load(dir_name+n+smg_name)
            sml=np.load(dir_name+n+sml_name)
            indexes=np.load(dir_name+n+smi_name)
        except:
            print("Incomplete data for:",n)
            continue
        df=df.loc[indexes]

        branch_length_species,branch_length_homology_species,distance,dist_p_s,dist_p_hs=create_tree_data(sptree,df)
        assert(len(branch_length_species)==len(df))
        assert(len(sml)==len(distance))

        for i in range(len(df)):
            index=indexes[i]
            row=df.loc[index]
            r={}
            r["species"]=row["species"]
            r["homology_species"]=row["homology_species"]
            r["gene_stable_id"]=row["gene_stable_id"]
            r["homology_gene_stable_id"]=row["homology_gene_stable_id"]
            r["label"]=label[row["homology_type"]]
            r["global_alignment_matrix"]=smg[i]
            r["local_alignment_matrix"]=sml[i]
            r["index_homology_dataset"]=index
            r["bls"]=branch_length_species[i]
            r["blhs"]=branch_length_homology_species[i]
            r["dis"]=distance[i]
            r["dps"]=dist_p_s[i]
            r["dphs"]=dist_p_hs[i]
            rows.append(r)
    return rows

def main():
    arg=sys.argv
    nfname=arg[-1]
    a_h,d_h=read_data_homology("data_homology",nfname)
    labels=dict(ortholog_one2one=1,
                other_paralog=0,
                non_homolog=2,
                ortholog_one2many=1,
                ortholog_many2many=1,
                within_species_paralog=0,
                gene_split=4)
    rows=prepare_features(a_h,d_h,"species_tree.tree",labels)
    with open("dataset","wb") as file:
        pickle.dump(rows,file)
    print("Dataset_Finalized")

if __name__=="__main__":
    main()
