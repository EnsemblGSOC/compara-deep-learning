import pandas as pd
import numpy as np 
import json
import gc
import pickle
from read_data import read_data_homology
from tree_data import create_tree_data

a_h,d_h=read_data_homology("data_homology")
d_h=list(d_h.keys())

smg_name="_synteny_matrices_global.npy"
sml_name="_synteny_matrices_local.npy"
smi_name="_indexes.npy"
dir_name="processed/synteny_matrices/"
label=dict(ortholog_one2one=0,other_paralog=1,ortholog_one2many=1,ortholog_many2many=1,within_species_paralog=0)

rows=[]
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
    branch_length_species,branch_length_homology_species,distance,dist_p_s,dist_p_hs=create_tree_data("species_tree.tree",df)
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

with open("dataset","wb") as file:
    pickle.dump(rows,file)

print("Data Saved Successfully:)")