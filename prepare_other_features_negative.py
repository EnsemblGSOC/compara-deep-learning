import pandas as pd
import numpy as np 
import json
import gc
import pickle
import progressbar
import sys
from tree_data import create_tree_data



smg_name="negative_dataset_synteny_matrices_global"
sml_name="negative_dataset_synteny_matrices_local"
smi_name="negative_dataset_indexes"
dir_name="processed/synteny_matrices/"
label=2

rows=[]
df=pd.read_hdf("negative_dataset.h5",key="ndf")
try:
    smg=np.load(dir_name+smg_name+".npy")
    sml=np.load(dir_name+sml_name+".npy")
    indexes=np.load(dir_name+smi_name+".npy")
except:
    print("Incomplete data")
    sys.exit(1)
df=df.loc[indexes]
for _,row in progressbar.progressbar(df.iterrows()):
    row["homology_species"]=row["homology_species"].lower()
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
    r["label"]=label
    r["global_alignment_matrix"]=smg[i]
    r["local_alignment_matrix"]=sml[i]
    r["index_homology_dataset"]=index
    r["bls"]=branch_length_species[i]
    r["blhs"]=branch_length_homology_species[i]
    r["dis"]=distance[i]
    r["dps"]=dist_p_s[i]
    r["dphs"]=dist_p_hs[i]
    rows.append(r)

with open("negative_dataset","wb") as file:
    pickle.dump(rows,file)

print("Data Saved Successfully:)")