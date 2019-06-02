import sys
import numpy as np
import pandas as pd
import json
import os
from get_data import get_data_homology,get_data_genome
from process_data import create_data_homology_ls
from read_get_gene_seq import read_gene_sequences
from create_synteny_matrix import synteny_matrix
from tree_data import create_tree_data
from prepare_train_data import train_data

if not os.path.exists("processed"):
    os.mkdir("processed")

arg=sys.argv
arg=arg[1:]

enable_break=0

if arg[-1]=="-test":
    enable_break=1

arg=arg[:-1]

if len(arg)!=5:
    print("No. of arguments more or less. Please check")
    sys.exit(1)

dir_g="data"
ld,ldg,a,d=get_data_genome(arg,dir_g)

dir_hom="data_homology"
a_h,d_h=get_data_homology(arg,dir_hom)

if arg[-1]=="-d":
    sys.exit(1)

print("Data Read")

n=3 #no. of numbers neighbors
save_after=0 #to save data after n steps

if enable_break==1:
    save_after=3
else:
    save_after=1000000

lsy=create_data_homology_ls(a_h,d_h,n,a,d,ld,ldg,save_after,enable_break)
print(len(lsy))

print("Neighbor Genes Found")

if enable_break==1:
    gene_sequences=read_gene_sequences(a_h[0][0:10],lsy,"geneseq","gene_sequences")
else:    
    gene_sequences=read_gene_sequences(a_h[0],lsy,"geneseq","gene_sequences")

if enable_break==1:
    synteny_matrices_global,synteny_matrices_local,indexes=synteny_matrix(gene_sequences,a_h[0][0:10],lsy,n)
else:    
    synteny_matrices_global,synteny_matrices_local,indexes=synteny_matrix(gene_sequences,a_h[0],lsy,n)

print("Synteny Matrices are created successfully\n",len(indexes),"\n",len(synteny_matrices_global))
np.save("processed/synteny_matrices_global",synteny_matrices_global)
np.save("processed/synteny_matrices_local",synteny_matrices_local)
np.save("processed/indexes",indexes)

with open("processed/gene_seq_updated.json","w") as file:
    json.dump(gene_sequences,file)

df=a_h[0].loc[indexes]
branch_length_species,branch_length_homology_species,distance,dist_p_s,dist_p_hs=create_tree_data("species_tree.tree",df)

train_synteny_matrices_global,train_synteny_matrices_local,train_branch_length_species,train_branch_length_homology_species,train_mean_gene_length,train_dist_p_s,train_dist_p_hs,train_distance,train_labels=train_data(indexes,synteny_matrices_global,synteny_matrices_local,df,branch_length_species,branch_length_homology_species,distance,dist_p_s,dist_p_hs,gene_sequences)

np.save("processed/train_synteny_matrices_global",train_synteny_matrices_global)
np.save("processed/train_synteny_matrices_local",train_synteny_matrices_local)
np.save("processed/train_branch_length_species",train_branch_length_species)
np.save("processed/train_branch_length_homology_species",train_branch_length_homology_species)
np.save("processed/train_mean_gene_length",train_mean_gene_length)
np.save("processed/train_dist_p_s",train_dist_p_s)
np.save("processed/train_dist_p_hs",train_dist_p_hs)
np.save("processed/train_distance",train_distance)
np.save("processed/train_labels",train_labels)

print("Data Saved Successfully to processed :)")
