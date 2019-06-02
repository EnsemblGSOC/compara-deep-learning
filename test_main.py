import sys
import numpy as np
import pandas as pd
import json
from get_data import get_data_homology,get_data_genome
from process_data import create_data_homology_ls
from read_get_gene_seq import read_gene_sequences
from create_synteny_matrix import synteny_matrix
from tree_data import create_tree_data
from create_train_data import train_data
from train import train

arg=sys.argv
arg=arg[1:]

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

n=2 #no. of numbers neighbors
save_after=3 #to save data after n steps
lsy=create_data_homology_ls(a_h,d_h,n,a,d,ld,ldg,save_after)
print(len(lsy))

print("Neighbor Genes Found")

gene_sequences=read_gene_sequences(a_h[0],lsy,"geneseq","gene_sequences")

synteny_matrices,indexes=synteny_matrix(gene_sequences,a_h[0][0:1000],lsy,n)
print("Synteny Matrices are created successfully\n",len(indexes),"\n",len(synteny_matrices))
np.save("synteny_matrices",synteny_matrices)
np.save("indexes",indexes)

with open("gene_seq_updated.json","w") as file:
    json.dump(gene_sequences,file)

df=a_h[0].loc[indexes]
branch_length_species,branch_length_homology_species,distance,dist_p_s,dist_p_hs=create_tree_data("species_tree.tree",df)

train_synteny_matrices,train_branch_length_species,train_branch_length_homology_species,train_mean_gene_length,train_dist_p_s,train_dist_p_hs,train_distance,train_labels=train_data(indexes,synteny_matrices,df,branch_length_species,branch_length_homology_species,distance,dist_p_s,dist_p_hs,gene_sequences)

train(train_synteny_matrices,train_branch_length_species,train_branch_length_homology_species,train_mean_gene_length,
            train_dist_p_s,train_dist_p_hs,train_distance,train_labels)
