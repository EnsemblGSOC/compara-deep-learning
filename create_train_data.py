import pandas as pd
import numpy as np
def create_branch_length_padding(bl):
    maxlen=0
    for x in bl:
        if len(x)>maxlen:
            maxlen=len(x)

    for x in bl:
        for i in range(len(x),maxlen):
            x.append(0)

def train_data(indexes,synteny_matrices,df,branch_length_species,branch_length_homology_species,distance,dist_p_s,dist_p_hs,gene_sequences):

    homology_type_counts=dict(df.homology_type.value_counts())
    homology_species_counts=dict(df.homology_species.value_counts())

    #limits for each species and homology type in the trainig data so that the dataset is balanced
    max_species_count=10000
    max_homology_type_count=25000
    for species in homology_species_counts:
        homology_species_counts[species]=0

    for homology_type in homology_type_counts:
        homology_type_counts[homology_type]=0

    labels=dict(ortholog_one2one=0,other_paralog=1,ortholog_one2many=2,ortholog_many2many=3,within_species_paralog=4)

    train_data_dataframe=pd.DataFrame()
    train_branch_length_species=[]
    train_branch_length_homology_species=[]
    train_distance=[]
    train_dist_p_s=[]
    train_dist_p_hs=[]
    train_labels=[]
    train_indexes=[]
    train_mean_gene_length=[]
    for i in range(len(indexes)):
        row=df.loc[indexes[i]]
        if homology_species_counts[row["homology_species"]]>=max_species_count and row["homology_type"]!="within_species_paralog":
            continue
        if homology_type_counts[row["homology_type"]]>=max_homology_type_count:
            continue

        train_branch_length_species.append(branch_length_species[i])
        train_branch_length_homology_species.append(branch_length_homology_species[i])
        train_distance.append(distance[i])
        train_dist_p_s.append(dist_p_s[i])
        train_dist_p_hs.append(dist_p_hs[i])
        train_data_dataframe=train_data_dataframe.append(row)
        train_indexes.append(i)
        train_labels.append(labels[row["homology_type"]])
        homology_species_counts[row["homology_species"]]+=1
        homology_type_counts[row["homology_type"]]+=1
        train_mean_gene_length.append((len(gene_sequences[row["gene_stable_id"]])+len(gene_sequences[row["homology_gene_stable_id"]]))/2)

    train_synteny_matrices=synteny_matrices[train_indexes]

    create_branch_length_padding(train_branch_length_species)
    train_branch_length_species=np.array(train_branch_length_species)
    create_branch_length_padding(train_branch_length_homology_species)
    train_branch_length_homology_species=np.array(train_branch_length_homology_species)

    #renormalize the train_mean_gene_length by (x-mean)/std
    train_mean_gene_length=(train_mean_gene_length-np.mean(train_mean_gene_length))/np.std(train_mean_gene_length)

    #create a random array of permutations to shuffle the indices
    shi=np.random.permutation(len(train_labels))

    train_branch_length_species=train_branch_length_species[shi]
    print(train_branch_length_species.shape)

    train_branch_length_homology_species=train_branch_length_homology_species[shi]
    print(train_branch_length_homology_species.shape)

    train_dist_p_s=np.array(train_dist_p_s)
    train_dist_p_s=train_dist_p_s[shi]
    print(train_dist_p_s.shape)

    train_dist_p_hs=np.array(train_dist_p_hs)
    train_dist_p_hs=train_dist_p_hs[shi]
    print(train_dist_p_hs.shape)

    train_synteny_matrices=train_synteny_matrices[shi]
    print(train_synteny_matrices.shape)

    train_indexes=np.array(train_indexes)
    train_indexes=train_indexes[shi]
    print(train_indexes.shape)

    train_labels=np.array(train_labels)
    train_labels=train_labels[shi]
    print(train_labels.shape)

    train_mean_gene_length=np.array(train_mean_gene_length)
    train_mean_gene_length=train_mean_gene_length[shi]
    print(train_mean_gene_length.shape)

    train_distance=np.array(train_distance)
    train_distance=train_distance[shi]
    train_distance=(train_distance-np.mean(train_distance))/np.std(train_distance)
    train_distance.shape

    return train_synteny_matrices,train_branch_length_species,train_branch_length_homology_species,train_mean_gene_length,train_dist_p_s,train_dist_p_hs,train_distance,train_labels
