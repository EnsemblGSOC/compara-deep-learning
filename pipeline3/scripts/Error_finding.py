import numpy as np
import pandas as pd
import json

Na_dict = {}
Void_dict = {}

species_list = np.loadtxt("config/species_list2.txt", dtype="str")
neighour_json_path = "/nfs/production/flicek/ensembl/compara/amarshall/Data_Storage/Species_Neighbours"

for species in species_list:
    # read in the samples for both the species and the candidate homologies
    with open( neighour_json_path + "/" + species + "_all_neighbours.json","r") as f:
        sample_neighbours = json.load(f)

    # with open( neighour_json_path + "/" + species + "_homology_sample_neighbour.json","r") as f:
    #     sample_homology_neighbours = json.load(f)
    
    print(species)

    gene_series = pd.Series(sample_neighbours.keys()).map(sample_neighbours)
    gene_series[gene_series.isna()] = gene_series[gene_series.isna()].apply(lambda x: [["Void","Void","Void"],["Void","Void","Void"]]) # fill the null values
    array1 = np.concatenate(gene_series.apply( np.array ).values ).reshape((-1,6))
    array1 = np.concatenate((pd.Series(sample_neighbours.keys()).values[:,np.newaxis], array1), axis=1) #.reshape(1,-1)

    Na_prop1 = float(((array1 == "NA").sum() / array1.ravel().shape))
    Void_prop1 = float(((array1 == "Void").sum() / array1.shape[0]))

    # gene_series = pd.Series(sample_homology_neighbours.keys()).map(sample_homology_neighbours)
    # gene_series[gene_series.isna()] = gene_series[gene_series.isna()].apply(lambda x: [["Void","Void","Void"],["Void","Void","Void"]]) # fill the null values
    # array2 = np.concatenate(gene_series.apply( np.array ).values ).reshape((-1,6))
    # array2 = np.concatenate((pd.Series(sample_homology_neighbours.keys()).values[:,np.newaxis], array2), axis=1) #.ravel()[np.newaxis,:]

    # Na_prop2 = float(((array2 == "NA").sum() / array2.ravel().shape))
    # Void_prop2 = float(((array2 == "Void").sum() / array2.shape[0]))

    Na_dict[species] = Na_prop1#,Na_prop2]
    Void_dict[species] = Void_prop1#,Void_prop2]

with open("./Na.json", "w") as f:
    json.dump(Na_dict,f)
    
with open("./Void.json", "w") as f:
    json.dump(Na_dict,f)