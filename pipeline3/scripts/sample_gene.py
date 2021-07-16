import pandas as pd
import numpy as np
import argparse
import os
import json

np.random.seed(42)

parser = argparse.ArgumentParser()
parser.add_argument("--homo", help="homology database path")
parser.add_argument(
    "--homo_paths",
    help="homology database paths file",
    default="config/homology_database_paths.txt",
)
parser.add_argument("--out", help="output directory")
parser.add_argument("--species", help="species name")
parser.add_argument("--n_samples", type=int, help="The number of samples to get", default=100)
parser.add_argument("--n_species", type = int, help="Number of species", default=10)
parser.add_argument(
    "--avail_species", type = str,
    help="Text file with list of available species",
    default="config/database_species_intersection.txt ",
)
parser.add_argument(
    "--sp_names", help="Path to sp name with dist matrix", default="scripts/sp_names"
)
parser.add_argument(
    "--dist_matrix", help="Path to distance matrix", default="scripts/dist_matrix"
)
parser.add_argument(
    "--neighbours",
    help="directory containing all neighbour pairs for each species",
)
args = parser.parse_args()
print(args)


# read in the distance matrix
sp_names = np.loadtxt(args.sp_names, dtype="str")
distance = np.loadtxt(args.dist_matrix, delimiter="\t")
distance_df = pd.DataFrame(distance, index=sp_names, columns=sp_names)

# read in the homology database
homology = pd.read_csv(args.homo, compression="gzip", delimiter="\t")

# read in the available species
print(args.avail_species, os.getcwd())
available_species = np.loadtxt("config/species_list.txt", dtype="str")

# Select samples only from available lists
intersection_species = list(set(sp_names) & set(homology.homology_species.drop_duplicates()) & set(available_species))

# reduce the species distance matrix to consider only species with available data
distance_df = distance_df[intersection_species].loc[intersection_species]

# reduce the homology database to only contain those species which are available
homology = homology[ homology.homology_species.isin(intersection_species)]


"""
Selection rule for sampling from different species based on distance.
Here, pick target species to  obtain homologous gene pairs from 
using equally spaced ranks from the distance matrix.

If dissatisfied, right your own rule to sample target species
"""
# look at species pairs based on the rank of their distance value
# get equally spaced rank increments
distance_increment = int(distance_df.shape[0] / args.n_species)
# sometimes there are less species than the desired number
if distance_df.shape[0] <= args.n_species:
    distance_increment = 1
print("Distance increment: ")
print(distance_increment)
# get the rank numbers
print("Number of species in the distance matrix")
print(distance_df.shape[0])
distance_ranks = np.arange(1, int(distance_df.shape[0]), distance_increment)
# get the species names
query_species = distance_df[args.species].sort_values()[distance_ranks].index.values

"""
Read in the homology database of the species we're running the script. 
We will then get samples from this database based on the targets above.
Eg, we read in human database, and get homologous pairs of genes beween 
genes and the species we selected above
"""


# map the homology types to the simpler ortho, para, gene split classes
homo_type_map = {
    i: j
    for i, j in zip(
        homology.homology_type.drop_duplicates(),
        ["ortholog", "ortholog", "ortholog", "paralog", "paralog", "gene_split"],
    )
}
homology["ortho_para"] = homology.homology_type.map(homo_type_map)
# sample equal numbers of orthologs and paralogs for each of the queary species
print(query_species)
num_ortho_sub_samples = int(args.n_samples / (args.n_species * 2)) # factor 1/2 as half data should be ortholog.
ortho_samples_df = (
    homology[homology.homology_species.isin(query_species) & (homology.ortho_para == "ortholog")]
    .groupby("homology_species")
    .sample(num_ortho_sub_samples, replace=True).drop_duplicates()
)

print(ortho_samples_df.shape[0])

num_para_samples = ortho_samples_df.shape[0] # get the same number of paralogs as orthologs
para_samples_df = (
    homology[homology.ortho_para == "paralog"].sample(num_para_samples, replace=True).drop_duplicates()
)

samples_df = pd.concat([ortho_samples_df,para_samples_df])


samples_df.to_csv(args.out + "/" + args.species + "_samples_pairs.csv")

"""
Now we get the neighbours of all the genes using the neighbours.json
file generated earlier in the pipeline
"""

# load the species map

print("loading neighbourhood genes")
with open(args.neighbours + "/" + args.species + "_all_neighbours.json") as f:
    neighbour_map = json.load(f)

neighbour_genes = dict(
    zip(samples_df.gene_stable_id, samples_df.gene_stable_id.map(neighbour_map))
)
print("writing pickle file: " + args.out + "/" + args.species + "_sample_neighbour.txt")
with open(args.out + "/" + args.species + "_sample_neighbour.json", "w") as json_file:
    json.dump(neighbour_genes, json_file)

# read in the list of paths to homology species
homo_paths = pd.Series(np.loadtxt(args.homo_paths, dtype="str"))

# function to get path to homology database of species
def homo_path_map(species):
    return homo_paths[homo_paths.str.contains(species)].values[0]


# loop through the homolog species to generate neighbours
super_homology_neighbour_map = {}
negs_main_species = []
negs_query_species = []
for species in pd.Series(query_species).append(pd.Series(args.species)):
    print("loading homology neighbourhood genes")
    with open(args.neighbours + "/" + species + "_all_neighbours.json") as f:
        homology_neighbour_map = json.load(f)
    super_homology_neighbour_map.update(homology_neighbour_map)
    # select negative samples
    """
    this part of the for loop deals with obtaining negative samples using negative log.
    First sample the genes from the main species. Then for those genes not in the homolgy database, 
    sample from the query species randomly because there's no risk of overlap.
    For the genes sampled from the main species that are in the homology database, use funky logic
    to ensure you don't pick homologs from the query species
    """
    temp_df = homology[homology.homology_species == species] # only look at the current species
    query_series = pd.Series(homology_neighbour_map.keys()) # get all genes from the query species
    negative_samples = pd.Series(list(neighbour_map.keys())).sample(int(num_ortho_sub_samples)) # sample them
    temp = negative_samples.to_frame("gene_stable_id")
    temp["species"] = species
    negs_main_species.append(temp)
    non_homo_samples_array = np.array(negative_samples) # initialise array for non-homologous pairs
    index = ~negative_samples.isin(temp_df.gene_stable_id) # indices of samples not in the homology db
    non_homo_samples_array[index] = query_series.sample(index.sum(), replace=True) # sample any genes for those indices

    def sample(gene):
        """
        Helper function to change the dataframe
        """
        bad_genes = temp_df[temp_df.gene_stable_id == gene].groupby("gene_stable_id")["homology_gene_stable_id"].apply(list)
        return list(query_series[~query_series.isin(bad_genes)].sample(1))[0]

    non_homo_samples_array[~index] = negative_samples[~index].apply(sample) # get samples for the rest of the positions
    temp = pd.DataFrame(non_homo_samples_array,columns=["non_homo_gene_stable_id"])
    temp["non_homo_species"] = species
    negs_query_species.append(temp)


# save the positive example homologous pairs
homology_neighbour_genes = dict(
    zip(samples_df.homology_gene_stable_id, samples_df.homology_gene_stable_id.map(super_homology_neighbour_map))
)

with open(
    args.out + "/" + args.species + "_homology_sample_neighbour.json", "w"
) as json_file:
    json.dump(homology_neighbour_genes, json_file)


# Deal with the negative examples
negs_main_species = pd.concat(negs_main_species)
negs_query_species = pd.concat(negs_query_species)

pd.concat([negs_main_species.reset_index(), negs_query_species.reset_index()], axis=1).to_csv(args.out + "/" + args.species + "_negative_samples_pairs.csv")


negative_neighbour_genes = dict(
    zip(negs_main_species.gene_stable_id, negs_main_species.gene_stable_id.map(neighbour_map))
)
with open(args.out + "/" + args.species + "_negative_sample_neighbour.json", "w") as json_file:
    json.dump(negative_neighbour_genes, json_file)

negative_homology_neighbour_genes = dict(
    zip(negs_query_species.non_homo_gene_stable_id, negs_query_species.non_homo_gene_stable_id.map(super_homology_neighbour_map))
)

with open(
    args.out + "/" + args.species + "_negative_homology_sample_neighbour.json", "w"
) as json_file:
    json.dump(negative_homology_neighbour_genes, json_file)


    
"""
ZOMBIE CODE. IGNORE
"""

# homology_neighbour_dict = {}
# for species in query_species:
#     # get only samples from that species
#     species_samples = samples_df[samples_df.homology_species == species]
#     neighbour_genes = dict(zip(species_samples.gene_stable_id, species_samples.gene_stable_id.map(neighbour_map)))
#     species_neighbour_dict.update(neighbour_genes)
#     # read in the query species neighbour list
#     with open(args.neighbours + "/" + species + "_all_neighbours.json") as f:
#         homology_neighbour_map = json.load(f)
#     # get the dictionary of each homolog gene and it's neighbours
#     homology_neighbour_genes = dict(
#         zip(
#             species_samples.homology_gene_stable_id,
#             species_samples.homology_gene_stable_id.map(homology_neighbour_map),
#         )
#     )
#     print("species: " + species)

#     # append to the homology_neighbours dict
#     homology_neighbour_dict.update(homology_neighbour_genes)

# with open(args.out + "/" + args.species + "_sample_neighbour.json", "w") as json_file:
#     json.dump(species_neighbour_dict, json_file)


# with open(
#     args.out + "/" + args.species + "_homology_sample_neighbour.json", "w"
# ) as json_file:
#     json.dump(homology_neighbour_dict, json_file)

# # get negative samples
# for species in query_species:
#     with open(args.neighbours + "/" + species + "_all_neighbours.json") as f:
#         homology_neighbour_map = json.load(f)
#     negative_samples = pd.Series(list(neighbour_map.keys())).sample(num_sub_samples) # select genes
#     negative_samples_array = np.array(negative_samples) # bad initialisation should be changed, to be more readable
#     query_series = pd.Series(negative_neighbour_map.keys())
#     negative_samples[~negative_samples.isin(homo.gene_stable_id)] = query_series.sample((~negative_samples.isin(homo.gene_stable_id))).values
#         def get_samples(gene_list):
#         non_list = pd.Series(gene_list)
#         return query_series[~query_series.isin(non_list)].sample(1)
#     negative_samples[~negative_samples.isin(homo.gene_stable_id)] = 
