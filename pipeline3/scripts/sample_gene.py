import pandas as pd
import numpy as np
import argparse
import os
import json

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
# get the rank numbers
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
num_sub_samples = int(args.n_samples / (args.n_species * 2))
samples_df = (
    homology[homology.homology_species.isin(query_species)]
    .groupby(["homology_species", "ortho_para"])
    .sample(num_sub_samples)
)

samples_df.to_csv(args.out + "/" + args.species + "_samples_pairs.csv")

"""
Now we get the neighbours of all the genes using the neighbours.json
file generated earlier in the pipeline
"""

# load the species map

print("loading neighbourhood genes")
with open(args.neighbours + "/" + args.species + "_all_neighbours.json") as f:
    neighbour_map = json.load(f)

# neighbour_genes = dict(
#     zip(samples_df.gene_stable_id, samples_df.gene_stable_id.map(neighbour_map))
# )
# print("writing pickle file: " + args.out + "/" + args.species + "_sample_neighbour.txt")
# with open(args.out + "/" + args.species + "_sample_neighbour.json", "w") as json_file:
#     json.dump(neighbour_genes, json_file)

# read in the list of paths to homology species
homo_paths = pd.Series(np.loadtxt(args.homo_paths, dtype="str"))

# function to get path to homology database of species
def homo_path_map(species):
    return homo_paths[homo_paths.str.contains(species)].values[0]


# loop throug the homolog species to generate neighbours
species_neighbour_dict = {}
homology_neighbour_dict = {}
for species in query_species:
    # get only samples from that species
    species_samples = samples_df[samples_df.homology_species == species]
    neighbour_genes = dict(zip(species_samples.gene_stable_id, species_samples.gene_stable_id.map(neighbour_map)))
    species_neighbour_dict.update(neighbour_genes)
    # read in the query species neighbour list
    with open(args.neighbours + "/" + species + "_all_neighbours.json") as f:
        homology_neighbour_map = json.load(f)
    # get the dictionary of each homolog gene and it's neighbours
    homology_neighbour_genes = dict(
        zip(
            species_samples.homology_gene_stable_id,
            species_samples.homology_gene_stable_id.map(homology_neighbour_map),
        )
    )
    print("species: " + species)

    # append to the homology_neighbours dict
    homology_neighbour_dict.update(homology_neighbour_genes)

with open(args.out + "/" + args.species + "_sample_neighbour.json", "w") as json_file:
    json.dump(species_neighbour_dict, json_file)


with open(
    args.out + "/" + args.species + "_homology_sample_neighbour.json", "w"
) as json_file:
    json.dump(homology_neighbour_dict, json_file)

