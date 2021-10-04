import pandas as pd
import numpy as np
import argparse
import os
import json

parser = argparse.ArgumentParser()
parser.add_argument("--out", help="output directory")
parser.add_argument("--species", help="name of species we're making comparisons to")
parser.add_argument(
    "--neighbours",
    help="directory containing all neighbour pairs for each species",
)
args = parser.parse_args()
print(args)


out='/nfs/production/flicek/ensembl/compara/amarshall/Data_Storage/Sample_Neighbours'
species='homo_sapiens'
neighbours='/nfs/production/flicek/ensembl/compara/amarshall/Data_Storage/Species_Neighbours'

# read in the neighbour samples
orthologs = pd.read_csv(out + "/" + species + "_orthologs_samples.csv")
paralogs = pd.read_csv(out + "/" + species + "_paralogs_samples.csv")
negatives = pd.read_csv(out + "/" + species + "_negatives_samples.csv")

# generate dictionary mapping all genes to neighbours from the species in ortho_samples
print("loading neighbourhood genes")
with open(neighbours + "/" + species + "_all_neighbours.json") as f:
    neighbour_map = json.load(f)

for comparison_species in orthologs.homology_species.drop_duplicates():
    with open(neighbours + "/" + comparison_species + "_all_neighbours.json") as f:
        neighbour_map.update(json.load(f)) # add more genes to the neighbour map

print(orthologs.gene_stable_id.map(neighbour_map))
# get neighbours of orthologs of the main species
with open(out + "/" + species + "_ortholog_sample_neighbour.json", "w") as json_file:
    json.dump(dict(zip(orthologs.gene_stable_id,orthologs.gene_stable_id.map(neighbour_map))), json_file)

# get neighbours of orthologs of the comparison species
with open(out + "/" + species + "_ortholog_homology_sample_neighbour.json", "w") as json_file:
    json.dump(dict(zip(orthologs.gene_stable_id,orthologs.gene_stable_id.map(neighbour_map))), json_file)

# get neighbours of paralogs of the main species
with open(out + "/" + species + "_paralog_sample_neighbour.json", "w") as json_file:
    json.dump(dict(zip(paralogs.gene_stable_id,paralogs.gene_stable_id.map(neighbour_map))), json_file)

# get neighbours of paralogs of the comparison species
with open(out + "/" + species + "_paralog_homology_sample_neighbour.json", "w") as json_file:
    json.dump(dict(zip(paralogs.gene_stable_id,paralogs.gene_stable_id.map(neighbour_map))), json_file)

# get neighbours of negatives of the main species
with open(out + "/" + species + "_negative_sample_neighbour.json", "w") as json_file:
    json.dump(dict(zip(negatives.gene_stable_id,negatives.gene_stable_id.map(neighbour_map))), json_file)

# get neighbours of negatives of the comparison species
with open(out + "/" + species + "_negative_homology_sample_neighbour.json", "w") as json_file:
    json.dump(dict(zip(negatives.gene_stable_id,negatives.gene_stable_id.map(neighbour_map))), json_file)