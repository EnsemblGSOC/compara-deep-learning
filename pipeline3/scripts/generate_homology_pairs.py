import pandas as pd
import numpy as np
import itertools
import argparse
import os
import json
import sys

np.random.seed(42)

parser = argparse.ArgumentParser()
parser.add_argument("--homo_db_path", help="homology database path")

parser.add_argument("--out", help="output directory")
parser.add_argument("--species", help="name of species we're making comparisons to")
parser.add_argument("--species_list", help="text file with comparison species in")
parser.add_argument("--n_samples", type=int, help="The number of pairs to get of each homology type, from each species comparison", default=100)
parser.add_argument("--hmm_file", type=str, help="The hmm family text file", default="./config/hmm_table.tsv.gz")

parser.add_argument(
    "--neighbours",
    help="directory containing all neighbour pairs for each species",
)
args = parser.parse_args()
print(args)
print(args.n_samples)









# read the homology database file for the current species or interest
df = pd.read_csv(args.homo_db_path, compression="gzip", delimiter="\t")
# read the hmm file
hmm = pd.read_csv(args.hmm_file, compression="gzip", delimiter="\t",names = ["family","gene_id","sequence_id", "species"])
# get the species pairs to get samples from
sp_list = pd.read_csv(args.species_list, delimiter=" ", names = ["species"])
# reduce the species list to only those in the sp_file
sp_list = sp_list[sp_list.species.isin(df.homology_species)]
print("comparison species list: " + str(sp_list))

# use the neighbourhood file to get the gene list
print("loading neighbourhood genes file: " + args.neighbours + "/" + args.species + "_all_neighbours.json")
with open(args.neighbours + "/" + args.species + "_all_neighbours.json") as f:
    neighbour_map = json.load(f)

print("neighbour map length: " + str(len(neighbour_map)))
species_genes = pd.Series(neighbour_map.keys())
print( "length of the species_genes" + str(len(species_genes)))


orthologs = []
paralogs = []
negatives = []


# iterate through all species comparisons
for comparison_species in sp_list["species"]:
    print("comparing to " + comparison_species)
    """
    generate the list of all pairs of genes in the same family 
    between the comparison and main species. Eg all gene
    family pairs in the human (main) and mouse (comparison 1), then
    family pairs in the human (main) and monkey (comparison 2)
    """
    # get the families of genes in reference and comparison species
    temp_hmm  = hmm[hmm.species.isin([args.species,comparison_species])]
    family_pairs = []
    # iterate through the families
    for family in temp_hmm.family.drop_duplicates():
        # get all genes sharing the family
        temp = temp_hmm[temp_hmm.family == family]
        # split this into two databases for the main species and the comparison
        main = temp[temp.species == args.species]
        comp= temp[temp.species == comparison_species]
        # skip if the family is not present in one of the species
        if (main.shape[0] == 0) | (comp.shape[0] == 0):
            continue
        # enumerate all pairs of gene between the species, in the same family
        pairs = pd.Series(itertools.product(list(main.gene_id),list(comp.gene_id)))
        # add it to the list
        family_pairs.append(pairs)    
    # concatenate this into a series of all gene pairs sharing a family
    family_pairs = pd.concat(family_pairs)
    # divide family into orthologs and not ortholog (ie, paralogs)
    # get all orthologous pairs
    ortholog_pairs = df[(df.homology_species == comparison_species) & (df.homology_type.str.contains("ortholog"))][["gene_stable_id","homology_gene_stable_id"]].apply(tuple, axis = 1)
    if len(ortholog_pairs) == 0:
        continue
    # paralogs are those genes which share a family, but are note orthologs
    paralog_pairs = family_pairs[~family_pairs.isin(ortholog_pairs)]
    # randomly sample from both sets
    print(ortholog_pairs)
    ortho_samples = ortholog_pairs.sample(args.n_samples)
    para_samples = paralog_pairs.sample(args.n_samples)
    # convert to dataframe for saving
    ortho_samples = pd.DataFrame(np.concatenate(ortho_samples.apply(np.array).values).reshape(-1,2), columns=["gene_stable_id","homology_gene_stable_id"] )
    para_samples = pd.DataFrame(np.concatenate(para_samples.apply(np.array).values).reshape(-1,2), columns=["gene_stable_id","homology_gene_stable_id"] )
    ortho_samples["species"] = args.species
    para_samples["species"]     = args.species
    ortho_samples["homology_species"] = comparison_species
    para_samples["homology_species"] = comparison_species
    orthologs.append(ortho_samples)
    paralogs.append(para_samples)


    """
    Now get genes which are not orthologs or paralogs. Do this by randomly choosing pairs between the two species,
    and removing any pairs which are in the family or orthology pair list
    """
    # get the comparison species gene list
    # use the neighbourhood file to get the gene list
    print("loading neighbourhood genes")
    with open(args.neighbours + "/" + comparison_species + "_all_neighbours.json") as f:
        comparison_map = json.load(f)
    comparison_species_genes = pd.Series(comparison_map.keys())
    print(comparison_species_genes)
    
    #randomly sample n-samples from both main and comparison, and make a combined series
    print("sampling negative genes")
    s1 = species_genes.sample(args.n_samples) 
    s2 = comparison_species_genes.sample(args.n_samples)
    negative_pairs = pd.Series(zip(s1,s2))
    # get boolean of the negative pairs in the ortho or family pairs
    bl = negative_pairs.isin(family_pairs) | negative_pairs.isin(ortholog_pairs)
    # retain only values not in either family of negative pairs
    negative_pairs = negative_pairs[~bl].drop_duplicates()
    print(negative_pairs)
    print(negative_pairs.shape)
    """
    remove pairs in family or ortholog pairs and drop duplicates
    resample remaining number of pairs n_samples - negative_pairs.shape[0]
    do this until n_samples - negative_pairs.shape[0] == 0 
    """
    while args.n_samples != negative_pairs.shape[0]:
        # resample to top up to the desired value
        s1 = species_genes.sample(n_samples - negative_pairs.shape[0]) 
        s2 = comparison_species_genes.sample(n_samples - negative_pairs.shape[0])
        negative_pairs = pd.concat([negative_pairs,resample], axis=0)
        # get boolean of the negative pairs in the ortho or family pairs
        bl = negative_pairs.isin(family_pairs) | negative_pairs.isin(ortholog_pairs)
        print(negative_pairs.shape)
        # retain only values not in either family of negative pairs
        negative_pairs = negative_pairs[~bl].drop_duplicates()

    negative_pairs = pd.DataFrame(np.concatenate(negative_pairs.apply(np.array).values).reshape(-1,2), columns=["gene_stable_id","homology_gene_stable_id"] )
    negative_pairs["species"] = args.species
    negative_pairs["homology_species"] = comparison_species
    negatives.append(negative_pairs)

if len(orthologs) == 0 or len(paralogs) == 0 or len(negatives) == 0:
    exit(0)


orthologs = pd.concat(orthologs, axis = 0)
paralogs = pd.concat(paralogs, axis = 0)
negatives = pd.concat(negatives, axis = 0)


orthologs.to_csv(args.out + "/" + args.species + "_orthologs_samples.csv")
paralogs.to_csv(args.out + "/" + args.species + "_paralogs_samples.csv")
negatives.to_csv(args.out + "/" + args.species + "_negatives_samples.csv")