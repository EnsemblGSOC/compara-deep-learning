"""
This script reads the hmm_table.tsv file and returns the all species along with their id.

"""
import pandas as pd
df = pd.read_csv("../config/hmm_table.tsv", delimiter="\t", names = ["family","gene_id","sequence_id", "species"])
uniques = df.iloc[df.species.drop_duplicates().index]
uniques["species_id"] = uniques.gene_id.str.split("G0").str[0] + "G0"
species_id_map = dict(zip(list(uniques.species),list(uniques.gene_id.str.split("G0").str[0] + "G0")))
# save the species id's
uniques[["species","species_id"]].to_csv("species_keys.txt", index=False, sep=" ")