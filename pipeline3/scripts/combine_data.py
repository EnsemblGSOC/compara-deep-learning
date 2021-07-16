import numpy as np
import pandas as pd
import glob
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument( 
    "--samples_dir",
    help="Directory containing all neighbourhood data",
)
parser.add_argument( 
    "--Pfam_dir",
    help="directory containing pfam alignments",
)

parser.add_argument( 
    "--synt",
    help="directory containing synteny matrices",
)

parser.add_argument( 
    "--neg_synt",
    help="directory containing negtive synteny matrices",
)

parser.add_argument( 
    "--out_dir",
    help="output directory for finalised dataset",
)
parser.add_argument( 
    "--date",
    help="Used to name the output file",
)
args = parser.parse_args()


paths = glob.glob( args.samples_dir + "/*negative_sample_neighbour.json")
# get the list of species with data available
species_list = [ path.split("/")[-1].split("_neg")[0] for path in paths ]

all_dataframes = []
all_matrices = []

files = pd.Series(["global1","global2","local1","local2","Pfam"])
for species in species_list:
    print(species)
    # read in the meta data frame
    df_pos = pd.read_csv( args.samples_dir + "/" + species + "_samples_pairs.csv")
    df_neg = pd.read_csv( args.samples_dir + "/" + species + "_negative_samples_pairs.csv")
    df_neg["ortho_para"] = "not"

    print(df_pos.shape)
    print(df_neg.shape)

    # read in the data
    pos_paths = args.synt + "/" + species + "_" + files + ".txt"
    # correct the paths for the pfam matrices
    pos_paths[pos_paths.str.contains("Pfam")] = pos_paths[pos_paths.str.contains("Pfam")].str.replace(args.synt,args.Pfam_dir )
    # if not os.path.isfile(pos_paths.iloc[-1]):
    #     print(species + " has missing pfam!")
    #     continue

    pos_data = []
    for path in pos_paths:
        pos_data.append(np.loadtxt(path).reshape(-1,7,7)[:,:,:,np.newaxis])

    pos_data = np.concatenate(pos_data, axis = -1)

    neg_paths = pos_paths.str.replace( species + "_", species + "_negative_" ).str.replace( "synteny_matrices", "negative_synteny_matrices" )
    neg_data = []
    for path in neg_paths:
        neg_data.append(np.loadtxt(path).reshape(-1,7,7)[:,:,:,np.newaxis])
    neg_data = np.concatenate(neg_data, axis = -1)

    print(pos_data.shape)
    print(neg_data.shape)

    #combine the positive and negative data matrices
    data = np.concatenate([pos_data,neg_data], axis=0)
    print("#Matrices: " + str(data.shape[0]))

    # coerce the positive and negative dataframes to have matching columns
    pos_keys = ["gene_stable_id","homology_gene_stable_id","species","homology_species","ortho_para"]
    pos_column_map = {"non_homo_species":"query_species","non_homo_gene_stable_id":"query_gene_stable_id"}
    df_pos = df_pos[pos_keys].rename(columns=pos_column_map)
    neg_keys = ["gene_stable_id","non_homo_gene_stable_id","species","non_homo_species","ortho_para"]
    neg_column_map = {"non_homo_species":"query_species","non_homo_gene_stable_id":"query_gene_stable_id"}
    df_neg = df_neg[neg_keys].rename(columns=neg_column_map)
    # combine the frame and append to the total data list
    dataframe = pd.concat([df_pos,df_neg])
    print("#Metadata Entries: " + str(dataframe.shape[0]))

    # check the data frames match
    if dataframe.shape[0] != data.shape[0]:
        raise ValueError("Different number of entries in Meta Dataframe and matrix entries. Major pipeline flaw!")
    
    all_matrices.append(data)
    all_dataframes.append(dataframe)

# combine all the data into a massive file and save
big_dataframe = pd.concat(all_dataframes)
big_data = np.concatenate(all_matrices, axis=0)

big_dataframe.to_csv(args.out_dir + "/" + args.date + "_big_final.csv")
np.save(args.out_dir + "/" + args.date + "_big_final.npy",big_data)