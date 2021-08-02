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
    df_neg["ortho_para"] = "not" # label the non homolog pairs
    # pfam matrix paths
    Pfam_path = args.out_dir + "/Pfam_matrices/" + species + "_Pfam.npy"
    negative_Pfam_path = args.out_dir + "/Pfam_matrices/" + species + "_negative_Pfam.npy"
    # global matrix paths
    negative_global1_path = args.out_dir + "/negative_synteny_matrices/" + species + "_negative_global1.npy"
    negative_global2_path = args.out_dir + "/negative_synteny_matrices/" + species + "_negative_global2.npy"
    global1_path = args.out_dir + "/synteny_matrices/" + species + "_global1.npy"
    global2_path = args.out_dir + "/synteny_matrices/" + species + "_global2.npy"
    # local matrix paths
    negative_local1_path = args.out_dir + "/negative_synteny_matrices/" + species + "_negative_local1.npy"
    negative_local2_path = args.out_dir + "/negative_synteny_matrices/" + species + "_negative_local2.npy"
    local1_path = args.out_dir + "/synteny_matrices/" + species + "_local1.npy"
    local2_path = args.out_dir + "/synteny_matrices/" + species + "_local2.npy"


    # load the positive data
    Pfam = np.load(Pfam_path)
    global1 = np.load(global1_path)
    global2 = np.load(global2_path)
    pos_global = np.stack([global1,global2], axis = -1) # stack the global alignments
    local1 = np.load(local1_path)
    local2 = np.load(local2_path)
    print(Pfam.shape)
    print(pos_global.shape)
    print(local1.shape, local2.shape)
    pos_data = np.concatenate((Pfam,pos_global,local1,local2), axis=-1)
    # load the negative data
    negative_Pfam = np.load(negative_Pfam_path)
    negative_global1 = np.load(negative_global1_path)
    negative_global2 = np.load(negative_global2_path)
    neg_global = np.stack([negative_global1,negative_global2], axis = -1) # stack the global alignments
    negative_local1 = np.load(negative_local1_path)
    negative_local2 = np.load(negative_local2_path)
    neg_data = np.concatenate((negative_Pfam,neg_global,negative_local1,negative_local2), axis=-1)
    # assert the data has the correct shape
    print(df_pos.shape)
    print(df_neg.shape)
    assert df_pos.shape[0] == pos_data.shape[0], "The metadata pos_dataframe and pos_data shapes do not match"
    assert df_neg.shape[0] == neg_data.shape[0], "The metadata neg_dataframe and neg_data shapes do not match"
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

    # # read in the data
    # pos_paths = args.synt + "/" + species + "_" + files + ".txt"
    # # correct the paths for the pfam matrices
    # pos_paths[pos_paths.str.contains("Pfam")] = pos_paths[pos_paths.str.contains("Pfam")].str.replace(args.synt,args.Pfam_dir )
    # # if not os.path.isfile(pos_paths.iloc[-1]):
    # #     print(species + " has missing pfam!")
    # #     continue

    # pos_data = []
    # for path in pos_paths:
    #     pos_data.append(np.loadtxt(path, dtype="str").reshape(-1,7,7)[:,:,:,np.newaxis])

    # pos_data = np.concatenate(pos_data, axis = -1)

    # neg_paths = pos_paths.str.replace( species + "_", species + "_negative_" ).str.replace( "synteny_matrices", "negative_synteny_matrices" )
    # neg_data = []
    # for path in neg_paths:
    #     neg_data.append(np.loadtxt(path, dtype="str").reshape(-1,7,7)[:,:,:,np.newaxis])
    # neg_data = np.concatenate(neg_data, axis = -1)

    # print(pos_data.shape)
    # print(neg_data.shape)

    # #combine the positive and negative data matrices
    # data = np.concatenate([pos_data,neg_data], axis=0)
    # print("#Matrices: " + str(data.shape[0]))

    # # coerce the positive and negative dataframes to have matching columns
    # pos_keys = ["gene_stable_id","homology_gene_stable_id","species","homology_species","ortho_para"]
    # pos_column_map = {"non_homo_species":"query_species","non_homo_gene_stable_id":"query_gene_stable_id"}
    # df_pos = df_pos[pos_keys].rename(columns=pos_column_map)
    # neg_keys = ["gene_stable_id","non_homo_gene_stable_id","species","non_homo_species","ortho_para"]
    # neg_column_map = {"non_homo_species":"query_species","non_homo_gene_stable_id":"query_gene_stable_id"}
    # df_neg = df_neg[neg_keys].rename(columns=neg_column_map)
    # # combine the frame and append to the total data list
    # dataframe = pd.concat([df_pos,df_neg])
    # print("#Metadata Entries: " + str(dataframe.shape[0]))

    # # check the data frames match
    # if dataframe.shape[0] != data.shape[0]:
    #     raise ValueError("Different number of entries in Meta Dataframe and matrix entries. Major pipeline flaw!")
    
    # all_matrices.append(data)
    # all_dataframes.append(dataframe)

# combine all the data into a massive file and save
big_dataframe = pd.concat(all_dataframes)
big_data = np.concatenate(all_matrices, axis=0)

big_dataframe.to_csv(args.out_dir + "/Final_data/" + args.date + "_big_final.csv")
np.save(args.out_dir + "/Final_data/" + args.date + "_big_final.npy",big_data)