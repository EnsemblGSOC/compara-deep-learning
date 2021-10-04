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
    "--synt_dir",
    help="directory containing synteny matrices",
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

"""
merge the synteny data
"""

# combine the paralog data
para_paths = glob.glob( args.synt_dir + "/*paralogs_synteny.npy")
# read in the data
para_list = [ np.load(path) for path in para_paths]
# combine it
para_synt_data = np.concatenate(para_list, axis = 0)

# combine the ortholog data
ortho_paths = glob.glob( args.synt_dir + "/*orthologs_synteny.npy")
# read in the data
ortho_list = [ np.load(path) for path in ortho_paths]
# combine it
ortho_synt_data = np.concatenate(ortho_list, axis = 0)

# combine the paralog data
negs_paths = glob.glob( args.synt_dir + "/*negatives_synteny.npy")
# read in the data
negs_list = [ np.load(path) for path in negs_paths]
# combine it
negs_synt_data = np.concatenate(negs_list, axis = 0)

""" 
Combine the PFAM data
"""
# combine the paralog data
para_paths = glob.glob( args.Pfam_dir + "/*paralogs_Pfam.npy")
# read in the data
para_list = [ np.load(path) for path in para_paths]
# combine it
para_Pfam_data = np.concatenate(para_list, axis = 0)

# combine the ortholog data
ortho_paths = glob.glob( args.Pfam_dir + "/*orthologs_Pfam.npy")
# read in the data
ortho_list = [ np.load(path) for path in ortho_paths]
# combine it
ortho_Pfam_data = np.concatenate(ortho_list, axis = 0)

# combine the paralog data
negs_paths = glob.glob( args.Pfam_dir + "/*negatives_Pfam.npy")
# read in the data
negs_list = [ np.load(path) for path in negs_paths]
# combine it
negs_Pfam_data = np.concatenate(negs_list, axis = 0)


"""
Combine the dataframes
"""


# combine the paralog data
df_ortho_paths = glob.glob( args.samples_dir + "/*_orthologs_samples.csv")
# read in the data
ortho_list = [ pd.read_csv(path) for path in df_ortho_paths]
# combine it
ortho_data = pd.concat(ortho_list)
ortho_data["homology_type"] = "ortho"

# combine the paralog data
df_para_paths = glob.glob( args.samples_dir + "/*_paralogs_samples.csv")
# read in the data
para_list = [ pd.read_csv(path) for path in df_para_paths]
# combine it
para_data = pd.concat(para_list, axis = 0)
para_data["homology_type"] = "para"

# combine the paralog data
df_negs_paths = glob.glob( args.samples_dir + "/*_negatives_samples.csv")
# read in the data
negs_list = [ pd.read_csv(path) for path in df_negs_paths]
# combine it
negs_data = pd.concat(negs_list, axis = 0)
negs_data["homology_type"] = "negs"



big_dataframe = pd.concat([para_data,ortho_data,negs_data])
big_synt_data = np.concatenate((para_synt_data,ortho_synt_data,negs_synt_data), axis=0)
big_Pfam_data = np.concatenate((para_Pfam_data,ortho_Pfam_data,negs_Pfam_data), axis=0)
big_data = np.concatenate((big_synt_data,big_Pfam_data), axis=-1)

big_dataframe.to_csv(args.out_dir + "/" + args.date + "_big_final.csv")
np.save(args.out_dir + "/" + args.date + "_big_final.npy",big_data)

print(big_dataframe.shape, big_data.shape)