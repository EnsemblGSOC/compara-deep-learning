# from skbio import Protein
# from skbio.alignment import local_pairwise_align_protein, global_pairwise_align_protein
from Bio import SeqIO
import json
import pandas as pd
import numpy as np
import progressbar
import edlib as ed
import argparse
from multiprocessing import Pool
from tqdm import tqdm


parser = argparse.ArgumentParser()
parser.add_argument( 
    "--species",
    help="Species data to process",
)
parser.add_argument( 
    "--samples_dir",
    help="Directory containing all neighbourhood data",
)
parser.add_argument( 
    "--Pfam_dir",
    help="directory containing pfam alignments",
)
parser.add_argument( 
    "--out_dir",
    help="output directory for syneny matrices",
)
parser.add_argument( 
    "--threads",
    help="output directory for syneny matrices",
    type=int
)
args = parser.parse_args()




samples_df = pd.read_csv(args.samples_dir + "/" + args.species + "_negative_samples_pairs.csv")

# read in the samples for both the species and the candidate homologies
with open( args.samples_dir + "/" + args.species + "_negative_sample_neighbour.json","r") as f:
    sample_neighbours = json.load(f)
    
with open( args.samples_dir + "/" + args.species + "_negative_homology_sample_neighbour.json","r") as f:
    sample_homology_neighbours = json.load(f)

gene_series = samples_df.gene_stable_id.map(sample_neighbours)
gene_series[gene_series.isna()] = gene_series[gene_series.isna()].apply(lambda x: [[0,0,0],[0,0,0]]) # fill the null values
array1 = np.concatenate(gene_series.apply( np.array ).values ).reshape((-1,6))
array1 = np.concatenate((samples_df.gene_stable_id.values[:,np.newaxis], array1), axis=1) #.reshape(1,-1)
array1 = array1[:,[1,2,3,0,4,5,6]] # put the main genes in the centre

gene_series = samples_df.non_homo_gene_stable_id.map(sample_homology_neighbours)
gene_series[gene_series.isna()] = gene_series[gene_series.isna()].apply(lambda x: [[0,0,0],[0,0,0]]) # fill the null values
array2 = np.concatenate(gene_series.apply( np.array ).values ).reshape((-1,6))
array2 = np.concatenate((samples_df.non_homo_gene_stable_id.values[:,np.newaxis], array2), axis=1) #.ravel()[np.newaxis,:]
array2 = array2[:,[1,2,3,0,4,5,6]] # put the main genes in the centre

# array1 = pd.Series(array1[0])

# array2 = pd.Series(array2[0])

keys = "qseqid qlen sseqid slen qstart qend sstart send".split(" ")

keys = ['qseqid', 'qlen', 'sseqid', 'slen', 'qstart', 'qend', 'sstart', 'send','evalue']
keys = ['qseqid', 'sseqid', 'sstart', 'send']

# pfam = pd.read_csv(args.Pfam_dir + "/" + args.species + ".outs.tsv.gz", delimiter="\t", names=keys, compression="gzip", usecols=[0,2,6,7])

# pfam["qseqid"] = pfam["qseqid"].str.split(".").str[0]

# pfam_dfs = []
# for query_species in samples_df.non_homo_species.drop_duplicates():
#     temp_pfam = pd.read_csv( args.Pfam_dir + "/" + query_species + ".outs.tsv.gz", delimiter="\t", names=keys, compression="gzip", usecols=[0,2,6,7])
#     temp_pfam["qseqid"] = temp_pfam["qseqid"].str.split(".").str[0]
#     pfam_dfs.append(temp_pfam)
#     # pfam_map.update(temp_pfam.groupby("qseqid")["sseqid"].apply(list).to_dict())
# print("Concetenating frames")

# pfam_df = pd.concat(pfam_dfs)

import pyranges as pr

# pfam_array = np.zeros(shape=(array1.shape[0],7,7))

pbar = tqdm(total=np.product(array1.shape), position=0,leave=True)
# for i,j,k in progressbar.progressbar(np.ndindex((array1.shape[0],7,7))):
def overlap(tup):
    pbar.update()
    i,j,k = tup # disagregate the tuple
    """
    Sometimes the genes aren't present in the pfam list. This shouldn't be possible as the queried genes and the
    and the Diamond alignments both come from the longest_pep_fasta file
    """
    if (type(array1[i,j]) != str) or (type(array2[i,k]) != str):
        return (float("NaN"),float("NaN"),float("NaN"))
    if (array1[i,j] == "NA") or (array2[i,k] == "NA"):
        return ("NA","NA","NA")
    species = samples_df.iloc[i,:].species
    homology_species = samples_df.iloc[i,:].non_homo_species
    print(species, homology_species)
    gene1 = array1[i,j]
    gene2 = array2[i,k]
    # read in the relavant files from the directory
    try:
        temp_df1 = pd.read_csv( args.Pfam_dir + "/" + args.species + "/" + gene1 + "_" + args.species + ".outs.tsv.gz", compression="gzip", delimiter="\t", names=keys,usecols=[0,2,6,7])
        temp_df2 = pd.read_csv( args.Pfam_dir + "/" + homology_species + "/" + gene2 + "_" + homology_species + ".outs.tsv.gz", compression="gzip", delimiter="\t",names=keys,usecols=[0,2,6,7] )
    except:
        # NP denotes No Pfam
        print("error")
        return ("NP","NP","NP")
    
    if (temp_df1.shape[0] == 0) or (temp_df2.shape[0] == 0):
        return (float("NaN"),float("NaN"),float("NaN"))
        
    else:
        """
        Pyranges requires the chromosome term, and cannot be changed.
        Mentally read Pfam domain, instead of chromosome!
        """   
        
        temp_df1 = temp_df1[["sseqid", "sstart", "send"]].rename(
            columns={"sseqid": "Chromosome", "sstart": "Start", "send": "End"}
        )
        temp_df2 = temp_df2[["sseqid", "sstart", "send"]].rename(
            columns={"sseqid": "Chromosome", "sstart": "Start", "send": "End"}
        )
        # print("df shapes: ")
        # print(temp_df1.shape,temp_df2.shape)
        ranges_1 = pr.PyRanges(temp_df1)
        ranges_2 = pr.PyRanges(temp_df2)
        # print(ranges_1, ranges_2)
        # Package them into a dictionary
        grs = {n: m for n, m in zip(["ranges_1", " ranges_2"], [ranges_1, ranges_2])}
        #  Check for sequence aligmnets overlapping on a Pfam domain
        overlaps = pr.count_overlaps(grs).df
        # Mark the overlaps
        overlaps["overlap_bool"] = (overlaps[["ranges_1", " ranges_2"]].sum(axis=1) > 1) + 0
        # Get the total number of actually overlapping Pfam domains
        overlapped_pfams = (overlaps.groupby("Chromosome").sum()["overlap_bool"] > 0).sum()
        # Get the total number of possible Pfam overlaps; ie total across both genes
        total_pfams = overlaps["Chromosome"].drop_duplicates().shape[0]
        Jaccard = overlapped_pfams / total_pfams
        # Get the coverage of the overlaps in both directions
        coverage1 = ranges_1.coverage(ranges_2).FractionOverlaps.sum()
        coverage2 = ranges_2.coverage(ranges_1).FractionOverlaps.sum()

        
        return (Jaccard, coverage1, coverage2)

pool = Pool(processes=args.threads)
pfam_list = pool.map( overlap, list(np.ndindex(array1.shape[0],7,7)))
pool.close()
# Separate the different features into lists
pfam_array1 = [ tup[0] for tup in pfam_list ]
pfam_array2 = [ tup[1] for tup in pfam_list ]
pfam_array3 = [ tup[2] for tup in pfam_list ]
# print(pfam_array)
pfam_array1= np.array(pfam_array1).reshape((array1.shape[0],7,7))
pfam_array2= np.array(pfam_array2).reshape((array1.shape[0],7,7))
pfam_array3= np.array(pfam_array3).reshape((array1.shape[0],7,7))
stacked_pfam = np.stack((pfam_array1,pfam_array2,pfam_array3), axis = -1)
print(stacked_pfam.shape)

np.save(args.out_dir + "/" + args.species + "_negative_Pfam.npy", stacked_pfam) # save the stacked pfam as an npy

# np.savetxt(args.out_dir + "/" + args.species + "_negative_Pfam1.txt", pfam_array1.reshape((array1.shape[0],-1)), fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_negative_Pfam2.txt", pfam_array2.reshape((array1.shape[0],-1)), fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_negative_Pfam3.txt", pfam_array3.reshape((array1.shape[0],-1)), fmt="%s")

#         return overlapped_pfams / total_pfams

# pool = Pool(processes=args.threads)
# pfam_array = pool.map( overlap, list(np.ndindex(array1.shape[0],7,7)))
# pool.close()
# # print(pfam_array)
# pfam_array = np.array(pfam_array).reshape((array1.shape[0],7,7))


# np.savetxt(args.out_dir + "/" + args.species + "_negative_Pfam.txt", pfam_array.reshape((array1.shape[0],-1)), fmt="%s")


# pfam_values = []
# print(array1, array2)
# for gene1, gene2 in progressbar.progressbar(zip(array1, array2)):
#     temp_df1 = pfam[pfam.qseqid == gene1]
#     temp_df2 = pfam_df[pfam_df.qseqid == gene2]
#     """
#     Sometimes the genes aren't present in the pfam list. This shouldn't be possible as the queried genes and the
#     and the Diamond alignments both come from the longest_pep_fasta file
#     """
    
#     if temp_df1.empty or temp_df2.empty:
#         pfam_values.append(float("NaN"))
#     else:
#         """
#         Pyranges requires the chromosome term, and cannot be changed.
#         Mentally read Pfam domain, instead of chromosome!
#         """
#         temp_df1 = temp_df1[["sseqid", "sstart", "send"]].rename(
#             columns={"sseqid": "Chromosome", "sstart": "Start", "send": "End"}
#         )
#         temp_df2 = temp_df2[["sseqid", "sstart", "send"]].rename(
#             columns={"sseqid": "Chromosome", "sstart": "Start", "send": "End"}
#         )
#         ranges_1 = pr.PyRanges(temp_df1)
#         ranges_2 = pr.PyRanges(temp_df2)
#         # Package them into a dictionary
#         grs = {n: m for n, m in zip(["ranges_1", " ranges_2"], [ranges_1, ranges_2])}
#         #  Check for sequence aligmnets overlapping on a Pfam domain
#         overlaps = pr.count_overlaps(grs).df
#         # Mark the overlaps
#         overlaps["overlap_bool"] = (overlaps[["ranges_1", " ranges_2"]].sum(axis=1) > 1) + 0
#         # Get the total number of actually overlapping Pfam domains
#         overlapped_pfams = (overlaps.groupby("Chromosome").sum()["overlap_bool"] > 0).sum()
#         # Get the total number of possible Pfam overlaps; ie total across both genes
#         total_pfams = overlaps["Chromosome"].drop_duplicates().shape[0]
#         print(overlapped_pfams / total_pfams)
#         pfam_values.append(overlapped_pfams / total_pfams)


# pfam_values = np.array(pfam_values).reshape(-1,7)

# np.savetxt(args.out_dir + "/" + args.species + "_Pfam.txt", pfam_values, fmt="%s")