from skbio import Protein
from skbio.alignment import local_pairwise_align_protein, global_pairwise_align_protein
from Bio import SeqIO
import json
import pandas as pd
import numpy as np
import progressbar
import edlib as ed
import argparse


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
    "--longest_fasta_dir",
    help="directory containing fastas with species longest gene sequences",
)
parser.add_argument( 
    "--out_dir",
    help="output directory for syneny matrices",
)
args = parser.parse_args()




samples_df = pd.read_csv(args.samples_dir + "/" + args.species + "_samples_pairs.csv")

# read in the samples for both the species and the candidate homologies
with open( args.samples_dir + "/" + args.species + "_sample_neighbour.json","r") as f:
    sample_neighbours = json.load(f)
    
with open( args.samples_dir + "/" + args.species + "_homology_sample_neighbour.json","r") as f:
    sample_homology_neighbours = json.load(f)

gene_series = samples_df.gene_stable_id.map(sample_neighbours)
gene_series[gene_series.isna()] = gene_series[gene_series.isna()].apply(lambda x: [[0,0,0],[0,0,0]]) # fill the null values
array1 = np.concatenate(gene_series.apply( np.array ).values ).reshape((-1,6))
array1 = np.concatenate((samples_df.gene_stable_id.values[:,np.newaxis], array1), axis=1).reshape(1,-1)

gene_series = samples_df.homology_gene_stable_id.map(sample_homology_neighbours)
gene_series[gene_series.isna()] = gene_series[gene_series.isna()].apply(lambda x: [[0,0,0],[0,0,0]]) # fill the null values
array2 = np.concatenate(gene_series.apply( np.array ).values ).reshape((-1,6))
array2 = np.concatenate((samples_df.homology_gene_stable_id.values[:,np.newaxis], array2), axis=1).ravel()

# read in the first fasta file
fasta_sequences = SeqIO.parse(open(args.longest_fasta_dir + "/" + args.species + "_longest_pep.fa"),'fasta')
seq_map = {}
for row in fasta_sequences:
    seq_map[row.id.split(".")[0]] = str(row.seq)

sequences1 = pd.Series(array1[0]).map(seq_map)

indices = np.repeat(samples_df.homology_species.values, 7)


for query_species in samples_df.homology_species.drop_duplicates():
    gene_series = pd.Series(array2[indices == query_species])

    # read in the fasta file
    fasta_sequences = SeqIO.parse(open(args.longest_fasta_dir + "/" + query_species + "_longest_pep.fa"),'fasta')
    seq_map = {}
    for row in fasta_sequences:
        seq_map[row.id.split(".")[0]] = str(row.seq)

    indices[indices == query_species] = gene_series.map(seq_map).values

sequences2 = pd.Series(indices)


seq_array1 = sequences1.values.reshape(-1,7)

seq_array2 = sequences2.values.reshape(-1,7)

global1 = np.zeros(shape=(seq_array1.shape[0],7,7))
global2 = np.zeros(shape=(seq_array1.shape[0],7,7))
local1 = np.zeros(shape=(seq_array1.shape[0],7,7))
local2 = np.zeros(shape=(seq_array1.shape[0],7,7))

for i,j,k in progressbar.progressbar(np.ndindex((seq_array1.shape[0],7,7))[:5]):
    # print(seq_array1[i,j])
    if (type(seq_array1[i,j]) != str) or (type(seq_array2[i,k]) != str):
        global1[i,j,k] = float("NaN")
        global2[i,j,k] = float("NaN")
        local1[i,j,k] = float("NaN")
        local2[i,j,k] = float("NaN")
    else:
        
        seq_array1[i,j] seq_array2[i,k]
        # norm = max(len(seq_array1[i,j]),len(seq_array2[i,k]))
        # global1[i,j,k] = ed.align(seq_array1[i,j], seq_array2[i,k], mode="NW", task="distance")["editDistance"] / norm
        # global2[i,j,k] = ed.align(seq_array1[i,j], seq_array2[i,k][::-1], mode="NW", task="distance")["editDistance"] / norm
        # _,result,_ = local_pairwise_align_protein(Protein(seq_array1[i,j]), Protein(seq_array2[i,k]))
        # local1[i,j,k] = result/norm
        # _,result,_ = local_pairwise_align_protein(Protein(seq_array1[i,j]), Protein(seq_array2[i,k])[::-1])
        # local2[i,j,k] = result/norm

# global1 = np.array(global1).reshape(seq_array1.shape[0],-1)
# global2 = np.array(global2).reshape(seq_array1.shape[0],-1)
# local1 = np.array(local1).reshape(seq_array1.shape[0],-1)
# local2 = np.array(local2).reshape(seq_array1.shape[0],-1)


# np.savetxt(args.out_dir + "/" + args.species + "_global1.txt", global1, fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_global2.txt", global2, fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_local1.txt", local1, fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_local2.txt", local2, fmt="%s")

"""
Extinct code from earlier work, but saved if needed
"""
## perform the pairwise alignments

# # perform local alignments and reciprocal alignment
# local1 = [] # forward alignments
# local2 = [] # reciprocal alignments

# global1 = [] # forward
# global2 = [] # reverse

# for i,j in progressbar.progressbar(zip(sequences1,sequences2)):
#     if (type(i) != str) or (type(j) != str) :
#         local1.append(float("NaN"))
#         local2.append(float("NaN"))
#         global1.append(float("NaN"))
#         global2.append(float("NaN"))
#     else:
#         seq1, seq2 = Protein(i), Protein(j)
#         norm = max(len(i),len(j))
#         result = ed.align(i,j, mode="NW", task="distance")
#         global1.append(result["editDistance"] / norm)
#         result = ed.align(i,j[::-1], mode="NW", task="distance")
#         global2.append(result["editDistance"] / norm)
#         # _,result,_ = local_pairwise_align_protein(seq1, seq2)
#         # local1.append(result/norm)
#         # _,result,_ = local_pairwise_align_protein(seq1, seq2[::-1])
#         # local2.append(result/norm)
#         # _,result,_ = global_pairwise_align_protein(seq1, seq2)
# #         global1.append(result/norm)
# #         _,result,_ = global_pairwise_align_protein(seq1, seq2)
# #         global2.append(result/norm)

# global1 = np.array(global1).reshape(-1,7)
# global2 = np.array(global2).reshape(-1,7)
# # local1 = np.array(local1).reshape(-1,7)
# # local2 = np.array(local2).reshape(-1,7)
