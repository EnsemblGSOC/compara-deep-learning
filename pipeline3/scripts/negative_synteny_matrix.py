from skbio.alignment import StripedSmithWaterman
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




samples_df = pd.read_csv(args.samples_dir + "/" + args.species + "_negative_samples_pairs.csv")

# read in the samples for both the species and the candidate homologies
with open( args.samples_dir + "/" + args.species + "_negative_sample_neighbour.json","r") as f:
    sample_neighbours = json.load(f)
    
with open( args.samples_dir + "/" + args.species + "_negative_homology_sample_neighbour.json","r") as f:
    sample_homology_neighbours = json.load(f)

gene_series = samples_df.gene_stable_id.map(sample_neighbours)
gene_series[gene_series.isna()] = gene_series[gene_series.isna()].apply(lambda x: [[0,0,0],[0,0,0]]) # fill the null values
array1 = np.concatenate(gene_series.apply( np.array ).values ).reshape((-1,6))
array1 = np.concatenate((samples_df.gene_stable_id.values[:,np.newaxis], array1), axis=1).reshape(1,-1)

gene_series = samples_df.non_homo_gene_stable_id.map(sample_homology_neighbours)
gene_series[gene_series.isna()] = gene_series[gene_series.isna()].apply(lambda x: [[0,0,0],[0,0,0]]) # fill the null values
array2 = np.concatenate(gene_series.apply( np.array ).values ).reshape((-1,6))
array2 = np.concatenate((samples_df.non_homo_gene_stable_id.values[:,np.newaxis], array2), axis=1).ravel()

# read in the first fasta file
fasta_sequences = SeqIO.parse(open(args.longest_fasta_dir + "/" + args.species + "_longest_pep.fa"),'fasta')
seq_map = {}
for row in fasta_sequences:
    seq_map[row.id.split(".")[0]] = str(row.seq)

sequences1 = pd.Series(array1[0]).map(seq_map)

indices = np.repeat(samples_df.non_homo_species.values, 7)


for query_species in samples_df.non_homo_species.drop_duplicates():
    gene_series = pd.Series(array2[indices == query_species])

    # read in the fasta file
    fasta_sequences = SeqIO.parse(open(args.longest_fasta_dir + "/" + query_species + "_longest_pep.fa"),'fasta')
    seq_map = {}
    for row in fasta_sequences:
        seq_map[row.id.split(".")[0]] = str(row.seq)

    indices[indices == query_species] = gene_series.map(seq_map).values

sequences2 = pd.Series(indices)

seq_array1 = sequences1.values.reshape(-1,7)[:,[1,2,3,0,4,5,6]] # Swap the orders to put the main genes centrally

seq_array2 = sequences2.values.reshape(-1,7)[:,[1,2,3,0,4,5,6]] # Swap the orders to put the main genes centrally

global1 = np.zeros(shape=(seq_array1.shape[0],7,7))
global2 = np.zeros(shape=(seq_array1.shape[0],7,7))
local1 = np.zeros(shape=(seq_array1.shape[0],7,7,3))
local2 = np.zeros(shape=(seq_array1.shape[0],7,7,3))

for i,j,k in progressbar.progressbar(np.ndindex((seq_array1.shape[0],7,7))):
    # print(seq_array1[i,j])
    if (type(seq_array1[i,j]) != str) or (type(seq_array2[i,k]) != str):
        global1[i,j,k] = float("NaN")
        global2[i,j,k] = float("NaN")
        local1[i,j,k] = float("NaN")
        local2[i,j,k] = float("NaN")
    else:
        len_gene1, len_gene2 = len(seq_array1[i,j]),len(seq_array2[i,k])
        norm = max(len_gene1, len_gene2)/2 # factor of half SSW score close to unity for idential sequences
        global1[i,j,k] = ed.align(seq_array1[i,j], seq_array2[i,k], mode="NW", task="distance")["editDistance"] / norm
        global2[i,j,k] = ed.align(seq_array1[i,j], seq_array2[i,k][::-1], mode="NW", task="distance")["editDistance"] / norm
        query = StripedSmithWaterman(seq_array1[i,j])
        # get the forward alignment
        forward_alignment = query(seq_array2[i,k])
        local1[i,j,k,0] = forward_alignment["optimal_alignment_score"]/norm # get the normalised Smith-Waterman alignment score
        alignment_length = (forward_alignment["target_end_optimal"] - forward_alignment["target_begin"])
        local1[i,j,k,1] = alignment_length/len_gene1 # coverage of gene1 computed by normalising alignment length
        local1[i,j,k,2] = alignment_length/len_gene2 # coverage of gene2 computed by normalising alignment length
        # repeat steps for alignment against the reversed sequences
        reverse_alignment = query(seq_array2[i,k][::-1]) # alignment for the reversed sequence
        local2[i,j,k,0] = reverse_alignment["optimal_alignment_score"]/norm # get the normalised Smith-Waterman alignment score
        alignment_length = (reverse_alignment["target_end_optimal"] - reverse_alignment["target_begin"])
        local2[i,j,k,1] = alignment_length/len_gene1 # coverage of gene1 computed by normalising alignment length
        local2[i,j,k,2] = alignment_length/len_gene2 # coverage of gene2 computed by normalising alignment length
        
np.save(args.out_dir + "/" + args.species + "_negative_global1.npy", global1)
np.save(args.out_dir + "/" + args.species + "_negative_global2.npy", global2)
np.save(args.out_dir + "/" + args.species + "_negative_local1.npy", local1)
np.save(args.out_dir + "/" + args.species + "_negative_local2.npy", local2)

# global1 = np.array(global1).reshape(seq_array1.shape[0],-1)
# global2 = np.array(global2).reshape(seq_array1.shape[0],-1)
# local1 = np.array(local1).reshape(seq_array1.shape[0],-1)
# local2 = np.array(local2).reshape(seq_array1.shape[0],-1)

# np.savetxt(args.out_dir + "/" + args.species + "_negative_global1.txt", global1, fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_negative_global2.txt", global2, fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_negative_local1.txt", local1, fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_negative_local2.txt", local2, fmt="%s")
# ## perform the pairwise alignments

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

# np.savetxt(args.out_dir + "/" + args.species + "_negative_global1.txt", global1, fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_negative_global2.txt", global2, fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_negative_local1.txt", local1, fmt="%s")
# np.savetxt(args.out_dir + "/" + args.species + "_negative_local2.txt", local2, fmt="%s")