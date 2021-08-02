from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument( 
    "--input_fasta",
    help="Species data to process",
)
parser.add_argument( 
    "--out_fasta",
    help="Directory containing all neighbourhood data",
)
args = parser.parse_args()

fasta_sequences = SeqIO.parse(open(args.input_fasta),'fasta') # read in the input fasta file
reved_seqs = []
   
for row in fasta_sequences:
    row.seq = row.seq[::-1] # reverse the peptide sequence
    reved_seqs.append(row)

SeqIO.write(reved_seqs, args.out_fasta, "fasta")  # write the reversed sequences to disc