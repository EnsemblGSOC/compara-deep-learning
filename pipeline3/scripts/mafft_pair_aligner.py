import pandas as pd
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
import os

df = pd.read_csv("sequence_pairs.csv")
"""
Mafft doesn't naturally like the amino acid U selenocysteine. The command line version has the option --anysymbol
see here: https://mafft.cbrc.jp/alignment/software/anysymbol.html

However the Bio Python version does not appear to support this symbol, causing a bug. I simply replace these occurences
with X, the unkown symbol, which mafft can handle. These occurrences are rare and should have minimal impact on performances.
But something better should be done than this hacky solution
"""
df.species_sequence = df.species_sequence.str.replace("U","X")
df.homolog_sequence = df.homolog_sequence.str.replace("U","X")

def align(a):
    record1 = SeqRecord(
    Seq(a[0]),
    id="species_seq",
    description=""
    )

    record2 = SeqRecord(
    Seq(a[1]),
    id="homo_seq",
    description="")
    output = [record1,record2]
    SeqIO.write(output, a[0][:5] + ".fa", "fasta")

    mafft_cline = MafftCommandline( input=a[0][:5] + ".fa")
    stdout, stderr = mafft_cline()
    aligned_pairs = str(stdout).replace("\n","").replace(">homo_","").split("seq")[1:]
    os.remove(a[0][:5] + ".fa")
    return aligned_pairs

aligned_sequences = df[["species_sequence","homolog_sequence"]].apply(align, axis = 1)
df["foo1"], df["foo2"] = [elem[0] for elem in aligned_sequences], [elem[1] for elem in aligned_sequences]

df.to_csv("aligned_sequence_pairs.csv")

"""
This serves as diagnostic copde
"""
# i, j = df.iloc[802][["species_sequence","homolog_sequence"]]

# record1 = SeqRecord(
#     Seq(i),
#     id="species_seq",
#     description=""
# )
# record2 = SeqRecord(
#     Seq(j),
#     id="homo_seq",
#     description="")

# output = [record1,record2]
# SeqIO.write(output, "temp.fa", "fasta")

# from Bio.Align.Applications import MafftCommandline
# mafft_cline = MafftCommandline( input="temp.fa", localpair=True, amino=True)
# stdout, stderr = mafft_cline()
# aligned_pairs = str(stdout).replace("\n","").replace(">homo_","").split("seq")[1:]
# print(aligned_pairs)