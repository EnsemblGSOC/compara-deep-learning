import sys
from Bio import SeqIO
import argparse
import gzip


def selectGeneID(description):
    tab = description.split()
    for elem in tab:
        if "gene:" in elem:
            return elem.split(":")[-1]
    raise Exception("No gene stable id in file")


def loadSequenceByGene(fileName):
    result = {}
    with gzip.open(fileName, "rt") as handle:
        for seq_record in SeqIO.parse(handle, "fasta"):
            gene_id = selectGeneID(seq_record.description)
            if not gene_id in result:
                result[gene_id] = {}
            result[gene_id][seq_record.id] = seq_record
        return result


def selectLongest(dico_gene):
    result = []
    for gene in dico_gene.keys():
        max = 0
        max_pep = None
        for pep in dico_gene[gene].keys():
            if max < len(dico_gene[gene][pep]):
                max_pep = dico_gene[gene][pep]
                max = len(max_pep)
        max_pep.id = gene
        max_pep.description = ""
        result.append(max_pep)
    return result


def main(args):
    dico_gene = loadSequenceByGene(args.seq)
    long_peptides = selectLongest(dico_gene)
    with open(args.out_fasta, "w") as output_handle:
        SeqIO.write(long_peptides, output_handle, "fasta")


parser = argparse.ArgumentParser()
parser.add_argument("--seq", help="sequence fasta file")
parser.add_argument("--out_fasta", help="out fasta sequence")
args = parser.parse_args()

main(args)
