import numpy as np
import pandas as pd
import json
import os
import sys
import pickle
from select_data import read_db_homology
from threads import Procerssrunner
from read_get_gene_seq import read_gene_sequences
from access_data_rest import update_rest, update_rest_protein
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def write_fasta(sequences, name):
    with open(name + ".fa", "w") as file:
        for seq in sequences:
            if sequences[seq] == "":
                continue
            record = SeqRecord(Seq(sequences[seq], IUPAC.protein), id=seq)
            SeqIO.write(record, file, "fasta")


def read_data_synteny(nop, name):
    smg = []
    sml = []
    indexes = []
    for i in range(nop):
        try:
            with open(
                "temp_" + name + "/thread_" + str(i + 1) + "_smg.temp", "rb"
            ) as file:
                smg = smg + pickle.load(file)
            with open(
                "temp_" + name + "/thread_" + str(i + 1) + "_sml.temp", "rb"
            ) as file:
                sml = sml + pickle.load(file)
            with open(
                "temp_" + name + "/thread_" + str(i + 1) + "_indexes.temp", "rb"
            ) as file:
                indexes = indexes + pickle.load(file)
        except Exception as e:
            print("Problem with thread", i + 1, "detected for", name, e)
            continue
    print(len(indexes))
    return smg, sml, indexes


def load_neighbor_genes():
    with open("processed/neighbor_genes.json", "r") as file:
        lsy = dict(json.load(file))
    print(len(lsy))
    print("Neighbor Genes Loaded")
    return lsy


def read_data_homology(dirname):
    lf = os.listdir(dirname)
    if len(lf) == 0:
        print("No Files in the Directory!!!!!!!")
        sys.exit(1)
    a_h = []
    d_h = []
    for x in lf:
        df, n = read_db_homology(dirname, x)
        n = n.split()[0]
        try:
            indexes = np.load("processed/" + n + "_selected_indexes.npy")
        except:
            print("Incomplete data for:", n)
        df = df.loc[indexes]
        a_h.append(df)
        d_h.append(n)
    return a_h, d_h


def main():
    arg = sys.argv
    nop = abs(int(arg[-1]))
    n = 3
    a_h, d_h = read_data_homology("downloads/homology_databases")
    print("Data Read")
    lsy = load_neighbor_genes()
    gene_sequences = read_gene_sequences(a_h, lsy, "downloads/cds", "gene_seq_positive")
    gene_sequences = update_rest(gene_sequences, "gene_seq_positive")
    print("Gene Sequences Loaded.")
    if not os.path.isdir("processed/synteny_matrices"):
        os.mkdir("processed/synteny_matrices")
    ndir = "processed/synteny_matrices/"
    nf1 = "synteny_matrices_global"
    nf2 = "synteny_matrices_local"
    nf3 = "indexes"
    for i in range(len(a_h)):
        df = a_h[i]
        part = len(df) // nop
        pr = Procerssrunner()
        pr.start_processes(nop, df, gene_sequences, lsy, part, n, d_h[i])
        smg, sml, indexes = read_data_synteny(nop, d_h[i])
        print(len(indexes))
        np.save(ndir + str(d_h[i]) + "_" + nf1, smg)
        np.save(ndir + str(d_h[i]) + "_" + nf2, sml)
        np.save(ndir + str(d_h[i]) + "_" + nf3, indexes)
        a_h[i] = df.loc[indexes]
    print("Synteny Matrices Created Successfully :)")
    protein_sequences = read_gene_sequences(
        a_h, lsy, "downloads/pep", "pro_seq_positive"
    )
    protein_sequences = update_rest_protein(protein_sequences, "pro_seq_positive")
    write_fasta(protein_sequences, "protein_seq_positive")
    print("Protein Sequences Loaded.")


if __name__ == "__main__":
    main()
