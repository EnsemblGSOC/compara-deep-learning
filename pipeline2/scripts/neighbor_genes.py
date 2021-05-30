import sys
import numpy as np
import os
import pickle
from select_data import read_db_homology
from process_data import create_data_homology_ls


def read_genome_maps():
    data = {}
    with open("genome_maps", "rb") as file:
        data = pickle.load(file)
    cmap = data["cmap"]
    cimap = data["cimap"]
    ld = data["ld"]
    ldg = data["ldg"]
    a = data["a"]
    d = data["d"]
    return a, d, ld, ldg, cmap, cimap


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
        except BaseException:
            print("Incomplete data for:", n)
        df = df.loc[indexes]
        a_h.append(df)
        d_h.append(n)
    return a_h, d_h


def main():
    a, d, ld, ldg, cmap, cimap = read_genome_maps()
    print("Genome Maps Loaded.")
    a_h, d_h = read_data_homology("homology_databases")
    print("Data Read.")
    n = 3
    _ = create_data_homology_ls(a_h, d_h, n, a, d, ld, ldg, cmap, cimap, 1)
    print("Neighbor Genes Found and Saved Successfully:)")


if __name__ == "__main__":
    main()
