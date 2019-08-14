import pandas as pd
import numpy as np
import os
import sys
from selector_sp import select, create_map_reverse


def read_db_homology(dir_name, filename):
    df = pd.read_csv(dir_name + "/" + filename, compression='gzip', sep='\t')
    n = filename.split(".")[0]
    n = n.split(" ")[0]
    return df, n


def get_selection_data():
    with open("dist_matrix", "r") as file:
        matrix = file.readlines()
    matrix = [x.split("\t") for x in matrix]
    matrix = [[float(y) for y in x] for x in matrix]
    matrix = np.array(matrix)
    with open("sp_names", "r") as file:
        dname = file.readlines()
    dname = [x.split("\n")[0] for x in dname]
    spnmap, nspmap = create_map_reverse(dname)
    return matrix, spnmap, nspmap


def read_select_data(dirname, matrix, spnmap, nspmap, nos):
    lf = os.listdir(dirname)
    if len(lf) == 0:
        print("No Files in the Directory!!!!!!!")
        sys.exit(1)
    if not os.path.isdir("processed"):
        os.mkdir("processed")
    for x in lf:
        df, n = read_db_homology(dirname, x)
        n = n.split(" ")[0]
        df = select(df, nos, matrix, spnmap, nspmap, n)
        (len(df) == nos)
        indexes = np.array(list(df.index.values))
        np.save("processed/" + n + "_selected_indexes", indexes)


def main():
    arg = sys.argv
    nos = int(arg[-1])
    matrix, spnmap, nspmap = get_selection_data()
    read_select_data("data_homology", matrix, spnmap, nspmap, nos)


if __name__ == "__main__":
    main()
