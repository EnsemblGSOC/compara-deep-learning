import pandas as pd
import numpy as np
import progressbar
import pyranges as pr
import os
import json
import sys
from prepare_synteny_matrix import read_data_homology, load_neighbor_genes
from process_data import create_map_list
from process_negative import read_database_txt


# def get_score_overlap(x, y, pfam_db, pfam_map):
#     df_1 = pfam_db.loc[pfam_map[x]]
#     df_2 = pfam_db.loc[pfam_map[y]]
#     l = list(df_2.sseqid)
#     c = 0
#     c_1 = 0
#     for _, row in df_1.iterrows():
#         if row.sseqid in l:  # check if the sseqid exists in the list
#             c_1 += 1
#             st = int(df_2[df_2["sseqid"] == row.sseqid].sstart)  # get the start
#             end = int(df_2[df_2["sseqid"] == row.sseqid].send)  # get the end
#             if (int(row.sstart) > st and int(row.sstart) < end) or (
#                 int(row.send) > st and int(row.sstart) < end
#             ):  # check if the sseqid is a ovelapping sseqid
#                 c += 1
#     return c / max(len(df_1), len(df_2))


def get_score_overlap(x, y, pfam_db, pfam_map):
    """
    Get's the Jaccard index of the SET of aligned Pfam domains of two genes.
    We piggyback on pyranges to check for the overlaps
    """
    # Get the set of Pfam domains for each gene
    df_1 = pfam_db.loc[pfam_map[x]]
    df_2 = pfam_db.loc[pfam_map[y]]
    # Put these dataframes into a pyranges friendly format.
    """
    Pyranges requires the chromosome term, and cannot be changed.
    Mentally read Pfam domain, instead of chromosome!
    """
    temp_df1 = df_1[["sseqid", "sstart", "send"]].rename(
        columns={"sseqid": "Chromosome", "sstart": "Start", "send": "End"}
    )
    temp_df2 = df_2[["sseqid", "sstart", "send"]].rename(
        columns={"sseqid": "Chromosome", "sstart": "Start", "send": "End"}
    )
    # Convert to the pyranges datatype
    ranges_1 = pr.PyRanges(temp_df1)
    ranges_2 = pr.PyRanges(temp_df2)
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
    return overlapped_pfams / total_pfams


def pfam_matrix(g1, g2, n, pfam_db, gmap, pfam_map):
    pm = np.zeros((n, n))
    for i in range(n):
        if g1[i] == "NULL_GENE":
            continue
        try:
            _ = gmap[g1[i]]
        except:
            continue
        for j in range(n):
            if g2[j] == "NULL_GENE":
                continue
            try:
                _ = gmap[g2[j]]
            except:
                continue
            pm[i][j] = get_score_overlap(g1[i], g2[j], pfam_db, pfam_map)
    return pm


def create_pfam_matrix(df, lsy, pfam_db, pfam_map):
    n = 3
    glist = list(pfam_db.gene_stable_id)
    gmap = create_map_list(glist)
    pg = []
    indexes = []
    for index, row in progressbar.progressbar(df.iterrows()):
        g1 = str(row["gene_stable_id"])
        g2 = str(row["homology_gene_stable_id"])
        x = []
        y = []
        try:
            _ = lsy[g1]
            _ = lsy[g2]
        except:
            continue
        for i in range(len(lsy[g1]["b"]) - 1, -1, -1):
            x.append(lsy[g1]["b"][i])
        x.append(g1)
        for k in lsy[g1]["f"]:
            x.append(k)

        for i in range(len(lsy[g2]["b"]) - 1, -1, -1):
            y.append(lsy[g2]["b"][i])
        y.append(g2)
        for k in lsy[g2]["f"]:
            y.append(k)

        assert len(x) == len(y)
        assert len(x) == (2 * n + 1)
        pmtemp = pfam_matrix(x, y, 2 * n + 1, pfam_db, gmap, pfam_map)
        pg.append(pmtemp)
        indexes.append(index)
    return np.array(pg), np.array(indexes)


def create_pfam_map(pfam_db):
    pfam_map = {}
    for index, row in pfam_db.iterrows():
        try:
            _ = pfam_map[row.gene_stable_id]
        except:
            pfam_map[row.gene_stable_id] = []
        pfam_map[row.gene_stable_id].append(index)
    return pfam_map


def main_positive():
    if not os.path.isdir("processed/pfam_matrices"):
        os.mkdir("processed/pfam_matrices")
    print("Reading Homology Data")
    a_h, d_h = read_data_homology("downloads/homology_databases")
    print("Loading neighbour genes")
    lsy = load_neighbor_genes()
    print("Reading pfam alignments")
    pfam_db = pd.read_hdf("Diamond/outputs/positive/pfam_db_positive.h5")
    print("Creating Pfam map")
    pfam_map = (
        pfam_db.reset_index().groupby("gene_stable_id")["index"].apply(list).to_dict()
    )
    ndir = "processed/pfam_matrices/"
    nf1 = "pfam_matrices"
    nf3 = "pfam_indexes"
    for i in progressbar.progressbar(range(len(a_h))):
        df = a_h[i]
        print(len(df))
        pfam_matrices, indexes = create_pfam_matrix(df, lsy, pfam_db, pfam_map)
        np.save(ndir + str(d_h[i]) + "_" + nf1, pfam_matrices)
        np.save(ndir + str(d_h[i]) + "_" + nf3, indexes)
        print(len(indexes))


def read_data_negative(arg):
    df = read_database_txt(arg[-1])
    name = os.path.basename(arg[-1]).split(".")[0]
    print("Working on " + name)
    ind = np.load("processed/synteny_matrices/" + name + "_indexes.npy")
    df = df.loc[ind]
    pfam_db = pd.read_hdf("Diamond/outputs/negative/pfam_db_negative.h5")
    pfam_map = create_pfam_map(pfam_db)
    with open("processed/neighbor_genes_negative.json", "r") as file:
        lsy = dict(json.load(file))
    return df, pfam_db, pfam_map, lsy, name


def main_negative(arg):
    print("Reading negative data")
    df, pfam_db, pfam_map, lsy, name = read_data_negative(arg)
    ndir = "processed/pfam_matrices/"
    nf1 = "pfam_matrices"
    nf3 = "pfam_indexes"
    print(len(df))
    print("Creating negative Pfam matrices")
    pfam_matrices, indexes = create_pfam_matrix(df, lsy, pfam_db, pfam_map)
    np.save(ndir + name + "_" + nf1, pfam_matrices)
    np.save(ndir + name + "_" + nf3, indexes)
    print(len(indexes))


def main():
    arg = sys.argv
    main_positive()
    # main_negative(arg)


if __name__ == "__main__":
    main()
