import pandas as pd
import numpy as np


def create_map_reverse(arr):
    m = {}
    rm = {}
    for i in range(len(arr)):
        m[arr[i]] = i
        rm[i] = arr[i]
    return m, rm


def get_data_prop(df, nspmap, sp, prop, nos):
    nos = int(nos * prop)
    sp = [nspmap[x] for x in sp]
    data = df[df["homology_species"].isin(sp)]
    if nos < len(data):
        random_indexes = np.random.permutation(len(data))
        data = data.loc[data.index.values[random_indexes[:nos]]]
        return data
    else:
        return data


def create_balanced_dataset_paralog(
        df,
        matrix,
        spnmap,
        nspmap,
        nos,
        hom_type,
        spname):
    df = df[df["homology_type"] == hom_type]
    dist = matrix[spnmap[spname]]
    dist_sort = np.argsort(dist)
    sp_far = dist_sort[-5:]
    sp_near = dist_sort[1:6]
    # get the records for the species which are far away
    df_far = get_data_prop(df, nspmap, sp_far, 0.2, nos)
    # get the records for the nearby species
    df_near = get_data_prop(df, nspmap, sp_near, 0.2, nos)
    df_dist = pd.concat([df_far, df_near])
    nos = nos - len(df_dist)
    df = df.drop(df_dist.index.values)
    random_ind = np.random.permutation(len(df))
    df_r = df.loc[df.index.values[random_ind[:nos]]]
    df = pd.concat([df_r, df_dist])
    return df


def select_data_goc(df, prop, nos):
    df = df[df["goc_score"] == 0.0]
    nos = int(prop * nos)
    rind = np.random.permutation(len(df))
    if len(df) < nos:
        return df
    else:
        return df.loc[df.index.values[rind[:nos]]]


def create_balanced_dataset_ortholog(
        df,
        matrix,
        spnmap,
        nspmap,
        nos,
        hom_type,
        spname):
    df = df[df["homology_type"] == hom_type]
    dist = matrix[spnmap[spname]]
    dist_sort = np.argsort(dist)
    sp_far = dist_sort[-5:]
    sp_near = dist_sort[1:6]
    # get the records for the species which are far away
    df_far = get_data_prop(df, nspmap, sp_far, 0.2, nos)
    # get the records for the nearby species
    df_near = get_data_prop(df, nspmap, sp_near, 0.2, nos)
    df_dist = pd.concat([df_far, df_near])
    len(df_dist)
    df = df.drop(df_dist.index.values)
    df_goc = select_data_goc(df, 0.1, nos)
    len(df_goc)
    df = df.drop(df_goc.index.values)
    nos = nos - len(df_dist) - len(df_goc)
    random_ind = np.random.permutation(len(df))
    df_r = df.loc[df.index.values[random_ind[:nos]]]
    df = pd.concat([df_r, df_goc, df_dist])
    return df


def select(df, nos, matrix, spnmap, nspmap, sp):
    nos_p = int((0.5 * nos) / 2)
    nos_o = int((0.5 * nos) / 3)

    # get the paralogy data
    df_p1 = create_balanced_dataset_paralog(
        df, matrix, spnmap, nspmap, nos_p, "within_species_paralog", sp)
    df_p2 = create_balanced_dataset_paralog(
        df, matrix, spnmap, nspmap, nos_p, "other_paralog", sp)
    # get the orthology data
    df_o1 = create_balanced_dataset_ortholog(
        df, matrix, spnmap, nspmap, nos_o, "ortholog_one2many", sp)
    df_o2 = create_balanced_dataset_ortholog(
        df, matrix, spnmap, nspmap, nos_o, "ortholog_many2many", sp)
    df_o3 = create_balanced_dataset_ortholog(
        df, matrix, spnmap, nspmap, nos_o, "ortholog_one2one", sp)

    # concatenate everything
    df = pd.concat([df_o1, df_o2, df_o3, df_p1, df_p2])
    return df
