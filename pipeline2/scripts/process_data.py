import numpy as np
import progressbar
from save_data import write_dict_json


def create_map_list(l):  # this function maps the indexes to values
    t = {}
    for i in range(len(l)):
        t[l[i]] = i
    return t


def create_chromosome_maps(a, n):
    cmap = []
    cimap = []
    for df in progressbar.progressbar(a):
        chmap = {}
        chindmap = {}
        for index, row in df.iterrows():
            g = row.gene_id
            try:
                _ = chmap[g]
            except BaseException:
                chmap[g] = str(row.Chr)
                if str(row.Chr) in chindmap:
                    chindmap[str(row.Chr)].append(index)
                else:
                    chindmap[str(row.Chr)] = []
                    chindmap[str(row.Chr)].append(index)
        cmap.append(chmap)
        cimap.append(chindmap)
    return cmap, cimap


def list_dict_genomes(a, n):
    lst = []
    ldt = []
    for x in a:
        ldgt = {}
        uc = list(x["gene_id"])
        for i, r in x.iterrows():
            ldgt[r.gene_id] = i
        lst.append(uc)
        ldt.append(ldgt)
    return lst, ldt


def get_nearest_neighbors(g, gs, n, a, d, ld, ldg, cmap, cimap):
    # print("Finding Neighbor Genes")
    ne = []  # list to store the backward genes
    nr = []  # list to store the forward genes
    # get the address of the corresponding species to which the gene belongs
    # whose neighbor has to be found
    gi = d[gs.capitalize()]
    sldf = a[gi]  # select the dataframe
    scmap = cmap[gi]  # select the correct chromosome map
    scimap = cimap[gi]  # select the correct index maps
    try:
        _ = ld[gi]  # see if the corresponding gene map exists
    except BaseException:
        # print("Length of Dataframes:{} \t Length of Loaded Genes:{} \t Length of Loaded Genomes Dictionaries:{}".format(len(a),len(ld),len(ldg)))
        return ne, nr
    sldg = ldg[gi]  # select the corresponding map
    if g not in sldg:  # if the gene is not present in the dataframe return empty lists
        # print(g,"\t",gs)
        return ne, nr
    i = sldg[g]  # find the index of the gene
    chromosome_id = scmap[g]  # get the chromosome no from the database.
    scimap = scimap[chromosome_id]  # select the correct chromosomes indexes
    sldf = sldf.loc[scimap]  # select only the same chromosome genes.
    # get the -n neighbors
    start = int(sldf.loc[i]["start"])  # get the start location of the gene
    flag = 0
    for j in range(n):
        if flag == 1:
            ne.append("NULL_GENE")
            continue
        itemp = 0
        # select the column
        end = list(sldf.end)
        end = np.array(end)
        assert len(end) == len(sldf)
        end = end - start  # subtract start from it so as to get relative position
        end_s = np.argsort(end)  # sort them by the order of distance
        if end[end_s[0]] >= 0:  # if all the genes end ahead of the one in consideration
            flag = 1  # increment the pointer
            ne.append("NULL_GENE")  # append the NULL_GENE value
            continue
        for k in end_s:  # iterate through the sorted array
            # find the first value that is negative and the next one is
            # positive to get the nearest gene
            if end[k] < 0 and end[k + 1] >= 0:
                itemp = k
                break
        itemp = scimap[itemp]
        ne.append(sldf.loc[itemp].gene_id)
        # make "start" the start location of the current gene
        start = int(sldf.loc[itemp].start)
        # print(start)
    # get the +n neighbors
    flag = 0
    end = int(sldf.loc[i].end)
    for j in range(n):
        if flag == 1:
            nr.append("NULL_GENE")
            continue
        itemp = 0
        start = list(sldf.start)
        start = np.array(start)
        start = start - end
        start_s = np.argsort(start)
        if start[start_s[-1]] < 0:
            flag = 1
            nr.append("NULL_GENE")
            continue
        for k in start_s:
            if start[k] > 0:
                itemp = k
                break
        itemp = scimap[itemp]
        nr.append(sldf.loc[itemp].gene_id)
        end = int(sldf.loc[itemp].end)
    return ne, nr


def create_data_homology_ls(a_h, d_h, n, a, d, ld, ldg, cmap, cimap, to_write):
    # dictionary which stores +/- n genes of the given gene by id. Each key is
    # a gene id which corresponds to the one in center.
    lsy = {}
    lsytemp = {}
    name = "neighbor_genes"
    for df in a_h:
        for _, row in progressbar.progressbar(df.iterrows()):
            x = row["gene_stable_id"]
            y = row["homology_gene_stable_id"]
            xs = row["species"]
            ys = row["homology_species"]
            try:
                _ = lsy[x]
            except BaseException:
                try:
                    # see if the species exist in genomic maps
                    _ = d[xs.capitalize()]
                    xl, xr = get_nearest_neighbors(x, xs, n, a, d, ld, ldg, cmap, cimap)
                    if (
                        len(xl) != 0
                    ):  # check if neighboring genes were successfully found
                        lsy[x] = dict(b=xl, f=xr)
                        lsytemp[x] = dict(b=xl, f=xr)
                except BaseException:
                    continue
            try:
                _ = lsy[y]
            except BaseException:
                try:
                    _ = d[ys.capitalize()]
                    yl, yr = get_nearest_neighbors(y, ys, n, a, d, ld, ldg, cmap, cimap)
                    if len(yl) != 0:
                        lsy[y] = dict(b=yl, f=yr)
                        lsytemp[y] = dict(b=yl, f=yr)
                except BaseException:
                    continue
    if to_write == 1:
        write_dict_json(name, "processed", lsy)
    return lsy
