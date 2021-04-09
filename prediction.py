import pandas as pd
import numpy as np
import tensorflow as tf
import sys
import progressbar
from neighbor_genes import read_genome_maps
from process_data import create_data_homology_ls
from read_get_gene_seq import read_gene_sequences
from access_data_rest import update_rest
from threads import Procerssrunner
from prepare_synteny_matrix import read_data_synteny
from tree_data import create_tree_data
from process_data import create_map_list


def read_database(fname):
    df = pd.read_csv(fname, sep="\t", header=None)
    label_dict = dict(ortholog_one2one=1,
                      other_paralog=0,
                      non_homolog=2,
                      ortholog_one2many=1,
                      ortholog_many2many=1,
                      within_species_paralog=0,
                      gene_split=4)
    label = []
    for _, row in df.iterrows():
        label.append(label_dict[row[7]])
    df = df.assign(label=label)
    df = df.drop(7, axis=1)
    df = df.drop(0, axis=1)
    df.columns = [
        "gene_stable_id",
        "species",
        "homology_gene_stable_id",
        "homology_species",
        "goc",
        "wga",
        "label"]
    return df


def select_data_by_length(df, st, end):
    try:
        if end < len(df):
            if st < end:
                df = df.loc[df.index.values[st:end]]
            else:
                raise ValueError()
    except BaseException:
        print("Making Predictions for the complete dataframe:)")
    print(len(df))
    return df


def create_synteny_features(a_h, d_h, n, a, d, ld, ldg, cmap, cimap, name):
    lsy = create_data_homology_ls(a_h, d_h, n, a, d, ld, ldg, cmap, cimap, 0)
    gene_sequences = read_gene_sequences(
        a_h, lsy, "geneseq", "prediction_" + name)
    gene_sequences = update_rest(gene_sequences, "prediction_" + name)
    print("Gene Sequences Loaded.")
    return lsy, gene_sequences


def threadmaker(nop, df, lsy, gene_sequences, n, name):
    part = len(df) // nop
    pr = Procerssrunner()
    pr.start_processes(nop, df, gene_sequences, lsy, part, n, name)
    smg, sml, indexes = read_data_synteny(nop, name)
    sml = np.array(sml)
    smg = np.array(smg)
    indexes = np.array(indexes)
    return sml, smg, indexes


def get_prediction(smg, sml, indexes, bls, blhs, dis, dps, dphs, model_name):
    preds = np.zeros((len(smg), 3))
    w = [0.86, 0.8, 0.06]
    for i in range(1, 4):
        try:
            model = tf.train.import_meta_graph(
                model_name + '_v' + str(i) + '/model.ckpt.meta')
        except BaseException:
            print("Something wrong with the model.")
            continue
        with tf.Session() as sess:
            try:
                model.restore(sess, model_name + '_v' + str(i) + "/model.ckpt")
                graph = tf.get_default_graph()
                synmgt, synmlt, blst, \
                    blhst, dpst, dphst, \
                    dist, lrt, yt = graph.get_collection("input_nodes")
                predictions = graph.get_tensor_by_name("Predictions/BiasAdd:0")
                print("Model Loaded Successfully :)")
            except BaseException:
                print(":(")
                sys.exit()

            fd = {synmgt: smg,
                  synmlt: sml,
                  blst: bls,
                  blhst: blhs,
                  dpst: dps.reshape((len(blhs), 1)),
                  dist: dis.reshape((len(blhs), 1)),
                  dphst: dphs.reshape((len(blhs), 1))}
            preds_t_1 = sess.run([predictions], feed_dict=fd)
            preds_t_1 = np.array(preds_t_1)[0]
            fd = {synmgt: smg.transpose((0, 2, 1, 3)),
                  synmlt: sml.transpose((0, 2, 1, 3)),
                  blst: blhs,
                  blhst: bls,
                  dpst: dphs.reshape((len(blhs), 1)),
                  dist: dis.reshape((len(blhs), 1)),
                  dphst: dps.reshape((len(blhs), 1))}
            preds_t_2 = sess.run([predictions], feed_dict=fd)
            preds_t_2 = np.array(preds_t_2)[0]
            preds = preds + w[i - 1] * (preds_t_1 + preds_t_2) / 2
        tf.reset_default_graph()
    preds = np.argmax(preds, axis=1)
    print(preds.shape)
    return preds


def write_preds(fname, model_name, name, preds, index_dict, df):
    print(
        "Writing predcitions to:",
        "prediction_" +
        fname +
        "_" +
        model_name +
        "_" +
        name +
        "_multiple.txt")
    with open("prediction_" + fname + "_" + model_name + "_" + name +
              "_multiple.txt", "w") as file:
        for index, row in progressbar.progressbar(df.iterrows()):
            file.write(str(row[0]))
            file.write("\t")
            file.write(row[1])
            file.write("\t")
            file.write(row[3])
            file.write("\t")
            file.write(str(row["label"]))
            file.write("\t")
            if index in index_dict:
                file.write(str(preds[index_dict[index]]))
                file.write("\t")
                if preds[index_dict[index]] == row["label"]:
                    file.write(str(1))
                else:
                    file.write(str(0))
            else:
                file.write("Error")
                file.write("\t")
                file.write("NaN")
            file.write("\n")


def main():
    arg = sys.argv
    fname = arg[-6]
    model_name = arg[-5]
    nop = int(arg[-4])
    st = int(arg[-3])
    end = int(arg[-2])
    name = arg[-1]
    df = read_database(fname)
    df = select_data_by_length(df, st, end)
    n = 3

    a, d, ld, ldg, cmap, cimap = read_genome_maps()  # read the genome mapd
    print("Genome Maps Loaded.")
    a_h = [df]
    d_h = ["prediction"]

    lsy, gene_sequences = create_synteny_features(
        a_h, d_h, n, a, d, ld, ldg, cmap, cimap, name)
    sml, smg, indexes = threadmaker(nop, df, lsy, gene_sequences, n, name)
    df_temp = df.loc[indexes]
    bls, blhs, dis, dps, dphs = create_tree_data("species_tree.tree", df_temp)
    index_dict = create_map_list(indexes)
    preds = get_prediction(
        smg,
        sml,
        indexes,
        bls,
        blhs,
        dis,
        dps,
        dphs,
        model_name)

    write_preds(fname, model_name, name, preds, index_dict, df)


if __name__ == "__main__":
    main()
