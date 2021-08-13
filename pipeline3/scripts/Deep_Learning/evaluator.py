import tensorflow as tf
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import glob
import argparse

from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Dense, Conv2D, Flatten, MaxPool2D, Dropout, Input, Concatenate, Average
from tensorflow.keras.losses import CategoricalCrossentropy
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.regularizers import l2
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument( 
    "--csv",
    help="path to csv with the species and gene pair labels",
    type=str
)
parser.add_argument( 
    "--npy",
    help="path to npy file with matrix data for gene pairs",
    type=str
)

parser.add_argument( 
    "--model_name",
    help="path to the model directory",
    type=str
)

parser.add_argument( 
    "--models_dir",
    help="path to the model directory",
    type=str
)

parser.add_argument( 
    "--out_dir",
    help="output directory to save to",
    type=str
)

args = parser.parse_args()

data = np.load(args.npy, mmap_mode="r") # eg "Jul-16-2021_big_final.npy", and leave the data on disc

dataframe = pd.read_csv(args.csv) # eg "Jul-16-2021_big_final.csv"
dataframe = dataframe[~(dataframe.ortho_para == "not")] # compare only homologs and paralogs to one another
# give label indexes to the homology categories
dataframe["homo_idx"] = dataframe.ortho_para.astype("category").cat.codes

# load the trained moedel
model = tf.keras.models.load_model(args.models_dir + "/" + args.model_name)

def npy_generator(index):
    """
    Dataset generator for reading np data directly on disc instead of loading into memory
    """
    # use the indices to iterate through the dataset
    for i in index:
        x = data[i]
        x = np.where( ~(x == "NP"), x, 0.)
        x = np.where( ~(x == "NA"), x, 0.)
        x = x.astype(float)
        x = np.nan_to_num(x,nan=0.0)
        y = np.array(dataframe["homo_idx"].loc[i])
        yield tf.convert_to_tensor(x), tf.convert_to_tensor(y)

accuracy_dict = {}

for species in dataframe.species.drop_duplicates().sort_values():

    idx = dataframe.loc[dataframe.species == species].index

    # evaluation_labels = dataframe.iloc[idx].homo_idx
    # set up two dataset objects for training the model
    dataset = tf.data.Dataset.from_generator(npy_generator, 
                                        output_signature=(tf.TensorSpec(shape=( 7, 7,11), dtype=tf.float32), 
                                        tf.TensorSpec(shape=(), dtype=tf.float32) ), 
                                        args=[idx] )
    dataset = dataset.batch(64)


    accuracy_dict[species] = model.evaluate(dataset, verbose = 2)[1]

accuracy = pd.Series(accuracy_dict)
fig, axs = plt.subplots(1,1,figsize=(15,7))
accuracy.plot.bar(ax=axs)
plt.xticks(rotation=90)
plt.title("Accuracy to distriguish homologs and between species paralogs", size= 20)
plt.ylabel("Accuracy", size= 20)
plt.xlabel("Species", size= 20)
plt.savefig(args.out_dir + "/" + args.model_name + ".png")

# direct species pairwise comparisons
pairs = dataframe[["species","homology_species"]].drop_duplicates().apply( lambda x: tuple(sorted([x[0], x[1]])), axis = 1)
pairs = set(pairs)
rows, cols = sorted(list(set([pair[0] for pair in pairs]))), sorted(list(set([pair[0] for pair in pairs])))
temp_array = np.zeros(shape=(len(rows),len(rows)))
df = pd.DataFrame(data=temp_array, index=rows, columns=cols) # set up an empty dataframe to make all pairwise comparisons


for i, j in pairs:
    idx = dataframe[((dataframe.species == i) & (dataframe.homology_species == j)) | 
                    ((dataframe.species == j) & (dataframe.homology_species == i)) ].index
#     idx = dataframe.loc[dataframe.species == species].index

    # evaluation_labels = dataframe.iloc[idx].homo_idx
    # set up two dataset objects for training the model
    dataset = tf.data.Dataset.from_generator(npy_generator, 
                                        output_signature=(tf.TensorSpec(shape=( 7, 7,11), dtype=tf.float32), 
                                        tf.TensorSpec(shape=(), dtype=tf.float32) ), 
                                        args=[idx] )
    dataset = dataset.batch(64)


    df[i][j] = model.evaluate(dataset)[1]

fig, axs = plt.subplots(1,1, figsize=(15,15))
sns.heatmap(df.values, xticklabels=rows, yticklabels=rows, ax = axs)
plt.savefig(args.out_dir + "/" + args.model_name + "_comparison_matrix.png")