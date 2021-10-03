# tf.debugging.set_log_device_placement(True)
# # Create some tensors
# a = tf.constant([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
# b = tf.constant([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
# c = tf.matmul(a, b)
# print(c)

"""
Evaluate the output from the distance_trainer.py script
"""


import tensorflow as tf
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import glob
import argparse

from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Dense, Conv2D, Flatten, MaxPool2D, Dropout, Input, Concatenate, Average
from tensorflow.keras.losses import SparseCategoricalCrossentropy
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.regularizers import l2

from models import get_model


from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession




def fix_gpu():
    config = ConfigProto()
    config.gpu_options.allow_growth = True
    session = InteractiveSession(config=config)


fix_gpu()

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
    help="the model to train",
    type=str
)

parser.add_argument( 
    "--out_dir",
    help="output directory to save to",
    type=str
)

args = parser.parse_args()

dataframe = pd.read_csv(args.csv) # eg "Jul-16-2021_big_final.csv"
# add the species distance to the dataframe
species = np.loadtxt("scripts/sp_names", dtype="str")
distance = pd.read_csv("scripts/dist_matrix", delimiter="\t", names=species).set_index(species)
dataframe["query_species"].fillna(dataframe["homology_species"], inplace=True)
distance[dataframe.species[0]].loc[dataframe.query_species[0]]
dataframe["sp_distance"] = dataframe[["species","query_species"]].apply(lambda x: distance[x.species].loc[x.query_species], axis =1 )

data = np.load(args.npy, mmap_mode="r") # eg "Jul-16-2021_big_final.npy", and leave the data on disc

dataframe = dataframe[~(dataframe.ortho_para == "not")] # compare only homologs and paralogs to one another
# give label indexes to the homology categories
dataframe["homo_idx"] = dataframe.ortho_para.astype("category").cat.codes

# load the model
model = tf.keras.models.load_model(args.out_dir + "/" + args.model_name)

# see the architecture
print(model.summary())



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
        x2 = np.array(dataframe["sp_distance"].loc[i])[np.newaxis]
        y = np.array(dataframe["homo_idx"].loc[i])
        yield {"input1":tf.convert_to_tensor(x),"input2":tf.convert_to_tensor(x2)}, tf.convert_to_tensor(y)



# distance_dict = {}

# for distance in dataframe.sp_distance.drop_duplicates().sort_values():

#     idx = dataframe.loc[dataframe.sp_distance == distance].index

#     # evaluation_labels = dataframe.iloc[idx].homo_idx
#     # set up two dataset objects for training the model
#     dataset = tf.data.Dataset.from_generator(npy_generator, 
#                                         output_signature=({"input1":tf.TensorSpec(shape=( 7, 7,11), dtype=tf.float32),
#                                         "input2":tf.TensorSpec(shape=(1), dtype=tf.float32)},
#                                         tf.TensorSpec(shape=(), dtype=tf.float32) ), 
#                                         args=[idx] )
#     dataset = dataset.batch(64)
#     distance_dict[distance] = model.evaluate(dataset, verbose = 2)[1]

# import json
# print("saving distance data")
# with open('distance.json', 'w') as f:
#     json.dump(distance_dict, f)

print("running the human comparisons")

homo_data =  dataframe[(dataframe.species == "homo_sapiens") | (dataframe.query_species == "homo_sapiens" )]

homo_distance_dict = {}

for distance in homo_data.sp_distance.drop_duplicates().sort_values():

    idx = homo_data.loc[homo_data.sp_distance == distance].index

    # evaluation_labels = homo_data.iloc[idx].homo_idx
    # set up two dataset objects for training the model
    dataset = tf.data.Dataset.from_generator(npy_generator, 
                                        output_signature=({"input1":tf.TensorSpec(shape=( 7, 7,11), dtype=tf.float32),
                                        "input2":tf.TensorSpec(shape=(1), dtype=tf.float32)},
                                        tf.TensorSpec(shape=(), dtype=tf.float32) ), 
                                        args=[idx] )
    dataset = dataset.batch(64)
    homo_distance_dict[distance] = model.evaluate(dataset, verbose = 2)[1]

import json

print("saving homo distance data")
with open('homo_distance.json', 'w') as f:
    json.dump(homo_distance_dict, f)