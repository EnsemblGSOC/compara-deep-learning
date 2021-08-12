# tf.debugging.set_log_device_placement(True)
# # Create some tensors
# a = tf.constant([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
# b = tf.constant([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
# c = tf.matmul(a, b)
# print(c)

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

data = np.load(args.npy, mmap_mode="r") # eg "Jul-16-2021_big_final.npy", and leave the data on disc

dataframe = dataframe[~(dataframe.ortho_para == "not")] # compare only homologs and paralogs to one another
# give label indexes to the homology categories
dataframe["homo_idx"] = dataframe.ortho_para.astype("category").cat.codes

# get the indices of the dataset to split
train_idx, test_idx = train_test_split(dataframe.index,test_size=0.2, random_state=42) # use 80% of the data for training
val_idx, test_idx = train_test_split(test_idx, test_size=0.5, random_state=42) # take 10% of all data for each of test and dev set

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


# set up two dataset objects for training the model
dataset = tf.data.Dataset.from_generator(npy_generator, 
                                        output_signature=(tf.TensorSpec(shape=( 7, 7,11), dtype=tf.float32), 
                                        tf.TensorSpec(shape=(), dtype=tf.float32) ), 
                                        args=[train_idx] )

val_dataset = tf.data.Dataset.from_generator(npy_generator, 
                                            output_signature=(tf.TensorSpec(shape=( 7, 7,11), dtype=tf.float32), 
                                                tf.TensorSpec(shape=(), dtype=tf.float32) ),
                                            args=[val_idx] )

# specify batch size for sgd, dropping the remaining values in the final batch
batch_size = 64
dataset = dataset.batch(batch_size, drop_remainder=True)
# make the dataset infinitley repeat so model can be train for arbitrary batch sizes
dataset = dataset.repeat()
# Batch the valiation set data generator which will be used once per epoch
val_dataset = val_dataset.batch(batch_size)

# setup the default training callbacks

# Learning rate decay scheduler
def scheduler(epoch, lr):
  if epoch < 10:
    return lr
  else:
    return lr * 1/2

callback = tf.keras.callbacks.LearningRateScheduler(scheduler)

stopper = EarlyStopping(monitor="val_loss", min_delta = 0.005, patience = 4)

# use this callback whilst training all models
callback_list = [callback,stopper]

# get the model
model = get_model(args.model_name)




# see the architecture
print(model.summary())





# handy function for plotting metrics
def plot_history(history, metrics_list, savepath):
    """
    plot the learning curve for the history object
    """
    for metric in metrics_list:
        plt.plot(history.history[metric], label = metric)
    plt.legend()
    plt.savefig(savepath)

# get the metrics
precis = tf.keras.metrics.Precision()
recall = tf.keras.metrics.Recall()
AUC = tf.keras.metrics.AUC()
categorical_loss = SparseCategoricalCrossentropy()
# compile the model
model.compile(optimizer="adam", loss=categorical_loss, metrics=["accuracy"])
# train the model using dataset generators, and tell the trainer how many examples to expect in an epoch
history = model.fit(dataset, 
                    epochs=100,
                    steps_per_epoch=train_idx.shape[0] // batch_size,  
                    validation_data=val_dataset,
                    callbacks=[callback, stopper],
                    verbose=2)

# save the model
model.save(args.out_dir + "/" + args.model_name)
# save the loss curves
plot_history(history, ["loss","val_loss"], args.out_dir + "/" + args.model_name + "/" + args.model_name + "_loss_curves.png")
# save the accuracy curves
plot_history(history, ["accuracy","val_accuracy"], args.out_dir + "/" + args.model_name + "/" + args.model_name + "_accuracy_curves.png")
