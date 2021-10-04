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

from models import get_model

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

data = np.load(args.npy) # eg "Jul-16-2021_big_final.npy"

dataframe = dataframe[~(dataframe.ortho_para == "not")]
data = data[dataframe.index]

# handy function for plotting metrics
def plot_history(history, metrics_list, savepath):
    """
    plot the learning curve for the history object
    """
    for metric in metrics_list:
        plt.plot(history.history[metric], label = metric)
    plt.legend()
    plt.savefig(savepath)


# Learning rate decay scheduler
def scheduler(epoch, lr):
  if epoch < 10:
    return lr
  else:
    return lr * 1/2

# setup the default

callback = tf.keras.callbacks.LearningRateScheduler(scheduler)
stopper = EarlyStopping(monitor="val_loss", min_delta = 0.005, patience = 4)

# use this callback whilst training all models
callback_list = [callback,stopper]
#clean the missing values from the data
data = np.where( ~(data == "NP"), data, 0.)
data = np.where( ~(data == "NA"), data, 0.)
data = data.astype(float)
data = np.nan_to_num(data,nan=0.0)
mean = data.mean(axis = (0,1,2))
std = data.std(axis = (0,1,2))
# scale the data using broadcasting
data = ((data - mean)/std)
# give label indexes to the homology categories
dataframe["homo_idx"] = dataframe.ortho_para.astype("category").cat.codes
# save the labels
labels = tf.keras.utils.to_categorical(dataframe.homo_idx)
# train test split
train_X, test_X, train_labels, test_labels = train_test_split(data, labels)


# get the model
model = get_model(args.model_name)


# see the architecture
print(model.summary())
# get the metrics
precis = tf.keras.metrics.Precision()
recall = tf.keras.metrics.Recall()
AUC = tf.keras.metrics.AUC()
categorical_loss = CategoricalCrossentropy()
# compile the model
model.compile(optimizer="adam", loss=categorical_loss, metrics=["accuracy"])
# train the model
history = model.fit(train_X, 
                    train_labels,
                    validation_split=0.1, 
                    epochs  = 100, # probably too many but the stopper call back should terminate training early
                    batch_size=256, # always have match size as power of 2 for efficient memory use
                    callbacks=callback_list)
# save the model
model.save(args.out_dir + "/" + args.model_name)
# save the loss curves
plot_history(history, ["loss","val_loss"], args.out_dir + "/" + args.model_name + "/" + args.model_name + "_loss_curves.png")
# save the accuracy curves
plot_history(history, ["accuracy","val_accuracy"], args.out_dir + "/" + args.model_name + "/" + args.model_name + "_accuracy_curves.png")
