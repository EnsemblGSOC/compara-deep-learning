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
        x2 = np.array(dataframe["sp_distance"].loc[i])[np.newaxis]
        y = np.array(dataframe["homo_idx"].loc[i])
        yield {"input1":tf.convert_to_tensor(x),"input2":tf.convert_to_tensor(x2)}, tf.convert_to_tensor(y)


# set up two dataset objects for training the model
dataset = tf.data.Dataset.from_generator(npy_generator, 
                                        output_signature=({"input1":tf.TensorSpec(shape=( 7, 7,11), dtype=tf.float32),
                                        "input2":tf.TensorSpec(shape=(1), dtype=tf.float32)},
                                        tf.TensorSpec(shape=(), dtype=tf.float32) ), 
                                        args=[train_idx] )

val_dataset = tf.data.Dataset.from_generator(npy_generator, 
                                        output_signature=({"input1":tf.TensorSpec(shape=( 7, 7,11), dtype=tf.float32),
                                        "input2":tf.TensorSpec(shape=(1), dtype=tf.float32)},
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

def split_feature_model():
    """
    This model attempts to use a separate convnet on each of the different features.
    This should allow each distinct data type to learn its own useful feature map.
    """
    # create an intermediary model which allows
    def Simple_model_filters(filters):
        """
        Like the simple model above, but with more MLP hidden units at the end
        """
        model = Sequential( [Conv2D(20,(2,2), activation="relu", strides=(1,1), padding="same", input_shape=(7,7,filters), name="conv2d_1"),
                    MaxPool2D(2,2, name="max_pooling2d_1"),
                        Conv2D(20,(2,2), activation="relu", strides=(1,1), padding="same", name="conv2d_2"),
                        Conv2D(20,(2,2), activation="relu", strides=(1,1), padding="same"),                     
                        MaxPool2D(2,2, name="max_pooling2d_2"),
                        Flatten(name="Flatten"),
                        Dense(100, activation="relu", name = "dense1", kernel_regularizer="l2"),
                        Dense(50, activation="relu", name = "dense2", kernel_regularizer="l2")])
                        
        return model

    # combine the intermediary model with the custom numbers of filters for each feature
    inputs = Input(shape = (7,7,11), name = "input1")
    input2 = Input(shape = (1,), name = "input2")
    h1a = Simple_model_filters(1)(tf.expand_dims(inputs[:,:,:,0], axis=-1))
    h1b = Simple_model_filters(2)(inputs[:,:,:,1:3])
    h2 = Simple_model_filters(2)(inputs[:,:,:,3:5])
    h3a = Simple_model_filters(1)(tf.expand_dims(inputs[:,:,:,5], axis=-1))
    h3b = Simple_model_filters(2)(inputs[:,:,:,6:8])
    h4a = Simple_model_filters(1)(tf.expand_dims(inputs[:,:,:,8], axis =-1))
    h4b = Simple_model_filters(2)(inputs[:,:,:,9:11])
    h =  Concatenate()([h1a,h1b,h2,h3a,h3b,h4a,h4b])
    h = Dense(40)(h)
    h = Concatenate(axis=-1)([h,input2])
    outputs= Dense(2, activation="softmax")(h)
    model = Model(inputs=[inputs, input2],outputs=outputs)
    return model

model = split_feature_model()





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
