import tensorflow as tf
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import glob
import argparse

from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import (
    Dense,
    Conv2D,
    Flatten,
    MaxPool2D,
    Dropout,
    Input,
    Concatenate,
    Average,
)
from tensorflow.keras.losses import CategoricalCrossentropy
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.regularizers import l2

# the number of categories to classify in the data
n = 2

def simple_model():
    """
    Simple model architecture that doesn't separate features
    """
    model = Sequential( [Conv2D(20,(4,4), activation="relu", strides=(1,1), padding="same", input_shape=(7,7,11)),
                        MaxPool2D(2,2),
                        Conv2D(20,(4,4), activation="relu", strides=(1,1), padding="same"),
                        Conv2D(20,(4,4), activation="relu", strides=(1,1), padding="same"),                     
                        MaxPool2D(2,2),
                        Flatten(),
                        Dense(20, activation="relu"),
                        Dense(20, activation="relu"),
                        Dense(n, activation="softmax")] 
    )
    return model

def ensemble_model():
    """
    Stacks 3 layers of the original model and concatenates their input.
    In theory, this should exploit different parameter values arising due to random initalisation.
    To do this, we have to make use of the functional API to build the model.
    """
    inputs = Input(shape = (7,7,11), name = "inputs")
    h1 = simple_model()(inputs)
    h2 = simple_model()(inputs)
    h3 = simple_model()(inputs)
    h =  Concatenate()([h1,h2,h3])
    h = Dense(20)(h)
    outputs= Dense(n, activation="softmax")(h)
    model = Model(inputs=inputs, outputs=outputs)
    return model

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
        model = Sequential( [Conv2D(100,(2,2), activation="relu", strides=(1,1), padding="same", input_shape=(7,7,filters), name="conv2d_1"),
                    MaxPool2D(2,2, name="max_pooling2d_1"),
                        Conv2D(100,(2,2), activation="relu", strides=(1,1), padding="same", name="conv2d_2"),
                        Conv2D(100,(2,2), activation="relu", strides=(1,1), padding="same"),                     
                        MaxPool2D(2,2, name="max_pooling2d_2"),
                        Flatten(name="Flatten"),
                        Dense(100, activation="relu", name = "dense1", kernel_regularizer="l2"),
                        Dense(50, activation="relu", name = "dense2", kernel_regularizer="l2")])
                        
        return model

    # combine the intermediary model with the custom numbers of filters for each feature
    inputs = Input(shape = (7,7,11), name = "inputs1")
    # inputs1a = Input(shape =  (7,7,1), tf.gather(inputs, [0], axis = -1), name = "inputs1a")
    # inputs1b = Input(shape =  (7,7,2),   tf.gather(inputs, [1,2], axis = -1), name = "inputs1b")
    # inputs2 = Input(shape =  (7,7,2),   tf.gather(inputs, [3,4], axis = -1), name = "inputs2")
    # inputs3a = Input(shape =  (7,7,1), tf.gather(inputs, [5], axis = -1), name = "inputs3a")
    # inputs3b = Input(shape =  (7,7,2),   tf.gather(inputs, [6,7], axis = -1), name = "inputs3b")
    # inputs4a = Input(shape =  (7,7,1),   tf.gather(inputs, [8], axis = -1), name = "inputs4a")
    # inputs4b = Input(shape =  (7,7,2),    tf.gather(inputs, [9,10], axis = -1), name = "inputs4b")

    h1a = Simple_model_filters(1)(tf.expand_dims(inputs[:,:,:,0], axis=-1))
    h1b = Simple_model_filters(2)(inputs[:,:,:,1:3])
    h2 = Simple_model_filters(2)(inputs[:,:,:,3:5])
    h3a = Simple_model_filters(1)(tf.expand_dims(inputs[:,:,:,5], axis=-1))
    h3b = Simple_model_filters(2)(inputs[:,:,:,6:8])
    h4a = Simple_model_filters(1)(tf.expand_dims(inputs[:,:,:,8], axis =-1))
    h4b = Simple_model_filters(2)(inputs[:,:,:,9:11])
    h =  Concatenate()([h1a,h1b,h2,h3a,h3b,h4a,h4b])
    h = Dense(40)(h)
    outputs= Dense(n, activation="softmax")(h)
    model = Model(inputs=inputs,outputs=outputs )
    return model

    
models = {"simple_model": simple_model(), "ensemble_model":ensemble_model(),"split_feature_model":split_feature_model()}

model = split_feature_model()

def get_model(model_name):
    return models[model_name]
