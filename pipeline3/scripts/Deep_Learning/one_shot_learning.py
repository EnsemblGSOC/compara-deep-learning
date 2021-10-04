import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import tensorflow as tf
import argparse

from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.preprocessing.text import Tokenizer
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, SimpleRNN, LSTM, Embedding, Bidirectional, Input
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.losses import CategoricalCrossentropy, SparseCategoricalCrossentropy

parser = argparse.ArgumentParser()
parser.add_argument( 
    "--csv",
    help="csv with the species and gene info",
)
parser.add_argument( 
    "--out_dir",
    help="output directory for syneny matrices",
)
args = parser.parse_args()


def scheduler(epoch, lr):
  if epoch < 5:
    return lr
  else:
    return lr * 1/2

callback = tf.keras.callbacks.LearningRateScheduler(scheduler)

stopper = EarlyStopping(monitor="val_loss", min_delta = 0.005, patience = 5)

from sklearn.model_selection import train_test_split

df = pd.read_csv(args.csv)
print("Original data frame shape" )
print(df.shape)

data = df[~df.sequence.isna() & ~df.query_sequence.isna()]

data["homo_idx"] = data.ortho_para.astype("category").cat.codes

print("Number of homo categories")
print(data["homo_idx"].drop_duplicates())
print("\n")
print("Data size")
print(data.shape)

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

# Create a Tokenizer object

tokenizer = Tokenizer(num_words=None, 
                      filters='!"#$%&()*+,-./:;<=>?@[\\]^_`{|}~\t\n',
                      lower=True,
                      split=None,
                      char_level=True,
                      document_count=0)

tokenizer.fit_on_texts(data.sequence)

tokenizer_config = tokenizer.get_config()
tokenizer_config.keys()

import json

print(tokenizer_config['word_index'])

word_counts = json.loads(tokenizer_config['word_counts'])

index_word = json.loads(tokenizer_config['index_word'])
word_index = json.loads(tokenizer_config['word_index'])

species_seq = tokenizer.texts_to_sequences(data.sequence)
homo_seq = tokenizer.texts_to_sequences(data.query_sequence)

lens = [ len(seq) for seq in species_seq]

max(lens)

lens = [ len(seq) for seq in homo_seq]

max(lens)

data.homo_idx.drop_duplicates()

species_padded_seq = pad_sequences(species_seq, 125, padding="post", truncating="post")
homo_padded_seq = pad_sequences(homo_seq, 125, padding="post", truncating="post")

labels = tf.keras.utils.to_categorical(data.homo_idx)

train_Species, test_Species, train_Homo, test_Homo, train_labels, test_labels = train_test_split(species_padded_seq,homo_padded_seq, labels)

def plot_history(history, metrics_list):
    """
    plot the learning curve for the history object
    """
    for metric in metrics_list:
        plt.plot(history.history[metric], label = metric)
    plt.legend()
    plt.savefig(str(metrics_list[0]) + ".png")
    plt.show()

# Define a recurrent Model
Recurrent = Sequential( [ Embedding(26, 26,  mask_zero=True),
                     Bidirectional(LSTM(units = 125, activation="tanh", return_sequences=True )),
                     Bidirectional(LSTM(units = 250, activation="tanh" ))
                     Dense(30)])
#                     Dense(28, activation="softmax")])

for layer in Recurrent.layers:
    print(layer.name)
    print(layer.input_shape)
    print(layer.output_shape)

Recurrent.summary()

input1 = Input(shape=(125))
input2 = Input(shape=(125))

encoding1 = Recurrent(input1)
encoding2 = Recurrent(input2)

difference = tf.math.subtract(
    encoding1,encoding2, name="difference"
)
absolute = tf.math.abs(
    difference, name="absolute"
)

output = Dense(3,activation="softmax")(absolute)


model = Model(inputs=[input1,input2], outputs=output)

model.summary()

categorical_loss = CategoricalCrossentropy()



precis = tf.keras.metrics.Precision()
recall = tf.keras.metrics.Recall()
AUC = tf.keras.metrics.AUC()
model.compile(optimizer="adam", loss=categorical_loss, metrics=["accuracy",AUC])


categorical_loss = CategoricalCrossentropy()

model.compile(optimizer="adam", loss=categorical_loss, metrics=["accuracy",AUC])
history = model.fit([train_Species,train_Homo], train_labels, validation_split=0.1, epochs  = 100, batch_size=128, callbacks=[stopper,callback])