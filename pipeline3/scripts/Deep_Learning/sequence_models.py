import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import tensorflow as tf

from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.preprocessing.text import Tokenizer
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, SimpleRNN, LSTM, Embedding, Bidirectional, Input
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.losses import CategoricalCrossentropy, SparseCategoricalCrossentropy


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

def scheduler(epoch, lr):
  if epoch < 10:
    return lr
  else:
    return lr * 1/2

callback = tf.keras.callbacks.LearningRateScheduler(scheduler)

stopper = EarlyStopping(monitor="val_loss", min_delta = 0.005, patience = 4)

from sklearn.model_selection import train_test_split

df = pd.read_csv("sequence_pairs.csv")

type(df.homolog_sequence.iloc[0])

data = df[~df.species_sequence.isna() & ~df.homolog_sequence.isna()]

data["homo_idx"] = data.ortho_para.astype("category").cat.codes

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

# Create a Tokenizer object

tokenizer = Tokenizer(num_words=None, 
                      filters='!"#$%&()*+,-./:;<=>?@[\\]^_`{|}~\t\n',
                      lower=True,
                      split=None,
                      char_level=True,
                      document_count=0)

tokenizer.fit_on_texts(data.species_sequence)

tokenizer_config = tokenizer.get_config()
tokenizer_config.keys()

tokenizer_config['word_counts']

import json

word_counts = json.loads(tokenizer_config['word_counts'])

tokenizer_config['word_index']

index_word = json.loads(tokenizer_config['index_word'])
word_index = json.loads(tokenizer_config['word_index'])

species_seq = tokenizer.texts_to_sequences(data.species_sequence)
homo_seq = tokenizer.texts_to_sequences(data.homolog_sequence)

lens = [ len(seq) for seq in species_seq]

max(lens)

lens = [ len(seq) for seq in homo_seq]

max(lens)

data.homo_idx.drop_duplicates()

species_padded_seq = pad_sequences(species_seq, 100, padding="post", truncating="post")
homo_padded_seq = pad_sequences(homo_seq, 100, padding="post", truncating="post")

labels = tf.keras.utils.to_categorical(data.homo_idx)

# use species as labels
labels = data.species_idx

train_Species, test_Species, train_Homo, test_Homo, train_labels, test_labels = train_test_split(species_padded_seq,homo_padded_seq, labels)

def plot_history(history, metrics_list):
    """
    plot the learning curve for the history object
    """
    for metric in metrics_list:
        plt.plot(history.history[metric], label = metric)
    plt.legend()
    plt.show()

embed = Embedding(23, 23, mask_zero=True)

train_X.shape

Recurrent = Sequential( [ Embedding(24, 10,  mask_zero=True),
                     Bidirectional(LSTM(units = 100, activation="tanh", return_sequences=True )),
                     Bidirectional(LSTM(units = 200, activation="tanh" )),
                     Dense(20)])
#                     Dense(28, activation="softmax")])

for layer in Recurrent.layers:
    print(layer.name)
    print(layer.input_shape)
    print(layer.output_shape)

Recurrent.summary()

input1 = Input(shape=(100))
input2 = Input(shape=(100))

encoding1 = Recurrent(input1)
encoding2 = Recurrent(input2)

difference = tf.math.subtract(
    encoding1,encoding2, name="difference"
)
absolute = tf.math.abs(
    difference, name="absolute"
)

output = Dense(2,activation="softmax")(absolute)


model = Model(inputs=[input1,input2], outputs=output)

model.summary()

categorical_loss = CategoricalCrossentropy()

model.compile(optimizer="adam", loss=categorical_loss, metrics=["accuracy",AUC])

precis = tf.keras.metrics.Precision()
recall = tf.keras.metrics.Recall()
AUC = tf.keras.metrics.AUC()

categorical_loss = CategoricalCrossentropy()

model.compile(optimizer="adam", loss=categorical_loss, metrics=["accuracy",AUC])

history = model.fit([train_Species,train_Homo], train_labels, validation_split=0.1, epochs  = 10, batch_size=128, callbacks=[stopper,callback])