# Thomas Pranzatelli
# 11-2-2017
# Must run generate_data.py first.
# Runs neural network on Karpas-422 data.

import numpy as np
import json

def read_in_data():
	input = np.load('input.npy')
	return input


input = read_in_data()
print input.shape

import tensorflow as tf
from keras.layers import Input,Dense,Lambda
from keras.models import Model
from keras import optimizers
from keras import backend as K
from keras import metrics
sgd = optimizers.SGD(lr=0.1)

encoding_dim = 200

input_img = Input(shape=(10000,))
encoded = Dense(encoding_dim,activation='sigmoid')(input_img)
decoded = Dense(10000,activation='relu')(encoded)

autoencoder = Model(input_img,decoded)
encoder = Model(input_img,encoded)

encoded_input = Input(shape=(encoding_dim,))
decoder_layer = autoencoder.layers[-1]
decoder = Model(encoded_input,decoder_layer(encoded_input))

autoencoder.compile(optimizer=sgd,loss='binary_crossentropy',metrics=['mae','acc'])

np.random.shuffle(input)
input_train = input[:190000]
input_test = input[190000:192000]

history = autoencoder.fit(input_train,input_train,
	epochs=20,batch_size=1000,shuffle=True,
	validation_data=(input_test,input_test))


encoded_imgs = encoder.predict(input_test)
decoded_imgs = decoder.predict(encoded_imgs)
print input_test[2].tolist()
print encoded_imgs[2].tolist()
print decoded_imgs[2].tolist()


from keras.layers import Conv1D,MaxPooling1D

sgd = optimizers.SGD(lr=0.1)
input_img = Input(shape=(10000,))
conv1 = Conv1D(5,10,strides=5,activation=None)(input_img)
pooling1 = MaxPooling1D(pool_size=5)(conv1)
conv2 = Conv1D(20,25,strides=20,activation=None)(pooling1)
pooling2 = MaxPooling1D(pool_size=20)(conv2)
dense1 = Dense(100,activation='sigmoid')(pooling2)
dense2 = Dense(100,activation='sigmoid')(dense1)
dense3 = Dense(1,activation='sigmoid')(dense2)

autoencoder = Model(input_img,dense3)

autoencoder.compile(optimizer=sgd,loss='binary_crossentropy',metrics=['mae','acc'])

np.random.shuffle(input)
input_train = input[:190000]
input_test = input[190000:192000]

history = autoencoder.fit(input_train,input_train,
	epochs=20,batch_size=1000,shuffle=True,
	validation_data=(input_test,input_test))


