from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import faulthandler
import os
import time
import numpy as np
import tensorflow as tf
from tensorflow import keras
from google import protobuf
print("Tensorflow version: ", tf.version.VERSION)
print("Protobuf version:", protobuf.__version__)

gpu_devices = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(gpu_devices[0], True)
tf.config.experimental.set_virtual_device_configuration(
            gpu_devices[0],
            [tf.config.experimental.VirtualDeviceConfiguration(
               memory_limit=2048)]) ## Crucial value, set lower than available GPU memory (note that Jetson shares GPU memory with CPU)

from tensorflow.python.compiler.tensorrt import trt_convert as trt

import numpy as np
import tensorflow as tf
print(tf.__version__)
print('GPU available:')
print(tf.config.list_physical_devices('GPU'))


from load_data import z_score_normalization, downsample_data, generate_overlapping_windows

# Downsample data
fs = 2000
downsampled_fs = 1250
overlapping = True
window_size = 0.0128

print("Generating windows...", end=" ")

# Overlapping window size
stride = 0.0064

nSamples=1000

# Separate the data into 12.8ms windows with 6.4ms overlapping
#Cast to 32bit precision because of optimization issues with tensorrt
X = tf.convert_to_tensor(generate_overlapping_windows(z_score_normalization(downsample_data(np.load("lfp_data.npy"), fs, downsampled_fs)), window_size, stride, downsampled_fs)[:nSamples, :, :], dtype=float)
print("X shape: ", X.shape)

print("Done!")

import tensorflow.keras as kr
print("Loading CNN model...", end=" ")
optimizer = kr.optimizers.Adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-07, amsgrad=False)
model = kr.models.load_model("./model", compile=False)
model.compile(loss="binary_crossentropy", optimizer=optimizer)
print("Done!")

nBurnIn=100
for i in range(nSamples):
	if(i==nBurnIn):
		startTime=time.time()
	model.predict(X[i:i+1,:,:], verbose=False)
endTime=time.time()
avgTime=(endTime-startTime)/(nSamples-nBurnIn)

#print("Detecting ripples...", end=" ")
#startTime=time.time()
#model.predict(X, verbose=False)
#endTime=time.time()
#avgTime=(endTime-startTime)/nSamples
print("Average time:",avgTime)
print("Done!")



