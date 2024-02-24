from __future__ import absolute_import, division, print_function, unicode_literals

import faulthandler
import os
import time
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.python.saved_model import tag_constants
from tensorflow.keras.preprocessing import image
from google import protobuf
print("Tensorflow version: ", tf.version.VERSION)
print("Protobuf version:", protobuf.__version__)

gpu_devices = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(gpu_devices[0], True)
tf.config.experimental.set_virtual_device_configuration(
            gpu_devices[0],
            [tf.config.experimental.VirtualDeviceConfiguration(
               memory_limit=2048)]) ## Crucial value, set lower than available GPU memory (note that Jetson shares GPU memory with CPU)

import sys
import os
from tensorflow.python.compiler.tensorrt import trt_convert as trt

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

params = trt.DEFAULT_TRT_CONVERSION_PARAMS
params = params._replace(max_workspace_size_bytes=(1<<30)) 
params = params._replace(precision_mode="FP32")
params = params._replace(maximum_cached_engines=10)

converter = trt.TrtGraphConverterV2(
   input_saved_model_dir="./model",
   conversion_params=params
)

trt_func = converter.convert()

MAX_BATCH_SIZE=1
def my_input_fn():
   #x = tf.constant(X[0:MAX_BATCH_SIZE, :, :])
   #yield [x]
   Inp1 = np.random.normal(size=(MAX_BATCH_SIZE, 16, 8)).astype(np.float32) 
   yield (Inp1, )

print("Building converter...", end=" ")

converter.build(input_fn=my_input_fn)

OUTPUT_SAVED_MODEL_DIR='./model_trt_test'
converter.save(output_saved_model_dir=OUTPUT_SAVED_MODEL_DIR)

#loaded = tf.saved_model.load('./converted_model/tftrt_saved_model/')
loaded = tf.saved_model.load('./model_trt_test')

print("The signature keys are: ",list(loaded.signatures.keys())) 
infer = loaded.signatures["serving_default"]
print(infer.structured_outputs)
faulthandler.enable()

nBurnIn=100
for i in range(nSamples):
	if(i==nBurnIn):
		startTime=time.time()
	infer(tf.constant(X[i:i+1, :, :]))['dense']
endTime=time.time()
avgTime=(endTime-startTime)/(nSamples-nBurnIn)

print("Averge time:",avgTime)
	


