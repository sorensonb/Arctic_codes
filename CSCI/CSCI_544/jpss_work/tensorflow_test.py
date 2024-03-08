#!/usr/bin/env python

"""
  NAME:
    final_project_complete.py
    
  PURPOSE:
    Run Python functions for training and running the DCGAN

    NOTE: BE SURE TO ACTIVATE THE pytf ENVIRONMENT BEFORE 
        RUNNING

  SYNTAX:
    python final_project_complete.py

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2022/02/27
      Written (modification of the Boundless colab tutorial
        from Tensorflow)
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import tensorflow as tf
from tensorflow.keras import layers
import time
import sys
import random

tf.debugging.set_log_device_placement(True)

print("GPU available:", tf.config.list_physical_devices('GPU'))
#print("GPU available:", tftest.is_gpu_available())
print("Built with CUDA:", tf.test.is_built_with_cuda())



a = tf.constant([[1.0,2.0,3.0], [4.0,5.0,6.0]])
b = tf.constant([[1.0,2.0], [3.0,4.0] ,[5.0, 6.0]])
c = tf.matmul(a,b)


print("SUCCESS")
