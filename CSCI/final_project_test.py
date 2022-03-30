#!/usr/bin/env python

"""
  NAME:
    
  PURPOSE:

    NOTE: the original example may be found at;
        https://www.tensorflow.org/hub/tutorials/boundless

    NOTE: BE SURE TO ACTIVATE THE pytf ENVIRONMENT BEFORE 
        RUNNING

  SYNTAX:
    ./bsorenson_tensorflow_test.py img_file

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2022/02/27
      Written (modification of the Boundless colab tutorial
        from Tensorflow)
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import PIL
from tensorflow.keras import layers
import time
import sys
from datetime import datetime
from netCDF4 import Dataset

from final_project_lib import *

#build_training_dataset('2006')
#sys.exit()
##!#
##!#date_str = ['200607060042',\
##!#            '200607060221',\
##!#            '200607060400',\
##!#            '200607060539',\
##!#            '200607060717',\
##!#            '200607060856',\
##!#            '200607061035',\
##!#            '200607061214',\
##!#            '200607061353',\
##!#            '200607061532',\
##!#            '200607061711',\
##!#            '200607061850',\
##!#            '200607062028',\
##!#            '200607062207',\
##!#            '200607062346']
##!#
##!##date_str = '200607060400'
##!##date_str = '201807260423'
##!##plot_resample_OMI(date_str, save = False)
##!#date_str = '20060728'
##!#dt_date_str = datetime.strptime(date_str,'%Y%m%d')
##!#filenames = glob.glob(dt_date_str.strftime('/home/bsorenson/data/OMI/H5_files/OMI-*%Ym%m%d*.he5'))
##!#for fname in filenames:
##!##for dstr in date_str:
##!#    dt_name = datetime.strptime(fname.strip().split('/')[-1][20:34],'%Ym%m%dt%H%M')
##!#    #print(dt_name.strftime('%Y%m%d%H%M'))
##!#    plot_resample_OMI_nomap(dt_name.strftime('%Y%m%d%H%M'), save = False)
##!#
##!#sys.exit()

print('reading training images')
data = Dataset('/home/bsorenson/CSCI/CSCI_543/final_project/train_dataset_2006.nc','r')

## Scale the AI data
## -----------------
data2 = contrast_stretch(data['AI'][:,:,:])
data2 = (data2 - 127.5) / 127.5

# Account for suspicious missing data?
# ------------------------------------

##!#(train_images, train_labels), (_, _) = tf.keras.datasets.mnist.load_data()
##!#
##!#train_images = train_images.reshape(train_images.shape[0], 28, 28, 1).astype('float32')
##!#train_images = (train_images - 127.5) / 127.5  # Normalize the images to [-1, 1]

BUFFER_SIZE = 60000
#BATCH_SIZE = 64
BATCH_SIZE = 256


##!## Batch and shuffle the data
##!#train_dataset = tf.data.Dataset.from_tensor_slices(train_images).shuffle(BUFFER_SIZE).batch(BATCH_SIZE)

print('making generator')
generator = make_generator_model()

noise = tf.random.normal([1, 100])
generated_image = generator(noise, training=False)

plt.imshow(generated_image[0, :, :, 0], cmap='gray')

discriminator = make_discriminator_model()
decision = discriminator(generated_image)
print (decision)

# This method returns a helper function to compute cross entropy loss
cross_entropy = tf.keras.losses.BinaryCrossentropy(from_logits=True)


generator_optimizer = tf.keras.optimizers.Adam(1e-4)
discriminator_optimizer = tf.keras.optimizers.Adam(1e-4)

##!#print('setting checkpoints')
##!#checkpoint_dir = './training_checkpoints'
##!#checkpoint_prefix = os.path.join(checkpoint_dir, "ckpt")
##!#checkpoint = tf.train.Checkpoint(generator_optimizer=generator_optimizer,
##!#                                 discriminator_optimizer=discriminator_optimizer,
##!#                                 generator=generator,
##!#                                 discriminator=discriminator)

EPOCHS = 10
noise_dim = 100
num_examples_to_generate = 16

# You will reuse this seed overtime (so it's easier)
# to visualize progress in the animated GIF)
seed = tf.random.normal([num_examples_to_generate, noise_dim])

# Notice the use of `tf.function`
# This annotation causes the function to be "compiled".
@tf.function
def train_step(images):
    noise = tf.random.normal([BATCH_SIZE, noise_dim])

    with tf.GradientTape() as gen_tape, tf.GradientTape() as disc_tape:
        generated_images = generator(noise, training=True)

        real_output = discriminator(images, training=True)
        fake_output = discriminator(generated_images, training=True)

        gen_loss = generator_loss(fake_output)
        disc_loss = discriminator_loss(real_output, fake_output)

    print(gen_loss, disc_loss)

    gradients_of_generator = gen_tape.gradient(gen_loss, generator.trainable_variables)
    gradients_of_discriminator = disc_tape.gradient(disc_loss, discriminator.trainable_variables)

    generator_optimizer.apply_gradients(zip(gradients_of_generator, generator.trainable_variables))
    discriminator_optimizer.apply_gradients(zip(gradients_of_discriminator, discriminator.trainable_variables))

def train(dataset, epochs):
    for epoch in range(epochs):
        start = time.time()
    
        loop_idx = np.arange(0,data2.shape[0], BATCH_SIZE)
        #loop_idx = np.arange(0,data['AI'].shape[0], BATCH_SIZE)
        for ii in range(len(loop_idx) - 1):
            print(loop_idx[ii], loop_idx[ii+1])
            train_step(np.array(data2[loop_idx[ii]:loop_idx[ii+1],:,:])[...,np.newaxis])
            #train_step(np.array(data['AI'][loop_idx[ii]:loop_idx[ii+1],:,:])[...,np.newaxis])
        ##!#for image_batch in dataset:
        ##!#  train_step(image_batch)
        #train_step(np.array(data['AI'][0:50,:,:])[...,np.newaxis])
        #train_step(np.array(data['AI'][50:100,:,:])[...,np.newaxis])
        #train_step(np.array(data['AI'][100:140,:,:])[...,np.newaxis])
      
        ##!## Produce images for the GIF as you go
        ##!#display.clear_output(wait=True)
        ##!#generate_and_save_images(generator,
        ##!#                         epoch + 1,
        ##!#                         seed)
      
        ##!## Save the model every 15 epochs
        ##!#if (epoch + 1) % 15 == 0:
        ##!#  checkpoint.save(file_prefix = checkpoint_prefix)
      
        print ('Time for epoch {} is {} sec'.format(epoch + 1, time.time()-start))
    
    ##!## Generate after the final epoch
    ##!#display.clear_output(wait=True)
    ##!#generate_and_save_images(generator,
    ##!#                         epochs,
    ##!#                         seed)

    print('done training')

#sys.exit()

data.close()
train(data2, EPOCHS)

