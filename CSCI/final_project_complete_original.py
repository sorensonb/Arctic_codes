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
#import PIL
from tensorflow.keras import layers
import time
import sys

##!#from final_project_lib import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# TensorFlow Model Functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def make_generator_model():
    model = tf.keras.Sequential()
    model.add(layers.Dense(7*7*256, use_bias=False, input_shape=(100,)))
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU())

    model.add(layers.Reshape((7, 7, 256)))
    assert model.output_shape == (None, 7, 7, 256)  # Note: None is the batch size

    model.add(layers.Conv2DTranspose(128, (5, 5), strides=(1, 1), padding='same', use_bias=False))
    assert model.output_shape == (None, 7, 7, 128)
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU())

    model.add(layers.Conv2DTranspose(64, (5, 5), strides=(2, 2), padding='same', use_bias=False))
    assert model.output_shape == (None, 14, 14, 64)
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU())

    model.add(layers.Conv2DTranspose(1, (5, 5), strides=(2, 2), padding='same', use_bias=False, activation='tanh'))
    assert model.output_shape == (None, 28, 28, 1)

    return model

def make_discriminator_model():
    model = tf.keras.Sequential()
    model.add(layers.Conv2D(64, (5, 5), strides=(2, 2), padding='same',
                                     input_shape=[28, 28, 1]))
    model.add(layers.LeakyReLU())
    model.add(layers.Dropout(0.3))

    model.add(layers.Conv2D(128, (5, 5), strides=(2, 2), padding='same'))
    model.add(layers.LeakyReLU())
    model.add(layers.Dropout(0.3))

    model.add(layers.Flatten())
    model.add(layers.Dense(1))

    return model

def discriminator_loss(real_output, fake_output):
    real_loss = cross_entropy(tf.ones_like(real_output), real_output)
    fake_loss = cross_entropy(tf.zeros_like(fake_output), fake_output)
    total_loss = real_loss + fake_loss
    return total_loss

def generator_loss(fake_output):
    return cross_entropy(tf.ones_like(fake_output), fake_output)

# Notice the use of `tf.function`
# This annotation causes the function to be "compiled".
@tf.function
#def train_step(images, epoch_gen_loss, epoch_disc_loss):
def train_step(images):
    noise = tf.random.normal([BATCH_SIZE, noise_dim])

    with tf.GradientTape() as gen_tape, tf.GradientTape() as disc_tape:
        
        #gen_loss_out  = tf.TensorArray(dtype = tf.float32, size = 1, dynamic_size = False)
        #disc_loss_out = tf.TensorArray(dtype = tf.float32, size = 1, dynamic_size = False)

        generated_images = generator(noise, training=True)

        real_output = discriminator(images, training=True)
        fake_output = discriminator(generated_images, training=True)

        gen_loss = generator_loss(fake_output)
        disc_loss = discriminator_loss(real_output, fake_output)

        ##!#epoch_gen_loss.update_state(gen_loss)
        ##!#epoch_disc_loss.update_state(disc_loss)

    gradients_of_generator = gen_tape.gradient(gen_loss, generator.trainable_variables)
    gradients_of_discriminator = disc_tape.gradient(disc_loss, discriminator.trainable_variables)

    generator_optimizer.apply_gradients(zip(gradients_of_generator, generator.trainable_variables))
    discriminator_optimizer.apply_gradients(zip(gradients_of_discriminator, discriminator.trainable_variables))

def train(dataset, epochs):
    ##!#loop_idx = np.arange(0,data2['AI'].shape[0], BATCH_SIZE)

    ##!#gen_losses = np.zeros((epochs, len(loop_idx) - 1))
    ##!#disc_losses = np.zeros((epochs, len(loop_idx) - 1))

    for epoch in range(epochs):
        start = time.time()

        epoch_gen_loss  = tf.keras.metrics.Mean()
        epoch_disc_loss = tf.keras.metrics.Mean()
    

        ##!##loop_idx = np.arange(0,data['AI'].shape[0], BATCH_SIZE)
        ##!#for ii in range(len(loop_idx) - 1):
        ##!#    local_arr = tf.convert_to_tensor(np.array(data2['AI'][loop_idx[ii]:loop_idx[ii+1],:,:])[...,np.newaxis])
        ##!#    print(loop_idx[ii], loop_idx[ii+1], local_arr.shape)
        ##!#    #gen_loss_out, disc_loss_out = train_step(np.array(data2['AI'][loop_idx[ii]:loop_idx[ii+1],:,:])[...,np.newaxis], ii)
        ##!#    train_step(local_arr, \
        ##!#        epoch_gen_loss, epoch_disc_loss)
        ##!#    print(float(epoch_gen_loss.result()), float(epoch_disc_loss.result()))

        ##!#    gen_losses[epoch, ii]  = float(epoch_gen_loss.result())
        ##!#    disc_losses[epoch, ii] = float(epoch_disc_loss.result())

        for image_batch in dataset:
          train_step(image_batch)
      
        # Produce images for the GIF as you go
        #display.clear_output(wait=True)
        plt.close('all')
        generate_and_save_images(generator,
                                 epoch + 1,
                                 seed)
      
        # Save the model every 15 epochs
        if (epoch + 1) % 15 == 0:
          checkpoint.save(file_prefix = checkpoint_prefix)
      
        print ('Time for epoch {} is {} sec'.format(epoch + 1, time.time()-start))
    
    # Generate after the final epoch
    #display.clear_output(wait=True)
    generate_and_save_images(generator,
                             epochs,
                             seed)

    print('done training')
    #return gen_losses, disc_losses


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Visualization Functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def generate_and_save_images(model, epoch, test_input):
  # Notice `training` is set to False.
  # This is so all layers run in inference mode (batchnorm).
  predictions = model(test_input, training=False)

  fig = plt.figure(figsize=(4, 4))

  for i in range(predictions.shape[0]):
      plt.subplot(4, 4, i+1)
      plt.imshow(predictions[i, :, :, 0] * 127.5 + 127.5, cmap='gray')
      plt.axis('off')

  plt.savefig('image_at_epoch_{:04d}_original.png'.format(epoch))
  #plt.show()

def save_gen_image(model, epoch, test_input):

  # Notice `training` is set to False.
  # This is so all layers run in inference mode (batchnorm).
  predictions = model(test_input, training=False)

  fig = plt.figure(figsize=(4, 4))
  ax1 = fig.add_subplot(1,1,1)
  ax1.imshow(predictions[0,:,:,0] * 127.5 + 127.5, cmap = 'gray')
  ax1.axis('off')

  outname = 'gen_image_epoch_{:04d}_original.png'.format(epoch)
  plt.savefig(outname)
  print("Saved image", outname)

def save_training_image(dataset, idx, mask = None):

  if(mask is not None):
    plot_data = (dataset[idx,:,:,0] * 127.5 + 127.5) * mask
  else:
    plot_data = (dataset[idx,:,:,0] * 127.5 + 127.5)

  fig = plt.figure(figsize=(4, 4))
  ax1 = fig.add_subplot(1,1,1)
  ax1.imshow(plot_data, cmap = 'gray')
  ax1.axis('off')

  if(mask is not None):
    outname = 'train_image_{:05d}_mask.png'.format(idx)
    plt.savefig(outname)
    print("Saved image", outname)
  else:
    outname = 'train_image_{:05d}.png'.format(idx)
    plt.savefig(outname)
    print("Saved image", outname)

def save_complete_image(mask_img, gen_img, mask, iter_num = 1, save = False):

    # Make completed image
    # --------------------
    #complete_img = mask_img * mask + (1. - mask) * gen_img
    complete_img = mask_img[:,:,0] * mask + (1. - mask) * gen_img[0,:,:,0]

    fig = plt.figure(figsize = (4,4))
    ax1 = fig.add_subplot(1,1,1)
    ax1.imshow(complete_img, cmap = 'gray')
    ax1.axis('off')

    if(save):
        outname = 'complete_image_test_v{:03d}.png'.format(iter_num)
        fig.savefig(outname)
        print("Saved image", outname)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Main code
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

print('reading training images')

(train_images, train_labels), (_, _) = tf.keras.datasets.mnist.load_data()

train_images = train_images.reshape(train_images.shape[0], 28, 28, 1).astype('float32')
train_images = (train_images - 127.5) / 127.5  # Normalize the images to [-1, 1]

BUFFER_SIZE = 60000
#BATCH_SIZE = 64
BATCH_SIZE = 256


# Batch and shuffle the data
train_dataset = tf.data.Dataset.from_tensor_slices(train_images).shuffle(BUFFER_SIZE).batch(BATCH_SIZE)

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


beta1=0.9
beta2=0.999
learn_rate = 1e-4
eps = 1e-8
generator_optimizer = tf.keras.optimizers.Adam(learn_rate)
discriminator_optimizer = tf.keras.optimizers.Adam(learn_rate)

print('setting checkpoints')
checkpoint_dir = './training_checkpoints'
checkpoint_prefix = os.path.join(checkpoint_dir, "ckpt")
checkpoint = tf.train.Checkpoint(generator_optimizer=generator_optimizer,
                                 discriminator_optimizer=discriminator_optimizer,
                                 generator=generator,
                                 discriminator=discriminator)

EPOCHS = 10
noise_dim = 100
num_examples_to_generate = 16

# You will reuse this seed overtime (so it's easier)
# to visualize progress in the animated GIF)
seed = tf.random.normal([num_examples_to_generate, noise_dim])


EPOCHS = 80
train(train_dataset, EPOCHS)
#gen_losses, disc_losses = train(data2, EPOCHS)
#plot_losses(gen_losses, disc_losses, EPOCHS, save = True)
sys.exit()
#gen_losses, disc_losses = train(data2, EPOCHS)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Image completion stuff
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

lam = tf.constant(0.1)


# Set up image mask
# -----------------
np_mask = np.ones((28,28))
np_mask[,:13:17] = 0.0
mask = tf.constant(np_mask, dtype = tf.float32)
#self.mask = tf.placeholder(tf.float32, [None] + self.image_shape, name='mask')

idx = 6789

# Save the original image and mask
# --------------------------------
save_training_image(train_images, idx, mask = None)
save_training_image(train_images, idx, mask = np_mask)

# Set up starting seeds 
# ---------------------
zhats = tf.Variable(np.random.uniform(-1, 1, size=(1, noise_dim)))

learn_rate = 2e-4
iter_nums = 1600
v_num = 0.0
losses = np.zeros(iter_nums)
for ii in range(iter_nums):
    
    with tf.GradientTape() as gen_tape:
        gen_tape.watch(zhats)
    
        # Set up loss calculations
        # ------------------------
        complete_loss = tf.reduce_sum(tf.compat.v1.layers.flatten(\
            tf.abs(tf.multiply(mask[np.newaxis,...,np.newaxis], \
            generator(zhats, training = False)) - tf.multiply(\
            mask[np.newaxis,...,np.newaxis], \
            tf.Variable(train_images[idx,:,::]))))) + \
            lam * generator_loss(\
            discriminator(generator(zhats, training = False), training = False))

        ##!#contextual_loss = tf.reduce_sum(\
        ##!#    tf.compat.v1.layers.flatten(tf.abs(tf.multiply(\
        ##!#    mask[np.newaxis,...,np.newaxis], generator(seed, training = False)) - \
        ##!#    tf.multiply(mask[np.newaxis,...,np.newaxis], \
        ##!#    tf.Variable(train_images[1234,:,:,:])))))
        ##!#
        ##!## Calculate generator loss with this seed
        ##!#generated_images = generator(zhats, training = False)
        ##!#fake_output = discriminator(generated_images, training = False)
        ##!#gen_loss = generator_loss(fake_output)
        ##!#perceptual_loss = gen_loss
        
        ##!## Calculate total loss
        ##!#complete_loss = contextual_loss + lam * perceptual_loss
    
    grad_complete_loss = gen_tape.gradient(complete_loss, zhats)
    #grad_complete_loss = tf.gradients(complete_loss, zhats)

    losses[ii] = complete_loss.numpy()

    ##!#option1 = True
    ##!#if(option1):
    v_prev = np.copy(v_num)
    v_num = beta1 * v_prev - learn_rate * grad_complete_loss[0]
    zhats = zhats.assign_add(\
        tf.expand_dims(-beta1 * v_prev + (1.0 + beta1) * v_num, axis = 0))
    zhats = tf.Variable(np.clip(zhats, -1, 1))
    ##!#else:
    ##!#    m_prev = np.copy(m_num)
    ##!#    v_prev = np.copy(v_num)
    ##!#    m = beta1 * m_prev + (1 - beta1) * grad_complete_loss[0]
    ##!#    v = beta2 * v_prev + (1 - beta2) * np.multiply(g[0], g[0])
    ##!#    m_hat = m / (1 - beta1 ** (i + 1))
    ##!#    v_hat = v / (1 - beta2 ** (i + 1))
    ##!#    zhats += - np.true_divide(learn_rate * m_hat, (np.sqrt(v_hat) + eps))
    ##!#    zhats = np.clip(zhats, -1, 1)

    if(ii % 100 == 0):
        zhat_gen_img = generator(zhats, training = False) 
        save_complete_image(train_images[idx], \
            zhat_gen_img, np_mask, iter_num = ii, save = True)
        

#self.contextual_loss = tf.reduce_sum(
#    tf.contrib.layers.flatten(
#        tf.abs(tf.mul(self.mask, self.G) - tf.mul(self.mask, self.images))), 1)
#self.perceptual_loss = self.g_loss
#self.complete_loss = self.contextual_loss + self.lam*self.perceptual_loss
#self.grad_complete_loss = tf.gradients(self.complete_loss, self.z)

# Define mask (extracted from the row anomaly mask)
# -------------------------------------------------

# Use projected gradient descent with minibatches and momemntum to project z
# to be in [-1, 1]
# ---------------------------------------------------------------------------
#for idx in xrange(0, batch_idxs):
#    batch_images = ...
#    batch_mask = np.resize(mask, [self.batch_size] + self.image_shape)
#    zhats = np.random.uniform(-1, 1, size=(self.batch_size, self.z_dim))
#
#    v = 0
#    for i in xrange(config.nIter):
#        fd = {
#            self.z: zhats,
#            self.mask: batch_mask,
#            self.images: batch_images,
#        }
#        run = [self.complete_loss, self.grad_complete_loss, self.G]
#        loss, g, G_imgs = self.sess.run(run, feed_dict=fd)
#
#        v_prev = np.copy(v)
#        v = config.momentum*v - config.lr*g[0]
#        zhats += -config.momentum * v_prev + (1+config.momentum)*v
#        zhats = np.clip(zhats, -1, 1)




