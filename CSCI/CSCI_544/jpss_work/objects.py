"""
  NAME:
    objects.py

  PURPOSE:
    Hold various objects to simplify the implementation of the DCGAN

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2022/05/05:
        Written

"""

dat_name = ''
beta1=0.9
learn_rate = 1e-4
epochs = 10
noise_dim = 100
num_examples_to_generate = 16
buffer_size = 60000   # size of training dataset
#BATCH_SIZE = 64
batch_size = 256
image_size = 28
complete_learn_rate = 2e-4

max_AI = None
min_AI = None
num_test_image = 1000

seed = None
generator = None
discriminator = None
train_images = None
train_lat = None
train_lon = None
train_time = None
test_images = None
test_lat = None
test_lon = None
test_time = None

generator_optimizer = None
discriminator_optimizer = None
checkpoint = None
checkpoint_prefix = None

np_mask = None
process_img = None
gen_img = None
complete_img = None
losses = None
min_loss = 20

errors = None
combined_errors = None
