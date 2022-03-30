#!/usr/bin/env python

"""
  NAME:
    bsorenson_tensorflow_test.py
    
  PURPOSE:
    Test the "Boundless" image extension tutorial from Tensorflow
    locally using the Tensorflow installed on my Anaconda. 
    The code from each cell in the tutorial was copied here, with
    minor edits made by me for using my own images and saving 
    the results.

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

import tensorflow as tf
import tensorflow_hub as hub
from io import BytesIO
from PIL import Image as PilImage
import numpy as np
from matplotlib import pyplot as plt
from six.moves.urllib.request import urlopen
import sys

def read_image(filename):
    fd = None
    if(filename.startswith('http')):
      fd = urlopen(filename)
    else:
      fd = tf.io.gfile.GFile(filename, 'rb')

    pil_image = PilImage.open(fd)
    width, height = pil_image.size
    # crop to make the image square
    pil_image = pil_image.crop((0, 0, height, height))
    pil_image = pil_image.resize((257,257),PilImage.ANTIALIAS)
    image_unscaled = np.array(pil_image)
    image_np = np.expand_dims(
        image_unscaled.astype(np.float32) / 255., axis=0)
    return image_np

def visualize_output_comparison(img_original, img_masked, img_filled,\
        img_name, save = True):
    plt.figure(figsize=(8,5))
    plt.subplot(131)
    plt.imshow((np.squeeze(img_original)))
    plt.title("Original", fontsize=12)
    plt.axis('off')
    plt.subplot(132)
    plt.imshow((np.squeeze(img_masked)))
    plt.title("Masked", fontsize=12)
    plt.axis('off')
    plt.subplot(133)
    plt.imshow((np.squeeze(img_filled)))
    plt.title("Generated", fontsize=12)
    plt.axis('off')
    if(save):
        name_split = img_name.strip().split('/')[-1].split('.')
        outname = name_split[0] + '_reconstruct_threequart.png'
        plt.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else: 
        plt.show()

# Check arguments
# ---------------
if(len(sys.argv) != 2):
    print("SYNTAX: ./bsorenson_tensorflow_test.py img_name")
    print("")
    print("     Note: Be sure to activate the pytf environment")
    sys.exit()

# Get the input image from the command line
# -----------------------------------------
wikimedia = sys.argv[1]

##!#wikimedia = "https://upload.wikimedia.org/wikipedia/commons/thumb/3/31/Nusfjord_road%2C_2010_09.jpg/800px-Nusfjord_road%2C_2010_09.jpg"
##!## wikimedia = "https://upload.wikimedia.org/wikipedia/commons/thumb/4/47/Beech_forest_M%C3%A1tra_in_winter.jpg/640px-Beech_forest_M%C3%A1tra_in_winter.jpg"
##!## wikimedia = "https://upload.wikimedia.org/wikipedia/commons/thumb/b/b2/Marmolada_Sunset.jpg/640px-Marmolada_Sunset.jpg"
##!## wikimedia = "https://upload.wikimedia.org/wikipedia/commons/thumb/9/9d/Aegina_sunset.jpg/640px-Aegina_sunset.jpg"

input_img = read_image(wikimedia)

# Load the model
# --------------

#@title Model Selection { display-mode: "form" }
model_name = 'Boundless Three Quarters' # @param ['Boundless Half', 'Boundless Quarter', 'Boundless Three Quarters']
model_handle_map = {
    'Boundless Half' : 'https://tfhub.dev/google/boundless/half/1',
    'Boundless Quarter' : 'https://tfhub.dev/google/boundless/quarter/1', 
    'Boundless Three Quarters' : 'https://tfhub.dev/google/boundless/three_quarter/1'
}

model_handle = model_handle_map[model_name]
print("Loading model {} ({})".format(model_name, model_handle))
model = hub.load(model_handle)

result = model.signatures['default'](tf.constant(input_img))
generated_image =  result['default']
masked_image = result['masked_image']

save = False
visualize_output_comparison(input_img, masked_image, generated_image, \
    wikimedia, save = save)


