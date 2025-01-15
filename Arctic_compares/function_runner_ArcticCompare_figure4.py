#!/usr/bin/env python

"""

  SYNTAX FOR PAPER: ./function_runner_ArcticCompare_figure4.py 
      neuralnet_output/test_calc_out_noland105_201908100308.hdf5

"""

import Arctic_compare_lib
from Arctic_compare_lib import *
import random
#from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

# Compare the plume locations between OMI and MODIS data
# ------------------------------------------------------
plot_compare_NN_output_overlay_v2(sys.argv[1], auto_zoom = True, \
    label_xloc = 83.6, label_yloc = 76.0, save = True)
sys.exit()
