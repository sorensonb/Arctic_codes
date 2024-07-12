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

from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras.layers import *
from tensorflow.keras import Input
from tensorflow.keras import Model
import time
import sys
import random
import h5py
from sklearn.model_selection import train_test_split
import json


#tf.debugging.set_log_device_placement(True)

print("GPU available:", tf.config.list_physical_devices('GPU'))
#print("GPU available:", tftest.is_gpu_available())
print("Built with CUDA:", tf.test.is_built_with_cuda())

l_load_model = False
l_save_data  = False
if(len(sys.argv) < 2):
    print("SYNTAX: ./tensorflow_ai_test.py simulation_name [load] [save]")
    print("        [load]: if added, will load the model specified by simulation_name")
    print("        [save]: if added, will generate output files for specs given below")
    sys.exit()
elif(len(sys.argv) > 2):
    #if(sys.argv[2] == 'load'):
    if('load' in sys.argv[2:]):
        l_load_model = True
    if('save' in sys.argv[2:]):
        l_save_data = True

print(sys.argv, sys.argv[2:])
print(l_load_model, l_save_data)

sim_name = sys.argv[1]


file_list = [
    '200607240029',\
    '200607240208',\
    '200607240347',\
    '200607240526',\
    '200607242016',\
    '200607242155',\
    '200607242334',\
    '200607250112',\
    '200607250251',\
    '200607250430',\
    '200607252238',\
    '200607260017',\
    '200607260156',\
    '200607260335',\
    '200607260513',\
    '200607260652',\
    '200607260831',\
    '200607262003',\
    '200607262142',\
    '200607262321',\
    '200607270100',\
    '200607270239',\
    '200607270418',\
    '200607270557',\
    '200607270736',\
    '200607270914',\
    '200607272047',\
    '200607272226',\
    '200804221841',\
    '200804222159',\
    '201408110046',\
    '201408110404',\
    '201408111853',\
    '201408112032',\
    '201408112211',\
    '201408120308',\
    '201408121758',\
    '201408121937',\
    '201408122115',\
    '201408122254',\
    '201506271220',\
    '201506271359',\
    '201506271538',\
    '201506271717',\
    '201506271856',\
    '201507061353',\
    '201507061532',\
    '201507061711',\
    '201507061850',\
    '201507062028',\
    '201507062207',\
    '201507071615',\
    '201507071754',\
    #NOT HERE'201507071933',\
    '201507081837',\
    '201507082016',\
    '201507082155',\
    '201507090113',\
    '201507090748',\
    '201507090927',\
    '201507091245',\
    '201507091424',\
    '201507091603',\
    '201507100514',\
    '201507100653',\
    #NO COD'201708160014',\
    '201708161146',\
    '201708161325',\
    #NOT HERE'201708161504',\
    '201708161643',\
    '201708161821',\
    '201708162000',\
    #NO ALB'201708162139',\
    '201708171050',\
    '201708171229',\
    '201708171547',\
    '201708171726',\
    '201708172222',\
    '201708181133',\
    '201708181312',\
    '201708181451',\
    '201708181630',\
    '201708181809',\
    '201708181948',\
    '201708191038',\
    '201708191217',\
    '201708191355',\
    '201708191534',\
    #NO ALB'201708191713',\
    '201708191852',\
    '201807040005',\
    '201807040144',\
    '201807040322',\
    '201807041633',\
    '201807041812',\
    '201807041951',\
    '201807042130',\
    '201807042309',\
    '201807050048',\
    '201807050227',\
    '201807051538',\
    '201807051717',\
    '201807051856',\
    '201807052034',\
    '201807052213',\
    '201807052352',\
    '201807210047',\
    '201807211358',\
    '201807211537',\
    '201807211716',\
    '201807211855',\
    '201807212034',\
    '201808100200',\
    '201808140135',\
    '201808141804',\
    '201808141942',\
    '201808142121',\
    '201808142300',\
    '201808260158',\
    '201808260337',\
    '201808260655',\
    '201908100129',\
    '201908100308',\
    '201908101936',\
    '201908102115',\
    '201908102254',\
    '201908110033',\
    '201908110212',\
    '201908110351',\
    '201908110708',\
    '201908111523',\
    '201908111702',\
    '201908111841',\
    '201908112019',\
    '201908112158',\
    '201908112337',\
]

aer_file_list = ['/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_' + \
    date + '.hdf5' for date in file_list]

print('Num aerosol swaths', len(aer_file_list))

# input variables
# - OMI AI
# - OMI SZA
# - NSIDC ICE
# - MODIS COD
# - MODIS cloud top temp?

# output variable:
# - CERES SWF

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Prep the neural network
#
#
# noland3: 5 hidden layers
#           8,8,8,8,4 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               100 epochs
#               32 batch size
#               ReLU activation hidden
#               Linear activation out
#           Ending MAE: 5.09  
#           INCLUDES LAND DATA
#
# noland4: 5 hidden layers
#           8,8,8,8,4 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               64 batch size
#               ReLU activation hidden
#               Linear activation out
#           Ending MAE: 5.12
#           INCLUDES LAND DATA
#           Job ID = 92323
#
# noland5: 5 hidden layers
#           8,8,8,8,4 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.97
#           INCLUDES LAND DATA
#           Job ID = 92333
#
# noland6: 5 hidden layers
#           8,8,8,8,4 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               sigmoid activation hidden
#               Linear activation out
#           Ending MAE: 5.38
#           INCLUDES LAND DATA
#           Job ID = 92363
#
# noland7: 5 hidden layers
#           8,8,8,8,4 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.84
#           INCLUDES LAND DATA
#           Job ID = 92400
#
# noland8: 5 hidden layers
#           8,16,16,16,8 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.63
#           INCLUDES LAND DATA
#           Job ID = 92426
#
# noland9: 5 hidden layers
#           16,16,16,16,16 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.59
#           INCLUDES LAND DATA
#           Job ID = 92458
#
# noland10: 7 hidden layers
#           8,8,8,8,8,8,8 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.78
#           INCLUDES LAND DATA
#           Job ID = 92506
#
# noland11: 9 hidden layers
#           5,5,5,5,5,5,5,5,5 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.998
#           INCLUDES LAND DATA
#           Job ID = 93278
#
# noland12: 5 hidden layers
#           16,32,32,32,16 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.51
#           INCLUDES LAND DATA
#           Job ID = 93283
#
# noland13: 5 hidden layers
#           16,16,32,16,16 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: ????
#           INCLUDES LAND DATA
#           Job ID = 93315
#
# noland14: 7 hidden layers
#           8,16,16,32,16,16,8 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.49
#           INCLUDES LAND DATA
#           Job ID = 93409
#
# noland15: 8 hidden layers
#           8,16,16,32,32,16,16,8 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.54
#           INCLUDES LAND DATA
#           Job ID = 93436
#
# noland16: 10 hidden layers
#           4,8,8,8,8,8,8,8,8,4 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.77
#           INCLUDES LAND DATA
#           Job ID = 93449
#
# noland17: 12 hidden layers
#           4,8,8,8,8,8,8,8,8,8,8,4 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.82
#           INCLUDES LAND DATA
#           Job ID = 93469
#
# noland18: 14 hidden layers
#           4,8,8,8,8,8,8,8,8,8,8,8,8,4 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.79
#           INCLUDES LAND DATA
#           Job ID = 93593
#
# noland19: 9 hidden layers
#           8,12,16,24,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.47
#           INCLUDES LAND DATA
#           Job ID = 93607
#
# noland20: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.51
#           INCLUDES LAND DATA
#           Job ID = 93622
#
# noland21: 6 hidden layers
#           8,16,8,6,4,2 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               200 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.87
#           INCLUDES LAND DATA
#           Job ID = 94505
#
# noland22: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               250 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.44
#           INCLUDES LAND DATA
#           Job ID = 94516 ****** GOOD? ******
#
# noland23: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               250 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 4.05
#           INCLUDES LAND DATA, CH7
#           Job ID = 94522 ******* BETTER *******
#
# noland24: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               250 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 2.69
#           INCLUDES LAND DATA, CH7, CH1
#           Job ID = 94529 **** OUTPUT IS TRASH *******
#           Problem with new processed data?
#
# noland25: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 through 20180708 for
#               250 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: ?????
#           INCLUDES LAND DATA, CH7
#           Job ID = 94612  ****** OUTPUT IS TRASH *******
#           Problem with new processed data?
#
# noland26: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 and 20180705 for
#               250 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 2.69
#           INCLUDES LAND DATA, CH7, CH1
#           ORIGINAL DATA            
#           Job ID = 94699 
#
# noland27: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 thru 20180708 for
#               250 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.99
#           INCLUDES LAND DATA, CH7
#           FIXED DATA (CLDPRES data still don't work)           
#           Job ID = 94832 
#
# noland28: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 thru 20180708 for
#               350 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.65
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES data still don't work)           
#           Job ID = 94842 
#
# noland29: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 thru 20180708 for
#               350 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.69
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           Job ID = 95056 
#
# noland30: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 thru 20180708 for
#               350 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.64
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025, slight tweaks to data selection)
#           Job ID = 96551   ** Worse, for some reason **
#           * CHANGED TO COMP_DATA FILES FOR THIS ITERATION
#
# nolandXX: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 thru 20180708 for
#               350 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: ????
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (Retrying noland30 to make sure that the system is trained on all data)
#           Job ID = ?????   ** Worse, for some reason **
#
#   DISCOVERED THAT ALL THE WORKING SIMULATIONS ARE TRAINED ONLY ON
#   DATA FROM 20180704 AND 20180705. 
#
# noland32: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 thru 20180705 for
#               350 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.54
#           INCLUDES LAND DATA, CH7, VZA (from original_prep_data)
#           Retrying with a larger batch size to see if this does anything
#           Job ID = 96721
#   
# noland33: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 thru 20180705 for
#               500 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.51
#           INCLUDES LAND DATA, CH7, VZA (from original_prep_data)
#           Retrying with a larger batch size to see if this does anything
#           Job ID = 96994
#   
# noland34: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 thru 20180705 for
#               500 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.48
#           INCLUDES LAND DATA, CH7, VZA (from original_prep_data)
#           Retrying with a larger batch size to see if this does anything
#           Job ID = 97003
#   
# noland35: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 thru 20180705 for
#               500 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.51
#           INCLUDES LAND DATA, CH7, VZA (from original_prep_data & for saving)
#           Retrying with a larger batch size to see if this does anything
#           Job ID = 97168
#   
# noland36: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 thru 20180705 for
#               500 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 7.49
#           INCLUDES LAND DATA, CH7, VZA (from fixed comp_data & for saving)
#           After fixing the colocated data. Retrying without NAN CTP data
#           Job ID = 97333
#   
# noland37: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 thru 20180705 for
#               350 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: ????
#           INCLUDES LAND DATA, CH7, VZA (from fixed original_comp_data & for saving)
#           After fixing the colocated data. Retrying without NAN CTP data
#           Job ID = 97343
#   
# noland38: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 thru 20180708 for
#               350 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.57
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           Trying to replicate noland29. using 'original_prep_data'
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 99750  *** LOOKS GOOD LIKE noland38
#
# Reprocessed all the data using the old
#     write code from the end of the arctic comp code (no out function)
#
# noland39: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180704 thru 20180705 for
#               350 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.36
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland38: switching
#               data to 'comp_data' from 'original_prep_data'
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 99926  *** FROM INITIAL OUTPUT, LOOKS VERY GOOD *****
#
# noland40: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 thru 20180708 for
#               350 epochs
#               32 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.79
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland39:
#               Processing data from 20180701 to 20180708
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 106976  ***** FROM INITIAL OUTPUT, LOOKS VERY GOOD *****
#
# noland41: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 thru 20180708 for
#               350 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.73
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland40:
#               larger batch size (128 instead of 32)
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 107907  ***** FROM INITIAL OUTPUT, BETTER THAN noland40 *****
#           Most of the 10 calculated values are within 10-15 W/m2 of original
#
# noland42: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 thru 20180708 for
#               500 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.65
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland40:
#               more epochs (500 instead of 350)
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 107992  ***** BETTER THAN NOLAND41 BY EPOCH 350 *****
#
# noland43: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20180701 thru 20180708 for
#               700 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.70
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland40:
#               more epochs (700 instead of 500)
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 108257  
#
# #######
# Colocated data for 200607
# ######
#
# noland44: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20060721 thru 20060730 for
#               500 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.45
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland42:
#               using 200607 data instead of 20180701-08 data
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 108962  ***** BETTER THAN NOLAND41 BY EPOCH 350 *****
#
# noland45: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20060721 thru 20060730 for
#               500 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.64
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland44:
#               using 200607 data AND 20180701-08 data
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 109108  
#
# #######
# Colocated 200804 and 201408 data
# #######
#
# noland46: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20060721 thru 20060730 for
#               500 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.85
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland45:
#               using 200607, 200804, 201408, and 2018070*
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 110446  
#
# #######
# Colocated 2015 through 2019
# #######
#
# noland47: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on 20060721 thru 20060730 for
#               500 epochs
#               256 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 3.64
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland45:
#               using 2006 - 2015 coloc data AND batch size of 256
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 112580  
#
# #######
# Modified arctic comp colocation code to include CERES albedo
#
# ALSO, ran colocation code on ALL (essentially) prep files. The previous
# data, which don't contain albedo, are housed in coloc_data_before_ceres_alb
# #######

# noland48: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on ALL data
#               100 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 2.90
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland45:
#               Adding CERES albedo to the thing
#               using ALL colocation data
#               Batch size of 128 and 100 epochs
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 124908  
#           **** DOWN TO 3.19 BY EPOCH 5 **** IS THIS TOO GOOD?
#
# noland49: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on ALL data
#               100 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 2.91
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland48:
#               All the aerosol-containing swaths and a testing aerosol-free 
#               swathare removed 
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 27399  
#           **** Down to 3.22 by Epoch 5 ****
#
# noland50: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on ALL data
#               100 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 2.89
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland49:
#               Dropped testing percentage to 0.1
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 27899  
#
# noland51: 15 hidden layers
#           8,12,16,24,32,64,64,64,64,64,32,24,16,12,8 nodes in each layer
#           Trained on ALL data
#               100 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 2.89
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland50:
#               Added 4 more hidden layers of 64 nodes each
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 27953  
#           Bad results, high bias in center
#
# noland52: 13 hidden layers
#           8,12,16,24,32,64,128,64,32,24,16,12,8 nodes in each layer
#           Trained on ALL data
#               100 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 2.86
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#           ONLY DIFFERENCE BETWEEN THIS AND noland50:
#               Added 2 more hidden layers of 128 and 64 nodes
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 28053  
#
# noland53: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on ALL data
#               100 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: 2.89
#           INCLUDES LAND DATA, CH7, VZA
#           FIXED DATA (CLDPRES NAN set to 1025)
#
#           SHOULD BE NO DIFFERENCES BETWEEN THIS AND noland50
#               Am just saving the scaling values to a json file, which
#               will be loaded in when the user wants to calculate output
#               for all the aerosol files.
#               
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 43157
#  
# noland54: 11 hidden layers
#           8,12,16,24,32,64,32,24,16,12,8 nodes in each layer
#           Trained on ALL data
#               100 epochs
#               128 batch size
#               Leaky ReLU activation hidden
#               Linear activation out
#           Ending MAE: ????
#           INCLUDES LAND DATA, CH7, VZA
#
#           Only differences between this and noland53
#               Removed the code to set 0-value and NaN-value CTP to 1025.
#               Looked at the data and discovered that there are no 0-value
#               or NaN-value CTPs in the input data.              
#
#           CLDPRES NANs are accounted for in 'select data'
#           Job ID = 43699
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Inputs
# - OMI SZA
# - OMI VZA
# - NSIDC ICE
# - MODIS COD
# - MODIS CH7
# - MODIS CLD TOP TEMP/PRESS
# - CERES ALB
#### - MODIS CH1
# Future inputs
# - OMI ALB
# - OMI AZM  NO! BAD DATA
# - OMI GPQF

batch_size = 128
epochs = 100

#minlat = 70.
min_ai = -2.0
max_ai = 1.0
minlat = 65.    # 2024/01/10: changed to 65.
min_swf = 0.
max_swf = 3000.
max_cod = 70.
min_ice = 0.
max_ice = 500.

if(l_load_model):

    # Load the model
    # --------------
    model = tf.keras.models.load_model('test_model_' + sim_name + '.keras')
    model.summary()

    # Load the min_max_dict values from the run
    # -----------------------------------------
    min_max_dict_name = 'min_max_dict_' + sim_name +  '.json'
    with open(min_max_dict_name, 'r') as fin:
        min_max_dict = json.load(fin)

else:

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Prep the data
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    #data_dir = '/home/blake.sorenson/OMI/arctic_comp/original_prep_data/'
    data_dir = '/home/blake.sorenson/OMI/arctic_comp/comp_data/'
    
    # Scale all input variables to 0 - 1
    
    #data_path = '/home/blake.sorenson/OMI/arctic_comp/original_prep_data/'
    data_path = data_dir
    ##!#fdates = [
    ##!#        '201807040005',\
    ##!#        '201807040144',\
    ##!#        '201807040322',\
    ##!#        '201807041633',\
    ##!#        '201807041812',\
    ##!#        '201807041951',\
    ##!#        '201807042130',\
    ##!#        '201807042309',\
    ##!#        '201807050048',\
    ##!#        '201807050227',\
    ##!#        '201807051538',\
    ##!#        '201807051717',\
    ##!#        '201807051856',\
    ##!#        '201807052034',\
    ##!#        '201807052213',\
    ##!#        '201807052352'
    ##!#    ]
    ##!#
    ##!#files = [data_path + 'colocated_subset_' + fdd + '.hdf5' for fdd in fdates]
    
    files = glob('/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated*.hdf5')
    #files = glob('/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_200607*.hdf5')
    
    # NEW FOR NOLAND49: Remove the aerosol swaths from the training dataset
    # NEW FOR NOLAND49: Remove a sample clear-sky testing swath
    # ---------------------------------------------------------------------
    print("Before removing, length of files = ", len(files))
    for fname in aer_file_list:
        if(fname in files):
            files.remove(fname)
            print('REMOVING', fname)
    
    clear_swaths = ['201807082244']
    files.remove('/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_201807082244.hdf5')
    print('REMOVING 201807082244')
    
    print("After removing, length of files = ", len(files))
    
    #files = ['/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_201807052213.hdf5']
    
    # Figure out the total size to insert the data
    # ---------------------------------------------
    #minlat = 70.
    min_ai = -2.0
    max_ai = 1.0
    minlat = 65.    # 2024/01/10: changed to 65.
    min_swf = 0.
    max_swf = 3000.
    max_cod = 70.
    min_ice = 0.
    max_ice = 500.
        
    
    # NOTE: added the "==254" section so that it temporarily removes land data.
    #       Want to test the system over only ocean and/or ice to start with
    def select_data(data, min_ai, max_ai, minlat, min_swf, max_swf, max_cod):
        local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
        return np.ma.masked_where( (local_data < min_ai) | \
                                   (local_data > max_ai) | \
                                   (data['omi_lat'][:,:] < minlat) | \
                                   (data['modis_ch1'][:,:] > 1.0)  | \
                                   (data['modis_ch7'][:,:] > 1.0)  | \
                                   ( (data['ceres_swf'][:,:] < min_swf) | \
                                   (data['ceres_swf'][:,:] > max_swf) ) | \
                                   (np.isnan(data['ceres_swf'][:,:]) == True) | \
                                   (data['ceres_alb'][:,:] == -999) | \
                                   (data['modis_cod'][:,:] == -999) | \
                                   (data['modis_cod'][:,:] > max_cod) | \
                                   (np.isnan(data['modis_cod'][:,:]) == True) | \
                                   (data['modis_cld_top_pres'][:,:] < 0) | \
                                   (data['modis_cld_top_pres'][:,:] > 1025) | \
                                   (np.isnan(data['modis_cld_top_pres'][:,:]) == True) | \
                                   (data['nsidc_ice'][:,:] == -999) | \
                                   (data['nsidc_ice'][:,:] == 251) | \
                                   #(data['nsidc_ice'][:,:] == 254) | \
                                   (data['nsidc_ice'][:,:] == 253), \
                                   local_data\
        )
        
    
    total_size = 0
    for ff in files:
        data = h5py.File(ff,'r')
    
        # See if CERES albedo is in the file
        # ----------------------------------
        if('ceres_alb' not in data.keys()):
            print("ERROR: ",ff," DOES NOT CONTAIN CERES ALBEDO")
        else:
            local_data = select_data(data, min_ai, max_ai, minlat, min_swf, \
                max_swf, max_cod)
    
            #print(ff, np.max(local_data))
            ##!#local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
            ##!#local_data = np.ma.masked_where( (abs(local_data) > 12) | \
            ##!#                                 (data['omi_lat'][:,:] < minlat) | \
            ##!#                                 ( (data['ceres_swf'][:,:] < 0) | \
            ##!#                                 (data['ceres_swf'][:,:] > 3000) ) | \
            ##!#                                 (data['modis_cod'][:,:] == -999) | \
            ##!#                                 (data['modis_cod'][:,:] > 70) | \
            ##!#                                 (np.isnan(data['modis_cod'][:,:]) == True) | \
            ##!#                                 (data['nsidc_ice'][:,:] == -999) | \
            ##!#                                 (data['nsidc_ice'][:,:] == 251) | \
            ##!#                                 (data['nsidc_ice'][:,:] == 253), \
            ##!#                                 local_data\
            ##!#)
    
            #local_data = np.ma.masked_where(np.isnan(local_data))
    
            local_size = local_data.compressed().shape[0]
            total_size += local_size
            print(ff, local_size, total_size)
            ##!#if(local_size > 0):
            ##!#    print('Overall AI max:', \
            ##!#        np.max(np.ma.masked_invalid(data['omi_uvai_pert'])), \
            ##!#        np.max(local_data))
            ##!#else:
            ##!#    print("NO DATA")
    
        data.close()
    
    combined_data = {}
    combined_data['omi_uvai_pert'] = np.full(total_size, np.nan)
    #combined_data['omi_uvai_raw']  = np.full(total_size, np.nan)
    #combined_data['modis_cld']     = np.full(total_size, np.nan)
    combined_data['modis_cld_top_pres']     = np.full(total_size, np.nan)
    combined_data['modis_cod']     = np.full(total_size, np.nan)
    combined_data['modis_ch1']     = np.full(total_size, np.nan)
    combined_data['modis_ch7']     = np.full(total_size, np.nan)
    combined_data['ceres_swf']     = np.full(total_size, np.nan)
    combined_data['ceres_alb']     = np.full(total_size, np.nan)
    combined_data['omi_sza']       = np.full(total_size, np.nan)
    combined_data['omi_vza']       = np.full(total_size, np.nan)
    #combined_data['omi_azm']       = np.full(total_size, np.nan)
    combined_data['omi_lat']       = np.full(total_size, np.nan)
    combined_data['nsidc_ice']     = np.full(total_size, np.nan)
    
    print("Loading data")
    
    # Loop back over the files and insert the data into the structure
    # ---------------------------------------------------------------
    total_size = 0
    beg_idx = 0
    end_idx = 0
    for ff in files:
    
        data = h5py.File(ff,'r')
    
    
        if('ceres_alb' in data.keys()):
            local_data = select_data(data, min_ai, max_ai, minlat, min_swf, \
                max_swf, max_cod)
    
    
            # NOTE: Changed the omi variable here from "pert" to "raw" on 20230623.
            #       This move should allow for coloc data to be read after 2020
            ##!#local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
    
            ##!#local_data = np.ma.masked_where( (abs(local_data) > 12) | \
            ##!#                                 (data['omi_lat'][:,:] < minlat) | \
            ##!#                                 ( (data['ceres_swf'][:,:] < 0) | \
            ##!#                                 (data['ceres_swf'][:,:] > 3000) ) | \
            ##!#                                 (data['modis_cod'][:,:] == -999) | \
            ##!#                                 (data['modis_cod'][:,:] > 70) | \
            ##!#                                 (np.isnan(data['modis_cod'][:,:]) == True) | \
            ##!#                                 (data['nsidc_ice'][:,:] == -999) | \
            ##!#                                 (data['nsidc_ice'][:,:] == 251) | \
            ##!#                                 (data['nsidc_ice'][:,:] == 253), 
            ##!#                                 local_data\
            ##!#)
    
            #local_data = np.ma.masked_invalid(local_data)
            #local_data = np.ma.masked_where(np.isnan(local_data))
            local_size = local_data.compressed().shape[0]
    
            beg_idx = end_idx
            end_idx = beg_idx + local_size
    
            for tkey in combined_data.keys():
                combined_data[tkey][beg_idx:end_idx] = \
                    data[tkey][~local_data.mask]
    
            print(beg_idx, end_idx, local_size)
            #print(local_size, 'AT THIS STEP', np.nanmax(combined_data['omi_uvai_pert']))
            total_size += local_size
    
        data.close()
    
    print(combined_data['omi_uvai_pert'].shape)
    
    # Remove any values that still have missing AI for some reason
    good_idxs = np.where(np.isnan(combined_data['omi_uvai_pert']) == False)
    for tkey in combined_data.keys():
        combined_data[tkey] = combined_data[tkey][good_idxs]
 
 
    # 2024-07-12: Testing how often the '0'-value CTPs show up
    # --------------------------------------------------------
    print(len(good_idxs[0]))
    zero_idxs = np.where(combined_data['modis_cld_top_pres'] == 0.)
    print('NUMBER OF ZERO CTP', len(zero_idxs[0]))  
    print('Unique CTPs', np.unique(combined_data['modis_cld_top_pres']))  

    # 2024-07-12: Testing how often the NaN-value CTPs show up
    # --------------------------------------------------------
    nan_idxs = np.where(np.isnan(combined_data['modis_cld_top_pres']) == True)
    print('NUMBER OF NAN CTP', len(nan_idxs[0]))  

    # Figure out what values of COD are found when the CTP is 1025
    # ------------------------------------------------------------
    """
    check_idxs = np.where(combined_data['modis_cld_top_pres'] == 1025.)
    check_cod_vals = combined_data['modis_cod'][check_idxs]
    print('COD value ranges for 1025 CTP', np.min(check_cod_vals), np.max(check_cod_vals))  
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist(check_cod_vals, bins = 100)
    ax.set_xlabel('COD')
    ax.set_title('MODIS COD when MODIS CTP == 1025 mb')
    fig.tight_layout()
    fig.savefig('modis_cod_ctp1025_check.png', dpi = 200)
    print("Saved image")
    """
    #sys.exit()
 
    """ 
    # Fix the cloud top pressure so that '0' values (which I assume mean 'clear-sky'?)
    # are at 1025 rather than 0. This will hopefully give continuity from low clouds
    # to no clouds
    combined_data['modis_cld_top_pres'] = \
        np.where(combined_data['modis_cld_top_pres'] == 0., 1025., \
        combined_data['modis_cld_top_pres'])
    # TESTING: Make all NAN values to be "clear-sky"
    combined_data['modis_cld_top_pres'] = \
        np.where(np.isnan(combined_data['modis_cld_top_pres'][:]) == True, 1025., \
        combined_data['modis_cld_top_pres'][:])
    """ 
    
    print(np.min(combined_data['modis_cld_top_pres']), np.max(combined_data['modis_cld_top_pres']))
    
    print("WHERE AI IS NAN", np.where(np.isnan(combined_data['omi_uvai_pert'])))
    
    print('AI ',np.min(combined_data['omi_uvai_pert']), np.max(combined_data['omi_uvai_pert']))
    print('SWF',np.min(combined_data['ceres_swf']), np.max(combined_data['ceres_swf']))
    print('COD',np.min(combined_data['modis_cod']), np.max(combined_data['modis_cod']))
    print('CH1',np.min(combined_data['modis_ch1']), np.max(combined_data['modis_ch1']))
    print('CH7',np.min(combined_data['modis_ch7']), np.max(combined_data['modis_ch7']))
    print('CTP',np.min(combined_data['modis_cld_top_pres']), np.max(combined_data['modis_cld_top_pres']))
    print('SZA',np.min(combined_data['omi_sza']), np.max(combined_data['omi_sza']))
    print('VZA',np.min(combined_data['omi_vza']), np.max(combined_data['omi_vza']))
    print('ALB',np.min(combined_data['ceres_alb']), np.max(combined_data['ceres_alb']))
    print('ICE',np.min(combined_data['nsidc_ice']), np.max(combined_data['nsidc_ice']))
    
    combined_data['nsidc_ice'][:] = \
        np.where(combined_data['nsidc_ice'][:] == 254., 101., combined_data['nsidc_ice'][:])
    
    min_max_dict = {}
    
    key_variables = ['omi_uvai_pert', 'omi_sza', 'omi_vza', 'modis_cod', 'modis_cld_top_pres', 'nsidc_ice', \
        'modis_ch1', 'modis_ch7','ceres_alb','ceres_swf']
    
    for key in key_variables:
        min_max_dict[key] = {}
        min_max_dict[key]['min'] = np.min(combined_data[key])
        min_max_dict[key]['max'] = np.max(combined_data[key])
    
        drange = min_max_dict[key]['max'] - min_max_dict[key]['min']
        combined_data[key] = ((combined_data[key] - min_max_dict[key]['min']) / drange) * 100.
        #combined_data[key] = ((combined_data[key] - min_max_dict[key]['min']) / drange) * 100.
   
    # Save the min_max_dict values to a json file for later loading
    # -------------------------------------------------------------
    min_max_dict_name = 'min_max_dict_' + sim_name + '.json'

    with open(min_max_dict_name,'w') as fout:
        json.dump(min_max_dict, fout, indent = 4)

 
    print('AI ',np.min(combined_data['omi_uvai_pert']), np.max(combined_data['omi_uvai_pert']))
    print('SWF',np.min(combined_data['ceres_swf']), np.max(combined_data['ceres_swf']))
    print('COD',np.min(combined_data['modis_cod']), np.max(combined_data['modis_cod']))
    print('CH1',np.min(combined_data['modis_ch1']), np.max(combined_data['modis_ch1']))
    print('CH7',np.min(combined_data['modis_ch7']), np.max(combined_data['modis_ch7']))
    print('CTP',np.min(combined_data['modis_cld_top_pres']), np.max(combined_data['modis_cld_top_pres']))
    print('SZA',np.min(combined_data['omi_sza']), np.max(combined_data['omi_sza']))
    print('VZA',np.min(combined_data['omi_vza']), np.max(combined_data['omi_vza']))
    print('ALB',np.min(combined_data['ceres_alb']), np.max(combined_data['ceres_alb']))
    print('ICE',np.min(combined_data['nsidc_ice']), np.max(combined_data['nsidc_ice']))
    
        

    pcnt_test = 0.10
    num_test = int(combined_data['omi_uvai_pert'].shape[0] * pcnt_test)
    num_train = combined_data['omi_uvai_pert'].shape[0] - num_test
    ranges = np.arange(0, combined_data['omi_uvai_pert'].shape[0])
    
    train_idxs, test_idxs = train_test_split(ranges, test_size = num_test)
    
    print(train_idxs.shape, test_idxs.shape)
    
    # Input format: OMI SZA, NSIDC ICE, MODIS COD
    x_train = np.array([combined_data['omi_sza'][train_idxs], \
                        combined_data['omi_vza'][train_idxs], \
                        combined_data['nsidc_ice'][train_idxs], \
                        combined_data['modis_cod'][train_idxs], \
                        combined_data['modis_ch7'][train_idxs], \
                        combined_data['modis_cld_top_pres'][train_idxs],\
                        combined_data['ceres_alb'][train_idxs]\
                      ]).T
    y_train = combined_data['ceres_swf'][train_idxs]
    
    x_val   = np.array([combined_data['omi_sza'][test_idxs], \
                        combined_data['omi_vza'][test_idxs], \
                        combined_data['nsidc_ice'][test_idxs], \
                        combined_data['modis_cod'][test_idxs], \
                        combined_data['modis_ch7'][test_idxs], \
                        combined_data['modis_cld_top_pres'][test_idxs],\
                        combined_data['ceres_alb'][test_idxs]\
                      ]).T
    y_val   = combined_data['ceres_swf'][test_idxs]
    
    print(np.max(y_train), np.max(y_val))

    print("EPOCHS:", epochs)
    print("BATCH SIZE:", batch_size)
    
    input_data = [combined_data['omi_sza'][train_idxs], \
                        combined_data['omi_vza'][train_idxs], \
                        combined_data['nsidc_ice'][train_idxs], \
                        combined_data['modis_cod'][train_idxs], \
                        #combined_data['modis_ch1'][train_idxs], \
                        combined_data['modis_ch7'][train_idxs], \
                        combined_data['modis_cld_top_pres'][train_idxs],\
                        combined_data['ceres_alb'][train_idxs]]
    
    print('INPUT DATA SHAPE', np.array(input_data).shape)


    # BUILD THE MODEL HERE

    x1 = Input(shape = (1,))
    x2 = Input(shape = (1,))
    x3 = Input(shape = (1,))
    x4 = Input(shape = (1,))
    x5 = Input(shape = (1,))
    x6 = Input(shape = (1,))
    x7 = Input(shape = (1,))

    # ADD DROPOUT?

    #input_layer = concatenate([x1, x2, x3, x4, x5], name = 'input')
    #input_layer = concatenate([x1, x2, x3, x4, x5, x6], name = 'input')
    input_layer = concatenate([x1, x2, x3, x4, x5, x6, x7], name = 'input')
    hidden_layer1  = Dense(units = 8, activation = 'leaky_relu', name = 'hidden1')(input_layer)
    hidden_layer2  = Dense(units = 12, activation = 'leaky_relu', name = 'hidden2')(hidden_layer1)
    hidden_layer3  = Dense(units = 16, activation = 'leaky_relu', name = 'hidden3')(hidden_layer2)
    hidden_layer4  = Dense(units = 24, activation = 'leaky_relu', name = 'hidden4')(hidden_layer3)
    hidden_layer5  = Dense(units = 32, activation = 'leaky_relu', name = 'hidden5')(hidden_layer4)
    hidden_layer6  = Dense(units = 64, activation = 'leaky_relu', name = 'hidden6')(hidden_layer5)
    hidden_layer7  = Dense(units = 32, activation = 'leaky_relu', name = 'hidden7')(hidden_layer6)
    hidden_layer8  = Dense(units = 24, activation = 'leaky_relu', name = 'hidden8')(hidden_layer7)
    hidden_layer9  = Dense(units = 16, activation = 'leaky_relu', name = 'hidden9')(hidden_layer8)
    hidden_layer10 = Dense(units = 12, activation = 'leaky_relu', name = 'hidden10')(hidden_layer9)
    hidden_layer11 = Dense(units = 8, activation = 'leaky_relu', name = 'hidden11')(hidden_layer10)
    prediction = Dense(1, activation = 'linear')(hidden_layer11)

    #model = Model(inputs = [x1, x2, x3, x4, x5], outputs = prediction)
    model = Model(inputs = [x1, x2, x3, x4, x5, x6, x7], outputs = prediction)
    #model = Model(inputs = [x1, x2, x3], outputs = prediction)





    model.compile(loss = 'mean_squared_error', optimizer = 'adam', metrics = ['mae'])
    model.summary()
    losses = model.fit(input_data, \
                        combined_data['ceres_swf'][train_idxs], \
                        epochs=epochs, batch_size = batch_size, verbose = 2)
    #model.fit([inp1, inp2, inp3], sumd, epochs=300, batch_size = 32, verbose = 2)

    #model.save('test_model_allsfc.keras')
    model.save('test_model_' + sim_name + '.keras')


    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Test the neural network
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #calc = model.predict(x_train[0:10,:])
    ##!#test_out = np.array([model.predict([np.expand_dims(np.array(combined_data['omi_sza'][test_idxs][xx]), 0), \
    ##!#                          np.expand_dims(np.array(combined_data['omi_vza'][test_idxs][xx]), 0), \
    ##!#                          np.expand_dims(np.array(combined_data['nsidc_ice'][test_idxs][xx]), 0), \
    ##!#                          np.expand_dims(np.array(combined_data['modis_cod'][test_idxs][xx]), 0), \
    ##!#                          #np.expand_dims(np.array(combined_data['modis_ch1'][test_idxs][xx]), 0), \
    ##!#                          np.expand_dims(np.array(combined_data['modis_ch7'][test_idxs][xx]), 0), \
    ##!#                          np.expand_dims(np.array(combined_data['modis_cld_top_pres'][test_idxs][xx]), 0)]) \
    ##!#                    for xx in range(10)]).squeeze()
    
    test_out = np.array([model.predict([combined_data['omi_sza'][test_idxs][:10], \
                                        combined_data['omi_vza'][test_idxs][:10], \
                                        combined_data['nsidc_ice'][test_idxs][:10], \
                                        combined_data['modis_cod'][test_idxs][:10], \
                                        #combined_data['modis_ch1'][test_idxs][xx]), 0), \
                                        combined_data['modis_ch7'][test_idxs][:10], \
                                        combined_data['modis_cld_top_pres'][test_idxs][:10], \
                                        combined_data['ceres_alb'][test_idxs][:10]])]).squeeze()
    
    
    print(  ((combined_data['ceres_swf'][test_idxs][0:10] / 100) * drange) + min_max_dict['ceres_swf']['min'])
    print(  ((test_out / 100) * drange) + min_max_dict['ceres_swf']['min'])



    """
    input_shape = 4
    model = tf.keras.Sequential([
            tf.keras.layers.Dense(units = 8, activation = 'relu', 
                                  input_shape = (input_shape,)),
            tf.keras.layers.Dense(units = 16, activation = 'relu'), 
            tf.keras.layers.Dense(units = 8, activation = 'relu'), 
            tf.keras.layers.Dense(units = 1, activation = 'linear')
            #tf.keras.layers.Dense(units = 32, activation = 'relu'), 
    
        ])
    
    model.compile(optimizer = 'adam',\
                  loss = 'mae')
                  #loss = 'mae')
    
    model.summary()
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Train the neural network
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    losses = model.fit(x_train, y_train, \
                       validation_data = (x_val, y_val), \
                       batch_size = 128, \
                       epochs = 100,
            )
    """

if(l_save_data):

    ##!#print(x_train[0:10,:])
    ##!#print(    ((y_val[0:10].T / 100) * drange) + min_max_dict['ceres_swf']['min'])
    ##!#print(    ((calc / 100) * drange) + min_max_dict['ceres_swf']['min'])
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Test saving the neural network
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    # HERE: Read in a single one of the OMI colocation files, scale the input data
    #       based on the scaling from above, and then use the model to calculate
    #       the predicted SWF values at each gridpoint in the OMI file. Want
    #       to then plot the results and see if there is any spatial information
    #       that makes it through the model
    #       
    #       ALSO: Must remember that this is for ice/ocean values ONLY. Land data
    #       are not included right now.
    #test_time = '201807051856'
    #test_time = '201807052034'
    #test_time = '201807052213'
    #test_time = '201807040819'
    #test_time = '201807081926'
    #test_time = '201807082244'
    #test_time = '201807082105'
    #test_time = '201807032226'
    #test_time = '201807040144'
    #test_time = '200607240029'
    #test_time = '200607270418'
    #infile = '/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_' + \
    #    test_time + '.hdf5'
    #infile = '/home/blake.sorenson/OMI/arctic_comp/original_prep_data/colocated_subset_' + \
    #infile = data_dir + 'colocated_subset_' + test_time + '.hdf5'

    #test_time = '20060724'
    #test_time = '20180705'
    test_time = '20170816'

    if(len(test_time) == 12): 
        files = glob('/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_' + test_time + '.hdf5')
    elif(len(test_time) < 12):
        files = glob('/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_' + test_time + '*.hdf5')
    else:
        print("ERROR: Invalid time value")
        sys.exit()  

    aer_file_list = ['/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_201807052213.hdf5', \
                     '/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_201807082244.hdf5']
    #aer_file_list = ['/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_201807082244.hdf5']
 
    #for infile in files:
    for infile in aer_file_list:

        print("Completing:", infile)

        local_time = infile.strip().split('/')[-1].split('_')[-1].split('.')[0]
 
        # Set up input data
        data = h5py.File(infile, 'r')
        
        # Set up output data
        outname = 'test_calc_out_' + sim_name + '_' + local_time + '.hdf5'
        dset = h5py.File(outname, 'w')
        
        calc_swf = np.full(data['ceres_swf'].shape, np.nan)

        input_shape = calc_swf.shape
   
        # Prep the input data here so the entire array can be passed to the model
        # -----------------------------------------------------------------------
        data_mask = np.ma.masked_where( \
                    (data['omi_lat'][:,:] < minlat) | \
                    (data['modis_cod'][:,:] == -999) | \
                    (data['modis_cod'][:,:] > max_cod) | \
                    (data['modis_ch1'][:,:] > 1.0) | \
                    (data['modis_ch7'][:,:] > 1.0) | \
                    #(data['modis_ch1'][:,:] <= 1.0) | \
                    #(data['modis_ch7'][:,:] <= 1.0) | \
                    (np.isnan(data['modis_cod'][:,:]) == True) | \
                    #(local_cldpres != -999) & \
                    #(local_cldpres <= 1025.) & \
                    #(np.isnan(data['modis_cod'][ii,jj]) == False) & \
                    (data['ceres_alb'][:,:] == -999) | \
                    (data['modis_cld_top_pres'][:,:] == -999) | \
                    (data['modis_cld_top_pres'][:,:] > 1025) | \
                    (np.isnan(data['modis_cld_top_pres'][:,:]) == True) | \
                    (data['nsidc_ice'][:,:] == -999) | \
                    (data['nsidc_ice'][:,:] == 251) | \
                    (data['nsidc_ice'][:,:] == 253)  , data['omi_lat'])

        """
        local_cldpres = np.where(data['modis_cld_top_pres'][:,:] == 0., 1025., \
                                    data['modis_cld_top_pres'][:,:])
        local_cldpres = np.where(np.isnan(local_cldpres) == True, 1025., \
                                    local_cldpres)
        """
        local_ice     = np.where(data['nsidc_ice'][:,:] == 254.,101., \
                                    data['nsidc_ice'])
 
        new_sza = ( ( data['omi_sza'][:,:] - min_max_dict['omi_sza']['min']) / \
                    ( min_max_dict['omi_sza']['max'] - min_max_dict['omi_sza']['min'] )) * 100.
        new_vza = ( ( data['omi_vza'][:,:] - min_max_dict['omi_vza']['min']) / \
                    ( min_max_dict['omi_vza']['max'] - min_max_dict['omi_vza']['min'] )) * 100.
        new_ice = ( ( local_ice - min_max_dict['nsidc_ice']['min']) / \
                    ( min_max_dict['nsidc_ice']['max'] - min_max_dict['nsidc_ice']['min'] )) * 100.
        new_cod = ( ( data['modis_cod'][:,:] - min_max_dict['modis_cod']['min']) / \
                    ( min_max_dict['modis_cod']['max'] - min_max_dict['modis_cod']['min'] )) * 100.
        #new_ch1 = ( ( data['modis_ch1'][:,:] - min_max_dict['modis_ch1']['min']) / \
        #            ( min_max_dict['modis_ch1']['max'] - min_max_dict['modis_ch1']['min'] )) * 100.
        new_ch7 = ( ( data['modis_ch7'][:,:] - min_max_dict['modis_ch7']['min']) / \
                    ( min_max_dict['modis_ch7']['max'] - min_max_dict['modis_ch7']['min'] )) * 100.
        new_cpr = ( ( local_cldpres - min_max_dict['modis_cld_top_pres']['min']) / \
                    ( min_max_dict['modis_cld_top_pres']['max'] - min_max_dict['modis_cld_top_pres']['min'] )) * 100.
        new_alb = ( ( data['ceres_alb'][:,:] - min_max_dict['ceres_alb']['min']) / \
                    ( min_max_dict['ceres_alb']['max'] - min_max_dict['ceres_alb']['min'] )) * 100.

        # Set the bad pixels to valid numbers for now to use in the calculation,
        # but mask these indices after the calculation
        # ----------------------------------------------------------------------
        new_sza[data_mask.mask] = 50.
        new_vza[data_mask.mask] = 50.
        new_ice[data_mask.mask] = 50.
        new_cod[data_mask.mask] = 50.
        new_ch7[data_mask.mask] = 50.
        new_cpr[data_mask.mask] = 50.
        new_alb[data_mask.mask] = 50.

        new_sza = new_sza.flatten()
        new_vza = new_vza.flatten()
        new_ice = new_ice.flatten()
        new_cod = new_cod.flatten()
        new_ch7 = new_ch7.flatten()
        new_cpr = new_cpr.flatten()
        new_alb = new_alb.flatten()

        print("INPUT SZA", new_sza)

        calc_swf = np.array([model.predict([new_sza, \
                                            new_vza, \
                                            new_ice, \
                                            new_cod, \
                                            #combined_data['modis_ch1'][test_idxs][xx]), 0), \
                                            new_ch7, \
                                            new_cpr, \
                                            new_alb])]).squeeze()


        calc_swf[:] = ((calc_swf[:] / 100) * \
                (min_max_dict['ceres_swf']['max'] - min_max_dict['ceres_swf']['min'])) + \
                min_max_dict['ceres_swf']['min']

        calc_swf = np.reshape(calc_swf, input_shape)
        calc_swf = np.where(data_mask.mask == True, np.nan, calc_swf)

        print('OUTPUT DATA', calc_swf)

        #calc_swf[data_mask.mask] == -999.
        #calc_swf = np.ma.masked_where(data_mask.mask == True, calc_swf)

        """
        for ii in range(calc_swf.shape[0]):
            for jj in range(calc_swf.shape[1]):
                
                ## Prep the MODIS data
                ## -------------------
                #if(np.isnan(data['modis_cld_top_pres'][ii,jj]) == True):
                #    local_cldpres = 1025.
                #elif(data['modis_cld_top_pres'][ii,jj] == 0.):
                #    local_cldpres = 1025.
                #else:
                #    local_cldpres = data['modis_cld_top_pres'][ii,jj]


                # Check the value here
                # --------------------  
                if( (data['omi_lat'][ii,jj] >= minlat) & \
                    (data['modis_cod'][ii,jj] != -999) & \
                    (data['modis_cod'][ii,jj] < max_cod) & \
                    (np.isnan(data['modis_cod'][ii,jj]) == False) & \
                    #(local_cldpres != -999) & \
                    #(local_cldpres <= 1025.) & \
                    #(np.isnan(data['modis_cod'][ii,jj]) == False) & \
                    (data['modis_cld_top_pres'][ii,jj] != -999) & \
                    (data['modis_cld_top_pres'][ii,jj] <= 1025) & \
                    (np.isnan(data['modis_cld_top_pres'][ii,jj]) == False) & \
                    (data['nsidc_ice'][ii,jj] != -999) & \
                    (data['nsidc_ice'][ii,jj] != 251) & \
                    (data['nsidc_ice'][ii,jj] != 253)):
        
                    if(np.isnan(data['modis_cld_top_pres'][ii,jj]) == True):
                        local_cldpres = 1025.
                    elif(data['modis_cld_top_pres'][ii,jj] == 0.):
                        local_cldpres = 1025.
                    else:
                        local_cldpres = data['modis_cld_top_pres'][ii,jj]

                    if(data['nsidc_ice'][ii,jj] == 254.):
                        local_ice = 101.
                    else:
                        local_ice = data['nsidc_ice'][ii,jj]

                    # Scale the input values here
                    # ---------------------------
                    new_sza = ( ( data['omi_sza'][ii,jj] - min_max_dict['omi_sza']['min']) / \
                                ( min_max_dict['omi_sza']['max'] - min_max_dict['omi_sza']['min'] )) * 100.
                    new_vza = ( ( data['omi_vza'][ii,jj] - min_max_dict['omi_vza']['min']) / \
                                ( min_max_dict['omi_vza']['max'] - min_max_dict['omi_vza']['min'] )) * 100.
                    new_ice = ( ( local_ice - min_max_dict['nsidc_ice']['min']) / \
                                ( min_max_dict['nsidc_ice']['max'] - min_max_dict['nsidc_ice']['min'] )) * 100.
                    new_cod = ( ( data['modis_cod'][ii,jj] - min_max_dict['modis_cod']['min']) / \
                                ( min_max_dict['modis_cod']['max'] - min_max_dict['modis_cod']['min'] )) * 100.
                    #new_ch1 = ( ( data['modis_ch1'][ii,jj] - min_max_dict['modis_ch1']['min']) / \
                    #            ( min_max_dict['modis_ch1']['max'] - min_max_dict['modis_ch1']['min'] )) * 100.
                    new_ch7 = ( ( data['modis_ch7'][ii,jj] - min_max_dict['modis_ch7']['min']) / \
                                ( min_max_dict['modis_ch7']['max'] - min_max_dict['modis_ch7']['min'] )) * 100.
                    new_cpr = ( ( local_cldpres - min_max_dict['modis_cld_top_pres']['min']) / \
                                ( min_max_dict['modis_cld_top_pres']['max'] - min_max_dict['modis_cld_top_pres']['min'] )) * 100.
        
                    # Use the model to estimate the SWF
                    # ---------------------------------
                    calc_swf[ii,jj] = np.array([model.predict([np.expand_dims(np.array(new_sza), 0), \
                                                               np.expand_dims(np.array(new_vza), 0), \
                                                               np.expand_dims(np.array(new_ice), 0), \
                                                               np.expand_dims(np.array(new_cod), 0), \
                                                               #np.expand_dims(np.array(new_ch1), 0), \
                                                               np.expand_dims(np.array(new_ch7), 0), \
                                                               np.expand_dims(np.array(new_cpr), 0)]) \
                                        ]).squeeze()
                    calc_swf[ii,jj] = ((calc_swf[ii,jj] / 100) * \
                            (min_max_dict['ceres_swf']['max'] - min_max_dict['ceres_swf']['min'])) + \
                            min_max_dict['ceres_swf']['min']
        
                    #print(ii, jj, data['ceres_swf'][ii,jj], calc_swf[ii,jj])
        """
        
        # Add the calculated values to the output file
        # --------------------------------------------
        dset.create_dataset('omi_lon', data = data['omi_lon'][:,:])
        dset.create_dataset('omi_lat', data = data['omi_lat'][:,:])
        dset.create_dataset('calc_swf', data = calc_swf[:,:])
        dset.create_dataset('ceres_swf', data = data['ceres_swf'][:,:])
        dset.create_dataset('omi_sza', data = data['omi_sza'][:,:])
        dset.create_dataset('omi_uvai_pert', data = data['omi_uvai_pert'][:,:])
        dset.create_dataset('modis_cld_top_pres', data = data['modis_cld_top_pres'][:,:])
        dset.create_dataset('modis_cod', data = data['modis_cod'][:,:])
        dset.create_dataset('nsidc_ice', data = data['nsidc_ice'][:,:])
        dset.close()

print("SUCCESS")
