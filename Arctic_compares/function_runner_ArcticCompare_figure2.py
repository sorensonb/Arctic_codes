#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *
import random
#from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

#make_gif('comp_images_20180705/', 'calc_swf_comp_20180705.gif')
date_str = '201807052213'
plot_compare_OMI_MODIS_v2(date_str, 7, \
    omi_dtype = 'ltc3', minlat = 65., zoom = True, save = False)
sys.exit()
