#!/usr/bin/env python
"""
Moves NSIDC bin files from the download ZIP location to a permanent
storage location.

"""

import numpy as np
import glob
import os
import sys

# base_path is the location of all the subdirectories containig
# each file
base_path = '/home/bsorenson/Research/Ice_analysis/data/temp_dir'

for root,dirs,files in os.walk(base_path):
    for filename in files:
        if(filename[-3:]=='bin'):
            #print(root+'/'+filename)
            mover = root+'/'+filename
            cmnd = 'cp '+mover+' /home/bsorenson/data/NSIDC/pre_2001/'
            print(cmnd)
            os.system(cmnd)
