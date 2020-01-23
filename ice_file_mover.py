#!/usr/bin/env python
"""


"""

import numpy as np
import glob
import os
import sys

base_path = '/home/bsorenson/HighLatitudeStudy/Ice_Analysis/data2/5000000437666'

for root,dirs,files in os.walk(base_path):
    for filename in files:
        if(filename[-3:]=='bin'):
            #print(root+'/'+filename)
            mover = root+'/'+filename
            cmnd = 'cp '+mover+' /home/bsorenson/HighLatitudeStudy/Ice_Analysis/data2/data_files/'
            print(cmnd)
            os.system(cmnd)
