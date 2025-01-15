#!/usr/bin/env python

"""


"""

import numpy as np
import sys

sys.path.append('/home/bsorenson')
from sounding_lib import *
from python_lib import *
from cn2_lib import *
from compare_cn2 import *
from datetime import datetime

date_str = '2018050505'
#date_str = '2019050403'
#plot_dual_vert_tempdiffs(date_str, save = False)
#sys.exit()

in_data = [file_dict[date_str]['radio_file'], \
    file_dict[date_str]['model_file'], file_dict[date_str]['thermo_file']]

# Use readsounding to read the current sounding file into a dictionary
radio = readsounding(in_data[0])

# Read the thermosonde data
thermo_scn2 = read_temp_diffs(file_dict[date_str]['radio_file_orig'], in_data[2])

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

good_indices = np.where(np.diff(radio['TEMP'])!=0.)
radio_diffs = np.diff(radio['TEMP'])[good_indices]/np.diff(radio['ALT'])[good_indices]
#if(calc_abs == True):
#    radio_diffs = abs(radio_diffs)
ax.plot(thermo_scn2['closetime_thermo'][thermo_scn2['masked_indices']],\
    thermo_scn2['tempdiff'][thermo_scn2['masked_indices']],label='Thermosonde [horizontal]')
ax.plot(radio['UTCTime'][:-1][good_indices],radio_diffs,label='Radiosonde [vertical]')

plt.show()

