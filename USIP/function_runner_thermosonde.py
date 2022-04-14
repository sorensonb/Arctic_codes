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

in_data = ['2018050505/2018050505/180505_051104_CKN_GRAW.txt',\
    '2018050505/2018050505/HRRR_2018050505_ANALYSIS_CKN_ARL']

##!#model_file = '2018050505/2018050505/HRRR_2018050505_ANALYSIS_CKN_ARL'
##!#thermo_file = '2018050505/2018050505/original_data/18_05_05_05_11_03.GRAW.tempdiff.1Hz'
##!#radio_file =  '2018050505/2018050505/180505_051104_CKN_GRAW.txt'
##!#model_file = '2018050505/2018050505/HRRR_2018050505_ANALYSIS_CKN_ARL'
##!#thermo_file = '2018050505/2018050505/original_data/18_05_05_05_11_03.GRAW.tempdiff.1Hz'
##!#radio_file =  '2018050505/2018050505/180505_051104_CKN_GRAW.txt'

##!#compare_cn2(radio_file,model_file,thermo_file,mlat=None,mlon=None)
##!#
##!#sys.exit()

date_str = '2018050505'
#date_str = '2019050403'

in_data = [file_dict[date_str]['radio_file'], \
    file_dict[date_str]['model_file']]

fig = plt.figure(figsize = (14,4))
gs = fig.add_gridspec(1,3)
#gs = fig.add_gridspec(1,3, hspace = 0.3)
ax1  = fig.add_subplot(gs[0,1])   # true color    
ax2  = fig.add_subplot(gs[0,2])   # true color 
plot_sounding_figure(in_data, fig = fig, skew = gs[0,0], \
    save = False)

thermo_scn2 = read_temp_diffs(file_dict[date_str]['radio_file_orig'], file_dict[date_str]['thermo_file'])
plot_tempdiffs(thermo_scn2, ax = ax1, save = False)


#plot_cn2_figure(in_data, ax = ax2, raw = False, old = False, \
plot_cn2_figure(date_str, ax = ax2, raw = False, old = False, \
        save = False, no_3km = False)
#print(thermo_scn2['CN2T'].shape, thermo_scn2['ALT'].shape)
#ax2.plot(np.log10(thermo_scn2['CN2T']),thermo_scn2['ALT']/1000. ,label='Thermosonde') 

#ax2.plot(np.log10(thermo_scn2['CN2']),(data['ALT']/1000.),color=colors[ii],label=olabel)
plt.show()
#plot_subplot_label(gs[0,0], '(a)', loc = 'upper_right')
