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

plot_monthly_cn2_boxplots_multi('201705 12','201705 00','201707 00','201702 00',save = True)

##!#date_str = '201705'
##!#time = '12'
##!##data = calc_monthly_climatology(date_str, time, location = 'BIS')
##!#fig = plt.figure(figsize = (8, 8))
##!#ax1 = fig.add_subplot(2,2,1)
##!#ax2 = fig.add_subplot(2,2,2)
##!#ax3 = fig.add_subplot(2,2,3)
##!#ax4 = fig.add_subplot(2,2,4)
##!#plot_monthly_cn2_boxplots('201705', '12', location = 'BIS', ax = ax1)
##!#plot_monthly_cn2_boxplots('201705', '00', location = 'BIS', ax = ax2)
##!#plot_monthly_cn2_boxplots('201711', '12', location = 'BIS', ax = ax3)
##!#plot_monthly_cn2_boxplots('201711', '00', location = 'BIS', ax = ax4)
##!#
##!#plot_subplot_label(ax1, '(a)', location = 'upper_right')
##!#plot_subplot_label(ax2, '(b)', location = 'upper_right')
##!#plot_subplot_label(ax3, '(c)', location = 'upper_right')
##!#plot_subplot_label(ax4, '(d)', location = 'upper_right')
##!#
##!#plt.suptitle('Monthly C$_{n}^{2}$ Estimated Climatology\nBismarck, ND')
##!#ax1.set_title('May 2017\nGenerated with 00Z Soundings ', fontsize = 10)
##!#ax2.set_title('May 2017\nGenerated with 12Z Soundings ', fontsize = 10)
##!#ax3.set_title('Nov. 2017\nGenerated with 00Z Soundings', fontsize = 10)
##!#ax4.set_title('Nov. 2017\nGenerated with 12Z Soundings', fontsize = 10)
##!#
##!#fig.tight_layout()
##!#plt.show()

sys.exit()

date_str = '2018050505'
plot_combined_figure(date_str, save = True)
date_str = '2019050403'
plot_combined_figure(date_str, save = True)

sys.exit()


in_data = [file_dict[date_str]['radio_file'], \
    file_dict[date_str]['model_file'], file_dict[date_str]['thermo_file']]

# Use readsounding to read the current sounding file into a dictionary
radio = readsounding(in_data[0])

# Read the thermosonde data
thermo_scn2 = read_temp_diffs(file_dict[date_str]['radio_file_orig'], in_data[2])

sys.exit()

plot_dual_vert_tempdiffs(date_str, save = True)

date_str = '201805 12'
plot_synthetic_figure(date_str, save = False)

sys.exit()
#read_synthetic_data('201805 00', 'BIS')
synth_data = '201905 12'
pvar = 'tempdiff'
plot_synthetic(synth_data, pvar, location = 'BIS', ax = None, save = False)

sys.exit()

#in_data = ['2018050505/2018050505/180505_051104_CKN_GRAW.txt',\
#    '2018050505/2018050505/HRRR_2018050505_ANALYSIS_CKN_ARL']

##!#model_file = '2018050505/2018050505/HRRR_2018050505_ANALYSIS_CKN_ARL'
##!#thermo_file = '2018050505/2018050505/original_data/18_05_05_05_11_03.GRAW.tempdiff.1Hz'
##!#radio_file =  '2018050505/2018050505/180505_051104_CKN_GRAW.txt'
##!#model_file = '2018050505/2018050505/HRRR_2018050505_ANALYSIS_CKN_ARL'
##!#thermo_file = '2018050505/2018050505/original_data/18_05_05_05_11_03.GRAW.tempdiff.1Hz'
##!#radio_file =  '2018050505/2018050505/180505_051104_CKN_GRAW.txt'

##!#compare_cn2(radio_file,model_file,thermo_file,mlat=None,mlon=None)
##!#
##!#sys.exit()


plot_calibration_curves_both(save = True)

sys.exit()

date_strs = ['2018050505', '2019050403']

for date_str in date_strs:
    plot_combined_figure(date_str, save = True)

##!#in_data = [file_dict[date_str]['radio_file'], \
##!#    file_dict[date_str]['model_file']]
##!#
##!#fig = plt.figure(figsize = (14,4))
##!#gs = fig.add_gridspec(1,3)
##!##gs = fig.add_gridspec(1,3, hspace = 0.3)
##!#ax1  = fig.add_subplot(gs[0,1])   # true color    
##!#ax2  = fig.add_subplot(gs[0,2])   # true color 
##!#plot_sounding_figure(in_data, fig = fig, skew = gs[0,0], \
##!#    save = False)
##!#
##!#thermo_scn2 = read_temp_diffs(file_dict[date_str]['radio_file_orig'], file_dict[date_str]['thermo_file'])
##!#plot_tempdiffs(thermo_scn2, ax = ax1, save = False)
##!#
##!#ax1_lims = ax1.get_ylim()
##!#print(ax1_lims, ax1_lims[1])
##!##plot_cn2_figure(in_data, ax = ax2, raw = False, old = False, \
##!#plot_cn2_figure(date_str, ax = ax2, raw = False, old = False, \
##!#        save = False, no_3km = False, ymax = ax1_lims[1])
##!##print(thermo_scn2['CN2T'].shape, thermo_scn2['ALT'].shape)
##!##ax2.plot(np.log10(thermo_scn2['CN2T']),thermo_scn2['ALT']/1000. ,label='Thermosonde') 
##!#
##!##ax2.plot(np.log10(thermo_scn2['CN2']),(data['ALT']/1000.),color=colors[ii],label=olabel)
##!#outname = 
##!#plt.show()
##!##plot_subplot_label(gs[0,0], '(a)', loc = 'upper_right')
