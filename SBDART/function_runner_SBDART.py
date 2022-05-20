#!/usr/bin/env python

"""


"""

import importlib, SBDART_Lib
from SBDART_Lib import *
import sys

satellite = 'modis_ch31'
#atms_file = ''

atms_file = '/home/bsorenson/Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
##!#process_SBDART_multi_plume_height_thick(atms_file = atms_file, save = True)
##!#sys.exit()
##!#process_SBDART_multi_lower_tmps(atms_file = '', save = True)


# This makes a re-creation of the SBDART GOES/MODIS figure from the paper
##!#atms_file = '/home/bsorenson/Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
##!## Run num 1
##!#z_maxs = np.array([5.])
##!#z_mins = np.array([0.])
##!#add_wv_mix = np.full((4, len(z_maxs)), np.nan)
##!#add_wv_mix[0,:] = 0. 
##!#add_wv_mix[1,:] = 2. 
##!#add_wv_mix[2,:] = 4. 
##!#add_wv_mix[3,:] = 8. 
##!#
##!#data1 = run_sbdart('goes17_ch08', calc_radiance = True, atms_file = atms_file, \
##!#    z_mins = z_mins, z_maxs = z_maxs, add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
##!#data2 = run_sbdart('goes17_ch09', calc_radiance = True, atms_file = atms_file, \
##!#    z_mins = z_mins, z_maxs = z_maxs, add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
##!#data3 = run_sbdart('goes17_ch10', calc_radiance = True, atms_file = atms_file, \
##!#    z_mins = z_mins, z_maxs = z_maxs, add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
##!#data4 = run_sbdart('modis_ch31', calc_radiance = True, atms_file = atms_file, \
##!#    z_mins = z_mins, z_maxs = z_maxs, add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
##!#
##!#plt.close('all')
##!#fig = plt.figure(figsize = (9, 7))
##!#ax1 = fig.add_subplot(2,2,1)
##!#ax2 = fig.add_subplot(2,2,2)
##!#ax3 = fig.add_subplot(2,2,3)
##!#ax4 = fig.add_subplot(2,2,4)
##!#
##!#plot_bright_vza(data1, pax = ax1)
##!#plot_bright_vza(data2, pax = ax2)
##!#plot_bright_vza(data3, pax = ax3)
##!#plot_bright_vza(data4, pax = ax4)
##!#fig.tight_layout()
##!#plt.show()
##!#
##!#sys.exit()

satellites = ['modis_ch31','goes17_ch08','goes17_ch09','goes17_ch10','goes17_ch13']

for satellite in satellites:

    # ADD MIX RUNS
    
    # Run num 1
    z_maxs = np.arange(1.0, 11., 1.00)
    z_mins = np.full(z_maxs.shape, 0.)
    add_wv_mix = np.full((3, len(z_maxs)), np.nan)
    add_wv_mix[0,:] = 0
    add_wv_mix[1,:] = 2
    add_wv_mix[2,:] = 4
    
    data1 = run_sbdart(satellite, calc_radiance = True, atms_file = atms_file, \
        z_mins = z_mins, z_maxs = z_maxs, add_wv_mix = add_wv_mix)
    
    # Run num 2
    z_maxs = np.arange(2, 11.)
    z_mins = np.arange(0, 9.)
    add_wv_mix = np.full((3, len(z_maxs)), np.nan)
    add_wv_mix[0,:] = 0
    add_wv_mix[1,:] = 2
    add_wv_mix[2,:] = 4
    
    data2 = run_sbdart(satellite, calc_radiance = True, atms_file = atms_file, \
        z_mins = z_mins, z_maxs = z_maxs, add_wv_mix = add_wv_mix)
    
    
    plot_SBDART_dual_height_thickness(data1, data2, save = True)
    plot_SBDART_atmos(data1, ptype = 'thicktop', save = True)
    plot_SBDART_atmos(data2, ptype = 'height', save = True)
    
    # RHUM runs
    
    # Run num 1
    z_maxs = np.arange(1.0, 11., 1.00)
    z_mins = np.full(z_maxs.shape, 0.)
    set_rh = np.full((3, len(z_maxs)), np.nan)
    set_rh[0,:] = None
    set_rh[1,:] = 40
    set_rh[2,:] = 80
    
    data3 = run_sbdart(satellite, calc_radiance = True, atms_file = atms_file, \
        z_mins = z_mins, z_maxs = z_maxs, set_rh = set_rh)
    
    # Run num 2
    z_maxs = np.arange(2, 11.)
    z_mins = np.arange(0, 9.)
    ##!#set_rh = np.full(len(z_maxs), 4.)
    set_rh = np.full((3, len(z_maxs)), np.nan)
    set_rh[0,:] = None
    set_rh[1,:] = 40
    set_rh[2,:] = 80
    
    data4 = run_sbdart(satellite, calc_radiance = True, atms_file = atms_file, \
        z_mins = z_mins, z_maxs = z_maxs, set_rh = set_rh)
    
    
    plot_SBDART_dual_height_thickness(data3, data4, save = True)
    plot_SBDART_atmos(data1, ptype = 'thicktop', save = True)
    plot_SBDART_atmos(data2, ptype = 'height', save = True)
    plot_SBDART_atmos(data3, ptype = 'thicktop', save = True)
    plot_SBDART_atmos(data4, ptype = 'height', save = True)
