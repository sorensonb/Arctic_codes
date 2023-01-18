#!/usr/bin/env python

"""


"""

import importlib, SBDART_Lib
from SBDART_Lib import *
import sys

atms_file = home_dir + '/Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
#process_SBDART_multi_sat_vza_tsfc(atms_file = atms_file, save = True)
#process_SBDART_multi_sat_vza_wv(atms_file = atms_file, save = True)
#process_SBDART_multi_sat_vza_co2(atms_file = atms_file, save = True)
#process_SBDART_multi_sat_vza_ch4(atms_file = atms_file, save = True)
#process_SBDART_multi_sat_vza_no2(atms_file = atms_file, save = True)
#process_SBDART_multi_sat_vza_co(atms_file = atms_file, save = True)
#process_SBDART_multi_sat_vza_nh3(atms_file = atms_file, save = True)


date_str = '20210722212119'
data1_smoke, data1_clear1, data1_clear1_wv, data1_clear1_wv2, data1_clear1_sktp = \
    process_SBDART_multi_sat_vza_CrIS(date_str = date_str, row_str = 'ml', \
        goes_sat = 'goes16', plot_obs = False, save = False)

sys.exit()

satellite = 'modis_ch31'
atms_file = ''
atms_file = '/home/bsorenson/Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
bght_tmps = process_SBDART_single_sat_wv_tmp_vza(satellite = 'modis_ch31',\
    atms_file = atms_file, save = False)

sys.exit()

process_SBDART_multi_plume_height_thick(atms_file = atms_file, save = True)
sys.exit()
#process_SBDART_multi_lower_tmps(atms_file = atms_file, save = False)
#sys.exit()


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
