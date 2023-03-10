#!/usr/bin/env python

"""
  To access the variables in the returned dictionary:
    Latitude (2D)
      >>> trop_data['lat']

    Longitude (2D)
      >>> trop_data['lon']

    Time (1D, in datetime format)
      >>> trop_data['time']

    UVAI (2D)
      >>> trop_data['AI']

"""

from TROPOMI_Lib import *

#date_str = '201807051819'
date_str = '201807052142'
#date_str = '201908110044'
#plot_compare_OMI_TROPOMI(date_str, minlat = 65., slope = 'linear', vmin = -2, vmax = 7, save = False)
convert_TROPOMI_to_HDF5(date_str, save_path = home_dir + '/Research/TROPOMI/')

sys.exit()

date_str = '20190506'
#date_str = '20190625'
#download_TROPOMI_match_OMI(date_str, \
#    save_path = home_dir + '/Research/TROPOMI')
plot_TROPOMI_row_avg(date_str, plot_swath = False, minlat = 65., \
    save = False)

sys.exit()

#date_str = '201807052142'
#date_str = '201908110044'
plot_compare_OMI_TROPOMI(date_str, minlat = 65., slope = 'linear', vmin = -2, vmax = 7, save = True)
#plot_TROPOMI_figure(date_str, minlat = 65., vmin = None, vmax = None, \
#        circle_bound = True, ptitle = '', zoom = True, \
#        save = False)
#trop_data = read_TROPOMI(date_str)
sys.exit()
#download_TROPOMI_file(date_str)
convert_TROPOMI_to_HDF5(date_str)
#

#filename = 'S5P_OFFL_L2__AER_AI_20190811T224359_20190812T002529_09471_01_010302_20190817T221032.nc'
