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

if(home_dir + '/Research/OMI/' not in sys.path):
    sys.path.append(home_dir + '/Research/OMI/')

from OMILib import readOMI_swath_shawn, plotOMI_single_swath, \
    readOMI_swath_hdf

#date_str = '201807052142'
#date_str = '201908110044'
date_str = '202108012225'
plot_compare_OMI_TROPOMI(date_str, minlat = 65., slope = 'linear', vmin = -2, vmax = 7, save = False)
#plot_TROPOMI_figure(date_str, minlat = 65., vmin = None, vmax = None, \
#        circle_bound = True, ptitle = '', zoom = True, \
#        save = False)
#trop_data = read_TROPOMI(date_str)
sys.exit()


#date_str = '20180705'
date_str = '20210801'
automate_TROPOMI_preprocess(date_str, download = True, images = True, \
    process = True, remove_bad_OMI = False)

sys.exit()


date_strs = [
    '202108010207', \
    '202108010349', \
    '202108010530', \
    '202108011539', \
    '202108011902', \
    '202108012044', \
    '202108012225', \
]

minlat = 65.
for date_str in date_strs:

    omi_date = trop_to_omi_dict[date_str]

    OMI_base = readOMI_swath_hdf(omi_date, 'control', latmin = minlat, \
        skiprows = [52])

    TROP_data = read_TROPOMI(date_str, minlat = minlat)
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 5))
    ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(1,2,2, projection = mapcrs)

    plotOMI_single_swath(ax1, OMI_base, pvar = 'UVAI', \
        title = 'OMI UVAI Raw', circle_bound = True, gridlines = False, \
        label = 'UVAI', vmin = -2, vmax = 6)

    plot_TROPOMI(TROP_data, ax = ax2, minlat = minlat, vmin = -2, vmax = 6, \
        circle_bound = True, zoom = True, save = False)
    ax2.set_title('TROPOMI UVAI')
    #plot_compare_OMI_TROPOMI(date_str, minlat = 65., slope = 'linear', \
    #    vmin = -2, vmax = 7, save = False)
#convert_TROPOMI_to_HDF5(date_str, save_path = home_dir + '/Research/TROPOMI/')

    plt.show()

sys.exit()

date_strs = ['201807051819',
            '201807052142',
            '201908110044',]
#date_str = '201807051819'
#date_str = '201807052142'
#date_str = '201908110044'
date_str = '20190506'
#date_str = '20190625'
#download_TROPOMI_match_OMI(date_str, \
#    save_path = home_dir + '/Research/TROPOMI')
plot_TROPOMI_row_avg(date_str, plot_swath = False, minlat = 65., \
    save = False)

sys.exit()

#download_TROPOMI_file(date_str)
convert_TROPOMI_to_HDF5(date_str)
#

#filename = 'S5P_OFFL_L2__AER_AI_20190811T224359_20190812T002529_09471_01_010302_20190817T221032.nc'
