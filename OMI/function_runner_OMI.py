#!/usr/bin/env python

"""


"""

import OMILib
from OMILib import *
import sys

# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------


###minlat = 65
###plt.close('all')
###fig1 = plt.figure(figsize = (6,6))
###mapcrs = ccrs.NorthPolarStereo()
###ax0 = fig1.add_subplot(1,1,1, projection = mapcrs)
###ax0.coastlines(resolution = '50m')
###ax0.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
###ax0.set_boundary(circle, transform = ax0.transAxes)
###plot_lat_circles(ax0, [80])
###plot_arctic_regions(ax0)
###ax0.gridlines()
###fig1.tight_layout()
####fig1.savefig('arctic_lat_circles_80.png',dpi=300)
###plt.show()

# pvars = 'glint', 'earthshine', 'sunshine', 'stray'
#plotOMI_single_flags('201204101833', pvar = 'glint', minlat = 55, zoom = True)

#plot_combined_fort_out('20170818', min_lat = 80., vtype = 'areas', save = False)


#date_str = '201807050048'
#date_str = '201807050227'
#date_str = '201807050406'
#date_str = '201807050545'
#date_str = '201807050723'
#date_str = '201807050902'
#date_str = '201807051041'
#date_str = '201807051220'
#date_str = '201807051359'
#date_str = '201807051538'
#date_str = '201807051717'
date_str = '201807051856'
#date_str = '201807052034'
#date_str = '201807052213'
#date_str = '201807052352'
minlat = 65.
OMI_base  = readOMI_swath_shawn(date_str, latmin = minlat)
sys.exit()
#write_shawn_to_HDF5(OMI_base)

#OMI_data = readOMI_swath_hdf(date_str, 'control', only_sea_ice = False, \
#    only_ice = False, no_ice = False, latmin = 65, skiprows = [52])
#date_str = '201605152104'
#date_str = '201605162009'


date_strs = ['200804221841',  # GOOD
             '200804222020',  # GOOD
             '200804222159',  # GOOD
             '201605151925',  # MEDIOCRE
             '201605152104',  # MEDIOCRE
             '201605152243',  # MEDIOCRE
             '201605162148',  # MEDIOCRE
             '201807051856',  # GOOD
             '201908102115',  # GOOD
             '201908102254',  # GOOD
             '201908110033',  # GOOD
             '201908110351',  # GOOD
             '200607260017',  # GOOD
             '200607252238',  # GOOD
             '200607260156',  # GOOD
             '200607260335',  # GOOD
             '200607260513',  # GOOD
             '201708161504',  # GOOD
             '201708161643',  # GOOD
             '201708161821',  # GOOD
             '201708171408',  # GOOD
             '201708171547',  # GOOD
             '201708171726',  # GOOD
             '201708171905',  # GOOD
             '201708172043',  # GOOD
             '201708181312',  # GOOD
             '201708181451',  # GOOD
             '201708181630',  # GOOD
             '201708181809',  # GOOD
             '201708181948',  # GOOD
             '201708191355',  # GOOD
             '201708191534',  # GOOD
             '201708191713' ] # GOOD


#date_str = '200804222020' # GOOD
#date_str = '200804222159' # GOOD
#date_str = '201605151925' # MEDIOCRE
#date_str = '201605152104' # MEDIOCRE
#date_str = '201605152243' # MEDIOCRE
#date_str = '201605162148' # MEDIOCRE
#date_str = '201807051856' # GOOD
#date_str = '201908102115' # GOOD
#date_str = '201908102254' # GOOD
#date_str = '201908110033' # GOOD
#date_str = '201908110351' # GOOD
#date_str = '200607260017' # GOOD
#date_str = '200607252238' # GOOD
#date_str = '200607260156' # GOOD
#date_str = '200607260335' # GOOD
#date_str = '200607260513' # GOOD
#date_str = '201708161504' # GOOD
#date_str = '201708161643' # GOOD
#date_str = '201708161821' # GOOD
#date_str = '201708171408' # GOOD
#date_str = '201708171547' # GOOD
#date_str = '201708171726' # GOOD
#date_str = '201708171905' # GOOD
#date_str = '201708172043' # GOOD
#date_str = '201708181312' # GOOD
#date_str = '201708181451' # GOOD
#date_str = '201708181630' # GOOD
#date_str = '201708181809' # GOOD
#date_str = '201708181948' # GOOD
#date_str = '201708191355' # GOOD
#date_str = '201708191534' # GOOD
#date_str = '201708191713' # GOOD

sys.path.append('/home/bsorenson/Research/MODIS/obs_smoke_forcing/')
from MODISLib import *


for date_str in date_strs[:5]:

    OMI_base = readOMI_swath_shawn(date_str, latmin = 65., \
        shawn_path = '/home/bsorenson/data/OMI/shawn_files/')

    CERES_date_str = np.min(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]).strftime('%Y%m%d%H')

    modis_list = download_MODIS_swath(CERES_date_str, \
            dest_dir = '/home/bsorenson/data/MODIS/Aqua/', download = False)

    print(date_str)
    print('    OMI - ', date_str)
    print('  CERES - ', CERES_date_str)
    print('  MODIS - ', *modis_list)
    print('  NSIDC - ', CERES_date_str[:10])
    
    #print('    ', np.min(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]), np.max(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]))

#plotOMI_single_swath_figure(date_str, dtype = 'shawn',  \
#        only_sea_ice = False, minlat = 65., skiprows = None, \
#        lat_circles = None, save = False, zoom = False, \
#        circle_bound = True, ax = None, \
#        shawn_path = '/home/bsorenson/data/OMI/shawn_files/')

sys.exit()

minlat = 65.
OMI_base = date_str
#OMI_base = '201908110351'
write_shawn_to_HDF5(OMI_base, save_path = '/home/bsorenson/Research/Arctic_compares/comp_data/20180705/', minlat = 65., \
    shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
sys.exit()

plot_compare_OMI_CERES_MODIS_NSIDC('201908110125', '7', \
    omi_dtype = 'shawn', minlat = 65., zoom = True, save = False)
sys.exit()

sys.exit()


#plot_compare_OMI_CERES_MODIS_NSIDC('201808241435', '7', \



#plotOMI_single_swath_figure(date_str, dtype = 'shawn',  \
#        only_sea_ice = False, minlat = 65., skiprows = None, \
#        lat_circles = None, save = False, zoom = True)
#
sys.exit()

plotOMI_daily_control_shawn('20170818', resolution = 1.0,
    shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
sys.exit()

plot_Arctic_row_coverage_compare(date_str = '20180726', save = False)

minlat = 65.
OMI_VBS0   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VBS0_2005_2020.nc', minlat = minlat)
plotOMI_NCDF_Climo_SpringSummer(OMI_VBS0,start_idx=0,end_idx=96,minlat=65,\
                   save=True)

sys.exit()
plot_combined_fort_out('20190811', min_lat = 70., vtype = 'areas', max_lat = 80., save = True)

# NOTE: for plotting the bias between OMI rows, use this line with the
#       CSCI netCDF data
#plt.plot(np.nanmean(np.nanmean(netdata['AI'], axis = 0), axis = 0))
##!#infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_100.txt'
##!#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = True)
##!#infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_150.txt'
##!#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = True)
##!#infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_200.txt'
##!#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = True)
##!#infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_250.txt'
##!#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = True)

sys.exit()

OMI_VBS0   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VBS0_2005_2020.nc', minlat = minlat)
OMI_VJZ211 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VJZ211_2005_2020.nc', minlat = minlat)
OMI_VSJ4   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', minlat = minlat)
#plotOMI_Compare_ClimoTrend_summer(OMI_VBS0,OMI_VJZ211,OMI_VSJ4,\
#        trend_type = 'standard', minlat=minlat,save=False)
plotOMI_Compare_ClimoTrend_all(OMI_VBS0,OMI_VJZ211, OMI_VSJ4,\
        trend_type = 'standard', minlat=minlat,save = False)

sys.exit()

plot_row_bias(save = True)
sys.exit()
plot_row_anomaly_combined(date_str = '201807260244', dtype = 'control', \
        minlat = 65., save = True)

sys.exit()


#OMI_data, CERES_data = plot_OMI_CERES_trend_compare_summer(minlat=72,\
#        ceres_type = 'sw', trend_type = 'standard', save=False)
##!#infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_100.txt'
##!#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = False)
#plot_OMI_fort_out_peaks(infile, minlat = 70., vtype = 'areas', save = False)

plotOMI_single_swath_multiple(dtype = 'control',  \
        only_sea_ice = False, minlat = 65., save = True)




##bad_row_file = 'row_anomaly_dates_20050401_20201001.txt'
##xtrack_file = 'row_anomaly_xtrack_dates_20050401_20201001.txt'
##plot_bad_row_table(bad_row_file, xtrack_file = xtrack_file, ax = None, \
##        save = False)

sys.exit()

#date_str = '200804222159'
#plotOMI_single_multipanel(date_str, only_sea_ice = False, minlat = 65., \
#        quad_panel = True, save = True)
#plotOMI_single_swath_multiple('22222222', dtype = 'control',  \
#    only_sea_ice = False, minlat = 65., save = True)

plot_OMI_fort_out_func(infile,\
     min_lat = 70., vtype = 'areas', save = False)

sys.exit()

date_str = ['200509270134',\
            '200509270313',\
            '200509270451',\
            '200509270630',\
            '200509270809',\
            '200509270948',\
            '200509271127',\
            '200509271306',\
            '200509271445',\
            '200509271624',\
            '200509271802',\
            '200509271941',\
            '200509272120',\
            '200509272259']

###
###    
###plotOMI_single_ground(date_str, only_sea_ice = False, minlat = 65., \
###    zoom = True, multi_panel = False, save = True)

plotOMI_single_swath_figure(date_str, dtype = 'control',  \
        only_sea_ice = False, minlat = 65., skiprows = [11,15, 16, 17, 18], lat_circles = None, save = False)

#plot_combined_figure1(save = True)
#plot_MODIS_temporary('202107222110', zoom = True, save = True)
#plot_MODIS_temporary_4panel('202108052125', zoom = True, composite = True, show_smoke = True, save = True)
#compare_MODIS_3panel('202107222110',31,1,5,zoom=True,save=True,\
#        plot_ASOS_loc = False, show_smoke = True)
#plot_combined_figure3(save = True)
#plot_figureS1(save=True, composite = True)
#plot_combined_figureS2(save=True)
