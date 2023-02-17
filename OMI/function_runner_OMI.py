#!/usr/bin/env python

"""


"""

from matplotlib.lines import Line2D
import OMILib
from OMILib import *
import sys

#date_str = '20190625'
date_str = '20190506'
#swath_ais, swath_lat, swath_xtk, swath_xrw, day_avgs, row_avgs = \
plot_OMI_row_avg(date_str, plot_swath = False, minlat = 65., save = True)

sys.exit()

##!#date_strs = [
##!#             '201904290123',
##!#             '201904290302',
##!#             '201904290440',
##!#             '201904290619',
##!#             '201904290758',
##!#             '201904290937',
##!#             '201904291116',
##!#             '201904291255',
##!#             '201904291434',
##!#             '201904291613',
##!#             '201904291752',
##!#             '201904291930',
##!#             '201904292109',
##!#             '201904292248',
##!#             #'201907170041',
##!#             #'201907170220',
##!#             #'201907170359',
##!#             #'201907170537',
##!#             #'201907170716',
##!#             #'201907170855',
##!#             #'201907171034',
##!#             #'201907171213',
##!#             #'201907171352',
##!#             #'201907171531',
##!#             #'201907171710',
##!#             #'201907171848',
##!#             #'201907172027',
##!#             #'201907172206',
##!#             #'201907172345',
##!#            ]

plot_combined_fort_out('20190811', min_lat = 70., vtype = 'areas', max_lat = 80., save = True)
plot_combined_fort_out('20170818', min_lat = 80., vtype = 'areas', save = True)
sys.exit()

# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------
minlat = 65.

OMI_VBS0   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VBS0_2005_2020.nc', minlat = minlat)
OMI_VJZ211 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VJZ211_2005_2020.nc', minlat = minlat)
OMI_VSJ4   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', minlat = minlat)
#plotOMI_Compare_ClimoTrend_summer(OMI_VBS0,OMI_VJZ211,OMI_VSJ4,\
#        trend_type = 'standard', minlat=minlat,save=False)
plotOMI_Compare_ClimoTrend_all(OMI_VBS0,OMI_VJZ211, OMI_VSJ4,\
        trend_type = 'standard', minlat=minlat,save = True)

sys.exit()

date_str = '200804222159'
plotOMI_single_multipanel(date_str, only_sea_ice = False, minlat = 65., \
       quad_panel = True, save = True)
#plotOMI_single_swath_multiple('22222222', dtype = 'control',  \
#    only_sea_ice = False, minlat = 65., save = True)
sys.exit()

plot_row_bias(dataset = 'normal', save = True)
sys.exit()

plotOMI_single_swath_multiple_v3(dtype = 'control',  \
        only_sea_ice = False, minlat = 65., save = True)
#plotOMI_single_swath_multiple_v4(dtype = 'control',  \
#        only_sea_ice = False, minlat = 65., save = True)

sys.exit()

plot_row_anomaly_combined(dtype = 'control', minlat = 65., save = True)
sys.exit()

minlat = 65.
OMI_VBS0   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VBS0_2005_2020.nc', minlat = minlat)
plotOMI_NCDF_Climo_SpringSummer(OMI_VBS0,start_idx=0,end_idx=96,minlat=65,\
                   save = True)
sys.exit()

min_lat = 70.
infile = home_dir + '/Research/OMI/shawn_analysis/count_analysis/' + \
    'omi_vsj4_areas_2005_2020_100.txt'

#calc_aerosol_event_dates(infile, minlat = min_lat, vtype = 'areas', \
#    max_lat = None)
#
#sys.exit()

vtype = 'areas'
omi_fort_dict = read_OMI_fort_out(infile, min_lat, vtype = vtype)
#omi_fort_dict_80 = read_OMI_fort_out(infile, min_lat, vtype = vtype)

fig1 = plt.figure(figsize = (8,5)) 
ax1 = fig1.add_subplot(1,1,1)
#ax2 = fig1.add_subplot(1,2,2)

loop_dates = omi_fort_dict['daily_dt_dates'][np.where( \
    (omi_fort_dict['daily_dt_dates'] >= datetime(2019,4,1)) & \
    (omi_fort_dict['daily_dt_dates'] <= datetime(2019,9,30)))]

plot_dates = datetime(year = 2000, month = 4, day = 1) + \
    (loop_dates - datetime(year = 2019, month = 4, day = 1))
year = 2019
ax1.plot(plot_dates,\
    omi_fort_dict['daily_data'][np.where( \
    (omi_fort_dict['daily_dt_dates'] >= datetime(year,4,1)) & \
    (omi_fort_dict['daily_dt_dates'] <= datetime(year,9,30)))], label = '2019')


# Plot the time series with the peaks as 'x's
# -------------------------------------------
#plot_OMI_fort_out_time_series(omi_fort_dict, ax1, minlat = min_lat)
#custom_lines = [Line2D([0], [0])]

print(ax1.get_xlim())

ax1.legend()
ax1.grid()
#ax1.legend(custom_lines, ['2019'], loc = 2)

exponent = 5
label_str = 'Area of high AI [10$^{' + \
    str(exponent) + '}$ km$^{2}$]'
ax1.set_ylabel(label_str, weight = 'bold')

ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax1.set_title('High AI ' + vtype + '\nPerturbed AI threshold of '+str(omi_fort_dict['ai_thresh'])+\
    '\nNorth of '+str(int(min_lat))+'$^{o}$N')
fig1.savefig('omi_daily_areas_minlat70_2019.png', dpi = 300)
plt.show()
sys.exit()







month_idx =  4
trend_type = 'standard'
inputfile = '/home/bsorenson/Research/OMI/JZ_analysis/climo_analysis/omi_vsj4_climo_2005_2020.txt'
OMI_data = readOMI_NCDF('/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', minlat = minlat)

#trends, pvals = calcOMI_grid_trend(OMI_data, month_idx, trend_type, minlat)
plotOMI_MonthTrend(OMI_data,month_idx=4,save=False,\
    trend_type='standard',minlat=65.,return_trend=False, colorbar = True, \
    title = '', label = '', colorbar_label_size = 14, pax = None)

sys.exit()

plotOMI_Compare_ClimoTrend_all(OMI_data1,OMI_data2,OMI_data3,\
    trend_type = 'standard', minlat=65,save=False)

plotOMI_Type_Trend_all(trend_type = 'standard', \
    minlat=65,save=True)

sys.exit()

start_date = 200504
end_date = 202009
inputfile = '/home/bsorenson/Research/OMI/JZ_analysis/climo_analysis/omi_jz211_climo_2005_2020.txt'
OMI_data = readOMI(inputfile,start_date,end_date,key=None)
filename = '/home/bsorenson/Research/OMI/sorenson_et_al_data/omi_ai_screened_2005_2020.nc'
writeOMI_toNCDF(OMI_data,filename, minlat = 65.)

inputfile = '/home/bsorenson/Research/OMI/shawn_analysis/climo_analysis/omi_vsj4_climo_2005_2020.txt'
OMI_data = readOMI(inputfile,start_date,end_date,key=None)
filename = '/home/bsorenson/Research/OMI/sorenson_et_al_data/omi_ai_perturbed_2005_2020.nc'
writeOMI_toNCDF(OMI_data,filename, minlat = 65.)
sys.exit()

#plotOMI_Type_Climo_all(trend_type = 'standard', \
#    minlat=65,save=True)

#minlat = 65.
#
#data_types = ['VBS1','VJZ2','VJZ28','VJZ29','VJZ211','VSJ2','VSJ4']    
#
## Read in all the data
#OMI_VBS1 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[0] + '_2005_2019.nc', minlat = minlat)
#OMI_VJZ2 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[1] + '_2005_2019.nc', minlat = minlat)
#OMI_VJZ28= readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[2] + '_2005_2019.nc', minlat = minlat)
#OMI_VJZ29= readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[3] + '_2005_2019.nc', minlat = minlat)
#OMI_VJZ211 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[4] + '_2005_2019.nc', minlat = minlat)
#OMI_VSJ2 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[5] + '_2005_2019.nc', minlat = minlat)
#OMI_VSJ4 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[6] + '_2005_2020.nc', minlat = minlat)
#sys.exit()
plotOMI_Type_MonthClimo_all(minlat=65, save = True, \
        save_dir = '/home/bsorenson/Research/OMI/monthly_images/combined_climo/')

sys.exit()


minlat = 65.
base_dir = '/home/bsorenson/Research/OMI/'
image_dir = base_dir + 'monthly_images/'
data_types = ['VJZ2','VJZ28','VJZ29','VJZ211']    
for dtype in data_types:
    print('\n\nNOW PROCESSING',dtype,'\n\n')

    local_dir = image_dir + dtype + '/'

    # Change to Arctic codes directory
    os.chdir(local_dir)

    print(os.system('pwd'))

    # Read in the data
    OMI_data = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
        'omi_ai_' + dtype + '_2005_2019.nc', minlat = minlat)

    for ii in range(len(OMI_data['DATES'])):
        plotOMI_NCDF_SingleMonth(OMI_data,ii,minlat = minlat, pax = None, \
            save=True)

sys.exit()


infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_100.txt'
#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = False)
plot_OMI_fort_out_peaks(infile, minlat = 75., vtype = 'areas', save = False)

sys.exit()

#plot_row_bias(dataset = 'normal', save = True)
#sys.exit()

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
#OMI_base  = readOMI_swath_shawn(date_str, latmin = minlat)
#sys.exit()
#write_shawn_to_HDF5(OMI_base)

#OMI_data = readOMI_swath_hdf(date_str, 'control', only_sea_ice = False, \
#    only_ice = False, no_ice = False, latmin = 65, skiprows = [52])
#date_str = '201605152104'
#date_str = '201605162009'

date_strs = ['200607240029', # GOOD
             '200607240208', # GOOD
             '200607240347', # GOOD
             '200607240526', # GOOD
             '200607240705', # GOOD
             '200607240844', # GOOD
             '200607242155', # GOOD
             '200607242334', # GOOD
             '200607250112', # GOOD
             '200607250251', # GOOD
             '200607250609', # GOOD
             '200607250748', # GOOD?
             '200607252238', # GOOD
             '200607260017', # GOOD
             '200607260156', # GOOD
             '200607260335', # GOOD
             '200607260513', # GOOD?
             '200607260831', # GOOD
             '200607262142', # GOOD
             '200607270100', # GOOD
             '200607270239', # GOOD?
             '200607270418', # GOOD?
             '200607270557', # GOOD?
             '200607270736', # GOOD?
             '200607270914',
             '200607271053',
             '200607271232',
             '200607271411',
             '200607271550',
             '200607271729',
             '200607271908',
             '200607272047', 
             '200607272226'] # GOOD



#date_strs = ['200804221841',  # GOOD
#             '200804222020',  # GOOD
#             '200804222159',  # GOOD
#date_strs = ['201207260009', 
#             '201207260148', 
#             '201207260327', 
#             '201207260506', 
#             '201207260645', 
#             '201207260824', 
#             '201207261003', 
#             '201207261142', 
#             '201207261320', 
#             '201207261459', 
#             '201207261638', 
#             '201207261817', 
#             '201207261956', 
#             '201207262135', 
#             '201207262314']
#date_strs = ['201207250105', 
#             '201207250244',
#             '201207250423', 
#             '201207250602', 
#             '201207250741', 
#             '201207250920', 
#             '201207251058', 
#             '201207251237', 
#             '201207251416', 
#             '201207251555', 
#             '201207251734', 
#             '201207251913', # good? 
#             '201207252052', # good?
#             '201207252231'] # good?
#date_strs = ['201308020123',
#             '201308020302',
#             '201308020441',
#             '201308020620',
#             '201308020759',
#             '201308020938',
#             '201308021117',
#             '201308021256',
#             '201308021434',
#             '201308021613',
#             '201308021752',
#             '201308021931',
#             '201308022110',
#             '201308022249']
#             '201605151925',  # MEDIOCRE
#             '201605152104',  # MEDIOCRE
#             '201605152243',  # MEDIOCRE
#             '201605162148',  # MEDIOCRE
#             '201807051856',  # GOOD
#             '201908102115',  # GOOD
#             '201908102254',  # GOOD
#             '201908110033',  # GOOD
#             '201908110351',  # GOOD
#             '200607260017',  # GOOD
#             '200607252238',  # GOOD
#             '200607260156',  # GOOD
#             '200607260335',  # GOOD
#             '200607260513',  # GOOD
#             '201708161504',  # GOOD
#             '201708161643',  # GOOD
#             '201708161821',  # GOOD
#             '201708171408',  # GOOD
#             '201708171547',  # GOOD
#             '201708171726',  # GOOD
#             '201708171905',  # GOOD
#             '201708172043',  # GOOD
#             '201708181312',  # GOOD
#             '201708181451',  # GOOD
#             '201708181630',  # GOOD
#             '201708181809',  # GOOD
#             '201708181948',  # GOOD
#             '201708191355',  # GOOD
#             '201708191534',  # GOOD
#             '201708191713' ] # GOOD


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

#sys.path.append('/home/bsorenson/Research/MODIS/obs_smoke_forcing/')
#from MODISLib import *

date_str = '200607270100'
plotOMI_single_swath_figure(date_str, dtype = 'shawn',  \
        only_sea_ice = False, minlat = 65., skiprows = None, \
        lat_circles = None, save = False, zoom = False, \
        circle_bound = True, ax = None, \
        shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
sys.exit()
#for date_str in date_strs[:5]:
for date_str in date_strs:

    plotOMI_single_swath_figure(date_str, dtype = 'shawn',  \
            only_sea_ice = False, minlat = 65., skiprows = None, \
            lat_circles = None, save = False, zoom = False, \
            circle_bound = True, ax = None, \
            shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
##!#    OMI_base = readOMI_swath_shawn(date_str, latmin = 65., \
##!#        shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
##!#
##!#    CERES_date_str = np.min(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]).strftime('%Y%m%d%H')
##!#
##!#    modis_list = download_MODIS_swath(CERES_date_str, \
##!#            dest_dir = '/home/bsorenson/data/MODIS/Aqua/', download = False)
##!#
##!#    print(date_str)
##!#    print('    OMI - ', date_str)
##!#    print('  CERES - ', CERES_date_str)
##!#    print('  MODIS - ', *modis_list)
##!#    print('  NSIDC - ', CERES_date_str[:10])
##!#    
##!#    #print('    ', np.min(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]), np.max(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]))


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




#OMI_data, CERES_data = plot_OMI_CERES_trend_compare_summer(minlat=72,\
#        ceres_type = 'sw', trend_type = 'standard', save=False)




##bad_row_file = 'row_anomaly_dates_20050401_20201001.txt'
##xtrack_file = 'row_anomaly_xtrack_dates_20050401_20201001.txt'
##plot_bad_row_table(bad_row_file, xtrack_file = xtrack_file, ax = None, \
##        save = False)

sys.exit()

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
