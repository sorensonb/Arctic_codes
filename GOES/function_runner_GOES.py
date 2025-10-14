#!/usr/bin/env python

"""


"""

from GOESLib import *
import sys

date_str = '202404081745'
ch1 = 2
ch2 = 13
#region = 'indiana'
#region = 'missouri'
region = 'missouri_bootheel'
#region = 'arkansas'
sat = 'goes16'
begin_date = '202404081200'
end_date = '202404082330'

#lats_points = np.array([region_dict[region]['point_coords'][str(pi)][0] for pi in range(1,6)])
#lons_points = np.array([region_dict[region]['point_coords'][str(pi)][1] for pi in range(1,6)])
#
#GOES_dict_points = read_GOES_time_series_auto(begin_date, end_date, \
#    channels = [2, 13], dlat = list(lats_points), \
#    dlon = list(lons_points), \
#    sat = sat, region = region)

GOES_dict_reg = \
        read_GOES_time_series_auto_regional(begin_date, end_date, \
        channels = [ch1, ch2], save_dir = './', \
        sat = sat, \
        minlat = region_dict[region]['minlat_data'], \
        maxlat = region_dict[region]['maxlat_data'], \
        minlon = region_dict[region]['minlon_data'], \
        maxlon = region_dict[region]['maxlon_data'], \
        min_max_use = ['max', 'max'])


plot_GOES_eclipse_comp(date_str, ch1, ch2, region, \
    GOES_dict_reg, sat = sat, plot_asos = True, \
    asos_site = 'TKX', \
    GOES_dict_points = None, \
    #GOES_dict_points = GOES_dict_points, \
    plot_point_BTs = False, save = False)
sys.exit()

begin_dt_date = datetime(2024,4,8,12,0)
end_dt_date   = datetime(2024,4,8,23,00)

local_dt_date = begin_dt_date

#date_str = '202404082100'
while(local_dt_date <= end_dt_date):
  
    print(local_dt_date) 
    date_str = local_dt_date.strftime('%Y%m%d%H%M')

    if( (local_dt_date != datetime(2024,4,8,16,30) ) & \
        (local_dt_date != datetime(2024,4,8,16,35) ) ):

        plot_GOES_eclipse_comp(date_str, ch1, ch2, region, \
            GOES_dict_reg, sat = sat, plot_asos = True, \
            asos_site = 'TKX', \
            GOES_dict_points = None, \
            #GOES_dict_points = GOES_dict_points, \
            plot_point_BTs = False, save = True)
   
    if( (local_dt_date >= datetime(2024,4,8,15,50) ) & \
        (local_dt_date <= datetime(2024,4,8,20,55) ) ):
        local_dt_date = local_dt_date + timedelta(minutes = 5)
    else:
        local_dt_date = local_dt_date + timedelta(minutes = 30)

sys.exit()




goes_file = 'goes16_cross_data_asos_202107221201_202107230231.nc'
GOES_dict1 = read_GOES_time_series_NCDF(goes_file)

asos_file = home_dir + '/Research/MODIS/obs_smoke_forcing/asos_data_20210722_4.csv'

df = pd.read_csv(asos_file)
df['tmpc'] = pd.to_numeric(df['tmpc'], errors = 'coerce').values + 273.15
#time_df = df.set_index('valid')
df['valid'] = pd.to_datetime(df['valid'])

stations = ['O05','AAT']

ch_idx = 2

fig = plt.figure(figsize = (8, 4))
ax1 = fig.add_subplot(1,1,1)

for jj in range(GOES_dict1['data'].shape[2]):

    ax1.plot(GOES_dict1['dt_dates'], GOES_dict1['data'][:,ch_idx,jj], \
        label = 'GOES16 TIR BT @ ' + stations[jj])

    df_sub = df[df['station'] == stations[jj]]

    good_idxs = np.where( (df_sub['valid'] >= GOES_dict1['dt_dates'][0]) & (df_sub['valid'] < GOES_dict1['dt_dates'][-1]))
    
    tmps  = df_sub['tmpc'].values[good_idxs]
    times = df_sub['valid'].values[good_idxs]
  
    # Get rid of masked values
    # ------------------------
    good_idxs = np.where(~np.isnan(tmps))
 
    ax1.plot(times[good_idxs], tmps[good_idxs], linestyle = ':', \
        label = stations[jj] + ' 2-m temp')
 
ax1.grid()
ax1.legend()
ax1.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
ax1.tick_params(axis="x", labelsize = 9)
ax1.set_ylabel('Temperature [K]')
ax1.set_title('GOES-16 vs ASOS\n20210722')

fig.tight_layout()

outname = 'goes16_asos_comp_20210722.png'
fig.savefig(outname, dpi = 200)
print("Saved image", outname)

plt.show()


sys.exit()

frame_folder = home_dir + '/Research/GOES/goes_eclipse_analysis/missouri_bootheel/'
gif_name = 'goes16_eclipse_comp_20240408_missouri_bootheel_asosTKX.gif'
make_gif(frame_folder, gif_name, duration = 250)

sys.exit()

#plot_GOES_eclipse_comp(date_str, ch1, ch2, region, \
#    GOES_dict_reg, sat = sat, plot_asos = False, \
#    GOES_dict_points = None, plot_point_BTs = False)

sys.exit()

#begin_date = '202404081200'
begin_date = '202107210300'
end_date   = '202107220230'
auto_GOES_download(begin_date, end_date, 30, sat = 'goes17', channels = [2,6,13])
#auto_GOES_download(begin_date, end_date, 30, sat = 'goes17', channels = [2,6,13])
sys.exit()



#plot_GOES_eclipse_comp(date_str, ch1, ch2, region, \
#    GOES_dict_reg, sat = sat, plot_asos = False, \
#    GOES_dict_points = None, plot_point_BTs = False)
plot_GOES_eclipse_comp(date_str, ch1, ch2, region, \
    GOES_dict_reg, sat = sat, plot_asos = False, \
    GOES_dict_points = GOES_dict_points, plot_point_BTs = False)



sys.exit()


begin_date = '202404081200'
end_date   = '202404082330'

# Coordinates for zoom area in AR area:
#minlat = 34.0
#maxlat = 36.0
#minlon = -94.5
#maxlon = -91.5

# Coordinates for wide area in TX/OK/AR area:
#minlat = 33.0
#maxlat = 36.0
#minlon = -95.0
#maxlon = -91.0

# Coordinates for zoom area in MO/AR/TN tri-state area:
#minlat = 36.0
#maxlat = 38.5
#minlon = -91.5
#maxlon = -89.0


kmaw_asos_tmps = [\
    8.00, \
    9.00, \
    10.00, \
    11.00, \
    13.00, \
    15.00, \
    16.00, \
    18.00, \
    19.00, \
    20.00, \
    21.00, \
    22.00, \
    23.00, \
    25.00, \
    25.00, \
    25.00, \
    25.00, \
    26.00, \
    26.00, \
    26.00, \
    26.00, \
    25.00, \
    23.00, \
    24.00, \
    26.00, \
    27.00, \
    27.00, \
    27.00, \
    28.00, \
    27.00, \
    28.00, \
    27.00, \
    27.00, \
    27.00, \
    27.00, \
    26.00, \
    25.00, \
]

kmaw_asos_tmps = np.array(kmaw_asos_tmps) + 273.15

kmaw_asos_times = [
    '2024-04-08 11:55',
    '2024-04-08 12:15',
    '2024-04-08 12:35',
    '2024-04-08 12:55',
    '2024-04-08 13:15',
    '2024-04-08 13:35',
    '2024-04-08 13:55',
    '2024-04-08 14:15',
    '2024-04-08 14:35',
    '2024-04-08 14:55',
    '2024-04-08 15:15',
    '2024-04-08 15:35',
    '2024-04-08 15:55',
    '2024-04-08 16:15',
    '2024-04-08 16:35',
    '2024-04-08 16:55',
    '2024-04-08 17:15',
    '2024-04-08 17:35',
    '2024-04-08 17:55',
    '2024-04-08 18:15',
    '2024-04-08 18:35',
    '2024-04-08 18:55',
    '2024-04-08 19:15',
    '2024-04-08 19:35',
    '2024-04-08 19:55',
    '2024-04-08 20:15',
    '2024-04-08 20:35',
    '2024-04-08 20:55',
    '2024-04-08 21:15',
    '2024-04-08 21:35',
    '2024-04-08 21:55',
    '2024-04-08 22:15',
    '2024-04-08 22:35',
    '2024-04-08 22:55',
    '2024-04-08 23:15',
    '2024-04-08 23:35',
    '2024-04-08 23:55']

kmaw_asos_times = [datetime.strptime(ttime, '%Y-%m-%d %H:%M') \
    for ttime in kmaw_asos_times]

## Coordinates for wide area in MO/AR/TN tri-state area:
#minlat = 35.0
#maxlat = 38.0
#minlon = -91.5
#maxlon = -88.5

# Coordinates for zoom area in IN tri-state area:
minlat = 38.5
maxlat = 41.0
minlon = -87.5
maxlon = -85.0

# Coordinates for wide area in IN/OH/KY tri-state area:
#minlat = 37.0
#maxlat = 41.0
#minlon = -87.0
#maxlon = -83.0

# Coordinates for area around points in IN/OH/KY tri-state area:
#minlat = 38.5
#maxlat = 39.5
#minlon = -86.0
#maxlon = -84.2

GOES_dict_reg = read_GOES_time_series_auto_regional(begin_date, end_date, \
        channels = [2, 13], save_dir = './', \
        sat = 'goes16', minlat = minlat, maxlat = maxlat, \
        minlon = minlon, maxlon = maxlon, \
        min_max_use = ['min', 'max'])


GOES_dict_reg['ptype'] = 'reg'

base_dir = '/home/bsorenson/Research/GOES/'

fig = plt.figure(figsize = (8, 4))
ax = fig.add_subplot(1,1,1)
ax.plot(GOES_dict_reg['dt_dates'], GOES_dict_reg['data'][:,1,0], label = 'Region Max GOES16')
#ax.plot(kmaw_asos_times, kmaw_asos_tmps, label = 'KMAW ASOS 2-m')
ax.grid()
ax.legend()
ax.set_ylabel('Temperature [K]')
ax.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
ax.set_title('GOES-16 vs ASOS Eclipse Comparison - Indiana')
plt.show()

#write_GOES_time_series_NCDF(GOES_dict_reg, save_dir = './')

sys.exit()





date_str = '202404081700'
#date_str = '202108052130'
#date_str = '202107202100'
plot_GOES_satpy_2panel(date_str, 2, 6, \
    zoom = False, save_dir = './', sat = 'goes17', save = False)

sys.exit()





#begin_date = '202404081200'
begin_date = '202108052100'
end_date   = '202108052200'
auto_GOES_download(begin_date, end_date, 30, sat = 'goes17', channels = [2,6,13])
#auto_GOES_download(begin_date, end_date, 30, sat = 'goes17', channels = [2,6,13])
sys.exit()



#goes_file = 'goes_cross_data_asos_202107201201_202107210331.nc'
#goes_file = 'goes16_cross_data_asos_202404081201_202404082301_v2.nc'
goes_file = 'goes16_cross_data_asos_202404081201_202404082301_v3.nc'
GOES_dict1 = read_GOES_time_series_NCDF(goes_file)


plot_GOES_time_series_points(GOES_dict1, time_idx = 25, \
    ch_idx = 0, save_dir = './', save = False)

sys.exit()



# Make a file containing GOES values for the ASOS sites
# These points are for the 2024 eclipse
# 
# - Batesville, IN ASOS:      39.3445, -85.2583
# - My location, Brookeville: 39.4681, -84.9475
# - other point 1           : 39.2551, -84.4600
# - other point 2           : 38.5965, -85.2676
# - other point 3           : 39.0130, -85.7749
#
# -----------------------------------------------------
lats_asos = np.array([39.3445, 39.4681, 39.2551, 38.5965, 39.0130])
lons_asos = np.array([-85.2583, -84.9475, -84.4600, -85.2676, -85.7749])

begin_date = '202404081200'
end_date   = '202404082330'
GOES_dict1_asos = read_GOES_time_series_auto(begin_date, end_date, \
    channels = [2, 13], dlat = list(lats_asos), \
    dlon = list(lons_asos), \
    sat = 'goes16')
    ####channels = [2, 6, 13, 8, 9, 10], dlat = list(interp_lats_low), \
    ####dlon = list(interp_lons_low))

GOES_dict1_asos['ptype'] = 'asos'

base_dir = '/home/bsorenson/Research/GOES/'

write_GOES_time_series_NCDF(GOES_dict1_asos, save_dir = './')

sys.exit()




begin_dt_date = datetime(2024,4,8,17,0)
end_dt_date   = datetime(2024,4,8,23,30)

local_dt_date = begin_dt_date

#date_str = '202404082100'
while(local_dt_date <= end_dt_date):
   
    print(local_dt_date) 
    date_str = local_dt_date.strftime('%Y%m%d%H%M')
    plot_GOES_satpy_2panel(date_str, 2, 13, \
        zoom = False, save_dir = './', sat = 'goes16', save = True)

    local_dt_date = local_dt_date + timedelta(minutes = 30)


sys.exit()


begin_date = '202107201200'
end_date   = '202107210300'
auto_GOES_download(begin_date, end_date, 30, sat = 'goes17', channels = [2,6,13])
sys.exit()


sys.exit()

#goes_file = 'goes_cross_data_asos_202107201201_202107210331.nc'
#GOES_dict1 = read_GOES_time_series_NCDF(goes_file)
#
#goes_file = 'goes_cross_data_asos_202107201201_202107210331_v2.nc'
#GOES_dict2 = read_GOES_time_series_NCDF(goes_file)
#
#sys.exit()

#download_GOES_bucket('202107210000', sat = 'goes17', channels = [1,2,3,6,8,9,10,13])
#sys.exit()

plot_GOES_figure2_v2(date_str = '202107210000', \
    goes_ch1 = 'true_color', goes_ch2 = 6, goes_ch3 = 13, \
    goes_ch4 = 8, goes_ch5 = 9, goes_ch6 = 10, \
    ch_idx1 = 0, ch_idx2 = 1, ch_idx3 = 2,\
    #ttype1 = 'mid', ttype2 = 'mid', \
    #idx1 = 0, idx2 = 3, idx3 = 6,
    ttype1 = 'low', ttype2 = 'ml', \
    idx1 = 3, idx2 = 8, idx3 = 5,
    date_idx = 25, 
    show_smoke = False, composite = True, double_fig = False, \
    zoom = True, save = False)
sys.exit()

date_str = '202107222120'
plot_GOES_satpy_5panel(date_str, 2, 13, 8, 9, 10, \
    zoom = True, save_dir = './', sat = 'goes16', save = True)

sys.exit()


plot_GOES_ASOS_comp(date_str = '202107230900', \
    goes_ch1 = 2, goes_ch2 = 6, goes_ch3 = 13, \
    ch_idx1 = 0, ch_idx2 = 1, ch_idx3 = 2,\
    ttype1 = 'asos', \
    idx1 = 0, idx2 = 1, \
    date_idx = 25, 
    show_smoke = False, composite = True, double_fig = False, \
    zoom = True, save=True)

sys.exit()

##!#mapcrs = init_proj('202107222110')
##!#fig = plt.figure(figsize = (9, 5))
##!#ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
##!#ax2 = fig.add_subplot(1,2,2, projection = mapcrs)
##!#
##!#plot_GOES_satpy('202107202100', 13, sat = 'goes16', zoom=True,\
##!#    use_xy = True, vmin = 275, vmax = 320, save = False, ax = ax1)
##!#plot_GOES_satpy('202107202100', 13, sat = 'goes17', zoom=True,\
##!#    use_xy = True, vmin = 275, vmax = 320, save = False, ax = ax2)
##!#ax1.set_title('GOES16')
##!#ax2.set_title('GOES17')
##!#fig.tight_layout()
##!#fig.savefig('test_goes_comp_ch13.png', dpi = 200)
##!#sys.exit()

#begin_date = '202306122100'
#end_date   = '202306122200'

###
###goes_file = 'goes_cross_data_up_202107201201_202107210331.nc'
###GOES_dict = read_GOES_time_series_NCDF(goes_file)
###plot_GOES_time_series_mesh(GOES_dict, date_idx = 42, ch_idx1 = 0, \
###    ch_idx2 = 0)
###
###sys.exit()

##!#begin_date = datetime(2021,7,20,21)
##!#end_date   = datetime(2021,7,21,2)
##!#local_date = begin_date
##!#while(local_date <= end_date):
##!#    date_str = local_date.strftime('%Y%m%d%H%M')
##!#    plot_GOES_figure2_v2(date_str = date_str, \
##!#        goes_ch1 = 'true_color', goes_ch2 = 6, goes_ch3 = 13, \
##!#        goes_ch4 = 8, goes_ch5 = 9, goes_ch6 = 10, \
##!#        ch_idx1 = 0, ch_idx2 = 1, ch_idx3 = 2,\
##!#        ttype1 = 'mid', ttype2 = 'mid', \
##!#        #ttype1 = 'low', ttype2 = 'ml', \
##!#        idx1 = 0, idx2 = 3, idx3 = 6,
##!#        date_idx = 25, 
##!#        show_smoke = False, composite = True, double_fig = False, \
##!#        zoom = True, save = True)
##!#
##!#    local_date = local_date + timedelta(minutes = 30)
##!#
##!#sys.exit()

date_str = '202107202130'

begin_date = '202107221200'
end_date   = '202107230400'
#begin_date  = '202107201200'
#end_date    = '202107210300'
begin_date2 = '202107211200'
end_date2   = '202107220300'
#end_date   = '202107220300'
save_dir = \
    '/home/bsorenson/Research/GOES/time_series_points/points_cross_section/20210720/'

# Prep the points - mid-plume
num_points = 29
upper_lat = 40.850577
upper_lon = -121.177811
lower_lat = 40.339418
lower_lon = -120.426204

# Prep the points - upper-plume
num_points_up = 32
upper_lat_up = 40.586326
upper_lon_up = -120.059622
lower_lat_up = 41.082106
lower_lon_up = -120.795154

##!## Prep the points - mid-lower-plume
##!#num_points_ml = 28
##!#upper_lat_ml = 40.600837
##!#upper_lon_ml = -121.319833
##!#lower_lat_ml = 40.184126
##!#lower_lon_ml = -120.721026
##!###!#
##!## Prep the points - lower-plume
##!#num_points_low = 24
##!#upper_lat_low = 40.417340
##!#upper_lon_low = -121.345775
##!#lower_lat_low = 40.007783
##!#lower_lon_low = -120.794903
##!#
##!##upper_lat = 40.75052
##!##upper_lon = -121.040965
##!##lower_lat = 40.339418
##!##lower_lon = -120.426204
##!#
##!## Select the interpolated lats and lons between the end points
interp_lats = np.linspace(lower_lat, upper_lat, num_points)
interp_lons = np.linspace(lower_lon, upper_lon, num_points)
interp_lats_up = np.linspace(lower_lat_up, upper_lat_up, num_points_up)
interp_lons_up = np.linspace(lower_lon_up, upper_lon_up, num_points_up)
##!#interp_lats_ml = np.linspace(lower_lat_ml, upper_lat_ml, num_points_ml)
##!#interp_lons_ml = np.linspace(lower_lon_ml, upper_lon_ml, num_points_ml)
##!#interp_lats_low = np.linspace(lower_lat_low, upper_lat_low, num_points_low)
##!#interp_lons_low = np.linspace(lower_lon_low, upper_lon_low, num_points_low)
##!#
##!##goes_var, goes_lat, goes_lon  = \
##!##    get_GOES_data_lat_lon(date_str, list(interp_lats_ml), \
##!##    list(interp_lons_ml), 2, version = 1, verbose = True)

##!#fig = plt.figure()
##!#
##!#labelsize = 10
##!#var4, crs0, lons2, lats2, lat_lims2, lon_lims2, plabel4 = \
##!#    read_GOES_satpy('202107210000', 13)
##!#
##!#ax = fig.add_subplot(1,1,1, projection = crs0)
##!#plot_GOES_satpy(date_str, 13, ax = ax, var = var4, crs = crs0, \
##!#    lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
##!#    vmin = None, vmax = None, \
##!#    ptitle = '', plabel = plabel4, colorbar = True, labelsize = labelsize, \
##!#    zoom=True,save=False)
##!#
##!#for lon, lat in zip(interp_lons, interp_lats):
##!#    plot_point_on_map(ax, lat, lon,\
##!#        markersize = 6)
##!#
##!#for lon, lat in zip(interp_lons_up, interp_lats_up):
##!#    plot_point_on_map(ax, lat, lon,\
##!#        markersize = 6)
##!#    
##!#ax.set_extent([-122,-119,39.5,41.5,])
##!#ax.coastlines()
##!#ax.add_feature(cfeature.STATES)
##!#ax.add_feature(cfeature.RIVERS)
##!#plt.show()
##!#
##!#sys.exit()

# Make a file containing GOES values for the ASOS sites
# -----------------------------------------------------
lats_asos = np.array([40.2824, 41.4914])
lons_asos = np.array([-121.2412, -120.5644])

GOES_dict1_low = read_GOES_time_series_auto(begin_date, end_date, \
    channels = [2, 6, 13], dlat = list(lats_asos), \
    dlon = list(lons_asos), \
    sat = 'goes16')
    ####channels = [2, 6, 13, 8, 9, 10], dlat = list(interp_lats_low), \
    ####dlon = list(interp_lons_low))

##!#GOES_dict2_low = read_GOES_time_series_auto(begin_date2, end_date2, \
##!#    channels = [2, 6, 13], dlat = list(interp_lats_low), \
##!#    dlon = list(interp_lons_low))
##!#GOES_dict1_ml = read_GOES_time_series_auto(begin_date, end_date, \
##!#    channels = [2, 6, 13, 8, 9, 10], dlat = list(interp_lats_ml), \
##!#    dlon = list(interp_lons_ml))
##!#GOES_dict2_ml = read_GOES_time_series_auto(begin_date2, end_date2, \
##!#    channels = [2, 6, 13], dlat = list(interp_lats_ml), \
##!#    dlon = list(interp_lons_ml))
##!#GOES_dict1_mid = read_GOES_time_series_auto(begin_date, end_date, \
##!#    channels = [2, 6, 13], dlat = list(interp_lats), \
##!#    dlon = list(interp_lons))
##!#GOES_dict2_mid = read_GOES_time_series_auto(begin_date2, end_date2, \
##!#    channels = [2, 6, 13], dlat = list(interp_lats), \
##!#    dlon = list(interp_lons))
##!#GOES_dict1_up = read_GOES_time_series_auto(begin_date, end_date, \
##!#    channels = [2, 6, 13], dlat = list(interp_lats_up), \
##!#    dlon = list(interp_lons_up))
##!#GOES_dict2_up = read_GOES_time_series_auto(begin_date2, end_date2, \
##!#    channels = [2, 6, 13], dlat = list(interp_lats_up), \
##!#    dlon = list(interp_lons_up))

GOES_dict1_low['ptype'] = 'asos'
##!#GOES_dict2_low['ptype'] = 'low'
##!#GOES_dict1_ml['ptype'] = 'ml'
##!#GOES_dict2_ml['ptype'] = 'ml'
##!#GOES_dict1_mid['ptype'] = 'mid'
##!#GOES_dict2_mid['ptype'] = 'mid'
##!#GOES_dict1_up['ptype'] = 'up'
##!#GOES_dict2_up['ptype'] = 'up'

base_dir = '/home/bsorenson/Research/GOES/'
##!#GOES_dict0_low  = read_GOES_time_series_NCDF(base_dir + \
##!#    'goes_cross_data_low_202107131201_202107140231.nc')
##!#GOES_dict1_low  = read_GOES_time_series_NCDF(base_dir + \
##!#    'goes_cross_data_low_202107201201_202107210231.nc')
##!#GOES_dict2_low  = read_GOES_time_series_NCDF(base_dir + \
##!#    'goes_cross_data_low_202107211201_202107220231.nc')
##!#GOES_dict0_ml   = read_GOES_time_series_NCDF(base_dir + \
##!#    'goes_cross_data_ml_202107131201_202107140231.nc')
##!#GOES_dict1_ml   = read_GOES_time_series_NCDF(base_dir + \
##!#    'goes_cross_data_ml_202107201201_202107210231.nc')
##!#GOES_dict2_ml   = read_GOES_time_series_NCDF(base_dir + \
##!#    'goes_cross_data_ml_202107211201_202107220231.nc')
##!#GOES_dict1_mid  = read_GOES_time_series_NCDF(base_dir + \
##!#    'goes_cross_data_mid_202107131201_202107140231.nc')
    #'goes_cross_data_mid_202107201201_202107210231.nc')
##!#GOES_dict2_mid  = read_GOES_time_series_NCDF(base_dir + \
##!#    'goes_cross_data_mid_202107211201_202107220231.nc')
##!#GOES_dict1_up   = read_GOES_time_series_NCDF(base_dir + \
##!#    'goes_cross_data_up_202107131201_202107140231.nc')
    #'goes_cross_data_up_202107201201_202107210231.nc')
##!#GOES_dict2_up   = read_GOES_time_series_NCDF(base_dir + \
##!#    'goes_cross_data_up_202107211201_202107220231.nc')

write_GOES_time_series_NCDF(GOES_dict1_low, save_dir = './')
##!#write_GOES_time_series_NCDF(GOES_dict2_low, save_dir = './')
##!#write_GOES_time_series_NCDF(GOES_dict1_ml,  save_dir = './')
##!#write_GOES_time_series_NCDF(GOES_dict2_ml,  save_dir = './')
##!#write_GOES_time_series_NCDF(GOES_dict1_mid, save_dir = './')
##!#write_GOES_time_series_NCDF(GOES_dict2_mid, save_dir = './')
##!#write_GOES_time_series_NCDF(GOES_dict1_up,  save_dir = './')
##!#write_GOES_time_series_NCDF(GOES_dict2_up,  save_dir = './')

##!#plot_GOES_time_series_points_auto(GOES_dict, 0, \
##!#        save_dir = save_dir + 'ch2/')
##!#plot_GOES_time_series_points_auto(GOES_dict, 1, \
##!#        save_dir = save_dir + 'ch6/')
##!#plot_GOES_time_series_points_auto(GOES_dict, 2, \
##!#        save_dir = save_dir + 'ch13/')
##!#plot_GOES_time_series_points_auto(GOES_dict, 3, \
##!#        save_dir = save_dir + 'ch8/')
##!#plot_GOES_time_series_points_auto(GOES_dict, 4, \
##!#        save_dir = save_dir + 'ch9/')
##!#plot_GOES_time_series_points_auto(GOES_dict, 5, \
##!#        save_dir = save_dir + 'ch10/')
sys.exit()




begin_date = '202306122100'
end_date   = '202306122200'
var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel0 = \
    read_GOES_satpy(begin_date, 13, sat = 'goes18', zoom = False)
#plot_GOES_satpy_2panel(begin_date, 2, 13, \
#    zoom = False, save_dir = './', sat = 'goes18', save = False)
sys.exit()
plot_GOES_satpy_6panel(date_str, 'true_color', 6, 13, 8, 9, 10, \
    zoom = True, save_dir = './', sat = 'goes17', save = True)

sys.exit()

#!#begin_date = '202107131200'
#!#end_date   = '202107140300'


#

#work_data = '/home/bsorenson/Research/MODIS/obs_smoke_forcing/asos_data_20210722_4.csv'
#df = pd.read_csv(work_data)
#
#sys.exit()

#begin_date = '202107202000'
#end_date   = '202107210330'
#plot_GOES_6panel_auto(begin_date, end_date, ch1 = 'true_color', \
#    save_dir = '/home/bsorenson/Research/GOES/six_panel/goes17_sixpanel_v4/20210720/', \
#    save = True)
#
#sys.exit()

plot_GOES_figure2(save=True, add_wv_time = False)
sys.exit()
##!##plot_GOES_satpy_6panel(date_str, 'true_color', 6, 8, 9, 10, 13, \
##!##    zoom = True, save_dir = './', save = False)
##!##sys.exit()
##!#download_GOES_bucket('202107201700', sat = 'goes17', channels = [1,3])
##!#download_GOES_bucket('202107202300', sat = 'goes17', channels = [1,3])
##!#download_GOES_bucket('202107210300', sat = 'goes17', channels = [1,3])
##!#download_GOES_bucket('202107211300', sat = 'goes17', channels = [1,3])

##!#date_str = '202107210300'
##!#plot_GOES_satpy_6panel(date_str, 'true_color', 6, 13, 8, 9, 10, \
##!#    zoom = True, save_dir = './', save = False)
##!#sys.exit()

plot_GOES_time_series_channel_comp_2loc(GOES_dict0_ml, GOES_dict1_ml, \
    0, 1, 18, ch_idx3 = 2, \
    date_idx = 20, save_dir = './', save = False)

sys.exit()

plot_GOES_time_series_channel_comp(GOES_dict1_ml, 0, 1, 27, 22, \
    ch_idx3 = 2, date_idx = 15, save = False)

plot_GOES_time_series_channel_comp(GOES_dict1_ml, 0, 1, 27, 22, \
    ch_idx3 = 2, date_idx = 15, save = True)
plot_GOES_time_series_channel_comp_2loc(GOES_dict1_low, GOES_dict1_ml, \
    0, 1, 3, 8, 27, 22, ch_idx3 = 2, date_idx = 23, save = True)
plot_GOES_time_series_mesh(GOES_dict, date_idx = 26, ch_idx1 = 0, \
    ch_idx2 = 1)
plot_GOES_cross_channels(GOES_dict, time_idx = 25)

#plot_GOES_time_series(GOES_dict, save = False)
#auto_GOES_download(begin_date, end_date, 30)

sys.exit()

#date_str = ['202107210000'] 
#date_str = ['202107202126'] 
date_str = ['202107201200',\
            '202107201500',\
            '202107201800',\
            '202107202100',\
            '202107202126',\
            '202107210000', \
            '202107210300',\
            '202107210600',\
            '202107210900',\
            '202107211200',\
            '202107211500',\
            '202107211800',\
            '202107212100',\
            '202107220000',\
            '202107220300'] 
#channel = 2
##plot_GOES_satpy(date_str, channel, ax = None, zoom=True,save=False)
#plot_GOES_satpy_6panel(date_str[4], 2, 6, 13, 8, 9, 10, zoom = True, save = False)


for dstr in date_str:
#dstr = date_str[2]
    plot_GOES_satpy_6panel(dstr, 2, 5, 13, 7, 9, 10, zoom = True, save = True)
