#!/usr/bin/env python

"""


"""
import importlib, NEXRADLib
from NEXRADLib import *
import sys

## Every 5 minutes from 202107202100 to 202107210330
#begin_date = '202107210215'
#end_date = '202107210230'
#auto_NEXRAD_download(begin_date, end_date, 15, 'KBBX')
#sys.exit()

radar = 'KRGX'
azimuth = 310.0      # 281.
angle_idx = 0
range_min = 100
range_max = 160
#radar = 'KBBX'
#azimuth = 29.       # 25, 37
#angle_idx = 4
#range_min = 65
#range_max = 120

date_str = '202107202100'
variable = 'cross_correlation_ratio'

"""
202107210000
- KBBX
  - azimuth = 37.
  - angle_idx = 4
  - range_min = 70
  - range_max = 120
- KRGX
  - azimuth = 281
  - angle_idx = 0
  - range_min = 120
  - range_max = 190
202107210230
- KBBX
  - azimuth = 29.
  - angle_idx = 4
  - range_min = 65
  - range_max = 120
- KRGX
  - azimuth = 287.
  - angle_idx = 0
  - range_min = 120
  - range_max = 190
202107210400
- KBBX
  - azimuth = 29.
  - angle_idx = 4
  - range_min = 65
  - range_max = 120
- KRGX
  - azimuth = 300.
  - angle_idx = 0
  - range_min = 100
  - range_max = 160


"""

#plot_NEXRAD_rhi_figure(date_str, radar, variable, azimuth = azimuth, \
#    save_dir = './', vmin = None, vmax = None,\
#    range_min = 60, \
#    mask_outside = True, zoom = True, save = False)

#np.where(  abs(NEXRAD_dict['radar'].get_azimuth(3) - azimuth) < 0.5)
#NEXRAD_dict['radar'].fields['cross_correlation_ratio']['data'].shape

# NEXRAD_dict['radar'].range['data']
# - the distance from the bin to the radar, along an azimuth
# - dimensioned to 1832
# NEXRAD_dict['radar'].azimuth['data']
# - 
# - dimensioned to 7920a
# NEXRAD_dict['radar'].sweep_number['data']
# - the index of each sweep (each spin of the radar)
# - dimensioned to 14
# NEXRAD_dict['radar'].get_azimuth(0)
# - returns the azimuths from a given sweep
# - dimensioned to 720
# - returns 720 values for the first 8 sweeps (0 - 7)
#   but returns 360 values for the last 6 (8 - 13)
# NEXRAD_dict['radar'].fields['cross_correlation_ratio']['data']
# - the actual correlation coeffieient data
# - dimensioned to 7920, 1832
#
# 1832 - each gate on the azimuth out from the radar
# 7920 - each individual azimuth for each sweep through the volume
#
#
# need an array dimensioned to
# (# sweeps, # gates)
# loop over the sweeps, find the azimuths that fall within the azimuth
# tolerance, extract 
#
# avg_rhi = np.full((NEXRAD_dict['radar'].nsweeps, NEXRAD_dict['radar'].ngates), np.nan)
# rng_rhi = np.full((NEXRAD_dict['radar'].nsweeps, NEXRAD_dict['radar'].ngates), np.nan)
# elv_rhi = np.full((NEXRAD_dict['radar'].nsweeps, NEXRAD_dict['radar'].ngates), np.nan)
#
# x_ranges = NEXRAD_dict['radar'].range['data']
# y_elevs  = 
#
# radar_height = NEXRAD_dict['radar'].altitude['data'][0] / 1e3
#
# az_idxs_swp = np.where( np.abs(NEXRAD_dict['radar'].get_azimuth(0) - check_az) < az_tol)
# - find the indices on the current sweep where the azimuths are within range
#   of the desired 
#
# for ii, sweep_slice in enumerate(NEXRAD_dict['radar'].iter_slice()):
#     start_idx = sweep_slice.start
#     end_idx   = start_idx+NEXRAD_dict['radar'].get_azimuth(ii).shape[0]
#     az_idxs_swp = np.where( np.abs(NEXRAD_dict['radar'].get_azimuth(ii) - check_az) < az_tol)[0]
#     avg_rhi[ii,:] = np.nanmean(NEXRAD_dict['radar'].fields['cross_correlation_ratio']['data'][\
#        start_idx:end_idx][az_idxs_swp,:], axis = 0)
#     elv_rhi[ii,:] = (np.nanmean(NEXRAD_dict['radar'].gate_altitude['data'][start_idx:end_idx,:][\
#        az_idxs_swp,:], axis = 0) / 1e3) - radar_height
#     rng_rhi[ii,:] = NEXRAD_dict['radar'].range['data'] / 1e3

radar = 'KRGX'
#azimuth = 310.0      # 281.
angle_idx = 0
range_min = 100
range_max = 190
#radar = 'KBBX'
#azimuth = 29.       # 25, 37
#angle_idx = 4
#range_min = 65
#range_max = 120

date_str = '202107210000'
variable = 'cross_correlation_ratio'
radar = 'KBBX'
azimuth = 38
plot_NEXRAD_GOES_12panel(date_str, 2, 6, 13, \
    36, 287, 315,
    variable = 'composite_reflectivity', \
    ax = None, ptitle = None, plabel = None, \
    vmin = -5, vmax = 90, \
    labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=True)
####plot_NEXRAD_GOES_5panel(date_str, 2, 6, 13, \
####    variable = 'composite_reflectivity', \
####    ax = None, ptitle = None, plabel = None, \
####    vmin = -5, vmax = 90, \
####    labelsize = 10, colorbar = True, counties = True, save_dir = './',\
####    alpha = 1.0, mask_outside = True, zoom=True, save=False)
#plot_NEXRAD_rhi_multipanel(date_str, radar, azimuth, \
#    angle_idx = 4, range_min = 45, range_max = 120, \
#    save = False, save_dir = './')
sys.exit()

radar = 'KBBX'
range_min = 45
range_max = 150
angle_idx = 4
#plot_NEXRAD_rhi_multipanel_auto_varyazm(date_str, radar, 15, \
#    45, delta_azm = 1, angle_idx = angle_idx, range_min = range_min, range_max = range_max, \
#    save = True)
plot_NEXRAD_rhi_multipanel_auto_varytime(radar, 32, \
    '202107210215', '202107210400', delta_time = 15, \
    angle_idx = angle_idx, range_min = range_min, range_max = range_max, \
    save_dir = './', save = True)

sys.exit()


#date_str = '202107210000'
date_str = '202107210330'
radar = 'KRGX'
variable = 'differential_reflectivity'
variable = 'cross_correlation_ratio'
variable = 'reflectivity'
variable = 'composite_reflectivity'
angle = 5 
#plot_NEXRAD_ppi_figure(date_str, radar, variable, angle = angle, \
#    save_dir = './', vmin = None, \
#    mask_outside = True, zoom = True, save = False)


radar = 'KRGX'
channel = 6
begin_date = datetime(2021,7,21,4)
#begin_date = datetime(2021,7,20,21)
end_date   = datetime(2021,7,22,0)
local_date = begin_date
while(local_date <= end_date):
    date_str = local_date.strftime('%Y%m%d%H%M')
    plot_NEXRAD_GOES_5panel(date_str, 2, 6, 13, \
        variable = 'composite_reflectivity', \
        ax = None, ptitle = None, plabel = None, \
        vmin = -5, vmax = 90, \
        labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=True)
    local_date = local_date + timedelta(minutes = 30)
sys.exit()
#plot_NEXRAD_GOES_3panel(date_str, 'KBBX', 'KRGX', variable, \
#    channel, ax = None, \
#    angle = angle, ptitle = None, plabel = None, vmin = -5, vmax = 90, \
#    labelsize = 10, colorbar = True, counties = True, save_dir = './',\
#    alpha = 1.0, mask_outside = True, zoom=False, save=False)
#plot_NEXRAD_multiradar_multiangle(date_str, 'KBBX', 'KRGX', 'reflectivity', \
#    ptitle = None, plabel = None, vmin = -5, vmax = 90, \
#    labelsize = 10, colorbar = True, counties = True, save_dir = './',\
#    alpha = 1.0, mask_outside = True, zoom=True, save=False)
plot_NEXRAD_multiangle(date_str, 'KRGX', 'reflectivity', \
    ptitle = None, plabel = None, vmin = -5, vmax = 90, \
    labelsize = 10, colorbar = True, counties = True, save_dir = './',\
    alpha = 1.0, mask_outside = True, zoom=True, save=False)
sys.exit()

#plot_NEXRAD_GOES_2panel(date_str, radar, variable, channel, ax = None, \
#    angle = 0, ptitle = None, plabel = None, vmin = -5, vmax = 90, \
#    labelsize = 10, colorbar = True, counties = True, save_dir = './',\
#    alpha = 1.0, mask_outside = True, zoom=True, save=False)

xidx = 123
yidx = 456
field = 'reflectivity'
angle = 0

#NEXRAD_dict = read_NEXRAD(date_str, radar, angle = angle)
#point_value = find_NEXRAD_value_point(NEXRAD_dict['radar'], field, xidx, yidx, angle)
#profile     = find_NEXRAD_profile_point(NEXRAD_dict['radar'], field, xidx, yidx)

"""
NEXRAD_dict = read_NEXRAD(date_str, radar, angle = angle)
minimum_sweep = np.min(NEXRAD_dict['radar'].sweep_number['data'])
for sweep in sorted(NEXRAD_dict['radar'].sweep_number['data']):
    sweep_slice = NEXRAD_dict['radar'].get_slice(sweep)
    z = NEXRAD_dict['radar'].get_field(sweep, field)
    z__dtype = z.dtype
    z_dtype = z.dtype
    lon = NEXRAD_dict['radar'].gate_longitude['data'][sweep_slice,:]
    lat = NEXRAD_dict['radar'].gate_latitude['data'][sweep_slice,:]
    time = NEXRAD_dict['radar'].time['data']
    ranges = NEXRAD_dict['radar'].range['data']
    az = NEXRAD_dict['radar'].azimuth['data'][sweep_slice]
    az_ids = np.argsort(az)
    az = az[az_ids]
    z = z[az_ids]
    lon = lon[az_ids]
    lat = lat[az_ids]
    time = time[az_ids]
    if(sweep == minimum_sweep):
            azimuth_final = az
            time_final = time
            lon_0 = copy.deepcopy(lon)
            lon_0[-1, :] = lon_0[0, :]
            lat_0 = copy.deepcopy(lat)
            lat_0[-1, :] = lat_0[0, :]
    else:
            z_interpolator = interp2d(ranges, az, z, kind="linear")
            z = z_interpolator(ranges, azimuth_final)
    if(sweep == minimum_sweep):
             z_stack = copy.deepcopy(z[np.newaxis, :, :])
    else:
            z_stack = np.concatenate([z_stack, z[np.newaxis, :, :]])
    test_compz = z_stack.max(axis=0).astype(z_dtype)

#test_compz = np.ma.masked_where(test_compz < -20, test_compz)

sweep = 0
sweep_slice = NEXRAD_dict['radar'].get_slice(sweep)
z = NEXRAD_dict['radar'].get_field(sweep, field)

lon = NEXRAD_dict['radar'].gate_longitude['data'][sweep_slice,:]
lat = NEXRAD_dict['radar'].gate_latitude['data'][sweep_slice,:]

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

ax1.pcolormesh(lon, lat, z, shading = 'auto')
ax2.pcolormesh(lon, lat, test_compz, shading = 'auto')

plt.show()
"""

sys.exit()
begin_date = datetime(2021,7,20,21)
end_date   = datetime(2021,7,21,2)
local_date = begin_date
while(local_date <= end_date):
    date_str = local_date.strftime('%Y%m%d%H%M')
    plot_NEXRAD_GOES_3panel(date_str, 'KBBX', 'KRGX', variable, \
        channel, ax = None, \
        angle = 0, ptitle = None, plabel = None, vmin = -5, vmax = 90, \
        labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=True)
    local_date = local_date + timedelta(minutes = 30)
sys.exit()


# Read in the NEXRAD data
# -----------------------
NEXRAD_dict = read_NEXRAD(date_str, radar, angle = angle)
mapcrs = init_proj(NEXRAD_dict)

# Plot the NEXRAD data
# --------------------
fig = plt.figure(figsize = (7, 3))
ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
ax2 = fig.add_subplot(1,2,2, projection = mapcrs)

plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', ax = ax1, angle = angle,\
    mask_outside = True, zoom = True, save = False) 
plot_NEXRAD_ppi(NEXRAD_dict, 'composite_reflectivity', ax = ax2, angle = angle,\
    mask_outside = True, zoom = True, save = False) 

plt.show()

sys.exit()

radar = 'KBBX'
channel = 6
begin_date = datetime(2021,7,20,21)
#end_date   = datetime(2021,7,21,2)
local_date = begin_date
#while(local_date <= end_date):
date_str = local_date.strftime('%Y%m%d%H%M')
plot_NEXRAD_GOES_4panel(date_str, radar, variable, channel, ax = None, \
    angle1 = 4, angle2 = 6, angle3 = 8, ptitle = None, plabel = None, \
    vmin = -5, vmax = 90, \
    labelsize = 10, colorbar = True, counties = True, save_dir = './',\
    alpha = 1.0, mask_outside = True, zoom=True, save=False)
#local_date = local_date + timedelta(minutes = 30)
sys.exit()

begin_date = '202107221430'
end_date   = '202107230330'

#plot_NEXRAD_ppi_auto(begin_date, end_date, 'KBBX', 'differential_reflectivity', \
#    save_dir = './', angle_idx = 4, zoom = True, save = True)
plot_NEXRAD_ppi_auto(begin_date, end_date, 'KRGX', 'reflectivity', \
    save_dir = './', angle_idx = 2, zoom = True, save = True)
sys.exit()


begin_date = '202107201200'
end_date   = '202107140300'
#begin_date  = '202107201200'
#end_date    = '202107210300'
begin_date2 = '202107211200'
end_date2   = '202107220300'
#end_date   = '202107220300'
save_dir = \
    '/home/bsorenson/Research/NEXRAD/time_series_points/points_cross_section/20210713/'

##!#plot_NEXRAD_6panel_auto(begin_date, end_date,\
##!#    save_dir = '/home/bsorenson/Research/NEXRAD/six_panel/goes17_sixpanel_v3/20210713/', \
##!#    save = True)
##!#
##!#sys.exit()
