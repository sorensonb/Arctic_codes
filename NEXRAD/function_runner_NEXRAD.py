#!/usr/bin/env python

"""


"""
import importlib, NEXRADLib
from NEXRADLib import *
import sys

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

radar = 'KBBX'
azimuth = 37.
plot_NEXRAD_rhi_figure(date_str, radar, variable, azimuth = azimuth, \
    save_dir = './', vmin = None, vmax = None,\
    mask_outside = True, zoom = True, save = True)


sys.exit()

begin_date = '202107221430'
end_date   = '202107230330'

#plot_NEXRAD_ppi_auto(begin_date, end_date, 'KBBX', 'differential_reflectivity', \
#    save_dir = './', angle_idx = 4, zoom = True, save = True)
plot_NEXRAD_ppi_auto(begin_date, end_date, 'KRGX', 'reflectivity', \
    save_dir = './', angle_idx = 2, zoom = True, save = True)
sys.exit()

# Every 5 minutes from 202107202100 to 202107210330
auto_NEXRAD_download(begin_date, end_date, 15, 'KRGX')
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
