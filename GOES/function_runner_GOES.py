#!/usr/bin/env python

"""


"""

from GOESLib import *
import sys

begin_date = '202107201200'
end_date   = '202107210300'
#end_date   = '202107220300'
save_dir = '/home/bsorenson/Research/GOES/time_series_points/points_cross_section/'

# Prep the points - mid-plume
num_points = 29
upper_lat = 40.850577
upper_lon = -121.177811
lower_lat = 40.339418
lower_lon = -120.426204


#upper_lat = 40.75052
#upper_lon = -121.040965
#lower_lat = 40.339418
#lower_lon = -120.426204

# Select the interpolated lats and lons between the end points
interp_lats = np.linspace(lower_lat, upper_lat, num_points)
interp_lons = np.linspace(lower_lon, upper_lon, num_points)


GOES_dict = read_GOES_time_series_auto(begin_date, end_date, \
    channels = [2, 6, 13], dlat = list(interp_lats), \
    dlon = list(interp_lons))

sys.exit()

plot_GOES_time_series_mesh(GOES_dict, date_idx = 26, ch_idx1 = 0, ch_idx2 = 1)
plot_GOES_cross_channels(GOES_dict, time_idx = 25)

plot_GOES_time_series_points_auto(GOES_dict, 0, \
        save_dir = save_dir + 'ch2/')
plot_GOES_time_series_points_auto(GOES_dict, 1, \
        save_dir = save_dir + 'ch6/')
plot_GOES_time_series_points_auto(GOES_dict, 2, \
        save_dir = save_dir + 'ch13/')
plot_GOES_time_series_points_auto(GOES_dict, 3, \
        save_dir = save_dir + 'ch8/')
plot_GOES_time_series_points_auto(GOES_dict, 4, \
        save_dir = save_dir + 'ch9/')
plot_GOES_time_series_points_auto(GOES_dict, 5, \
        save_dir = save_dir + 'ch10/')
#plot_GOES_time_series(GOES_dict, save = False)
#auto_GOES_download(begin_date, end_date, 30)
#plot_GOES_6panel_auto(begin_date, end_date,\
#    save_dir = '/home/bsorenson/Research/GOES/six_panel/goes17_sixpanel_v3/', \
#    save = True)

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
