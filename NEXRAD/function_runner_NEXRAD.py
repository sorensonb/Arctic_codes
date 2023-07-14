#!/usr/bin/env python

"""


"""
import importlib, NEXRADLib
from NEXRADLib import *
import sys

date_str = '202107210000'
radar = 'KBBX'
variable = 'differential_reflectivity'
variable = 'cross_correlation_ratio'
variable = 'reflectivity'
#variable = 'composite_reflectivity'
#angle = 8
#plot_NEXRAD_ppi_figure(date_str, radar, variable, angle = angle, \
#    save_dir = './', vmin = None, \
#    mask_outside = True, zoom = True, save = False)
#sys.exit()

radar = 'KRGX'
channel = 6
plot_NEXRAD_GOES_2panel(date_str, radar, variable, channel, ax = None, \
    angle = 2, ptitle = None, plabel = None, vmin = -5, vmax = 90, \
    labelsize = 10, colorbar = True, counties = True, save_dir = './',\
    alpha = 1.0, mask_outside = True, zoom=True, save=False)
sys.exit()
radar = 'KBBX'
azimuth = 37.
plot_NEXRAD_rhi_figure(date_str, radar, variable, azimuth = azimuth, \
    save_dir = './', vmin = None, vmax = None,\
    mask_outside = True, zoom = True, save = True)


sys.exit()
begin_date = '202107201200'
end_date   = '202107211200'
sys.exit()
plot_NEXRAD_ppi_auto(begin_date, end_date, 'KBBX', 'differential_reflectivity', \
    save_dir = './', angle_idx = 4, zoom = True, save = True)
plot_NEXRAD_ppi_auto(begin_date, end_date, 'KRGX', 'differential_reflectivity', \
    save_dir = './', angle_idx = 2, zoom = True, save = True)
#auto_NEXRAD_download(begin_date, end_date, 30, 'KRGX')

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
