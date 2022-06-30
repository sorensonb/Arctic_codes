#!/usr/bin/env python

"""


"""
import importlib, NEXRADLib
from NEXRADLib import *
import sys

begin_date = '202107210300'
end_date   = '202107220300'
sys.exit()
plot_NEXRAD_ppi_auto(begin_date, end_date, 'KBBX', 'reflectivity', \
    save_dir = './', angle_idx = 4, zoom = True, save = True)
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
