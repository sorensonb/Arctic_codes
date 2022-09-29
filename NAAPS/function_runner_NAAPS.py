#!/usr/bin/env python

"""


"""

from NAAPSLib import *

date_str = '20170816'
#date_str = '20080422'
var = 'smoke_conc_sfc'
plot_NAAPS_event(date_str, var, minlat = 65., vmin = None, vmax = 300, \
    plot_log = False, ptitle = '', zoom = True, save = False)
#base_date = datetime.strptime('2017081600','%Y%m%d%H')
#for ii in range(30):
#    new_date = base_date + timedelta(hours = ii * 6)
#    plot_NAAPS_figure(new_date.strftime('%Y%m%d%H'), var, minlat = 65., \
#        vmax = 20, plot_log = False, zoom = True, save = True)
