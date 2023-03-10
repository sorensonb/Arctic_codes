#!/usr/bin/env python

"""


"""

from NAAPSLib import *

#date_str = '20080422'
date_str = '20120615'
#date_str = '20170816'
#date_str = '20200824'
var = 'smoke_conc_sfc'
#date_str = '20080422'
#date_str = '20170816'

# beg: 2020082400
# end: 2020090800
base_date = datetime.strptime('2014072400','%Y%m%d%H')
for ii in range(80):
    new_date = base_date + timedelta(hours = ii * 6)
    plot_NAAPS_figure(new_date.strftime('%Y%m%d%H'), var, minlat = 60., \
        vmax = 50, plot_log = False, zoom = True, save = True)

sys.exit()


if(date_str == '20170816'):
    minlat = 78.
    min_smoke = 60
    max_smoke = 300
    vmax = 300
    lat_bounds = [78, 85]
    lon_bounds = [180, 340]
    min_ice = 70.
    alb_min = 0.2
elif(date_str == '20080422'):
    minlat = 60.
    min_smoke = 20
    max_smoke = 2e5
    vmax = 300
    lat_bounds = [65, 75]
    lon_bounds = [150, 210]
    min_ice = 80.
    vmin2 = 0.4
elif(date_str == '20120615'):
    minlat = 60.
    min_smoke = 15
    max_smoke = 50
    vmax = 50
    lat_bounds = [60, 75]
    lon_bounds = [305, 345]
    min_ice = 0.
    alb_min = 0.4
elif(date_str == '20200824'):
    minlat = 75.
    min_smoke = 10
    max_smoke = 150
    vmax = 150
    lat_bounds = [75, 90]
    lon_bounds = [0, 360]
    min_ice = 0.
    vmin2 = 0.2

second_arrays, second_labels, second_arrays_nosmoke = \
    plot_NAAPS_event_CERES(date_str, var, ceres_var = 'swf_clr', \
#ceres = plot_NAAPS_event_CERES(date_str, var, ceres_var = 'alb_clr', \
    satellite = 'All', minlat = minlat, vmin = None, vmax = vmax, \
    #satellite = 'All', minlat = 60., vmin = None, vmax = 300, \
    vmin2 = alb_min, vmax2 = 300, \
    plot_log = False, min_ice = min_ice, min_smoke = min_smoke, max_smoke = max_smoke, \
    ptitle = '', zoom = True, lat_bounds = lat_bounds, lon_bounds = lon_bounds, \
    save = True)
sys.exit()



plot_NAAPS_event(date_str, var, minlat = 60., vmin = None, vmax = 150, \
    plot_log = False, ptitle = '', zoom = True, save = False)
sys.exit()


# Bugged
#2013 466 501
#2014 488 500
#2015 452 460
#2016 392 477
#2017 249 505
#2018 408 488
#2019 452 452
#2020 464 464
#2021 465 484

# Fixed
#2013 466 489
#2014 488 485
#2015 452 451
#2016 392 476
#2017 249 391
#2018 408 452
#2019 452 243
#2020 464 356
#2021 465 391


for ii in range(second_arrays.shape[0]):
    u_stat, p_val = mannwhitneyu(second_arrays[ii,0], second_arrays[ii,1],\
        alternative = 'greater')
    print(ii, p_val)

#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.hist(second_arrays[4,1], label = 'After')
#ax.hist(second_arrays[4,0], label = 'Before')
#ax.legend()
plt.show()


