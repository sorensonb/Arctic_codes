#!/usr/bin/env python

"""


"""

from NAAPSLib import *

date_str = '20080422'
date_str = '20120615'
var = 'smoke_conc_sfc'
#date_str = '20080422'
#date_str = '20170816'

if(date_str == '20170816'):
    minlat = 78.
    min_smoke = 60
elif(date_str == '20080422'):
    minlat = 60.
    min_smoke = 40
elif(date_str == '20120615'):
    minlat = 60.
    min_smoke = 40

second_arrays, second_labels = plot_NAAPS_event_CERES(date_str, var, ceres_var = 'alb_clr', \
#ceres = plot_NAAPS_event_CERES(date_str, var, ceres_var = 'alb_clr', \
    satellite = 'All', minlat = minlat, vmin = None, vmax = 300, \
    #satellite = 'All', minlat = 60., vmin = None, vmax = 300, \
    vmin2 = 0.2, vmax2 = 0.7, \
    plot_log = False, min_ice = 80., min_smoke = min_smoke, max_smoke = 2e5, \
    ptitle = '', zoom = True, lat_bounds = [60, 75], lon_bounds = [270, 360], \
    save = False)
sys.exit()


plot_NAAPS_event(date_str, var, minlat = 60., vmin = None, vmax = 50, \
    plot_log = False, ptitle = '', zoom = True, save = True)
sys.exit()


base_date = datetime.strptime('2012061300','%Y%m%d%H')
for ii in range(30):
    new_date = base_date + timedelta(hours = ii * 6)
    plot_NAAPS_figure(new_date.strftime('%Y%m%d%H'), var, minlat = 60., \
        vmax = 8, plot_log = False, zoom = True, save = True)


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


