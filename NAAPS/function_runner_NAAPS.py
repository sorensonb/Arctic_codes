#!/usr/bin/env python

"""


"""

from NAAPSLib import *

date_str = '20170816'
#date_str = '20080422'
var = 'smoke_conc_sfc'
#plot_NAAPS_event(date_str, var, minlat = 65., vmin = None, vmax = 300, \
#    plot_log = False, ptitle = '', zoom = True, save = True)
second_arrays, second_labels = plot_NAAPS_event_CERES(date_str, var, ceres_var = 'alb_clr', \
    satellite = 'All', minlat = 78., vmin = None, vmax = 300, \
    vmin2 = 0.2, vmax2 = 0.7, \
    plot_log = False, min_ice = 80., min_smoke = 60, max_smoke = 2e5, \
    ptitle = '', zoom = True, save = True)

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


#base_date = datetime.strptime('2017081600','%Y%m%d%H')
#for ii in range(30):
#    new_date = base_date + timedelta(hours = ii * 6)
#    plot_NAAPS_figure(new_date.strftime('%Y%m%d%H'), var, minlat = 65., \
#        vmax = 20, plot_log = False, zoom = True, save = True)
