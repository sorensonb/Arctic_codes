#!/usr/bin/env python

"""


"""

from CrISLib import *

date_str_day   = '20210722212119'
date_str_night = '20210723093719'
press = 500.
pvar = 'temp'

wv_max = 1.4
wv_min = 0.0
sfc_tmp_max = 330
sfc_tmp_min = 280
#pvar = 'sfc_temp'

smoke_lat = 41.0000
smoke_lon = -120.4676
clear_lat1 = 41.3951
clear_lon1 = -121.0405
clear_lat2 = 40.6634
clear_lon2 = -120.1207

CrIS_level_day   = readCrIS_retrieval_level(date_str_day,   press)
CrIS_level_night = readCrIS_retrieval_level(date_str_night, press)

# Read the daytime profiles
CrIS_data_day_smoke  =  readCrIS_retrieval_profile(date_str_day, smoke_lat, smoke_lon)
CrIS_data_day_clear1 =  readCrIS_retrieval_profile(date_str_day, clear_lat1, clear_lon1)
CrIS_data_day_clear2 =  readCrIS_retrieval_profile(date_str_day, clear_lat2, clear_lon2)
# Read the nighttime profiles
CrIS_data_night_smoke  =  readCrIS_retrieval_profile(date_str_night, smoke_lat, smoke_lon)
CrIS_data_night_clear1 =  readCrIS_retrieval_profile(date_str_night, clear_lat1, clear_lon1)
CrIS_data_night_clear2 =  readCrIS_retrieval_profile(date_str_night, clear_lat2, clear_lon2)

fig = plt.figure(figsize = (8, 10))
mapcrs = init_proj('202107232155')
ax1 = fig.add_subplot(3,2,1, projection = mapcrs)
ax2 = fig.add_subplot(3,2,2, projection = mapcrs)
ax3 = fig.add_subplot(3,2,3, projection = mapcrs)
ax4 = fig.add_subplot(3,2,4, projection = mapcrs)
ax5 = fig.add_subplot(3,2,5)
ax6 = fig.add_subplot(3,2,6)
#ax1 = fig.add_subplot(3,2,1)
#ax2 = fig.add_subplot(3,2,4)
#ax3 = fig.add_subplot(3,2,2)
#ax4 = fig.add_subplot(3,2,5)
#ax5 = fig.add_subplot(3,2,3)
#ax6 = fig.add_subplot(3,2,6)

# Plot the spatial data
# ---------------------
msize = 150
point_size = 15
plot_CrIS_retrieval_level(CrIS_level_day, 'sfc_temp', ax = ax1, labelsize = 12, \
    labelticksize = 10, zoom = True, show_smoke = False, vmin = sfc_tmp_min, \
    vmax = sfc_tmp_max, save = False, markersize = msize)
plot_CrIS_retrieval_level(CrIS_level_night, 'sfc_temp', ax = ax2, labelsize = 12, \
    labelticksize = 10, zoom = True, show_smoke = False, vmin = sfc_tmp_min, \
    vmax = sfc_tmp_max, save = False, markersize = msize)
plot_CrIS_retrieval_level(CrIS_level_day, 'wv', ax = ax3, labelsize = 12, \
    labelticksize = 10, zoom = True, show_smoke = False, vmin = wv_min, \
    vmax = wv_max, save = False, markersize = msize)
plot_CrIS_retrieval_level(CrIS_level_night, 'wv', ax = ax4, labelsize = 12, \
    labelticksize = 10, zoom = True, show_smoke = False, vmin = wv_min, \
    vmax = wv_max, save = False, markersize = msize)

ax1.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax1.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:blue')
ax1.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax1.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:orange')
ax1.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax1.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:green')

ax2.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax2.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:blue')
ax2.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax2.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:orange')
ax2.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax2.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:green')

ax3.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax3.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:blue')
ax3.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax3.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:orange')
ax3.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax3.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:green')

ax4.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax4.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:blue')
ax4.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax4.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:orange')
ax4.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size + 2, marker='.',
        color = 'black', transform=datacrs)
ax4.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size, marker='.',
        transform=datacrs, color = 'tab:green')

#plot_CrIS_retrieval_profile(CrIS_data_day_smoke, pvar, smoke_lat, smoke_lon, ax = ax1)
#print(CrIS_data_day_smoke['sfc_press'])
#ax1.scatter(CrIS_data_day_smoke['sfc_temp'],CrIS_data_day_smoke['sfc_press'],  color = 'tab:blue')
#plot_CrIS_retrieval_profile(CrIS_data_day_clear1, pvar, clear_lat1, clear_lon1, ax = ax1)
#ax1.scatter(CrIS_data_day_clear1['sfc_temp'],CrIS_data_day_clear1['sfc_press'],  color = 'tab:orange')
#plot_CrIS_retrieval_profile(CrIS_data_day_clear2, pvar, clear_lat2, clear_lon2, ax = ax1)
#ax1.scatter(CrIS_data_day_clear2['sfc_temp'],CrIS_data_day_clear2['sfc_press'],  color = 'tab:green')

# Plot CrIS mixing ratio / relative humidity
plot_CrIS_retrieval_profile(CrIS_data_day_smoke, 'rh', smoke_lat, smoke_lon, ax = ax5)
plot_CrIS_retrieval_profile(CrIS_data_day_clear1, 'rh', clear_lat1, clear_lon1, ax = ax5)
plot_CrIS_retrieval_profile(CrIS_data_day_clear2, 'rh', clear_lat2, clear_lon2, ax = ax5)

plot_CrIS_retrieval_profile(CrIS_data_night_smoke, 'rh', smoke_lat, smoke_lon, ax = ax6)
plot_CrIS_retrieval_profile(CrIS_data_night_clear1, 'rh', clear_lat1, clear_lon1, ax = ax6)
plot_CrIS_retrieval_profile(CrIS_data_night_clear2, 'rh', clear_lat2, clear_lon2, ax = ax6)

ax5.axhline(500, linestyle = '--', color = 'k')
ax5.axhline(500, linestyle = '--', color = 'k')

plt.show()

sys.exit()


CrIS_data = readCrIS_retrieval_level(date_str, press)
plot_CrIS_retrieval_level(CrIS_data, pvar, ax = None, labelsize = 12, \
    labelticksize = 10, zoom = True, show_smoke = False, vmin = sfc_tmp_min, \
    vmax = sfc_tmp_max, save = False, markersize = 500)

sys.exit()

date_str = '202107231030'
CrIS_data = readCrIS_granule(date_str, satellite = 'JPSS', resolution = 'NSR')
plot_CrIS_granule_test(CrIS_data, zoom = True, save = False)
