#!/usr/bin/env python

"""


"""

from CrISLib import *

#date_str   = '20210722212119'
#date_str   = '20210723093719'
date_strs = ['20210720201735','20210722212119','20210723093719']
#row_str = 'mu'
row_strs = ['low','ml','md','mu']

#smoke_lat = 40.5672
#smoke_lon = -120.9731
#CrIS_prof = readCrIS_retrieval_profile(date_strs[1], smoke_lat, smoke_lon)
#sys.exit()


plot_CrIS_retrieval_combined(date_strs[1], press = 500., pvar = 'wv',\
    row_str = row_strs[1], plot_skin_temp = True, alpha = 1.0, dot_size = 300, save = False)
sys.exit()

for date_str in date_strs:
    for row_str in row_strs:
        plot_CrIS_retrieval_combined(date_str, press = 500., pvar = 'wv',\
            row_str = row_str, plot_skin_temp = False, alpha = 1.0, dot_size = 300, save = True)

sys.exit()

date_str_night = '20210723093719'
press = 500.
pvar = 'wv'

wv_max = 1.4
wv_min = 0.0
sfc_tmp_max = 330
sfc_tmp_min = 280
#pvar = 'sfc_temp'

# Middle Up
#smoke_lat = 41.2328
#smoke_lon = -120.1568
#clear_lat1 = 41.4879
#clear_lon1 = -120.5508
#clear_lat2 = 40.9871
#clear_lon2 = -119.7220

## Middle
smoke_lat = 41.0000
smoke_lon = -120.4676
clear_lat1 = 41.3951
clear_lon1 = -121.0405
clear_lat2 = 40.6634
clear_lon2 = -120.1207

# Middle Low
#smoke_lat = 40.5672
#smoke_lon = -120.9731
#clear_lat1 = 40.9128
#clear_lon1 = -121.3236
#clear_lat2 = 40.4173
#clear_lon2 = -120.5538

# Low
##smoke_lat = 40.1929
##smoke_lon = -121.2456
#smoke_lat = 40.3261
#smoke_lon = -121.0302
#clear_lat1 = 40.6438
#clear_lon1 = -121.5822
#clear_lat2 = 39.9617
#clear_lon2 = -120.4403

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

fig = plt.figure(figsize = (12, 6))
mapcrs = init_proj('202107232155')
ax1 = fig.add_subplot(2,4,1, projection = mapcrs)
ax2 = fig.add_subplot(2,4,5, projection = mapcrs)
ax3 = fig.add_subplot(2,4,2, projection = mapcrs)
ax4 = fig.add_subplot(2,4,6, projection = mapcrs)
ax5 = fig.add_subplot(2,4,3)
ax6 = fig.add_subplot(2,4,7)
ax7 = fig.add_subplot(2,4,4)
ax8 = fig.add_subplot(2,4,8)
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
#spatial_var = 'temp'
spatial_var = pvar
plot_CrIS_retrieval_level(CrIS_level_day, 'sfc_temp', ax = ax1, labelsize = 12, \
    labelticksize = 10, zoom = True, show_smoke = False, vmin = sfc_tmp_min, \
    vmax = sfc_tmp_max, save = False, markersize = msize)
plot_CrIS_retrieval_level(CrIS_level_night, 'sfc_temp', ax = ax2, labelsize = 12, \
    labelticksize = 10, zoom = True, show_smoke = False, vmin = sfc_tmp_min, \
    vmax = sfc_tmp_max, save = False, markersize = msize)
plot_CrIS_retrieval_level(CrIS_level_day, spatial_var, ax = ax3, labelsize = 12, \
    labelticksize = 10, zoom = True, show_smoke = False, vmin = wv_min, \
    vmax = wv_max, save = False, markersize = msize)
plot_CrIS_retrieval_level(CrIS_level_night, spatial_var, ax = ax4, labelsize = 12, \
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

# Plot CrIS mixing ratio / relative humidity
pvar = 'wv'
plot_CrIS_retrieval_profile(CrIS_data_day_smoke, pvar, smoke_lat, smoke_lon, ax = ax5)
plot_CrIS_retrieval_profile(CrIS_data_day_clear1, pvar, clear_lat1, clear_lon1, ax = ax5)
plot_CrIS_retrieval_profile(CrIS_data_day_clear2, pvar, clear_lat2, clear_lon2, ax = ax5)

plot_CrIS_retrieval_profile(CrIS_data_night_smoke, pvar, smoke_lat, smoke_lon, ax = ax6)
plot_CrIS_retrieval_profile(CrIS_data_night_clear1, pvar, clear_lat1, clear_lon1, ax = ax6)
plot_CrIS_retrieval_profile(CrIS_data_night_clear2, pvar, clear_lat2, clear_lon2, ax = ax6)

pvar = 'temp'
plot_CrIS_retrieval_profile(CrIS_data_day_smoke, pvar, smoke_lat, smoke_lon, ax = ax7)
plot_CrIS_retrieval_profile(CrIS_data_day_clear1, pvar, clear_lat1, clear_lon1, ax = ax7)
plot_CrIS_retrieval_profile(CrIS_data_day_clear2, pvar, clear_lat2, clear_lon2, ax = ax7)

plot_CrIS_retrieval_profile(CrIS_data_night_smoke, pvar, smoke_lat, smoke_lon, ax = ax8)
plot_CrIS_retrieval_profile(CrIS_data_night_clear1, pvar, clear_lat1, clear_lon1, ax = ax8)
plot_CrIS_retrieval_profile(CrIS_data_night_clear2, pvar, clear_lat2, clear_lon2, ax = ax8)

#ax7.scatter(CrIS_data_day_smoke['sfc_temp'],CrIS_data_day_smoke['sfc_press'],  color = 'tab:blue')
#ax7.scatter(CrIS_data_day_clear1['sfc_temp'],CrIS_data_day_clear1['sfc_press'],  color = 'tab:orange')
#ax7.scatter(CrIS_data_day_clear2['sfc_temp'],CrIS_data_day_clear2['sfc_press'],  color = 'tab:green')
#
#ax8.scatter(CrIS_data_night_smoke['sfc_temp'],CrIS_data_night_smoke['sfc_press'],  color = 'tab:blue')
#ax8.scatter(CrIS_data_night_clear1['sfc_temp'],CrIS_data_night_clear1['sfc_press'],  color = 'tab:orange')
#ax8.scatter(CrIS_data_night_clear2['sfc_temp'],CrIS_data_night_clear2['sfc_press'],  color = 'tab:green')

ax5.axhline(press, linestyle = '--', color = 'k')
ax6.axhline(press, linestyle = '--', color = 'k')
ax7.axhline(press, linestyle = '--', color = 'k')
ax8.axhline(press, linestyle = '--', color = 'k')

# Highlight regions where blue wv is higher than green wv with cyan
# Highlight regions where blue tp is lower than green tp with purple
higher_wv = np.where((CrIS_data_day_smoke['wv'] - CrIS_data_day_clear2['wv']) > 0.1)
lower_tmp = np.where((CrIS_data_day_clear2['temp'] - CrIS_data_day_smoke['temp']) > 0.5)

tmp_idxs = np.split(lower_tmp[0],np.where(np.squeeze(np.diff(lower_tmp)) != 1)[0]+1)

ax5.axhspan(CrIS_data_day_smoke['press'][higher_wv[0][0]], CrIS_data_day_smoke['press'][higher_wv[0][-1]], color = 'cyan', alpha = 0.5)
ax7.axhspan(CrIS_data_day_smoke['press'][higher_wv[0][0]], CrIS_data_day_smoke['press'][higher_wv[0][-1]], color = 'cyan', alpha = 0.5)

for ii in range(len(tmp_idxs)):
    ax5.axhspan(CrIS_data_day_smoke['press'][tmp_idxs[ii][0]], CrIS_data_day_smoke['press'][tmp_idxs[ii][-1]], color = 'purple', alpha = 0.5)
    ax7.axhspan(CrIS_data_day_smoke['press'][tmp_idxs[ii][0]], CrIS_data_day_smoke['press'][tmp_idxs[ii][-1]], color = 'purple', alpha = 0.5)

ax5.set_ylim(1000, 350)
ax6.set_ylim(1000, 350)
ax7.set_ylim(1000, 350)
ax8.set_ylim(1000, 350)

ax7.set_xlim(250, 315)
ax8.set_xlim(250, 315)

fig.tight_layout()

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
