#!/usr/bin/env python

"""


"""

from VIIRSLib import *

#lat_lims = [39., 41.5]
#lon_lims = [-105.5, -102.]
#lat_lims = [39.89, 41.4453]
#lon_lims = [-104.2245, -102.5303]

# GOOD:
#lat_lims = [40.5020, 41.9963]
#lon_lims = [-105.1769, -103.0956]
lat_lims = [40.2707, 41.9963]
lon_lims = [-105.3604, -103.0956]

work_lat_lims = [lat_lims[0] - 0.1, lat_lims[1] + 0.1]
work_lon_lims = [lon_lims[0] - 0.1, lon_lims[1] + 0.1]

date_str1 = '201709031930'  # Daytime GOOD MATCHES DENVER - aer-free
date_str2 = '201709030806'  # GOOD SMOKE-FREE NIGHTTIME
date_str3 = '201709041906'  # Daytime MATCHES HALF OF STUDY AREA
date_str4 = '201709040748'  # GOOD: SMOKY NIGHTTIME

#date_str1 = '201709040806'  # NO: WEST OF CHILE

satellite = 'SNPP'
param = 'MOD'

# Make the overall figure
# -----------------------
mapcrs = ccrs.LambertConformal(central_longitude = np.mean(lon_lims))
#fig = plt.figure(figsize = (7,6))
#ax1 = fig.add_subplot(1,1,1, projection = mapcrs)

##!#var1, crs1, lons1, lats1, plabel, xx, yy = read_VIIRS_satpy(date_str, str(channel), \
##!#        satellite = satellite, zoom = True, return_xy = True, \
##!#        add_time = 5)





plt.close('all')
fig = plt.figure(figsize = (8, 6.5))
##!#grid = ImageGrid(fig, 111, \
##!#    nrows_ncols = (2, 2), \
##!#    axes_pad = 0.20, \
##!#    share_all = True, \
##!#    cbar_location = 'right', \
##!#    cbar_mode = 'edge', \
##!#    cbar_size = '7%', \
##!#    cbar_pad = 0.15)
##!#ax1  = grid[0,0]  # control April 
##!#ax2  = grid[0,1]  # control May 
##!#ax3  = grid[1,0]  # control June 
##!#ax4  = grid[1,1]  # control July 
##!#sys.exit()
#ax1 = fig.add_subplot(1,1,1, projection = crs1)
ax1 = fig.add_subplot(2,2,2, projection = mapcrs)
ax2 = fig.add_subplot(2,2,1, projection = mapcrs)
ax3 = fig.add_subplot(2,2,4, projection = mapcrs)
ax4 = fig.add_subplot(2,2,3, projection = mapcrs)

# = = = = = = = = = = = = = = = = = 
#
# PLOT DAYTIME 1
#
# = = = = = = = = = = = = = = = = = 

##!#mask_lat = np.ma.masked_where((VIIRS_data['lat'] < lat_lims[0]) | \
##!#                              (VIIRS_data['lat'] > lat_lims[1]) | \
##!#                              (VIIRS_data['lon'] < lon_lims[0]) | \
##!#                              (VIIRS_data['lon'] > lon_lims[1]), \
##!#                              VIIRS_data['lat'])
##!#
##!#keep_long_idx = np.array([False in mask_lat.mask[ii,:] for ii in range(mask_lat.shape[0])])
##!#keep_short_idx = np.array([False in mask_lat.mask[:,ii] for ii in range(mask_lat.shape[1])])
##!#
##!#test_VIIRS = {}
##!#test_VIIRS['lat']  = VIIRS_data['lat'][keep_long_idx, :][:,keep_short_idx]
##!#test_VIIRS['lon']  = VIIRS_data['lon'][keep_long_idx, :][:,keep_short_idx]
##!#test_VIIRS['data'] = VIIRS_data['data'][keep_long_idx, :, :][:, keep_short_idx, :]
##!#test_VIIRS['cmap'] = VIIRS_data['cmap']
##!#test_VIIRS['vmax'] = VIIRS_data['vmax']
##!#test_VIIRS['dtype'] = VIIRS_data['dtype']
##!#test_VIIRS['satellite'] = VIIRS_data['satellite']
##!#test_VIIRS['filename'] = VIIRS_data['filename']
##!##test_VIIRS['col

print("READING DAYTIME FOR", date_str1)

channel = 'true_color'
dt_date_str = datetime.strptime(date_str1, '%Y%m%d%H%M')
param = 'MOD'
satellite = 'SNPP'
viirs_path = '/home/bsorenson/data/VIIRS/DNB/'
filename = glob(dt_date_str.strftime(viirs_path + satellite + '/*02' + param + '.A%Y%j.%H%M*.nc'))[0]

VIIRS_data = read_VIIRS_channel(date_str1, channel, zoom = False, swath = False, \
        satellite = 'SNPP', add_time = 5, lat_lims = work_lat_lims, lon_lims = work_lon_lims)

plot_VIIRS_granule(VIIRS_data, ax = ax1, labelsize = 12, \
    labelticksize = 10, zoom = False, show_smoke = False, vmin = None, \
    vmax = None, save = False, colorbar = False, show_title = False)

ax1.set_title(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
ax1.set_extent([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]], datacrs)  # DENVER AREA
ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
ax1.add_feature(cfeature.STATES)
plot_figure_text(ax1, 'True Color', xval = None, yval = None, \
    transform = None, color = 'red', fontsize = 11, \
    backgroundcolor = 'white', halign = 'right')

del(VIIRS_data)

#fig.tight_layout()

#outname = 'viirs_true_color_' + date_str1 + '.png'
#fig.savefig(outname, dpi = 200)
#print("Saved image", outname)

#plt.show()
#sys.exit()

# = = = = = = = = = = = = = = = = = 
#
# PLOT NIGHTTIME 1
#
# = = = = = = = = = = = = = = = = = 

print("READING NIGHTIME FOR", date_str2)
dt_date_str = datetime.strptime(date_str2, '%Y%m%d%H%M')
channel = 'DNB' # not actually necessary
param = 'DNB'
satellite = 'SNPP'
viirs_path = '/home/bsorenson/data/VIIRS/DNB/'
filename = glob(dt_date_str.strftime(viirs_path + satellite + '/*02' + param + '.A%Y%j.%H%M*.nc'))[0]
print("FILENAME")
VIIRS_data = read_VIIRS_channel(date_str2, channel, zoom = False, swath = False, \
        satellite = 'SNPP', add_time = 5, lat_lims = work_lat_lims, lon_lims = work_lon_lims)

plot_VIIRS_granule(VIIRS_data, ax = ax2, labelsize = 12, \
    labelticksize = 10, zoom = False, show_smoke = False, vmin = 0, \
    vmax = 15, save = False, colorbar = True, show_title = True)

ax2.set_title(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
ax2.set_extent([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]], datacrs)  # DENVER AREA
ax2.coastlines()
ax2.add_feature(cfeature.BORDERS)
ax2.add_feature(cfeature.STATES)
plot_figure_text(ax2, 'Day/Night Band', xval = None, yval = None, \
    transform = None, color = 'red', fontsize = 11, \
    backgroundcolor = 'white', halign = 'right')

del(VIIRS_data)

#outname = 'viirs_dnb_' + date_str2 + '.png'
#fig.savefig(outname, dpi = 200)
#print("Saved image", outname)

#plt.show()
#sys.exit()

# = = = = = = = = = = = = = = = = = 
#
# PLOT DAYTIME 2
#
# = = = = = = = = = = = = = = = = = 

channel = 'true_color'
dt_date_str = datetime.strptime(date_str3, '%Y%m%d%H%M')
param = 'MOD'
satellite = 'SNPP'
viirs_path = '/home/bsorenson/data/VIIRS/DNB/'
filename = glob(dt_date_str.strftime(viirs_path + satellite + '/*02' + param + '.A%Y%j.%H%M*.nc'))[0]

VIIRS_data = read_VIIRS_channel(date_str3, channel, zoom = False, swath = False, \
        satellite = 'SNPP', add_time = 12, lat_lims = work_lat_lims, lon_lims = work_lon_lims)

plot_VIIRS_granule(VIIRS_data, ax = ax3, labelsize = 12, \
    labelticksize = 10, zoom = False, show_smoke = False, vmin = None, \
    vmax = None, save = False, colorbar = False, show_title = False)

ax3.set_title(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
ax3.set_extent([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]], datacrs)  # DENVER AREA
ax3.coastlines()
ax3.add_feature(cfeature.BORDERS)
ax3.add_feature(cfeature.STATES)
plot_figure_text(ax3, 'True Color', xval = None, yval = None, \
    transform = None, color = 'red', fontsize = 11, \
    backgroundcolor = 'white', halign = 'right')

del(VIIRS_data)

#outname = 'viirs_true_color_' + date_str3 + '.png'
#fig.savefig(outname, dpi = 200)
#print("Saved image", outname)
#
##plt.show()
#sys.exit()

# = = = = = = = = = = = = = = = = = 
#
# PLOT NIGHTTIME 2
#
# = = = = = = = = = = = = = = = = = 

dt_date_str = datetime.strptime(date_str4, '%Y%m%d%H%M')
channel = 'DNB' # not actually necessary
param = 'DNB'
satellite = 'SNPP'
viirs_path = '/home/bsorenson/data/VIIRS/DNB/'
filename = glob(dt_date_str.strftime(viirs_path + satellite + '/*02' + param + '.A%Y%j.%H%M*.nc'))[0]

VIIRS_data = read_VIIRS_channel(date_str4, channel, zoom = False, swath = False, \
        satellite = 'SNPP', add_time = 5, lat_lims = work_lat_lims, lon_lims = work_lon_lims)

plot_VIIRS_granule(VIIRS_data, ax = ax4, labelsize = 12, \
    labelticksize = 10, zoom = False, show_smoke = False, vmin = 0, \
    vmax = 15, save = False, colorbar = True, show_title = True)

ax4.set_title(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
ax4.set_extent([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]], datacrs)  # DENVER AREA
ax4.coastlines()
ax4.add_feature(cfeature.BORDERS)
ax4.add_feature(cfeature.STATES)
plot_figure_text(ax4, 'Day/Night Band', xval = None, yval = None, \
    transform = None, color = 'red', fontsize = 11, \
    backgroundcolor = 'white', halign = 'right')

del(VIIRS_data)

plt.suptitle('Suomi-NPP VIIRS')

#outname = 'viirs_dnb_' + date_str4 + '.png'
#fig.savefig(outname, dpi = 200)
#print("Saved image", outname)

plot_subplot_label(ax2, '(a)', backgroundcolor = 'white', fontsize = 11, location = 'upper_left')
plot_subplot_label(ax1, '(b)', backgroundcolor = 'white', fontsize = 11, location = 'upper_left')
plot_subplot_label(ax4, '(c)', backgroundcolor = 'white', fontsize = 11, location = 'upper_left')
plot_subplot_label(ax3, '(d)', backgroundcolor = 'white', fontsize = 11, location = 'upper_left')




fig.tight_layout()
outname = 'viirs_4panel_denver_area_zoom_v3.png'
fig.savefig(outname, dpi = 200)
print('Saved image', outname)
plt.show()
sys.exit()


##!#plot_VIIRS_satpy(date_str, 'true_color', ax = ax1, var = var1, crs = crs1, \
##!##plot_VIIRS_satpy(date_str, 'true_color', ax = ax1, var = var1, crs = mapcrs, \
##!#    lons = lons1, lats = lats1, lat_lims = lat_lims, lon_lims = lon_lims, \
##!#    ptitle = '', plabel = plabel, \
##!#    labelsize = 10, zoom=True, save=False)


sys.exit()

#
# Nighttime
#   22 July (203): 1000
#   23 July (204): 0942
# Daytime
#   22 July (203): 2124
#   23 July (204): None yet
#

#date_str = '201709060718' NO. Panama
#date_str = '201709060812' #NO. Sri Lanka BAD DATA
#date_str = '201709061024' # CLOSE: Alberta to Nunavut
#date_str = '201708221030' # BAD
#date_str = '201708221006' # VERY CLOSE
#date_str = '201708221012' # GOOE
#date_str = '201708020942' # MATCHES, BUT BAD
#date_str = '201708151036' # CLOSE
#date_str = '201709111030' # CLOSE

#date_str = '201709021030' # BAD
#date_str = '201709021036' # BAD
#date_str = '201709031036' # BAD
#date_str = '201709031042' # BAD
#date_str = '201708291030' # BAD
#date_str = '201708291030' # BAD
#date_str = '201708271012' # BAD
#date_str = '201708271018' # MATCHES, BUT BAD IMAGE
#date_str = '201708211024' # BAD
#date_str = '201708211030' # GOOD? for Medfor OR showing aerosol
#date_str = '201708051030' # GOOD FOR SHOWING AEROSOL OVER WA


#date_str = '201709111036' # GOOD AER-FREE DAY FOR PACIFIC NW
#date_str = '201708151042' # GOOD AER-FREE DAY FOR BOISE AND SEATTLE
date_str = '201707091036' # GOOD: , but clouds
#date_str = '201707081048' # BAD
#date_str = '201707070930' # GOOD AER-FREE DAY FOR BOISE AND SEATTLE

# = = = = = = = = = = = = =
#
# GOOD PAIR for Pacific NW
#
# = = = = = = = = = = = = =
date_str1 = '201707070930' # GOOD AER-FREE DAY FOR BOISE AND SEATTLE
date_str2 = '201708051030' # GOOD FOR SHOWING AEROSOL OVER WA

#date_str1 = '201707091036' # GOOD: , but clouds
#date_str2 = '201709061030' # GOOD:  EXTREMELY SMOKY DAY

# = = = = = = = = = = = = = 
#
# Looking at Denver area
#
# = = = = = = = = = = = = = 

#date_str1 = '201709030818'  # NO: GALAPAGOS
#date_str1 = '201709030800'  # CLOSE: a bit far north
date_str1 = '201709030806'  # GOOD SMOKE-FREE NIGHTTIME
#date_str1 = '201709040806'  # NO: WEST OF CHILE
date_str1 = '201709040748'  # GOOD: SMOKY NIGHTTIME

#date_str1 = '201709031948'  # NO: ARCTIC
#date_str1 = '201709031924'  # ALMOST: ONE GRANULE TOO FAR SOUTH
date_str1 = '201709031930'  # GOOD MATCHES DENVER - aer-free

date_str1 = '201709041906'  # MATCHES HALF OF STUDY AREA
date_str1 = '201709041912'  # MATCHES HALF OF STUDY AREA

#plot_VIIRS_twopanel(date_str1, date_str2, 'DNB', satellite = 'SNPP', \
#        check_download = False, save = True, city = 'Yakima')
#        #check_download = False, save = False, city = None)
#        #check_download = False, save = True, city = 'TriCities')
#
#sys.exit()

lat_lims = [39., 41.5]
lon_lims = [-105.5, -102.]

#mapcrs = init_proj(None)
mapcrs = ccrs.LambertConformal(central_longitude = np.mean(lon_lims))
plt.close('all')
fig = plt.figure(figsize = (7,6))
ax1 = fig.add_subplot(1,1,1, projection = mapcrs)
#ax2 = fig.add_subplot(1,2,2, projection = mapcrs)


param = 'MOD'
#band = 'I01'
band = 'M05'

identify_VIIRS_file(date_str1, param, viirs_dir, satellite = 'SNPP')
#identify_VIIRS_file(date_str2, param, viirs_dir, satellite = 'SNPP')

dt_date_str = datetime.strptime(date_str1, '%Y%m%d%H%M')
satellite = 'SNPP'
viirs_path = '/home/bsorenson/data/VIIRS/DNB/'
filename = glob(dt_date_str.strftime(viirs_path + satellite + '/*02' + param + '.A%Y%j.%H%M*.nc'))[0]
#filename = glob(dt_date_str.strftime(viirs_path + satellite + '/*02MOD.A%Y%j.%H%M*.nc'))

print(dt_date_str)
print(filename)
#plot_VIIRS_figure(filename[0], band = 'DNB', vmax = None, zoom = False)

zoom = False
show_smoke = False
vmin = None
vmax = None


# Read the data for this granule
# ------------------------------
VIIRS_data = readVIIRS_granule(filename,\
    band = band, zoom = zoom)

# Plot the data for this granule
# ------------------------------
plot_VIIRS_granule(VIIRS_data, ax = ax1, zoom = zoom, \
        vmin = vmin, vmax = vmax, show_smoke = show_smoke)

#dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
ax1.set_title(dt_date_str.strftime('%Y-%m-%d %H:%M\nSmoke-free'))

ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
ax1.add_feature(cfeature.STATES)

#ax1.set_extent([-125., -113., 41., 49.5], datacrs)  # PACIFIC NW
#ax1.set_extent([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]], datacrs)  # DENVER AREA

fig.tight_layout()
plt.show()

"""


#date_str = '201709062012'
#date_str = '201709061024'
#date_str = '201709062012'
#date_str = '202107222124'
param = 'DNB'
identify_VIIRS_file(date_str2, param, viirs_dir, satellite = 'SNPP')
dt_date_str = datetime.strptime(date_str2, '%Y%m%d%H%M')
satellite = 'SNPP'
viirs_path = '/home/bsorenson/data/VIIRS/DNB/'
filename = glob(dt_date_str.strftime(viirs_path + satellite + '/*02DNB.A%Y%j.%H%M*.nc'))[0]
#filename = glob(dt_date_str.strftime(viirs_path + satellite + '/*02MOD.A%Y%j.%H%M*.nc'))

print(dt_date_str)
print(filename)
#plot_VIIRS_figure(filename[0], band = 'DNB', vmax = None, zoom = False)

band = 'DNB'
zoom = False
show_smoke = False
vmin = None
vmax = None


# Read the data for this granule
# ------------------------------
VIIRS_data = readVIIRS_granule(filename,\
    band = band, zoom = zoom)

# Plot the data for this granule
# ------------------------------
plot_VIIRS_granule(VIIRS_data, ax = ax2, zoom = zoom, \
        vmin = vmin, vmax = vmax, show_smoke = show_smoke)

dt_date_str = datetime.strptime(date_str2, '%Y%m%d%H%M')
ax2.set_title(dt_date_str.strftime('%Y-%m-%d %H:%M\nSmoky'))

ax2.coastlines()
ax2.add_feature(cfeature.BORDERS)
ax2.add_feature(cfeature.STATES)

ax2.set_extent([-125., -113., 41., 49.5], datacrs)

plt.suptitle('Suomi-NPP VIIRS DNB')
plot_subplot_label(ax1, '(a)', backgroundcolor = 'white')
plot_subplot_label(ax2, '(b)', backgroundcolor = 'white')

fig.tight_layout()




plt.show()
"""

sys.exit()

filename = glob('/home/bsorenson/data/VIIRS/DNB/JPSS/VNP02MOD.A2017249.1030*.nc')
plot_VIIRS_figure(filename, band = 'DNB', zoom = True)
#filename = glob('/home/bsorenson/data/VIIRS/DNB/JPSS/VJ102MOD.A2021204.2012*.nc')
#filename = glob('/home/bsorenson/data/VIIRS/DNB/JPSS/VJ102MOD.A2021204.1030*.nc')
#filename = glob('/home/bsorenson/data/VIIRS/DNB/JPSS/VJ102MOD.A2021204.0848*.nc')
filename = glob('/home/bsorenson/data/VIIRS/DNB/JPSS/VJ102DNB.A2021204.1030*.nc')
#filename = glob('/home/bsorenson/data/VIIRS/DNB/JPSS/VJ102MOD.A2021204.2154*.nc')
plot_VIIRS_sixpanel(satellite = 'SNPP', save = True)
sys.exit()
#plot_VIIRS_ninepanel(save = False)

##!##filename = glob('/home/bsorenson/data/VIIRS/DNB/VNP02MOD.A2021204.0942*.nc')
##!##filename = glob('/home/bsorenson/data/VIIRS/DNB/VNP02DNB.A2021204*')
##!##filename = glob('/home/bsorenson/data/VIIRS/DNB/VNP46A1.A2021204*')
##!##plot_VIIRS_DNB(filename, band = 'M05')
##!#plot_VIIRS_figure(filename, band = 'M11', zoom = True)
