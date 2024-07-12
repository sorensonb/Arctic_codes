#!/usr/bin/env python

"""


"""

from VIIRSLib import *

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

# = = = = = = = = =
#
# GOOD PAIR
#
# = = = = = = = = =
date_str1 = '201707070930' # GOOD AER-FREE DAY FOR BOISE AND SEATTLE
date_str2 = '201708051030' # GOOD FOR SHOWING AEROSOL OVER WA

#date_str1 = '201707091036' # GOOD: , but clouds
#date_str2 = '201709061030' # GOOD:  EXTREMELY SMOKY DAY

plot_VIIRS_twopanel(date_str1, date_str2, 'DNB', satellite = 'SNPP', \
        check_download = False, save = True, city = 'Yakima')
        #check_download = False, save = False, city = None)
        #check_download = False, save = True, city = 'TriCities')

sys.exit()
 


mapcrs = init_proj('202107232155')
plt.close('all')
fig = plt.figure(figsize = (11,5))
ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
ax2 = fig.add_subplot(1,2,2, projection = mapcrs)






param = 'DNB'
band = 'DNB'

identify_VIIRS_file(date_str1, param, viirs_dir, satellite = 'SNPP')
identify_VIIRS_file(date_str2, param, viirs_dir, satellite = 'SNPP')

dt_date_str = datetime.strptime(date_str1, '%Y%m%d%H%M')
satellite = 'SNPP'
viirs_path = '/home/bsorenson/data/VIIRS/DNB/'
filename = glob(dt_date_str.strftime(viirs_path + satellite + '/*02DNB.A%Y%j.%H%M*.nc'))[0]
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

dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
ax1.set_title(dt_date_str.strftime('%Y-%m-%d %H:%M\nSmoke-free'))

ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
ax1.add_feature(cfeature.STATES)

ax1.set_extent([-125., -113., 41., 49.5], datacrs)


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
