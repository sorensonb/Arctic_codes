#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *
import random
#from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

#make_gif('comp_images_20180705/', 'calc_swf_comp_20180705.gif')
date_str = '201807052213'
#plot_compare_OMI_MODIS_v2(date_str, 7, \
#    omi_dtype = 'ltc3', minlat = 65., zoom = True, save = False)
#sys.exit()

"""    

# PLOT IMAGE AND HISTOGRAM FROM RAW MODIS DATA

with open(json_time_database, 'r') as fin:
    file_date_dict = json.load(fin)
modis_date = file_date_dict[date_str]['MODIS'][0]
MODIS_ch7 = read_MODIS_channel(modis_date, 7, swath = True, \
       include_cloud_mask = True)

lat_lims = [64, 80]
lon_lims = [145,200]
mask_ch7 = np.ma.masked_where( (MODIS_ch7['data'][:,:] < 0) | \
                               (MODIS_ch7['data'][:,:] > 1), \
                               MODIS_ch7['data'])
mask_ch7 = np.ma.masked_where( ((MODIS_ch7['lon'][:,:] > 0 ) & (MODIS_ch7['lon'][:,:] < lon_lims[0])) | \
                               ((MODIS_ch7['lon'][:,:] < 0 ) & (MODIS_ch7['lon'][:,:] > (360 - lon_lims[0]))) | \
                               (MODIS_ch7['lat'][:,:] < lat_lims[0]) | \
                               (MODIS_ch7['lat'][:,:] > lat_lims[1]), \
                               mask_ch7)
#                               (MODIS_ch7['lon'][:,:] > lon_lims[1]) | \
#                               (MODIS_ch7['lat'][:,:] < lat_lims[0]) | \
#                               (MODIS_ch7['lat'][:,:] > lat_lims[1]), \
#                               mask_ch7)

mask_ch7_clr  = np.ma.masked_where(MODIS_ch7['cloud_mask'] != 3, mask_ch7)
mask_ch7_nclr = np.ma.masked_where(MODIS_ch7['cloud_mask'] == 3, mask_ch7)

mask_cld_clr  = np.ma.masked_where(MODIS_ch7['cloud_mask'] != 3, MODIS_ch7['cloud_mask'])
mask_cld_nclr = np.ma.masked_where(MODIS_ch7['cloud_mask'] == 3, MODIS_ch7['cloud_mask'])

mask_cld_nclr[mask_cld_nclr != 3] = 0


fig = plt.figure(figsize = (9, 5))
ax1 = fig.add_subplot(1,2,1, projection = ccrs.NorthPolarStereo())
ax3 = fig.add_subplot(1,2,2)

# Plot the data
# -------------
ax1.pcolormesh(MODIS_ch7['lon'][:,:], MODIS_ch7['lat'][:,:], mask_ch7, \
    cmap = 'Greys_r', transform = ccrs.PlateCarree(), shading = 'auto', \
    vmin = 0, vmax = 0.4)
ax1.coastlines()
ax1.set_extent([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]], ccrs.PlateCarree())

col_dict = {0: 'tab:blue',1: 'tab:red'}
labels = np.array(['Not clear','Clear'])
ccmm = mcolors.ListedColormap([col_dict[x] for x in col_dict.keys()])
bounds = np.array([-0.5,0.5,1.5])   
norm = mcolors.BoundaryNorm(bounds, ccmm.N)
cld_mask = ax1.pcolormesh(MODIS_ch7['lon'], MODIS_ch7['lat'], \
    mask_cld_nclr, cmap = ccmm, norm = norm, \
    transform = ccrs.PlateCarree(),
    shading = 'auto', alpha = 0.3) 
cld_mask = ax1.pcolormesh(MODIS_ch7['lon'], MODIS_ch7['lat'], \
    mask_cld_clr, cmap = ccmm, norm = norm, \
    transform = ccrs.PlateCarree(),
    shading = 'auto', alpha = 0.3) 
cbar = plt.colorbar(ScalarMappable(norm = norm, cmap = ccmm), \
    ticks = [0,1,2], ax = ax1, orientation='vertical', pad = 0.03, \
    fraction = 0.052)

# Plot the histograms
# -------------------
ax3.hist(mask_ch7_nclr.flatten().compressed(), bins = 50, alpha = 0.5, label = 'Not Clear')
ax3.hist(mask_ch7_clr.flatten().compressed(), bins = 50, alpha = 0.5, label = 'Clear')
ax3.set_yscale('log')
ax3.legend()

plt.show()
sys.exit()
"""    


# Load the data
# -------------
data = h5py.File(home_dir + \
    '/Research/Arctic_compares/comp_data/colocated_subset_' + \
    date_str + '.hdf5')
#home_dir = os.environ['HOME']
#sys.path.append(home_dir + '/Research/OMI')

# Plot the data
# -------------
mask_cld = np.ma.masked_where(data['modis_cld'][:,:] == -999., data['modis_cld'][:,:])
mask_ch7 = np.ma.masked_where(data['modis_ch7'][:,:] == -999., data['modis_ch7'][:,:])

mask_uvai = np.ma.masked_invalid(data['omi_uvai_pert'][:,:])
mask_cld_allsmoke = np.ma.masked_where(mask_uvai < 0.5, mask_cld)
mask_cld = np.ma.masked_where(mask_uvai < 0.5, mask_cld)


# Get the count of misclassified smoke pixels over land
# -----------------------------------------------------
mask_cld_smoke_land = np.ma.masked_where(\
    data['nsidc_ice'][:,:] != 254., mask_cld_allsmoke)
mask_cld_smoke_land_nocld = \
    np.ma.masked_where(mask_ch7 > 0.15, mask_cld_smoke_land)
mask_cld_smoke_land_miss  = np.ma.masked_where(\
    mask_cld_smoke_land_nocld == 3, \
    mask_cld_smoke_land_nocld)

# Get the count of misclassified smoke pixels over ocean/ice
# ----------------------------------------------------------
mask_cld_smoke_ocnice = \
    np.ma.masked_where(data['nsidc_ice'][:,:] > 100., \
    mask_cld_allsmoke)
mask_cld_smoke_ocnice_nocld = \
    np.ma.masked_where(mask_ch7 > 0.04, mask_cld_smoke_ocnice)
mask_cld_smoke_ocnice_miss = np.ma.masked_where(\
    mask_cld_smoke_ocnice_nocld == 3, \
    mask_cld_smoke_ocnice_nocld)

pcnt_misclass = \
    ((mask_cld_smoke_land_miss.compressed().shape[0] + \
     mask_cld_smoke_ocnice_miss.compressed().shape[0]) / \
     mask_cld_allsmoke.compressed().shape[0]) * 100.

print("Total # smokey (UVAI > 0.5) = ", \
    mask_cld_allsmoke.compressed().shape)
print("Total # misclass over land  = ", \
    mask_cld_smoke_land_miss.compressed().shape)
print("Total # misclass over ocean = ", \
    mask_cld_smoke_ocnice_miss.compressed().shape)
print("Percent misclassified = ", np.round(pcnt_misclass, 1))

#mask_cld = np.ma.masked_where(data['nsidc_ice'][:,:] > 252., mask_cld)
#mask_cld = np.ma.masked_where(mask_ch7 > 0.04, mask_cld)
   
#mask_cld = mask_cld_smoke_ocnice_miss

cld_bad_smoke_land = np.ma.masked_where( (data['nsidc_ice'] != 254.) & \
                                     (mask_uvai < 0.5) & \
                                     (mask_cld == 3), mask_cld)

#mask_ch7 = np.ma.masked_where(mask_uvai < 1, mask_ch7)
#mask_cld = np.ma.masked_where(data['nsidc_ice'][:,:] != 254., mask_cld)
    
col_dict = {0: 'tab:blue',1: 'tab:orange',2: 'tab:green',3: 'tab:red'}
labels = np.array(['Cloudy','Prob.\nCloudy','Prob.\nClear','Clear'])
ccmm = mcolors.ListedColormap([col_dict[x] for x in col_dict.keys()])
bounds = np.array([-0.5,0.5,1.5,2.5,3.5])   
norm = mcolors.BoundaryNorm(bounds, ccmm.N)



fig = plt.figure(figsize = (5, 8))
ax1 = fig.add_subplot(2,1,1, projection = ccrs.NorthPolarStereo())
ax2 = fig.add_subplot(2,1,2, projection = ccrs.NorthPolarStereo())

# Plot the cloud mask here
# ------------------------
cld_mask = ax1.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], \
    mask_cld, cmap = ccmm, norm = norm, \
    transform = ccrs.PlateCarree(),
    shading = 'auto', alpha = 1.0) 
cbar = fig.colorbar(ScalarMappable(norm = norm, cmap = ccmm), \
    ticks = [0,1,2,3], ax = ax1, orientation='vertical', pad = 0.03, \
    fraction = 0.052)
#ax1.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], \
#    cld_bad_smoke_land, cmap = 'Greys_r', vmin = 0, vmax = 100, \
#    alpha = 0.5, transform = ccrs.PlateCarree(), shading = 'auto')
cbar.ax.set_yticklabels(labels)

# Plot the MODIS 2.1 um here
# --------------------------
cld_mask = ax2.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], \
    mask_ch7, cmap = 'Greys_r', \
    transform = ccrs.PlateCarree(),
    shading = 'auto', alpha = 1.0) 
cld_mask = ax2.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], \
    mask_cld_smoke_ocnice_miss, cmap = ccmm, norm = norm, \
    transform = ccrs.PlateCarree(),
    shading = 'auto', alpha = 0.3) 
cld_mask = ax2.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], \
    mask_cld_smoke_land_miss, cmap = ccmm, norm = norm, \
    transform = ccrs.PlateCarree(),
    shading = 'auto', alpha = 0.3) 

#ax.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], mask_cld, \
#    transform = ccrs.PlateCarree(), shading = 'auto', cmap = 'jet')
ax1.coastlines()
ax2.coastlines()
#ax.set_extent([-180,180,65,90], ccrs.PlateCarree())
ax1.set_extent([ 145,200, 64, 80], datacrs)
ax2.set_extent([ 145,200, 64, 80], datacrs)
fig.tight_layout()
plt.show()

#data.close()
