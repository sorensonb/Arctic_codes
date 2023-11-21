#!/usr/bin/env python

"""


"""

from matplotlib.lines import Line2D
import OMILib
from OMILib import *
import sys

# Perform nonlinear histogram equalizing stuff
def hist_equal(data):
    # Determine the number of pixel values
    local_data = np.ma.masked_where(data < 0, data)
    local_data[(local_data.mask == True) | (local_data.data < 0)] = -9.
    local_data = local_data.data.astype(int)
    local_data = np.ma.masked_where(local_data < 0, local_data)
    work_data = local_data.compressed()
 
    pixel_vals = int(np.nanmax(work_data)) + 1 
     
    bins = np.arange(pixel_vals+1) 
    values = np.full((pixel_vals), -9.)

    # Calculate histogram
    for ii in range(len(values)):
        values[ii] = len(np.where(work_data == ii)[0])
    
    # Calculate the cumulative histogram
    cum_hist = np.cumsum(values)
  
    new_hist = ((pixel_vals - 1.)/cum_hist[-1]) * cum_hist
    new_bright = np.round(new_hist,0)

    new_data = np.copy(local_data)
    new_data = np.ma.masked_where(new_data < 0, new_data)

    # Loop over the new data and adjust using the new 
    # brightness values
    # -----------------------------------------------
    for ii in range(len(new_bright)):
        new_data[local_data == bins[ii]] = new_bright[ii]

    return new_data

infiles = [\
    '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2007m0728t0927-o16137_v003-2017m0721t040315.he5',
    '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2007m0728t1106-o16138_v003-2017m0721t040329.he5',
    '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2007m0728t1245-o16139_v003-2017m0721t040251.he5'
    ]

fig = plt.figure(figsize = (8, 6))
ax1 = fig.add_subplot(2,2,1, projection = ccrs.PlateCarree())
ax2 = fig.add_subplot(2,2,2, projection = ccrs.PlateCarree())
ax3 = fig.add_subplot(2,2,3, projection = ccrs.PlateCarree())
ax4 = fig.add_subplot(2,2,4, projection = ccrs.PlateCarree())

def plot_subplot_label(ax, label, xval = None, yval = None, transform = None, \
        color = 'black', backgroundcolor = None, fontsize = 14, \
        location = 'upper_left'):

    if(location == 'upper_left'):
        y_lim = 0.90
        #y_lim = 1.03
        x_lim = 0.05
    if(location == 'upper_upper_left'):
        y_lim = 1.05
        x_lim = 0.05
    elif(location == 'lower_left'):
        y_lim = 0.05
        x_lim = 0.05
    elif(location == 'upper_right'):
        y_lim = 0.90
        x_lim = 0.80
    elif(location == 'upper_upper_right'):
        y_lim = 1.03
        x_lim = 0.90
    elif(location == 'lower_right'):
        y_lim = 0.05
        x_lim = 0.90

    if(xval is None):
        xval = ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * x_lim
    if(yval is None):
        yval = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * y_lim
    print('Xval = ',xval, 'Yval = ',yval)

    if(transform is None):
        if(backgroundcolor is None):
            ax.text(xval,yval,label, \
                color=color, weight='bold', \
                fontsize=fontsize)
        else:
            ax.text(xval,yval,label, \
                color=color, weight='bold', \
                fontsize=fontsize, backgroundcolor = backgroundcolor)
    else:
        if(backgroundcolor is None):
            ax.text(xval,yval,label, \
                color=color, weight='bold', \
                transform = transform, fontsize=fontsize)
        else:
            ax.text(xval,yval,label, \
                color=color, weight='bold', \
                transform = transform, fontsize=fontsize, \
                backgroundcolor = backgroundcolor)

def plot_swath_on_map(infile, ax1, ax2, ax3, ax4, minlon = -60, maxlon = 50, \
        minlat = -20, maxlat = 40, equalize = True):

    data = h5py.File(infile)
   
    allrad =  data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/NormRadiance'][:,:,:]
    allrad = np.ma.masked_where(allrad <= 0, allrad)
    LAT    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
    LON    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
    AI     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
    data.close()

    rad354 = np.ma.masked_where( (LON < minlon) | (LON > maxlon)  | \
        (LAT < minlat) | (LAT > maxlat), allrad[:,:,0])
    rad388 = np.ma.masked_where( (LON < minlon) | (LON > maxlon)  | \
        (LAT < minlat) | (LAT > maxlat), allrad[:,:,1])
    rad500 = np.ma.masked_where( (LON < minlon) | (LON > maxlon)  | \
        (LAT < minlat) | (LAT > maxlat), allrad[:,:,2])
    AI     = np.ma.masked_where( (LON < minlon) | (LON > maxlon)  | \
        (LAT < minlat) | (LAT > maxlat), AI)
    AI = np.ma.masked_where(AI < -20, AI) 
    print(np.min(rad354), np.max(rad354))
    print(np.min(rad388), np.max(rad388))
    print(np.min(rad500), np.max(rad500))

    allmaxs = np.array([np.max(rad354), np.max(rad388), np.max(rad500)])

    red   = rad500*(255./np.max(allmaxs)) 
    green = rad388*(255./np.max(allmaxs)) 
    blue  = rad354*(255./np.max(allmaxs)) 
    #red   = rad354*(255./0.8) 
    #green = rad388*(255./0.8) 
    #blue  = rad500*(255./0.8) 

    if(equalize):
        red   = hist_equal(red)
        green = hist_equal(green)
        blue  = hist_equal(blue)
    
    # Create color scales for each RGB channel
    # ----------------------------------------
    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,110,160,210,240,255])
    red     = np.interp(red,old_val,ehn_val) / 255.
    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,110,160,200,230,240])
    green   = np.interp(green,old_val,ehn_val) / 255.
    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,100,150,210,240,255])
    blue    = np.interp(blue,old_val,ehn_val) / 255.
    
    # Combine the three RGB channels into 1 3-d array
    # -----------------------------------------------
    image = np.zeros((red.shape[0],red.shape[1],3))
    image[:,:,0] = red
    image[:,:,1] = green
    image[:,:,2] = blue
    
    # Convert the color values into a format usable for plotting
    # ----------------------------------------------------------
    colortuple = tuple(np.array([image[:,:,0].flatten(), \
        image[:,:,1].flatten(), image[:,:,2].flatten()]).transpose().tolist())
        
    datacrs = ccrs.PlateCarree()
    ax1.pcolormesh(LON, LAT, rad500, transform = datacrs, shading = 'auto')
    ax1.coastlines()
    ax1.set_extent([minlon, maxlon, minlat, maxlat], datacrs)

    ax2.pcolormesh(LON, LAT, rad388, transform = datacrs, shading = 'auto')
    ax2.coastlines()
    ax2.set_extent([minlon, maxlon, minlat, maxlat], datacrs)
    
    ax3.pcolormesh(LON, LAT, rad354, transform = datacrs, shading = 'auto')
    ax3.coastlines()
    ax3.set_extent([minlon, maxlon, minlat, maxlat], datacrs)

    ax4.pcolormesh(LON,LAT,\
        image[:,:,0], color= colortuple, \
        shading='auto', transform = datacrs) 
    ax4.coastlines()
    ax4.set_extent([minlon, maxlon, minlat, maxlat], datacrs)

equalize = True
for tfile in infiles:
    plot_swath_on_map(tfile, ax1, ax2, ax3, ax4, equalize = equalize)

ax1.set_title('500 nm Norm. Radiance')
ax2.set_title('388 nm Norm. Radiance')
ax3.set_title('354 nm Norm. Radiance')
if(equalize):
    ax4.set_title('False Color Image\n(R = 500 nm; G = 388 nm; B = 354 nm)\nHistogram Equalized')
else:
    ax4.set_title('False Color Image\n(R = 500 nm; G = 388 nm; B = 354 nm)')

plot_subplot_label(ax1, '(a)', location = 'upper_left')
plot_subplot_label(ax2, '(b)', location = 'upper_left')
plot_subplot_label(ax3, '(c)', location = 'upper_left')
plot_subplot_label(ax4, '(d)', location = 'upper_left')

fig.tight_layout()
plt.show()

sys.exit()


daily_VSJ4   = h5py.File('omi_shawn_daily_2005_2020.hdf5')
daily_VJZ211 = h5py.File('omi_VJZ211_daily_2005_2020.hdf5')

plot_monthly_AI_from_daily(daily_VSJ4, 86, min_AI = -20., max_AI = 20., \
    minlat = 65., maxlat = 90.)

sys.exit()

# Extract the dates
# -----------------

# Create datetime objects for each date
# -------------------------------------
#mod_dates = np.array([datetime.strptime(str(ttime), '%Y%m%d').strftime('%Y%m')
#    for ttime in daily_VSJ4['day_values']])

mod_dates = np.array([str(ttime)[:6] for ttime in daily_VSJ4['day_values']])
unique_months = np.unique(mod_dates)

#mod_dt_dates = np.array([ dtd.strftime('%Y%m') for dtd in dt_dates])

# Figure out which dates correspond to each month
# -----------------------------------------------

# Apply masking to the data according to lat/AI here
# --------------------------------------------------

# Calculate the monthly averages
# ------------------------------
monthly_VSJ4 = np.array([ np.nanmean(daily_VSJ4['grid_AI'][\
    np.where(mod_dates == t_unique_month)[0],:,:], axis = 0) \
    for t_unique_month in unique_months])
monthly_VJZ211 = np.array([ np.nanmean(daily_VJZ211['grid_AI'][\
    np.where(mod_dates == t_unique_month)[0],:,:], axis = 0) \
    for t_unique_month in unique_months])


def plot_temp_time(data1, data2, tidx, min_ai = -2., minlat = 65.):

    mask_data1 = np.ma.masked_where(daily_VSJ4['count_AI'][tidx,:,:] == 0, \
        daily_VSJ4['grid_AI'][tidx,:,:])
    mask_data2 = np.ma.masked_where(daily_VJZ211['count_AI'][tidx,:,:] == 0, \
        daily_VJZ211['grid_AI'][tidx,:,:])

    fig = plt.figure(figsize = (7, 4))

    ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(1,2,2, projection = mapcrs)

    ax1.pcolormesh(data1['lon_values'][:], data1['lat_values'][:], mask_data1, \
        transform = datacrs, shading = 'auto', cmap = 'jet', vmin = -2.0, \
        vmax = 4.0)
    ax1.coastlines()
    ax1.set_extent([-180,180,minlat, 90], datacrs)

    ax2.pcolormesh(data2['lon_values'][:], data2['lat_values'][:], mask_data2, \
        transform = datacrs, shading = 'auto', cmap = 'jet', vmin = -2.0, \
        vmax = 4.0)
    ax2.coastlines()
    ax2.set_extent([-180,180,minlat, 90], datacrs)

    plt.suptitle(str(data1['day_values'][tidx]))

    fig.tight_layout()

    plt.show() 

sys.exit()

#date_str = '201807052213'
date_str = '201507082016'
#download_OMI_single_HDF(date_str, dest_dir = h5_data_dir)
dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

filename = glob(dt_date_str.strftime(h5_data_dir + 'OMI-Aura_L2-OMAERUV_%Ym%m%dt%H%M-*.he5*'))[0]

data = h5py.File(filename)

lat = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
lon = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
cod = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/CloudOpticalDepth'][:,:]

cod = np.ma.masked_where(abs(cod) > 100, cod)

data.close()

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection = ccrs.NorthPolarStereo())
mesh = ax.pcolormesh(lon, lat, cod, vmax = 20, transform = datacrs, shading = 'auto')
cbar = plt.colorbar(mesh, ax = ax, label = 'Cloud Optical Depth')
ax.coastlines()
ax.set_extent([-180,180,65,90], datacrs)

ax.set_boundary(circle, transform=ax.transAxes)
ax.set_title('OMI Cloud Optical Depth\n' + date_str)

fig.tight_layout()
fig.savefig('omi_cod_' + date_str + '.png', dpi = 200)

plt.show()

sys.exit()

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection = ccrs.NorthPolarStereo())


sys.exit()

date_strs = [
    '202107040908',
    '202107050812',
    '202107061213',
    '202107102003',
    '202107110735',
    '202107120818',
    '202107131537',
    '202107140806',
    '202107151207',
    '202107291854',
    '202107301123',
    '202107311345',
    '202108010931',
    '202108021332',
    '202108031733',
    '202108040326',
    '202108051900',
    '202108062122',
    '202108071529',
    '202108081116',
    '202108090702',
    '202108101421',
    '202108111822',
    '202108120912',
    '202108131134',
    '202108140859',
    '202109110921',
    '202109121501',
    '202109130909',
]
#date_strs = [
#    '202207020850',
#    '202207030933',
#    '202207041652',
#    '202207051735',
#    '202207061004',
#    '202207071226',
#    '202207080813',
#    '202207091353',
#    '202207101436',
#    '202207110844',
#    ]

date_strs = [[
#    '20210704',
#    '20210705',
##    '20210706',
##    '20210710',
    '20210711',
    '20210712',
    '20210713',
##    '20210714',
##    '20210715',
    '20210729',
    '20210730',
    '20210731',
    '20210801',
    '20210802',
    '20210803',
    '20210804',
    '20210805',
    '20210806',
    '20210807',
    '20210808',
    '20210809',
    '20210810',
    '20210811',
    '20210812',
    '20210813',
##    '20210814',
    '20210911',
    '20210912',
    '20210913',
    ],
#date_strs = [
    [
    '20220702',
    '20220703',
    '20220704',
    '20220705',
    '20220706',
    '20220707',
    '20220708',
    '20220709',
    '20220710',
    '20220711',
    ]]

bad_rows = [52,21]

total_dict = {}
for date_str_list, bad_row in zip(date_strs, bad_rows):
    for date_str in date_str_list:
        print(date_str, bad_row)
        #plotOMI_single_swath_figure(date_str, dtype = 'control',  \
        #        only_sea_ice = False, minlat = 65., skiprows = None, \
        #        lat_circles = None, save = False, zoom = False, \
        #        circle_bound = True, ax = None, \
        #        shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
        good_list = download_identify_OMI_swaths(date_str, minlat = 70., min_AI = 2.0, \
            remove_bad = False, skiprows = [bad_row], omi_dtype = 'control', screen_SZA = False)
        total_dict[date_str] = good_list
            #remove_bad = False, skiprows = [52], omi_dtype = 'control', screen_SZA = False)

sys.exit()
    

date_strs = ['20210704','20210705','20210706','20210710','20210711',\
    '20210712','20210713','20210714','20210715','20210729','20210730',\
    '20210731','20210801','20210802','20210803','20210804','20210805',\
    '20210806','20210807','20210808','20210809','20210810','20210811',\
    '20210812','20210813','20210814','20210911','20210912','20210913',\
    '20220702','20220703','20220704','20220705','20220706','20220707',\
    '20220708','20220709','20220710','20220711']
total_dict = {}
for dstr in date_strs:
    bad_dict = identify_bad_rows(dstr, redownload_data = True, \
        identify_xtrack = False)

    total_dict[dstr] = bad_dict['bad_rows']
    print(dstr, bad_dict['bad_rows']) 

sys.exit()

date_str = '20180721'
good_list = download_identify_OMI_swaths(date_str, minlat = 70., min_AI = 2.0, \
    remove_bad = False, skiprows = [52], omi_dtype = 'ltc3', screen_SZA = True)
sys.exit()

date_str = '200607242155'
date_str = '200804210805'
date_str = '201507121137'
#date_str = '201807041812'
fig = plt.figure()
ax1 = fig.add_subplot(2,2,1, projection = mapcrs)
ax2 = fig.add_subplot(2,2,2, projection = mapcrs)
ax3 = fig.add_subplot(2,2,3, projection = mapcrs)
ax4 = fig.add_subplot(2,2,4, projection = mapcrs)

vmax = 1.0
plotOMI_single_swath_figure(date_str, dtype = 'shawn',  \
        only_sea_ice = False, minlat = 70., skiprows = [52], \
        lat_circles = None, save = False, zoom = False, \
        circle_bound = True, ax = ax1, \
        vmax = vmax, \
        shawn_path = '/home/bsorenson/data/OMI/shawn_files/ltc3/')
plotOMI_single_swath_figure(date_str, dtype = 'shawn',  \
        only_sea_ice = False, minlat = 70., skiprows = [52], \
        lat_circles = None, save = False, zoom = False, \
        circle_bound = True, ax = ax2, \
        vmax = vmax, \
        shawn_path = '/home/bsorenson/data/OMI/shawn_files/ltc3_old/')
plotOMI_single_swath_figure(date_str, dtype = 'shawn',  \
        only_sea_ice = False, minlat = 70., skiprows = [52], \
        lat_circles = None, save = False, zoom = False, \
        circle_bound = True, ax = ax3, \
        vmax = vmax, \
        shawn_path = '/home/bsorenson/data/OMI/shawn_files/ltc3_new/')
plotOMI_single_swath_figure(date_str, dtype = 'shawn',  \
        only_sea_ice = False, minlat = 70., skiprows = [52], \
        lat_circles = None, save = False, zoom = False, \
        circle_bound = True, ax = ax4, \
        vmax = vmax, \
        shawn_path = '/home/bsorenson/data/OMI/shawn_files/ltc4/')

ax1.set_title('Control')
ax2.set_title('LTC3 OLD')
ax3.set_title('LTC3 NEW')
ax4.set_title('LTC4')

plt.show()

sys.exit()

date_str = '201807052213'
minlat = 65.
OMI_base  = readOMI_swath_shawn(date_str, latmin = minlat,\
    shawn_path = home_dir + '/data/OMI/shawn_files/ltc3/')
dtype = 'ltc3'
OMI_hdf   = readOMI_swath_hdf(date_str, 'control', latmin = minlat)
#write_swath_to_HDF5(OMI_base, dtype, save_path = './', minlat = 65., \
#    shawn_path = home_dir + '/data/OMI/shawn_files/ltc3/', \
#    remove_empty_scans = True)

sys.exit()

#date_str = '20180721'
#good_list = download_identify_OMI_swaths(date_str, minlat = 70., min_AI = 2.0, \
#    remove_bad = False, skiprows = [52], omi_dtype = 'ltc3', screen_SZA = True)


date_strs = [\
    #'201206141245',
    #'201206141424',
    #'201206141603',
    #'201206150653',
    #'201206150832',
    #'201206151010',
    #'201206151149',
#    '201206151328',
    #'201206151507',
    #'201206151646',
    '201807210047',
    '201807210226',
    '201807211358',
    '201807211537',
    '201807211716',
    '201807211855',
    '201807212034',
    '201807212213',

    #'202108010117',
    #'202108010256',
    #'202108010435',
    #'202108010614',
    #'202108010931',  # BAD
    #'202108011110',  # BAD
    #'202108011428',
    #'202108011607',
    #'202108011925',
    #'202108012103',
    #'202108012242',


    ]

for date_str in date_strs:
    plotOMI_single_swath_figure(date_str, dtype = 'ltc3',  \
            only_sea_ice = False, minlat = 70., skiprows = [52], \
            lat_circles = None, save = False, zoom = False, \
            circle_bound = True, ax = None, \
            shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
sys.exit()
    

date_str = '20180705'
identify_OMI_HDF_swaths(date_str, minlat = 70., min_AI = 2.0, \
    remove_bad = False, skiprows = [52])

sys.exit()


date_strs = [\
        '202108010435',\
        '202108010614',\
        '202108010752',\
        '202108010931',\
        '202108011110',\
        '202108011249',\
        '202108011428',\
        '202108011607',\
        '202108011746',\
        '202108011925',\
        '202108012103',\
        '202108012242',\
    ]

date_str = '20210801'
base_date = datetime.strptime(date_str, '%Y%m%d')
quit_date = base_date + timedelta(days = 1)
local_date_str = base_date
while(local_date_str < quit_date):
    
#for date_str in date_strs:
    #dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    print(local_date_str)
    #print(dt_date_str, dt_date_str + timedelta(minutes = 99))

    date_str = local_date_str.strftime('%Y%m%d%H%m')
    download_OMI_file_HDF(date_str, dest_dir = h5_data_dir)

    local_date_str = local_date_str + timedelta(minutes = 99)

sys.exit()

#date_str = '200607270100'
#date_str = '201908101301'
date_str = '202108010117'
plotOMI_single_swath_figure(date_str, dtype = 'control',  \
        only_sea_ice = False, minlat = 50., skiprows = [52], \
        lat_circles = None, save = False, zoom = False, \
        circle_bound = True, ax = None, \
        shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
sys.exit()

OMI_base = '201908110033'
write_shawn_to_HDF5(OMI_base, save_path = './', minlat = 65., \
    shawn_path = home_dir + '/data/OMI/shawn_files/')

sys.exit()

plot_combined_fort_out('20190811', min_lat = 70., vtype = 'areas', max_lat = 80., save = True)
plot_combined_fort_out('20170818', min_lat = 80., vtype = 'areas', save = True)
sys.exit()

version = 'jz211'
month_idx = 4
plotOMI_drift_figure(version, month_idx, plot_counts = True, \
    plot_trend = False, deseasonalize = True, save = True)
sys.exit()

plotOMI_Compare_TrendUncertainty_all('VJZ211', 'VSJ4',\
    trend_type = 'standard', minlat=65,save=True)

sys.exit()

# For JZ211 vs SJ4
lat1, lon1 = 78., 46.
lat2, lon2 = 85., 46.

## For SJ4 vs SJ42
#lat1, lon1 = 78., 44.
#lat2, lon2 = 72., 44.

#plotOMI_month_avg_comp('VSJ4', 'VSJ43', 69, lat1, lon1, lat2, lon2)
plotOMI_month_avg_comp('VJZ211', 'VSJ4', 27, lat1, lon1, lat2, lon2)
    
sys.exit()


minlat = 65.

OMI_VBS0   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VBS0_2005_2020.nc', minlat = minlat)
OMI_VJZ211 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VBS01_2005_2020.nc', minlat = minlat)
OMI_VSJ4   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', minlat = minlat)
#    ai_trends, ai_pvals, ai_uncert = calcOMI_grid_trend(OMI_data, month_idx, trend_type, \
#plotOMI_Compare_ClimoTrend_summer(OMI_VBS0,OMI_VJZ211,OMI_VSJ4,\
#        trend_type = 'standard', minlat=minlat,save=False)
plotOMI_Compare_ClimoTrend_all(OMI_VBS0,OMI_VJZ211, OMI_VSJ4,\
        trend_type = 'standard', minlat=minlat,save = False)

sys.exit()

#start_date = 200504
#end_date = 202009
#inputfile = '/home/bsorenson/Research/OMI/JZ_analysis/climo_analysis/omi_bs01_climo_2005_2020.txt'
#OMI_data = readOMI(inputfile,start_date,end_date,key=None)
#filename = '/home/bsorenson/Research/OMI/omi_ai_VBS01_2005_2020.nc'
#writeOMI_toNCDF(OMI_data,filename, minlat = 65.)
#sys.exit()

#
#inputfile = '/home/bsorenson/Research/OMI/shawn_analysis/climo_analysis/omi_vsj43_climo_2005_2020.txt'
#OMI_data = readOMI(inputfile,start_date,end_date,key=None)
#filename = '/home/bsorenson/Research/OMI/omi_ai_VSJ43_2005_2020.nc'
#writeOMI_toNCDF(OMI_data,filename, minlat = 65.)
#sys.exit()

for ii in range(0, 6):
    plotOMI_TrendUncertainty_SingleMonth('VJZ211', 'VSJ4', ii,\
         minlat = 65., save = True)


month_idx = 4
minlat = 65.
#OMI_VBS0   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VBS0_2005_2020.nc', minlat = minlat)
#OMI_VJZ211 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VJZ211_2005_2020.nc', minlat = minlat)
OMI_VSJ4   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', minlat = minlat)
OMI_VSJ42  = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VSJ42_2005_2020.nc', minlat = minlat)

ai_trends, ai_pvals, ai_uncert = calcOMI_grid_trend(OMI_data, month_idx, trend_type, minlat)

sys.exit()



# Read the JZ211 and JZ212 data
minlat = 65.
OMI_VJZ211 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
    'omi_ai_VJZ211_2005_2020.nc', minlat = minlat)
OMI_VJZ212 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
    'omi_ai_VJZ212_2005_2020.nc', minlat = minlat)

# Compare JZ211 and JZ212 August trends
fig = plt.figure()
ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
ax2 = fig.add_subplot(1,2,2, projection = mapcrs)

month_idx = 5
trend_type = 'standard'
colorbar_label_size = 10
plotOMI_MonthTrend(OMI_VJZ211,month_idx=month_idx,trend_type=trend_type,label = ' ',\
    minlat=65.,title = 'VJZ211', pax = ax1, colorbar = False, \
    colorbar_label_size = colorbar_label_size)
plotOMI_MonthTrend(OMI_VJZ212,month_idx=month_idx,trend_type=trend_type,label = ' ',\
    minlat=65.,title = 'VJZ212', pax = ax2, colorbar = False, \
    colorbar_label_size = colorbar_label_size)
fig.tight_layout()
plt.show()
sys.exit()



#date_str = '20190625'
date_str = '20190506'
#swath_ais, swath_lat, swath_xtk, swath_xrw, day_avgs, row_avgs = \
plot_OMI_row_avg(date_str, plot_swath = False, minlat = 65., save = True)

sys.exit()

##!#date_strs = [
##!#             '201904290123',
##!#             '201904290302',
##!#             '201904290440',
##!#             '201904290619',
##!#             '201904290758',
##!#             '201904290937',
##!#             '201904291116',
##!#             '201904291255',
##!#             '201904291434',
##!#             '201904291613',
##!#             '201904291752',
##!#             '201904291930',
##!#             '201904292109',
##!#             '201904292248',
##!#             #'201907170041',
##!#             #'201907170220',
##!#             #'201907170359',
##!#             #'201907170537',
##!#             #'201907170716',
##!#             #'201907170855',
##!#             #'201907171034',
##!#             #'201907171213',
##!#             #'201907171352',
##!#             #'201907171531',
##!#             #'201907171710',
##!#             #'201907171848',
##!#             #'201907172027',
##!#             #'201907172206',
##!#             #'201907172345',
##!#            ]

# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------
date_str = '200804222159'
plotOMI_single_multipanel(date_str, only_sea_ice = False, minlat = 65., \
       quad_panel = True, save = True)
#plotOMI_single_swath_multiple('22222222', dtype = 'control',  \
#    only_sea_ice = False, minlat = 65., save = True)
sys.exit()

plot_row_bias(dataset = 'normal', save = True)
sys.exit()

plotOMI_single_swath_multiple_v3(dtype = 'control',  \
        only_sea_ice = False, minlat = 65., save = True)
#plotOMI_single_swath_multiple_v4(dtype = 'control',  \
#        only_sea_ice = False, minlat = 65., save = True)

sys.exit()

plot_row_anomaly_combined(dtype = 'control', minlat = 65., save = True)
sys.exit()

minlat = 65.
OMI_VBS0   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VBS0_2005_2020.nc', minlat = minlat)
plotOMI_NCDF_Climo_SpringSummer(OMI_VBS0,start_idx=0,end_idx=96,minlat=65,\
                   save = True)
sys.exit()

min_lat = 70.
infile = home_dir + '/Research/OMI/shawn_analysis/count_analysis/' + \
    'omi_vsj4_areas_2005_2020_100.txt'

#calc_aerosol_event_dates(infile, minlat = min_lat, vtype = 'areas', \
#    max_lat = None)
#
#sys.exit()

vtype = 'areas'
omi_fort_dict = read_OMI_fort_out(infile, min_lat, vtype = vtype)
#omi_fort_dict_80 = read_OMI_fort_out(infile, min_lat, vtype = vtype)

fig1 = plt.figure(figsize = (8,5)) 
ax1 = fig1.add_subplot(1,1,1)
#ax2 = fig1.add_subplot(1,2,2)

loop_dates = omi_fort_dict['daily_dt_dates'][np.where( \
    (omi_fort_dict['daily_dt_dates'] >= datetime(2019,4,1)) & \
    (omi_fort_dict['daily_dt_dates'] <= datetime(2019,9,30)))]

plot_dates = datetime(year = 2000, month = 4, day = 1) + \
    (loop_dates - datetime(year = 2019, month = 4, day = 1))
year = 2019
ax1.plot(plot_dates,\
    omi_fort_dict['daily_data'][np.where( \
    (omi_fort_dict['daily_dt_dates'] >= datetime(year,4,1)) & \
    (omi_fort_dict['daily_dt_dates'] <= datetime(year,9,30)))], label = '2019')


# Plot the time series with the peaks as 'x's
# -------------------------------------------
#plot_OMI_fort_out_time_series(omi_fort_dict, ax1, minlat = min_lat)
#custom_lines = [Line2D([0], [0])]

print(ax1.get_xlim())

ax1.legend()
ax1.grid()
#ax1.legend(custom_lines, ['2019'], loc = 2)

exponent = 5
label_str = 'Area of high AI [10$^{' + \
    str(exponent) + '}$ km$^{2}$]'
ax1.set_ylabel(label_str, weight = 'bold')

ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax1.set_title('High AI ' + vtype + '\nPerturbed AI threshold of '+str(omi_fort_dict['ai_thresh'])+\
    '\nNorth of '+str(int(min_lat))+'$^{o}$N')
fig1.savefig('omi_daily_areas_minlat70_2019.png', dpi = 300)
plt.show()
sys.exit()







month_idx =  4
trend_type = 'standard'
inputfile = '/home/bsorenson/Research/OMI/JZ_analysis/climo_analysis/omi_vsj4_climo_2005_2020.txt'
OMI_data = readOMI_NCDF('/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', minlat = minlat)

#trends, pvals = calcOMI_grid_trend(OMI_data, month_idx, trend_type, minlat)
plotOMI_MonthTrend(OMI_data,month_idx=4,save=False,\
    trend_type='standard',minlat=65.,return_trend=False, colorbar = True, \
    title = '', label = '', colorbar_label_size = 14, pax = None)

sys.exit()

plotOMI_Compare_ClimoTrend_all(OMI_data1,OMI_data2,OMI_data3,\
    trend_type = 'standard', minlat=65,save=False)

plotOMI_Type_Trend_all(trend_type = 'standard', \
    minlat=65,save=True)

sys.exit()

start_date = 200504
end_date = 202009
inputfile = '/home/bsorenson/Research/OMI/JZ_analysis/climo_analysis/omi_jz211_climo_2005_2020.txt'
OMI_data = readOMI(inputfile,start_date,end_date,key=None)
filename = '/home/bsorenson/Research/OMI/sorenson_et_al_data/omi_ai_screened_2005_2020.nc'
writeOMI_toNCDF(OMI_data,filename, minlat = 65.)

inputfile = '/home/bsorenson/Research/OMI/shawn_analysis/climo_analysis/omi_vsj4_climo_2005_2020.txt'
OMI_data = readOMI(inputfile,start_date,end_date,key=None)
filename = '/home/bsorenson/Research/OMI/sorenson_et_al_data/omi_ai_perturbed_2005_2020.nc'
writeOMI_toNCDF(OMI_data,filename, minlat = 65.)
sys.exit()

#plotOMI_Type_Climo_all(trend_type = 'standard', \
#    minlat=65,save=True)

#minlat = 65.
#
#data_types = ['VBS1','VJZ2','VJZ28','VJZ29','VJZ211','VSJ2','VSJ4']    
#
## Read in all the data
#OMI_VBS1 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[0] + '_2005_2019.nc', minlat = minlat)
#OMI_VJZ2 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[1] + '_2005_2019.nc', minlat = minlat)
#OMI_VJZ28= readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[2] + '_2005_2019.nc', minlat = minlat)
#OMI_VJZ29= readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[3] + '_2005_2019.nc', minlat = minlat)
#OMI_VJZ211 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[4] + '_2005_2019.nc', minlat = minlat)
#OMI_VSJ2 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[5] + '_2005_2019.nc', minlat = minlat)
#OMI_VSJ4 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
#    'omi_ai_' + data_types[6] + '_2005_2020.nc', minlat = minlat)
#sys.exit()
plotOMI_Type_MonthClimo_all(minlat=65, save = True, \
        save_dir = '/home/bsorenson/Research/OMI/monthly_images/combined_climo/')

sys.exit()


minlat = 65.
base_dir = '/home/bsorenson/Research/OMI/'
image_dir = base_dir + 'monthly_images/'
data_types = ['VJZ2','VJZ28','VJZ29','VJZ211']    
for dtype in data_types:
    print('\n\nNOW PROCESSING',dtype,'\n\n')

    local_dir = image_dir + dtype + '/'

    # Change to Arctic codes directory
    os.chdir(local_dir)

    print(os.system('pwd'))

    # Read in the data
    OMI_data = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
        'omi_ai_' + dtype + '_2005_2019.nc', minlat = minlat)

    for ii in range(len(OMI_data['DATES'])):
        plotOMI_NCDF_SingleMonth(OMI_data,ii,minlat = minlat, pax = None, \
            save=True)

sys.exit()


infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_100.txt'
#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = False)
plot_OMI_fort_out_peaks(infile, minlat = 75., vtype = 'areas', save = False)

sys.exit()

#plot_row_bias(dataset = 'normal', save = True)
#sys.exit()

###minlat = 65
###plt.close('all')
###fig1 = plt.figure(figsize = (6,6))
###mapcrs = ccrs.NorthPolarStereo()
###ax0 = fig1.add_subplot(1,1,1, projection = mapcrs)
###ax0.coastlines(resolution = '50m')
###ax0.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
###ax0.set_boundary(circle, transform = ax0.transAxes)
###plot_lat_circles(ax0, [80])
###plot_arctic_regions(ax0)
###ax0.gridlines()
###fig1.tight_layout()
####fig1.savefig('arctic_lat_circles_80.png',dpi=300)
###plt.show()

# pvars = 'glint', 'earthshine', 'sunshine', 'stray'
#plotOMI_single_flags('201204101833', pvar = 'glint', minlat = 55, zoom = True)



#date_str = '201807050048'
#date_str = '201807050227'
#date_str = '201807050406'
#date_str = '201807050545'
#date_str = '201807050723'
#date_str = '201807050902'
#date_str = '201807051041'
#date_str = '201807051220'
#date_str = '201807051359'
#date_str = '201807051538'
#date_str = '201807051717'
date_str = '201807051856'
#date_str = '201807052034'
#date_str = '201807052213'
#date_str = '201807052352'
minlat = 65.
#OMI_base  = readOMI_swath_shawn(date_str, latmin = minlat)
#sys.exit()
#write_shawn_to_HDF5(OMI_base)

#OMI_data = readOMI_swath_hdf(date_str, 'control', only_sea_ice = False, \
#    only_ice = False, no_ice = False, latmin = 65, skiprows = [52])
#date_str = '201605152104'
#date_str = '201605162009'

date_strs = ['200607240029', # GOOD
             '200607240208', # GOOD
             '200607240347', # GOOD
             '200607240526', # GOOD
             '200607240705', # GOOD
             '200607240844', # GOOD
             '200607242155', # GOOD
             '200607242334', # GOOD
             '200607250112', # GOOD
             '200607250251', # GOOD
             '200607250609', # GOOD
             '200607250748', # GOOD?
             '200607252238', # GOOD
             '200607260017', # GOOD
             '200607260156', # GOOD
             '200607260335', # GOOD
             '200607260513', # GOOD?
             '200607260831', # GOOD
             '200607262142', # GOOD
             '200607270100', # GOOD
             '200607270239', # GOOD?
             '200607270418', # GOOD?
             '200607270557', # GOOD?
             '200607270736', # GOOD?
             '200607270914',
             '200607271053',
             '200607271232',
             '200607271411',
             '200607271550',
             '200607271729',
             '200607271908',
             '200607272047', 
             '200607272226'] # GOOD



#date_strs = ['200804221841',  # GOOD
#             '200804222020',  # GOOD
#             '200804222159',  # GOOD
#date_strs = ['201207260009', 
#             '201207260148', 
#             '201207260327', 
#             '201207260506', 
#             '201207260645', 
#             '201207260824', 
#             '201207261003', 
#             '201207261142', 
#             '201207261320', 
#             '201207261459', 
#             '201207261638', 
#             '201207261817', 
#             '201207261956', 
#             '201207262135', 
#             '201207262314']
#date_strs = ['201207250105', 
#             '201207250244',
#             '201207250423', 
#             '201207250602', 
#             '201207250741', 
#             '201207250920', 
#             '201207251058', 
#             '201207251237', 
#             '201207251416', 
#             '201207251555', 
#             '201207251734', 
#             '201207251913', # good? 
#             '201207252052', # good?
#             '201207252231'] # good?
#date_strs = ['201308020123',
#             '201308020302',
#             '201308020441',
#             '201308020620',
#             '201308020759',
#             '201308020938',
#             '201308021117',
#             '201308021256',
#             '201308021434',
#             '201308021613',
#             '201308021752',
#             '201308021931',
#             '201308022110',
#             '201308022249']
#             '201605151925',  # MEDIOCRE
#             '201605152104',  # MEDIOCRE
#             '201605152243',  # MEDIOCRE
#             '201605162148',  # MEDIOCRE
#             '201807051856',  # GOOD
#             '201908102115',  # GOOD
#             '201908102254',  # GOOD
#             '201908110033',  # GOOD
#             '201908110351',  # GOOD
#             '200607260017',  # GOOD
#             '200607252238',  # GOOD
#             '200607260156',  # GOOD
#             '200607260335',  # GOOD
#             '200607260513',  # GOOD
#             '201708161504',  # GOOD
#             '201708161643',  # GOOD
#             '201708161821',  # GOOD
#             '201708171408',  # GOOD
#             '201708171547',  # GOOD
#             '201708171726',  # GOOD
#             '201708171905',  # GOOD
#             '201708172043',  # GOOD
#             '201708181312',  # GOOD
#             '201708181451',  # GOOD
#             '201708181630',  # GOOD
#             '201708181809',  # GOOD
#             '201708181948',  # GOOD
#             '201708191355',  # GOOD
#             '201708191534',  # GOOD
#             '201708191713' ] # GOOD


#date_str = '200804222020' # GOOD
#date_str = '200804222159' # GOOD
#date_str = '201605151925' # MEDIOCRE
#date_str = '201605152104' # MEDIOCRE
#date_str = '201605152243' # MEDIOCRE
#date_str = '201605162148' # MEDIOCRE
#date_str = '201807051856' # GOOD
#date_str = '201908102115' # GOOD
#date_str = '201908102254' # GOOD
#date_str = '201908110033' # GOOD
#date_str = '201908110351' # GOOD
#date_str = '200607260017' # GOOD
#date_str = '200607252238' # GOOD
#date_str = '200607260156' # GOOD
#date_str = '200607260335' # GOOD
#date_str = '200607260513' # GOOD
#date_str = '201708161504' # GOOD
#date_str = '201708161643' # GOOD
#date_str = '201708161821' # GOOD
#date_str = '201708171408' # GOOD
#date_str = '201708171547' # GOOD
#date_str = '201708171726' # GOOD
#date_str = '201708171905' # GOOD
#date_str = '201708172043' # GOOD
#date_str = '201708181312' # GOOD
#date_str = '201708181451' # GOOD
#date_str = '201708181630' # GOOD
#date_str = '201708181809' # GOOD
#date_str = '201708181948' # GOOD
#date_str = '201708191355' # GOOD
#date_str = '201708191534' # GOOD
#date_str = '201708191713' # GOOD

#sys.path.append('/home/bsorenson/Research/MODIS/obs_smoke_forcing/')
#from MODISLib import *

#for date_str in date_strs[:5]:
for date_str in date_strs:

    plotOMI_single_swath_figure(date_str, dtype = 'shawn',  \
            only_sea_ice = False, minlat = 65., skiprows = None, \
            lat_circles = None, save = False, zoom = False, \
            circle_bound = True, ax = None, \
            shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
##!#    OMI_base = readOMI_swath_shawn(date_str, latmin = 65., \
##!#        shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
##!#
##!#    CERES_date_str = np.min(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]).strftime('%Y%m%d%H')
##!#
##!#    modis_list = download_MODIS_swath(CERES_date_str, \
##!#            dest_dir = '/home/bsorenson/data/MODIS/Aqua/', download = False)
##!#
##!#    print(date_str)
##!#    print('    OMI - ', date_str)
##!#    print('  CERES - ', CERES_date_str)
##!#    print('  MODIS - ', *modis_list)
##!#    print('  NSIDC - ', CERES_date_str[:10])
##!#    
##!#    #print('    ', np.min(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]), np.max(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]))


sys.exit()

minlat = 65.
OMI_base = date_str
#OMI_base = '201908110351'
write_shawn_to_HDF5(OMI_base, save_path = '/home/bsorenson/Research/Arctic_compares/comp_data/20180705/', minlat = 65., \
    shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
sys.exit()

plot_compare_OMI_CERES_MODIS_NSIDC('201908110125', '7', \
    omi_dtype = 'shawn', minlat = 65., zoom = True, save = False)
sys.exit()

sys.exit()


#plot_compare_OMI_CERES_MODIS_NSIDC('201808241435', '7', \



#plotOMI_single_swath_figure(date_str, dtype = 'shawn',  \
#        only_sea_ice = False, minlat = 65., skiprows = None, \
#        lat_circles = None, save = False, zoom = True)
#
sys.exit()

plotOMI_daily_control_shawn('20170818', resolution = 1.0,
    shawn_path = '/home/bsorenson/data/OMI/shawn_files/')
sys.exit()

plot_Arctic_row_coverage_compare(date_str = '20180726', save = False)


# NOTE: for plotting the bias between OMI rows, use this line with the
#       CSCI netCDF data
#plt.plot(np.nanmean(np.nanmean(netdata['AI'], axis = 0), axis = 0))
##!#infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_100.txt'
##!#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = True)
##!#infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_150.txt'
##!#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = True)
##!#infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_200.txt'
##!#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = True)
##!#infile = '/home/bsorenson/Research/OMI/shawn_analysis/count_analysis/omi_vsj4_areas_2005_2020_250.txt'
##!#plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = True)

sys.exit()




#OMI_data, CERES_data = plot_OMI_CERES_trend_compare_summer(minlat=72,\
#        ceres_type = 'sw', trend_type = 'standard', save=False)




##bad_row_file = 'row_anomaly_dates_20050401_20201001.txt'
##xtrack_file = 'row_anomaly_xtrack_dates_20050401_20201001.txt'
##plot_bad_row_table(bad_row_file, xtrack_file = xtrack_file, ax = None, \
##        save = False)

sys.exit()

plot_OMI_fort_out_func(infile,\
     min_lat = 70., vtype = 'areas', save = False)

sys.exit()

date_str = ['200509270134',\
            '200509270313',\
            '200509270451',\
            '200509270630',\
            '200509270809',\
            '200509270948',\
            '200509271127',\
            '200509271306',\
            '200509271445',\
            '200509271624',\
            '200509271802',\
            '200509271941',\
            '200509272120',\
            '200509272259']

###
###    
###plotOMI_single_ground(date_str, only_sea_ice = False, minlat = 65., \
###    zoom = True, multi_panel = False, save = True)

plotOMI_single_swath_figure(date_str, dtype = 'control',  \
        only_sea_ice = False, minlat = 65., skiprows = [11,15, 16, 17, 18], lat_circles = None, save = False)

#plot_combined_figure1(save = True)
#plot_MODIS_temporary('202107222110', zoom = True, save = True)
#plot_MODIS_temporary_4panel('202108052125', zoom = True, composite = True, show_smoke = True, save = True)
#compare_MODIS_3panel('202107222110',31,1,5,zoom=True,save=True,\
#        plot_ASOS_loc = False, show_smoke = True)
#plot_combined_figure3(save = True)
#plot_figureS1(save=True, composite = True)
#plot_combined_figureS2(save=True)
