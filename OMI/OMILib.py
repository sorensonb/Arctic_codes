#!/usr/bin/env python
"""
  NAME:

  PURPOSE:
    Plot trends in OMI-measured aerosol index across the northern hemisphere.
    The average AI values are read from files produced by 
    Average_AerosolIndexCalculator.pro (included in ai_trend_codes_YYMMDD.tar.gz)

  PYTHON VERSION:
    2.6.6

  MODULES:
    - Custom module AILib
    - Matplotlib
    - mpl_toolkits.basemap
    - datetime
    - sys
    - numpy
    
"""

import numpy as np
import sys
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
#from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.path as mpath
import matplotlib.colors as color
from matplotlib.colors import rgb2hex,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase
import matplotlib.gridspec as gridspec
from scipy import stats
from netCDF4 import Dataset
import gzip
import h5py
import subprocess
from scipy.stats import pearsonr,spearmanr
from sklearn.linear_model import HuberRegressor
from sklearn.preprocessing import StandardScaler
import pandas as pd

def get_ice_flags(value):
    return int(format(value,"016b")[-15:-8],2)

def get_xtrack_flags(value, flag_idx = None):
    #return int(format(XTRACK[1000,yyy],"08b")[-3:],2))
    bit_vals = np.full(5,np.nan)
    bit_vals[0] = int(format(value,"08b")[-3:],2)
    bit_vals[1] = int(format(value,"08b")[-4],2)
    bit_vals[2] = int(format(value,"08b")[-5],2)
    bit_vals[3] = int(format(value,"08b")[-6],2)
    bit_vals[4] = int(format(value,"08b")[-7],2)
    if(flag_idx == None):
        return bit_vals
    else:
        return bit_vals[flag_idx]

# Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

# Set up mapping variables 
datacrs = ccrs.PlateCarree() 
colormap = plt.cm.jet
mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

label_dict = {
    'V003': 'Only XTrack == 0',
    'VJZ2': 'No Snow-free Land',
    'VJZ28': 'No Snow-free Land, Only Pure Good Rows',
    'VJZ282': 'No Snow-free Land, Only Pure Good Rows',
    'VJZ29': 'Include Snow-free Land, Only Pure Good Rows',
    'VJZ211': 'Include Snow-free Land, Rows 55 - 60',
    'VJZ4': 'XTrack == 0, not 4',
    'VJZ5': 'AI >= 0',
    'VBS0': 'No Bad Row Screening',
    'VBS1': 'Bad Row Screening Only',
    'VBS2': 'Only Rows 1-22',
    'VSJ2': 'Perturbation Analysis',
    'VSJ4': 'Perturbation Analysis'
}

var_dict = {
    'SurfaceAlbedo':           {'min': 0.0,  'max': 1.0},\
    'Reflectivity':            {'min': 0.0,  'max': 1.0},\
    'NormRadiance':            {'min': 0.0,  'max': 0.2},\
    'CloudFraction':           {'min': None, 'max': None},\
    'UVAerosolIndex':          {'min': -2.0, 'max': 3.0 },\
    'PixelQualityFlags':       {'min': None, 'max': None},\
    'MeasurementQualityFlags': {'min': None, 'max': None},\
    'FinalAlgorithmFlags':     {'min':    0, 'max':    8},\
    'ViewingZenithAngle':      {'min': 0.0,  'max': 180.},\
    'SolarZenithAngle':        {'min': None, 'max': None},\
    'RelativeZenithAngle':     {'min': None, 'max': None},\
    'GroundPixelQualityFlags': {'min': None, 'max': None},\
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Miscellaneous functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def plot_subplot_label(ax, label, xval = None, yval = None, transform = None, \
        color = 'black', backgroundcolor = None, fontsize = 14, \
        location = 'upper_left'):

    if(location == 'upper_left'):
        y_lim = 0.90
        x_lim = 0.05
    elif(location == 'lower_left'):
        y_lim = 0.05
        x_lim = 0.05
    elif(location == 'upper_right'):
        y_lim = 0.90
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

def plot_theil_sen_trend(pax, xdata, ydata, color, linestyle):
    res = stats.theilslopes(ydata, xdata, 0.90)
    print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))

    # Then, plot the trend line on the figure
    pax.plot(xdata, res[1] + res[0] * xdata, \
        color='k', linewidth = 2.5, linestyle = linestyle)
    # Then, plot the trend line on the figure
    pax.plot(xdata, res[1] + res[0] * xdata, \
        color=color, linestyle = linestyle)
    
def plot_lin_regress_trend(pax, xdata, ydata, color, linestyle):
    # First, calculate the trend
    zdata = np.polyfit(xdata, ydata, 1)

    print("Lin Regress: {0}x + {1}".format(*zdata))

    # Then, plot the trend line on the figure
    pax.plot(np.unique(xdata), np.poly1d(zdata)(np.unique(xdata)), \
        color='k', linewidth = 2.5, linestyle = linestyle)
    pax.plot(np.unique(xdata), np.poly1d(zdata)(np.unique(xdata)), \
        color=color, linestyle = linestyle)

def plot_trend_line(pax, xdata, ydata, color='black', linestyle = '-', \
        slope = 'theil-sen'):

    if(slope == 'theil-sen'):
        plot_theil_sen_trend(pax, xdata, ydata, color, linestyle)
    elif(slope == 'both'):
        plot_theil_sen_trend(pax, xdata, ydata, color, linestyle)
        plot_lin_regress_trend(pax, xdata, ydata, color, linestyle)
    else:
        plot_lin_regress_trend(pax, xdata, ydata, color, linestyle)



# Class MidpointNormalize is used to center the colorbar on zero
class MidpointNormalize(Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

# Function onclick performs an operation at the location of a click on the
# map. 
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    #print ix, iy  
    mapLon, mapLat = m(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = [OMI_data[dictkey][date]['avg'] for date in                           \
        sorted(OMI_data[dictkey].keys())]
    dates = sorted(OMI_data[dictkey].keys())
    x_vals = np.arange(0,len(OMI_data[dictkey].keys()))

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals,avgs)
    print(dictkey, slope*len(x_vals))
    #print slope/len(OMI_data[dictkey].keys())
    regress_y = x_vals*slope+intercept

    #The slope
    S=0
    sm=0
    nx = len(avgs)
    num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
    Sn=np.zeros(num_d)
    for si in range(0,nx-1):
        for sj in range(si+1,nx):
            # Find the slope between the two points
            Sn[sm] = (avgs[si]-avgs[sj])/(si-sj) 
            sm=sm+1
        # Endfor
    # Endfor
    Snsorted=sorted(Sn)
    sm=int(num_d/2.)
    if(2*sm    == num_d):
        tsslope=0.5*(Snsorted[sm]+Snsorted[sm+1])
    if(2*sm+1 == num_d): 
        tsslope=Snsorted[sm+1]
    regress_ts = x_vals*tsslope+intercept

    # Convert the dates to datetime objects
    ddates = [datetime.strptime(date,'%Y%m') for date in dates] 

    label_dates = [ddate.strftime('%b %Y') for ddate in ddates]

    tickNumber = 12

    fig1 = plt.figure()
    plt.plot(avgs)
    plt.plot(x_vals,regress_y,linestyle='--',color='red',label='Control')
    plt.plot(x_vals,regress_ts,linestyle='--',color='blue',label='Thiel-Sen')
    plt.title(dictkey)
    plt.xticks(np.arange(0,len(avgs))[::-int(len(avgs)/tickNumber)],label_dates[::\
        -int(len(avgs)/tickNumber)],rotation=45)
    plt.ylabel('Average Aerosol Index')
    plt.legend()
    plt.show()

def onclick_climo(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = m(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([OMI_data[dictkey][date] for date in sorted(OMI_data[dictkey].keys())])
    ai_avg = np.average(avgs)
    print(dictkey, ai_avg)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Reading functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Written for old old data
   
def readOMI(inputfile,start_date,end_date,key=None):
    global OMI_data
    OMI_data = {}

    if(key is not None):
        OMI_data[key]={}

    if(inputfile.strip().split('/')[-1].split('.')[-1]=='gz'):
        f = gzip.open(inputfile,'rb')
    else:
        f = open(inputfile,'r')
    #with open(inputfile,'r') as f:
    # Skip the header line
    for line in f:
        templine = line.strip().split()
        if(len(templine)>1):
            if(len(templine) == 5):
                loc_key = str(templine[1])+'x'+str(templine[2])
                avg_idx = 3
                cnt_idx = 4
            else:
                loc_key = templine[1] 
                avg_idx = 2
                cnt_idx = 3
            if((int(templine[0])>=start_date) & (int(templine[0])<=end_date)):
                if(key is not None):
                    if(loc_key==key):
                        OMI_data[key][templine[0]] = {}
                        OMI_data[key][templine[0]]['avg']=float(templine[avg_idx])
                        OMI_data[key][templine[0]]['#_obs']=int(templine[cnt_idx])
                else:
                    # If the current lat/lon pair are not found in the dictionary's
                    # keys, then make a new subdictionary for it.
                    if(loc_key not in OMI_data.keys()):
                        OMI_data[loc_key] = {}
                    # If the current lat/lon pair are already in the dictionary's
                    # keys, then add the new data to the subdictionary
                    OMI_data[loc_key][templine[0]]={}
                    OMI_data[loc_key][templine[0]]['avg']=float(templine[avg_idx])
                    OMI_data[loc_key][templine[0]]['#_obs']=int(templine[cnt_idx])
    f.close()    

    return OMI_data

# NOTE: This only works for plotting 1 file time at a time. No multiple swaths
# dtype is either "control" or "JZ"
def readOMI_swath_hdf(plot_time, dtype, only_sea_ice = False, latmin = 65):

    # Extract date information to use when finding the data files
    year = plot_time[:4]
    date = plot_time[4:8]
    if(len(plot_time)==13):
        time = plot_time[9:]
    elif(len(plot_time)==12):
        time = plot_time[8:]
    else:
        time = ''

    # Set up a dictionary to hold the data
    # -------------------------------------
    OMI_swath = {}
    OMI_swath['date'] = plot_time
    OMI_swath['dtype'] = dtype

    base_path = '/home/bsorenson/data/OMI/H5_files/'
    total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+\
        year+'m'+date+'t'+time+'*.he5',shell=True).decode('utf-8').strip().split('\n')

    # Depending on what the user wants, 
    print(total_list[0])
    data = h5py.File(total_list[0],'r')

    LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
    LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
    XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
    UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
    GPQF   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags'][:,:]
    GPQF_decode = np.array([[get_ice_flags(val) for val in GPQFrow] for GPQFrow in GPQF])
    XTRK_decode = np.array([[get_xtrack_flags(val) for val in GPQFrow] for GPQFrow in GPQF])
    mask_UVAI = np.ma.masked_where((XTRACK != 0) | (UVAI < -2e5) | (LAT < latmin), UVAI)
    #UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/' + path_dict[variable] + variable][:,:]
    ##!#if(len(UVAI.shape) == 3):
    ##!#    # If 3 dimensions, assume that the smallest dimension is the wavelength
    ##!#    # dimension. Find that dimension and grab the first index of that
    ##!#    # dimension, which is the 354 nm channel.
    ##!#    min_dim = np.min(UVAI.shape)
    ##!#    if(UVAI.shape[2] == min_dim):
    ##!#        UVAI = UVAI[:,:,channel_idx]
    ##!#    elif(UVAI.shape[1] == min_dim):
    ##!#        UVAI = UVAI[:,channel_idx,:]
    ##!#    else:
    ##!#        UVAI = UVAI[channel_idx,:,:]

    OMI_swath['LAT'] = LAT
    OMI_swath['LON'] = LON
    OMI_swath['XTRACK'] = XTRACK
    OMI_swath['GPQF'] = GPQF_decode

    if(dtype == 'JZ'):
        AZM    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle'][:,:]

        mask_UVAI[~((((GPQF_decode >= 0) & (GPQF_decode <= 101)) | (GPQF_decode == 104)) & \
                      (AZM > 100))] = np.nan
        mask_UVAI = np.ma.masked_invalid(mask_UVAI)

        OMI_swath['AZM'] = AZM
    OMI_swath['UVAI'] = mask_UVAI

    data.close()

    return OMI_swath

def readOMI_swath_shawn(plot_time, latmin = 65.):

    # Read in the OMI data for this file, for matching purposes
    # ---------------------------------------------------------
    OMI_data = readOMI_swath_hdf(plot_time, 'control', latmin = latmin)

    # Open the shawn file
    # -------------------
    shawn_path = '/home/bsorenson/data/OMI/shawn_files/'
    df_shawn = pd.read_csv(shawn_path + plot_time, delim_whitespace = True, \
        header = None)

    # Extract the needed variables
    # ----------------------------
    s_lat    = pd.to_numeric(df_shawn[0], errors = 'coerce').values
    s_lon    = pd.to_numeric(df_shawn[1], errors = 'coerce').values
    s_airaw  = pd.to_numeric(df_shawn[2], errors = 'coerce').values
    s_aiclim = pd.to_numeric(df_shawn[3], errors = 'coerce').values
    s_aipert = pd.to_numeric(df_shawn[4], errors = 'coerce').values

    # Make nan arrays dimensioned to the OMI swath data and fill with
    # nans. Then, loop over the shawn data, and wherever each shawn
    # value has a matching lat/lon pair in the original data, insert
    # the shawn data there
    # ---------------------------------------------------------------
    #sfile_LAT  = np.full(OMI_data['LAT'].shape, np.nan)
    #sfile_LON  = np.full(OMI_data['LAT'].shape, np.nan)
    sfile_AIraw  = np.full(OMI_data['LAT'].shape, np.nan)
    sfile_AIclim = np.full(OMI_data['LAT'].shape, np.nan)
    sfile_AIpert = np.full(OMI_data['LAT'].shape, np.nan)

    for slat, slon, sair, saic, saip in \
            zip(s_lat, s_lon, s_airaw, s_aiclim, s_aipert):
        match_loc = np.where((slat == OMI_data['LAT']) & \
            (slon == OMI_data['LON']))
        #print(slat, slon, sair, OMI_data['LAT'][match_loc], \
        #    OMI_data['LON'][match_loc], OMI_data['UVAI'][match_loc], match_loc)
        #sfile_LAT[match_loc] = slat
        #sfile_LON[match_loc] = slon
        sfile_AIraw[match_loc] = sair
        sfile_AIclim[match_loc] = saic
        sfile_AIpert[match_loc] = saip

    OMI_swath = {}
    OMI_swath['date'] = plot_time
    OMI_swath['dtype'] = 'shawn'
    OMI_swath['LAT'] = OMI_data['LAT']
    OMI_swath['LON'] = OMI_data['LON']
    OMI_swath['UVAI_raw']   = np.ma.masked_invalid(sfile_AIraw)
    OMI_swath['UVAI_climo'] = np.ma.masked_invalid(sfile_AIclim)
    OMI_swath['UVAI_pert']  = np.ma.masked_invalid(sfile_AIpert)

    return OMI_swath

def readOMI_single_swath(plot_time,row_max,only_sea_ice = True,latmin=65,coccolith = False):
    n_p = 1440
    nl = 720
    lonmin = -180.
    lonmax = 180.
    #latmin = 60.
    latmax = 90.
    # Set up values for gridding the AI data
    lat_gridder = latmin * 4.
    
    lat_ranges = np.arange(latmin,90.,0.25)
    lon_ranges = np.arange(-180.,180.,0.25)

    if(coccolith):
        g_NRAD_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        g_NRAD_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        g_NRAD_500 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        g_REFL_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        g_REFL_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        g_SALB_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_NRAD_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_NRAD_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_NRAD_500 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_REFL_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_REFL_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_SALB_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    g_UVAI_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_UVAI_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    g_SZA = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_SZA = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

    algae = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

    print("only_sea_ice = ",only_sea_ice)

    year = plot_time[:4]
    date = plot_time[4:8]
    if(len(plot_time)==13):
        time = plot_time[9:]
    elif(len(plot_time)==12):
        time = plot_time[8:]
    else:
        time = ''
    base_path = '/home/bsorenson/data/OMI/H5_files/'
    total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+year+'m'+date+'t'+time+'*.he5',\
              shell=True).decode('utf-8').strip().split('\n')

    for fileI in range(len(total_list)):
        # read in data directly from HDF5 files
        data = h5py.File(total_list[fileI],'r')
        print(total_list[fileI])
    
        LAT     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
        LON     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
        UVAI    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
        XTRACK  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
        GPQF    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags']
        SZA     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SolarZenithAngle']
        if(coccolith):
            NRAD    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/NormRadiance']
            REFL    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity']
            SALB    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo']
    
        #albedo = ALBEDO[:,:,0]   
        #reflectance = REFLECTANCE[:,:,0]   
        counter = 0
        #AI = AI[:,:,0]   
        # Loop over the values and rows 
        #for i in range(0,int(CBA2)):
        #for i in range(albedo.shape[0]):
        for i in range(UVAI.shape[0]):
            for j in range(0,row_max):
                #if(j == 52):
                #    continue
                #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
                if(LAT[i,j] > latmin):
                    if((UVAI[i,j]>-2e5) & ((XTRACK[i,j] == 0) | (XTRACK[i,j] == 4))):
                        # 0       :  snow-free land
                        # 1 - 100 :  sea ice concentration (percent)
                        # 101     :  permanent ice (greenland, antarctica)
                        # 103     :  dry snow
                        # 104     :  ocean   
                        GPQF_decode = get_ice_flags(GPQF[i,j])
                        #if(only_sea_ice and not (((GPQF_decode == 101) & (SZA[i,j] < 60)))):
                        if(only_sea_ice and not (((GPQF_decode >= 1) & (GPQF_decode <= 101)) )):
                        #if(only_sea_ice and not (((GPQF_decode >= 0) & (GPQF_decode <= 101)) | (GPQF_decode == 104)  )):
                            continue
                    # Only plot if XTrack flag is met
                        # Print values to text file
                        counter+=1
                        index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                        index2 = int(np.floor(LON[i,j]*4 + 720.))
                        #index1 = int(np.floor(plotLAT[i,j]*4 + 360.))
                        #index2 = int(np.floor(plotLON[i,j]*4 + 720.))
                        
                        if(index1 < 0): index1 = 0
                        if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                        if(index2 < 0): index2 = 0                                                                                            
                        if(index2 > 1439): index2 = 1439
                   
                        #if(diff<0.2): 
                        #    UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
                        if(coccolith):
                            g_NRAD_354[index2, index1] = (g_NRAD_354[index2,index1]*count_NRAD_354[index2,index1] + NRAD[i,j,0])/(count_NRAD_354[index2,index1]+1)
                            g_NRAD_388[index2, index1] = (g_NRAD_388[index2,index1]*count_NRAD_388[index2,index1] + NRAD[i,j,1])/(count_NRAD_388[index2,index1]+1)
                            g_NRAD_500[index2, index1] = (g_NRAD_500[index2,index1]*count_NRAD_500[index2,index1] + NRAD[i,j,2])/(count_NRAD_500[index2,index1]+1)
                            g_REFL_354[index2, index1] = (g_REFL_354[index2,index1]*count_REFL_354[index2,index1] + REFL[i,j,0])/(count_REFL_354[index2,index1]+1)
                            g_REFL_388[index2, index1] = (g_REFL_388[index2,index1]*count_REFL_388[index2,index1] + REFL[i,j,1])/(count_REFL_388[index2,index1]+1)
                            g_SALB_354[index2, index1] = (g_SALB_354[index2,index1]*count_SALB_354[index2,index1] + SALB[i,j,0])/(count_SALB_354[index2,index1]+1)
                            count_NRAD_354[index2,index1] = count_NRAD_354[index2,index1] + 1
                            count_NRAD_388[index2,index1] = count_NRAD_388[index2,index1] + 1
                            count_NRAD_500[index2,index1] = count_NRAD_500[index2,index1] + 1
                            count_REFL_354[index2,index1] = count_REFL_354[index2,index1] + 1
                            count_REFL_388[index2,index1] = count_REFL_388[index2,index1] + 1
                            count_SALB_354[index2,index1] = count_SALB_354[index2,index1] + 1
                        g_UVAI_354[index2, index1] = (g_UVAI_354[index2,index1]*count_UVAI_354[index2,index1] + UVAI[i,j])/(count_UVAI_354[index2,index1]+1)
                        count_UVAI_354[index2,index1] = count_UVAI_354[index2,index1] + 1
                        g_SZA[index2, index1] = (g_SZA[index2,index1]*count_SZA[index2,index1] + SZA[i,j])/(count_SZA[index2,index1]+1)
                        count_SZA[index2,index1] = count_SZA[index2,index1] + 1

  
    if(coccolith): 
        # Apply algae screening to 500 nm normalized radiance data 
        print("Applying algae screening to 500 nm normalized radiance data")
        mask_rad500 = np.ma.masked_where(((count_NRAD_500 == 0)), g_NRAD_500)
        mask_rad500 = np.ma.masked_where(((g_SALB_354 > 0.09)), mask_rad500)
        mask_rad500 = np.ma.masked_where(((g_REFL_354 > 0.18)), mask_rad500)
        mask_rad500 = np.ma.masked_where(((g_REFL_354 < 0.09)), mask_rad500)
        mask_NRAD500 = np.ma.masked_where(((g_UVAI_354 < 0.6)), mask_rad500)

        # Apply algae screening to UVAI data 
        print("Applying algae screening to UVAI data")
        mask_uvai = np.ma.masked_where(((count_UVAI_354 == 0)), g_UVAI_354)
        mask_uvai = np.ma.masked_where(((g_SALB_354 > 0.09)), mask_uvai)
        mask_uvai = np.ma.masked_where(((g_REFL_354 > 0.18)), mask_uvai)
        mask_uvai = np.ma.masked_where(((g_REFL_354 < 0.09)), mask_uvai)
        mask_UVAI354 = np.ma.masked_where(((g_UVAI_354 < 0.6)), mask_uvai)

    OMI_single_dict = {}
    OMI_single_dict['lat'] = lat_ranges
    OMI_single_dict['lon'] = lon_ranges
    OMI_single_dict['AI'] = g_UVAI_354
    OMI_single_dict['AI_count'] = count_UVAI_354
    OMI_single_dict['SZA'] = g_SZA
    OMI_single_dict['SZA_count'] = count_SZA
    OMI_single_dict['date'] = plot_time
    OMI_single_dict['row_max'] = row_max
    if(coccolith):
        OMI_single_dict['NRAD500'] = g_NRAD_500
        OMI_single_dict['NRAD500_count'] = count_NRAD_500
        OMI_single_dict['AI_algae'] = mask_UVAI354
        OMI_single_dict['NRAD500_algae'] = mask_NRAD500

    return OMI_single_dict

def readOMI_NCDF(infile='/home/bsorenson/Research/OMI/omi_ai_V003_2005_2020.nc',\
                 start_date = 200504, end_date = 202009, calc_month = True, \
                 minlat=50):
    # Read in data to netCDF object
    in_data = Dataset(infile,'r')

    # Pull the version type from the filename
    version = infile.split('/')[-1].split('_')[2]
   
    # Set up dictionary to hold data
    OMI_data = {}

    # Set up date strings in the file
    base_date = datetime(year = 2005, month = 1, day = 1)
    str_dates = \
        np.array([(base_date + relativedelta(months=mi)).strftime('%Y%m') for mi in \
            in_data['MONTH']])

    # Use the file dates to select just the data within the user-specified
    # timeframe
    # --------------------------------------------------------------------
    int_str_dates = np.array([int(tdate) for tdate in str_dates])

    time_indices = np.where((int_str_dates >= start_date) & \
        (int_str_dates <= end_date))[0]

    #check_int_dates
    #check_int_dates = np.array([\
    #    (start_date + relativedelta(months=mi)).strftime('%Y:%m')\
    #     for mi in in_data['MONTH']])

    
    OMI_data['AI']       = in_data['AI'][time_indices,:,:]
    OMI_data['OB_COUNT'] = in_data['OB_COUNT'][time_indices,:,:]
    OMI_data['LAT']      = in_data['Latitude'][:,:]
    OMI_data['LON']      = in_data['Longitude'][:,:]
    OMI_data['MONTH']    = in_data['MONTH'][time_indices]
    OMI_data['DATES']    = str_dates[time_indices]
    OMI_data['VERSION']  = version

    if(calc_month == True):
        OMI_data = calcOMI_MonthClimo(OMI_data)

    # to add months to datetime object, do
    ###from dateutil.relativedelta import relativedelta
    ###datetime.datetime(year=2004,month=10,day=1) + relativedelta(months=1)
    
    in_data.close()
   
    return OMI_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Data writing functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def writeOMI_toNCDF(OMI_data,file_name,minlat=30):
    lat_ranges = np.arange(minlat,90.5,1.0)
    lon_ranges = np.arange(-180,180,1.0)
   
    # Create a new netCDF dataset to write to
    # --------------------------------------- 
    nc = Dataset(file_name,'w',format='NETCDF4')
  
    # Dimensions = lat, lon, time
    # Create the sizes of each dimension in the file. In this case,
    # the dimensions are "# of latitude", "# of longitude", and 
    # "# of times in the file". 
    # -------------------------------------------------------------
    testdict = OMI_data['70x5']
    testkeys = list(testdict.keys())
    num_lat = len(lat_ranges)
    num_lon = len(lon_ranges)
    num_time = len(testdict.keys())
    times = np.arange(num_time)
  
    # Instead of using a simple 'arange' function to define the 'times'
    # variable, actually calculate the number of months between each date
    # variable and a reference date, which is January 2005. This will have
    # no effect on any total OMI processes but will allow compatibility with
    # testing processes where just a few months between January 2005 and
    # July 2019 are selected. 
    base_date = datetime(year=2005,month=1,day=1) 
    times = np.array([(datetime.strptime(tmpx,'%Y%m').year - base_date.year) \
        * 12 + datetime.strptime(tmpx,'%Y%m').month - base_date.month \
        for tmpx in testkeys])
   
    # Use the dimension size variables to actually create dimensions in 
    # the file.
    # ----------------------------------------------------------------- 
    n_time = nc.createDimension('nmth',num_time)
    n_lat  = nc.createDimension('dlat',num_lat)
    n_lon  = nc.createDimension('dlon',num_lon)

    # Create variables for the three dimensions. Note that since these
    # are variables, they are still given 'dimensions' using 'nmth',
    # 'dlat', and 'dlon'. Latitude and longitude are each 2-d grids, 
    # so they are given 2 dimensions (dlat, dlon).
    # ----------------------------------------------------------------
    MONTH=nc.createVariable('MONTH','i2',('nmth'))
    MONTH.description='Months since January 2005'
    LAT=nc.createVariable('Latitude','i2',('dlat','dlon'))
    LAT.description='Latitude'
    LAT.units='Degrees'
    LON=nc.createVariable('Longitude','i2',('dlat','dlon'))
    LON.description='Longitude'
    LON.units='Degrees'
   
    # Create a variable for the AI data, dimensioned using all three
    # dimensions.
    # --------------------------------------------------------------  
    AI = nc.createVariable('AI','f4',('nmth','dlat','dlon'))
    AI.description='Monthly Averaged Aerosol Index (using JZ method)'
    OB_COUNT=nc.createVariable('OB_COUNT','i2',('nmth','dlat','dlon'))
    OB_COUNT.description='# of OMI AI measurements used in each monthly average'

    # Fill in dimension variables one-by-one.
    # NOTE: not sure if you can insert the entire data array into the
    # dimension (for example, doing MONTH = times), so this could be
    # something for you to try. Might make this faster if it works
    # ---------------------------------------------------------------
    for i in range(num_time):
        MONTH[i] = times[i]
    for i in range(num_lat):
        for j in range(num_lon):
            LAT[i,j]=lat_ranges[i]
            LON[i,j]=lon_ranges[j]

    # Fill in actual variables
    # I have a bunch of extra stuff in here for handling the dictionary
    # keys, which you will likely not need for your data
    # ------------------------------------------------------------------
    for i in range(num_lat):
        print(lat_ranges[i])
        for j in range(num_lon):
            dictkey = (str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j])))
            if(dictkey not in OMI_data):
                # Insert missing values for AI and count
                AI[:,i,j] = [-999.9 for m in range(num_time)]
                OB_COUNT[:,i,j] = [-99 for m in range(num_time)]
            else:
                for m in range(num_time):
                    timekey = testkeys[m]
                    if(timekey not in OMI_data[dictkey]):
                        AI[m,i,j] = -999.9
                        OB_COUNT[m,i,j] = -99
                    else:
                        AI[m,i,j] = OMI_data[dictkey][timekey]['avg']
                        OB_COUNT[m,i,j] = OMI_data[dictkey][timekey]['#_obs']

    # Close the netCDF file, which actually saves it.
    # -----------------------------------------------
    nc.close()

def write_da_to_NCDF(avgAI,counts,latmin,da_time):

    # Set up values for gridding the AI data
    lat_gridder = latmin * 4.
    
    lat_ranges = np.arange(latmin,90.1,0.25)
    lon_ranges = np.arange(-180,180.1,0.25)

    # Create a new netCDF dataset to write to
    # --------------------------------------- 
    outfile = '/home/bsorenson/Research/OMI/omi_ai_da_'+da_time.strftime("%Y%m%d%H") + '.nc'
    nc = Dataset(outfile,'w',format='NETCDF4')
  
    # Dimensions = lon, lat
    # Create the sizes of each dimension in the file. In this case,
    # the dimensions are "# of latitude" and "# of longitude"
    # -------------------------------------------------------------
    num_lon = len(lon_ranges)
    num_lat = len(lat_ranges)
    
    # Use the dimension size variables to actually create dimensions in 
    # the file.
    # ----------------------------------------------------------------- 
    n_lon  = nc.createDimension('nlon',num_lon)
    n_lat  = nc.createDimension('nlat',num_lat)

    # Create variables for the three dimensions. Note that since these
    # are variables, they are still given 'dimensions' using 'dlat', 
    # and 'dlon'. Latitude and longitude are each 2-d grids, so they are 
    # given 2 dimensions (dlat, dlon).
    # ------------------------------------------------------------------
    LON=nc.createVariable('lon','f4',('nlon','nlat'))
    LON.description='Longitude'
    LON.units='Degrees'
    LAT=nc.createVariable('lat','f4',('nlon','nlat'))
    LAT.description='Latitude'
    LAT.units='Degrees'

    # Create a variable for the AI data, dimensioned using both 
    # dimensions.
    # ---------------------------------------------------------  
    AI = nc.createVariable('AI','f4',('nlon','nlat'))
    AI.description='Quality controlled aerosol index data during the data assimilation window ('+\
        da_time.strftime("%Y%m%d%H") + ' +/- 3 hrs)'
    AI_COUNT=nc.createVariable('AI_COUNT','i2',('nlon','nlat'))
    AI_COUNT.description='# of OMI AI measurements used in each grid box'

    # Fill in dimension variables one-by-one.
    # NOTE: not sure if you can insert the entire data array into the
    # dimension (for example, doing :
    #      LAT, LON = np.meshgrid(lat_ranges,lon_ranges) ),
    # so this could be
    # something for you to try. Might make this faster if it works
    # ---------------------------------------------------------------
    for yi in range(num_lon):
        for xj in range(num_lat):
            LON[yi,xj]=lon_ranges[yi]
            LAT[yi,xj]=lat_ranges[xj]

    # Fill in actual variables
    # ------------------------
    for yi in range(num_lon):
        for xj in range(num_lat):
            AI[yi,xj] = avgAI[yi,xj]
            AI_COUNT[yi,xj] = counts[yi,xj]

    # Save, write, and close the netCDF file
    # --------------------------------------
    nc.close()

    print("Saved file ",outfile)  

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Calculation functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def calcOMI_grid_trend(OMI_data, month_idx, trend_type, minlat):
    version = OMI_data['VERSION']

    if(month_idx == None):
        month_idx = 0
        index_jumper = 1
    else:
        month_adder = '_month'
        if(version == 'V003'):
            index_jumper = 12
        else:
            index_jumper = 6 

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    # Make copy of OMI_data array
    print(OMI_data['DATES'][month_idx::index_jumper])
    local_data   = np.copy(OMI_data['AI'][month_idx::index_jumper,:,:])
    local_counts = np.copy(OMI_data['OB_COUNT'][month_idx::index_jumper,:,:])
    local_mask = np.ma.masked_where(local_counts == 0, local_data)
    local_mask = np.ma.masked_where((local_mask == -999.9) & \
        (OMI_data['LAT'] < minlat), local_mask)
    ai_trends = np.zeros(local_data.shape[1:])
    print(local_data.shape)

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    for i in range(0,len(np.arange(np.min(OMI_data['LAT']),90))):
        for j in range(0,len(lon_ranges)):
            # Check the current max and min
            x_vals = np.arange(0,len(local_mask[:,i,j]))
            # Find the slope of the line of best fit for the time series of
            # average data
            if(trend_type=='standard'): 
                slope, intercept, r_value, p_value, std_err = \
                    stats.linregress(x_vals,local_mask[:,i,j])
                ai_trends[i,j] = slope * len(x_vals)
            else:
                res = stats.theilslopes(local_mask[:,i,j], x_vals, 0.90)
                ai_trends[i,j] = res[0]*len(x_vals)

    ai_trends = np.ma.masked_where(OMI_data['LAT'] < minlat, ai_trends)

    print('in trend calc')
    for x, y in zip(OMI_data['LAT'][:,10], ai_trends[:,10]):
        print(x,y)

    return ai_trends

   
# This function assumes the data is being read from the netCDF file
# NOTE: Assume user is using new OMI climo file which starts in January
def calcOMI_MonthClimo(OMI_data):

    # Set up arrays to hold monthly climatologies
    month_climo = np.zeros((12,OMI_data['AI'].shape[1],OMI_data['AI'].shape[2]))

    # Mask the monthly averages
    local_data   = np.copy(OMI_data['AI'][:,:,:])
    local_counts = np.copy(OMI_data['OB_COUNT'][:,:,:])
    local_mask = np.ma.masked_where(local_counts == 0, local_data)
    local_mask = np.ma.masked_where(local_mask == -999.9, local_mask)
 
    # Calculate monthly climatologies
    for m_i in range(12):
        month_climo[m_i,:,:] = np.nanmean(local_mask[m_i::12,:,:],axis=0)
        print("Month: ",OMI_data['DATES'][m_i][4:]) 
    #month_climo[0,:,:]  = np.nanmean(local_mask[0::12,:,:],axis=0)  # January 
    #month_climo[1,:,:]  = np.nanmean(local_mask[1::12,:,:],axis=0)  # February
    #month_climo[2,:,:]  = np.nanmean(local_mask[2::12,:,:],axis=0)  # March 
    #month_climo[3,:,:]  = np.nanmean(local_mask[3::12,:,:],axis=0)  # April 
    #month_climo[4,:,:]  = np.nanmean(local_mask[4::12,:,:],axis=0)  # May 
    #month_climo[5,:,:]  = np.nanmean(local_mask[5::12,:,:],axis=0)  # June 
    #month_climo[6,:,:]  = np.nanmean(local_mask[6::12,:,:],axis=0)  # July 
    #month_climo[7,:,:]  = np.nanmean(local_mask[7::12,:,:],axis=0)  # August 
    #month_climo[8,:,:]  = np.nanmean(local_mask[8::12,:,:],axis=0)  # September 
    #month_climo[9,:,:]  = np.nanmean(local_mask[9::12,:,:],axis=0)  # October 
    #month_climo[10,:,:] = np.nanmean(local_mask[10::12,:,:],axis=0) # November
    #month_climo[11,:,:] = np.nanmean(local_mask[11::12,:,:],axis=0) # December

    # Insert data into dictionary
    OMI_data['MONTH_CLIMO'] = month_climo

    return OMI_data

# Generate 6-hr cleaned average files to use for DA or aerosol event frequency 
# tests. Start_date and end_date are of format "YYYYMMDDHH". Average files are
# saved in netCDF format.
def omi_da_gen(start_date,end_date):

    latmin = 60 
    
    # Set up values for gridding the AI data
    lat_gridder = latmin * 4.
    
    lat_ranges = np.arange(latmin,90.1,0.25)
    lon_ranges = np.arange(-180,180.1,0.25)

    # Convert the date objects to datetime objects
    dtime_start = datetime.strptime(start_date,"%Y%m%d%H")
    dtime_end   = datetime.strptime(end_date,"%Y%m%d%H")

    base_path = '/home/bsorenson/data/OMI/H5_files/'
    total_list = np.array(subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_*.he5',\
              shell=True).decode('utf-8').strip().split('\n'))

    dtime_list = np.array([datetime.strptime(tname.split('/')[-1][20:24] + \
        tname.split('/')[-1][25:29] + tname.split('/')[-1][30:34],"%Y%m%d%H%M") \
        for tname in total_list]) 

    ##!## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ##!## Get the number of 6 - hour time slots desired   
    ##!#time_diff = dtime_end - dtime_start
    ##!#num_slots = int(((time_diff.days * 24) + (time_diff.seconds) / 3600) / 6)
    ##!## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 
    temp_start = dtime_start
    while(temp_start < dtime_end): 
        temp_end = temp_start + timedelta(hours = 6)
    
        da_time = temp_start + timedelta(hours = 3)
        

        keep_files = total_list[np.where((dtime_list >= temp_start) & \
            (dtime_list <= temp_end))]
   
        print("Grabbing files for ",da_time)
        if(len(keep_files) > 0):
    
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
            # This is where the averaging and outputting will occur
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
            ##print(keep_files)

            # Set up blank grid arrays to hold the counts and the data
            UVAI = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
            count = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

            for infile in keep_files:
                print(infile)
                data = h5py.File(infile,'r')

                LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
                LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
                AI    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
                GPQF  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags'][:,:]
                XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
                AZM    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle'][:,:]
       
                GPQF_decode = np.reshape(np.array([get_ice_flags(tvar) for tvar in GPQF.flatten()]),GPQF.shape)
                # IGNORE ROW 51

                for yi in range(len(lon_ranges)-1):
                    for xj in range(len(lat_ranges)-1):
                        ai_in_range = AI[np.where(((LAT >= lat_ranges[xj]) & (LAT < lat_ranges[xj + 1])) & \
                                                  ((LON >= lon_ranges[yi]) & (LON < lon_ranges[yi + 1])) & \
                                                  ((XTRACK == 0) | (XTRACK == 4)) & \
                                                  (abs(AI) < 50) & \
                                                  (((GPQF_decode >= 0) & (GPQF_decode <= 101)) | (GPQF_decode == 104)) & \
                                                  (AZM > 100))]
                        if(len(ai_in_range) > 0):
                            UVAI[yi,xj] += sum(ai_in_range)  
                            count[yi,xj] += ai_in_range.size
                                                  
                data.close()
            # end file loop

            # Mask any areas where counts are zero.
            #mask_count = np.ma.masked_where(count == 0,count)
            #mask_UVAI  = np.ma.masked_where(count == 0,UVAI)

            # Divide the accumulated AI values by the accumulated counts
            avg_UVAI = np.copy(UVAI)
            avg_UVAI[count > 0] = UVAI[count > 0]/count[count > 0]
            avg_UVAI[count == 0] = -999.

            print("Saving data for time",da_time)
            # Write the data to a netCDF file
            write_da_to_NCDF(avg_UVAI,count,latmin,da_time)
     
        # end number of files check
        temp_start = temp_end 
    # end time loop
    #return UVAI,count,avg_UVAI,latmin,da_time

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Plotting functions (OLD)
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def plotOMI_MK(OMI_data,start_date,end_date,save=False,file_type='XR123',season='',minlat=30.):
    if(file_type=='NXAR'):
        title_flabel = 'No XTrack and All Rows'
        outname_flabel = 'noX_allRows'
    elif(file_type=='XAR'):
        title_flabel = 'XTrack and All Rows'
        outname_flabel = 'xtrack_allRows'
    elif(file_type=='XR123'):
        title_flabel = 'XTrack and Rows 1-23'
        outname_flabel = 'xtrack_rows1to23'
    lat_ranges = np.arange(minlat,90,1.0)
    #lat_ranges = np.arange(lowest_lat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False
    if(season=='spring'):
        spring = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    OMI_data[lkey].pop(tkey)

    # Find the lowest lat in the file
    #lowest_lat = float(sorted(OMI_data.keys())[0].split('x')[0])
     
    # Set up the polar stereographic map
    fig1 = plt.figure()
    ax = plt.subplot(111)
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
   
    # Define the colors for each range of P-Value
    # Ranges: >=0.4,0.4-0.3,0.3-0.2,0.2-0.15,0.15-0.1,0.1-0.05,<0.05
    # Seven colors
    # - White,yellow,orange-yellow,orange,red-orange,red,dark red
    color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
    pRED   = np.array([255.,255.,255.,255.,255.,150.,150.])
    pGREEN = np.array([255.,255.,150., 75.,  0.,  0.,  0.])
    pBLUE  = np.array([255.,  0.,  0.,  0.,  0.,  0.,150.])
    pRED   = pRED/256.
    pGREEN = pGREEN/256.
    pBLUE  = pBLUE/256.

    noRED = 0.
    noGREEN = 0.
    noBLUE = 0.
    
    # - White,blue-green1,orange-yellow,orange,red-orange,red,dark red
    nRED   = np.array([255.,100.,  0.,  0.,  0.,  0.,  0.])
    nGREEN = np.array([255.,255.,255.,255.,200.,100.,  0.])
    nBLUE  = np.array([255.,  0.,  0.,255.,255.,255.,150.])
    nRED   = nRED/256.
    nGREEN = nGREEN/256.
    nBLUE  = nBLUE/256.
 
    ##cmap = plt.get_cmap('bwr')
    ### Center the colorbar on zero
    ##norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
    ###norm = Normalize(vmin=vmin,vmax=vmax)
    ##mapper = ScalarMappable(norm=norm,cmap=cmap)

    # Set up initial values for the analysis regions
    canada_avg = 0.
    canada_counts = 0
    siberia_avg = 0.
    siberia_counts = 0 
    eus_avg = 0.
    eus_counts = 0
    wus_avg = 0.
    wus_counts = 0
    ea_avg = 0.
    ea_counts = 0
    wa_avg = 0.
    wa_counts = 0
    europe_avg = 0.
    europe_counts = 0
    
    # Loop over all the keys and print(the regression slopes 
    # Grab the averages for the key
    max_pval = -10.
    min_pval = 10.
    #for i in range(15,20):
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = (str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j])))
            #print(dictkey)
            if(dictkey=='48x-97'):
                #print("sorted(OMI_data[dictkey].keys())=",sorted(OMI_data[dictkey].keys()))
                min_date = sorted(OMI_data[dictkey].keys())[0]
                max_date = sorted(OMI_data[dictkey].keys())[-1]
    
            # If no data are present for the curent lat/lon box, fill it with
            # black
            if(dictkey not in OMI_data.keys()):
                colorRED = 0.
                colorGREEN = 0.
                colorBLUE = 0.
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(OMI_data[dictkey].keys())==0):
                print(dictkey,pval,color_index_value,'NO DATA')
                #print('Here 2')
                colorRED = 0.
                colorGREEN = 0.
                colorBLUE = 0.
            else:
                avgs = [OMI_data[dictkey][date]['avg'] for date in \
                    sorted(OMI_data[dictkey].keys())]
                if(len(avgs)<2):
                    #print('Here 3')
                    colorRED = 0.
                    colorGREEN = 0.
                    colorBLUE = 0.
                else:
                    # Check the current max and min
                    x_vals = np.arange(0,len(OMI_data[dictkey].keys()))
                    avgs = np.ma.masked_array([OMI_data[dictkey][date]['avg'] for date in sorted(OMI_data[dictkey].keys())])
                    #avgs = avgs[np.where(avgs.mask==False)[0]]
                    temp_dates = np.array(sorted(OMI_data[dictkey].keys()))
                    dates = temp_dates
                    #dates = temp_dates[np.where(avgs.mask==False)[0]]
                    x_vals = np.arange(0,len(dates))
    
                    nx=len(avgs)
                    num_d=int(nx*(nx-1)/2)  # The number of elements in d
                    Sn=np.zeros(num_d)
                    S=0
                    for ti in range(nx-1):
                        for tj in range(ti+1,nx):
                            S+=np.sign(avgs[tj]-avgs[ti])
                        # Endfor
                    # Endfor

                    # Find the unique values in the data
                    uniq = np.unique(avgs)
                    g = len(uniq)
                    if(nx==g):
                        # No ties
                        Vs = (nx*(nx-1.)*(2.*nx-5.))/18.
                    else:
                        tp = np.zeros(uniq.shape)
                        for si in range(g):
                            tp[si] = sum(avgs==uniq[si])
                        Vs = (nx*(nx-1.)*(2.*nx-5.)-np.sum(tp*(tp-1.)*(2.*tp-5.)))/18.
                    if (S > 0.): 
                        z=(S-1.)/np.sqrt(Vs)
                    elif (S < 0.): 
                        z=(S+1.)/np.sqrt(Vs)
                    else: 
                        z=0.
                    # Calculate the p value of the trend
                    pval=2*(1.-stats.norm.cdf(abs(z)))  # (two-side)
                    alpha=0.05
                    ##h = abs(z) > stats.norm.ppf(1.-alpha/2.)
                    ### Determine the trend type 
                    ##if(z>0) and h:
                    ##    trend='increasing'
                    ##elif(z<0) and h:
                    ##    trend='decreasing'
                    ##else:
                    ##    trend='no trend'

                    if(pval<min_pval):
                        min_pval=pval
                    elif(pval>max_pval):
                        max_pval=pval

                    #color = mapper.to_rgba(slope)
    #color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
                    color_index_value = 0
                    color_index = np.where(pval<color_range)[0]
                    if(len(color_index)>0):
                        color_index_value = color_index[-1]+1
                    if(S>0):
                        colorRED = pRED[color_index_value]
                        colorGREEN = pGREEN[color_index_value]
                        colorBLUE = pBLUE[color_index_value]
                    else:
                        colorRED = nRED[color_index_value]
                        colorGREEN = nGREEN[color_index_value]
                        colorBLUE = nBLUE[color_index_value]
    
                    if(np.isnan(pval) == False):
                        # Add value to analysis regions (if possible)
                        if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                            canada_avg+=pval
                            canada_counts+=1
                        elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                            siberia_avg+=pval
                            siberia_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                            eus_avg+=pval
                            eus_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                            wus_avg+=pval
                            wus_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                            ea_avg+=pval
                            ea_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                            wa_avg+=pval
                            wa_counts+=1
                        elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                            europe_avg+=pval
                            europe_counts+=1
    
                        print(dictkey,pval,color_index_value)
            # Find the y coordinates of the current LatxLon box
            y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
            # Find the x coordinates of the current LatxLon box
            x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
            # Convert x and y into map coordinates
            mx, my = m(x,y)
            mxy = zip(mx,my)
            pair1 = (mx[0],my[0])
            pair2 = (mx[1],my[1])
            pair3 = (mx[2],my[2])
            pair4 = (mx[3],my[3])
            #mxy = zip(mx,my)
            # Plot the box on the map using color
            colors = [colorRED,colorGREEN,colorBLUE]
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors)
            #poly = Polygon(mxy,facecolor=color2,edgecolor=color2)
            plt.gca().add_patch(poly)
    
    y = [0.,0.,0.,0.]
    # Find the x coordinates of the current LatxLon box
    x = [0.,0.,0.,0.]
    color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
    # Convert x and y into map coordinates
    mx, my = m(x,y)
    pair1 = (mx[0],my[0])
    pair2 = (mx[1],my[1])
    pair3 = (mx[2],my[2])
    pair4 = (mx[3],my[3])
    # Positive
    # Add legend for 0.3-0.2
    colors = [pRED[1],pGREEN[1],pBLUE[1]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.4-0.3')
    plt.gca().add_patch(poly)
    # Add legend for 0.3-0.2
    colors = [nRED[1],nGREEN[1],nBLUE[1]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.4-0.3')
    plt.gca().add_patch(poly)
    # Add legend for 0.2-0.15
    colors = [pRED[2],pGREEN[2],pBLUE[2]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.3-0.2')
    plt.gca().add_patch(poly)
    # Add legend for 0.15-0.1
    colors = [nRED[2],nGREEN[2],nBLUE[2]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.3-0.2')
    plt.gca().add_patch(poly)
    # Add legend for 0.1-0.05
    colors = [pRED[3],pGREEN[3],pBLUE[3]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.2-0.15')
    plt.gca().add_patch(poly)
    # Add legend for <0.05
    colors = [nRED[3],nGREEN[3],nBLUE[3]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.2-0.15')
    plt.gca().add_patch(poly)
    colors = [pRED[4],pGREEN[4],pBLUE[4]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.15-0.1')
    plt.gca().add_patch(poly)
    # Negative
    # Add legend for 0.2-0.15
    colors = [nRED[4],nGREEN[4],nBLUE[4]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.15-0.1')
    plt.gca().add_patch(poly)
    # Add legend for 0.15-0.1
    colors = [pRED[5],pGREEN[5],pBLUE[5]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.1-0.05')
    plt.gca().add_patch(poly)
    # Add legend for 0.1-0.05
    colors = [nRED[5],nGREEN[5],nBLUE[5]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.1-0.05')
    plt.gca().add_patch(poly)
    # Add legend for <0.05
    colors = [pRED[6],pGREEN[6],pBLUE[6]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ <0.05')
    plt.gca().add_patch(poly)
    colors = [nRED[6],nGREEN[6],nBLUE[6]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- <0.05')
    plt.gca().add_patch(poly)

    box = ax.get_position()
    ax.set_position([box.x0,box.y0+box.height*0.1,box.width,box.height*0.9])
    ax.legend(loc='lower center',prop={'size':5},bbox_to_anchor=(0.48,-0.1),ncol=6)

    minmax_diff = (max_pval-min_pval) 
    minmax_range = np.arange(min_pval,max_pval,minmax_diff/6.0)
    print("min pval = ",min_pval)
    print("max pval = ",max_pval)
    print("Range = ",minmax_range)

    if(canada_counts>0):
        canada_avg = canada_avg/canada_counts
        print("Canada       avg = ",canada_avg,"  counts = ",canada_counts)
    if(siberia_counts>0):
        siberia_avg = siberia_avg/siberia_counts
        print("Siberia      avg = ",siberia_avg,"  counts = ",siberia_counts)
    if(eus_counts>0):
        eus_avg = eus_avg/eus_counts
        print("Eastern US   avg = ",eus_avg,"  counts = ",eus_counts)
    if(wus_counts>0):
        wus_avg = wus_avg/wus_counts
        print("Western US   avg = ",wus_avg,"  counts = ",wus_counts)
    if(ea_counts>0):
        ea_avg = ea_avg/ea_counts
        print("Eastern Asia avg = ",ea_avg,"  counts = ",ea_counts)
    if(wa_counts>0):
        wa_avg = wa_avg/wa_counts
        print("Western Asia avg = ",wa_avg,"  counts = ",wa_counts)
    if(europe_counts>0):
        europe_avg = europe_avg/europe_counts
        print("Europe       avg = ",europe_avg,"  counts = ",europe_counts)

    #start_date = min_date.decode("utf-8")
    #end_date = max_date.decode("utf-8")
    title_string = 'OMI Average Ultraviolet Aerosol Index Trend P-Values\n'+\
                   start_date+' to '+end_date
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    elif(summer is True):
        title_string = title_string+'\nJune, July, August'
    elif(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    elif(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)
    # Set up the colorbar
    #cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    #cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    #color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
    
    plt.xticks(rotation=45,fontsize=6)
    #plt.legend(prop={'size':10})
    fig1.canvas.mpl_connect('button_press_event',onclick)
    if(save is True):
        if(spring is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_spring_newData.png'
        elif(summer is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_summer_newData.png'
        elif(autumn is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_autumn_newData.png'
        elif(winter is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_winter_newData.png'
        else:
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_whole_year_newData.png'
        
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

def plotOMI_Climo(OMI_data,start_date,end_date,save=False,trend_type='standard',file_type='XR123',season='',minlat=30.):
    if(file_type=='NXAR'):
        title_flabel = 'No XTrack and All Rows'
        outname_flabel = 'noX_allRows'
    elif(file_type=='XAR'):
        title_flabel = 'XTrack and All Rows'
        outname_flabel = 'xtrack_allRows'
    elif(file_type=='XR123'):
        title_flabel = 'XTrack and Rows 1-23'
        outname_flabel = 'xtrack_rows1to23'
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False

    # If only summer months are being analyzed, remove all data except 
    # in summer
    if(season=='spring'):
        spring = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    OMI_data[lkey].pop(tkey)


    # Find the lowest lat in the file

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    # Set up the polar stereographic map
    fig1 = plt.figure()
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,         \
                resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    

    cmap = plt.get_cmap('jet')
    #bcmap = plt.cm.set_cmap(cmap)
    if(summer is True):
        v_min = -0.350 # Summer values
        mid_val = 0
        v_max = 0.900
    else:
        v_min = -0.350  # Whole-year values
        mid_val = 0
        v_max = 0.900
    v_min = -1.000  # Whole-year values
    mid_val = 0
    v_max = 1.500
    if(minlat>30.):
        v_min = 0.0
        mid_val = 0
        v_max = 1.00

    # Center the colorbar on zero
    #norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
    norm = Normalize(vmin=v_min,vmax=v_max)
    mapper = ScalarMappable(norm=norm,cmap=cmap)

    # Set up initial values for the analysis regions
    canada_avg = 0.
    canada_counts = 0
    siberia_avg = 0.
    siberia_counts = 0 
    eus_avg = 0.
    eus_counts = 0
    wus_avg = 0.
    wus_counts = 0
    ea_avg = 0.
    ea_counts = 0
    wa_avg = 0.
    wa_counts = 0
    europe_avg = 0.
    europe_counts = 0

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    max_avg_uvai = -10.
    min_avg_uvai = 10.
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = (str(int(lat_ranges[i]))+'x'+                            \
                       str(int(lon_ranges[j])))

            # If no data are present for the curent lat/lon box, fill it with
            # black
            keylist = [ky for ky in OMI_data.keys()]
            if(dictkey not in OMI_data.keys()):
                color=(0,0,0,0)
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(OMI_data[dictkey].keys())==0):
                color=(0,0,0,0)
            else:
                avgs = np.array([OMI_data[dictkey][date]['avg'] for date in \
                    sorted(OMI_data[dictkey].keys())])
                counts = np.array([OMI_data[dictkey][date]['#_obs'] for date in \
                    sorted(OMI_data[dictkey].keys())])
                avg_uvai = sum(avgs*counts)/sum(counts)
                #avg_uvai = np.average(avgs)
                temp_dates = np.array(sorted(OMI_data[dictkey].keys()))
                # Check the current max and min
                #x_vals = np.arange(0,len(OMI_data[dictkey].keys()))
                # Find the slope of the line of best fit for the time series of
                # average data
                #slope, intercept, r_value, p_value, std_err = \
                #    stats.linregress(x_vals,avgs)
                #slope *= len(x_vals)

                if(avg_uvai>max_avg_uvai):
                    max_avg_uvai=avg_uvai
                elif(avg_uvai<min_avg_uvai):
                    min_avg_uvai=avg_uvai
        
                color = mapper.to_rgba(avg_uvai)

                # Add value to analysis regions (if possible)
                if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                    canada_avg+=avg_uvai
                    canada_counts+=1
                elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                    siberia_avg+=avg_uvai
                    siberia_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                    eus_avg+=avg_uvai
                    eus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                    wus_avg+=avg_uvai
                    wus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                    ea_avg+=avg_uvai
                    ea_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                    wa_avg+=avg_uvai
                    wa_counts+=1
                elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                    europe_avg+=avg_uvai
                    europe_counts+=1

            # Find the y coordinates of the current LatxLon box
            y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
            # Find the x coordinates of the current LatxLon box
            x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
            # Convert x and y into map coordinates
            mx, my = m(x,y)
            mxy = zip(mx,my)
            pair1 = (mx[0],my[0])
            pair2 = (mx[1],my[1])
            pair3 = (mx[2],my[2])
            pair4 = (mx[3],my[3])
            color2 = rgb2hex(color)
            #color2 = rgb2hex(color[:2]+color[3:])
            # Plot the box on the map using color
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=color2,         \
                edgecolor=color2)
            plt.gca().add_patch(poly)
    
    #start_date = str(temp_dates[0].decode("utf-8"))
    #end_date = str(temp_dates[-1].decode("utf-8"))
    minmax_diff = (max_avg_uvai-min_avg_uvai) 
    minmax_range = np.arange(min_avg_uvai,max_avg_uvai,minmax_diff/8.0)
    print("min avg_uvai = ",min_avg_uvai)
    print("max avg_uvai = ",max_avg_uvai)
    print("Range = ",minmax_range)

    if(canada_counts>0):
        canada_avg = canada_avg/canada_counts
        print("Canada       avg = ",canada_avg,"  counts = ",canada_counts)
    if(siberia_counts>0):
        siberia_avg = siberia_avg/siberia_counts
        print("Siberia      avg = ",siberia_avg,"  counts = ",siberia_counts)
    if(eus_counts>0):
        eus_avg = eus_avg/eus_counts
        print("Eastern US   avg = ",eus_avg,"  counts = ",eus_counts)
    if(wus_counts>0):
        wus_avg = wus_avg/wus_counts
        print("Western US   avg = ",wus_avg,"  counts = ",wus_counts)
    if(ea_counts>0):
        ea_avg = ea_avg/ea_counts
        print("Eastern Asia avg = ",ea_avg,"  counts = ",ea_counts)
    if(wa_counts>0):
        wa_avg = wa_avg/wa_counts
        print("Western Asia avg = ",wa_avg,"  counts = ",wa_counts)
    if(europe_counts>0):
        europe_avg = europe_avg/europe_counts
        print("Europe       avg = ",europe_avg,"  counts = ",europe_counts)

    start_date = str(start_date)
    end_date = str(end_date)
    title_string = 'OMI Ultraviolet Average Aerosol Index Climatology\n'+      \
        title_flabel+'\n'+start_date+' to '+end_date
        ##start_date+' to '+end_date+'\nXTrack and All Rows'
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    if(summer is True):
        title_string = title_string+'\nJune, July, August'
    if(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    if(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)
    # Set up the colorbar
    cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    plt.xticks(rotation=45,fontsize=6)
    plt.xlabel('Ultraviolet Aerosol Index',fontsize=6)
    fig1.canvas.mpl_connect('button_press_event',onclick_climo)
    if(save is True):
        if(spring is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_spring_'+outname_flabel+'_newData.png'
        elif(summer is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_summer_'+outname_flabel+'_newData.png'
        elif(autumn is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_autumn_'+outname_flabel+'_newData.png'
        elif(winter is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_winter_'+outname_flabel+'_newData.png'
        else:
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_whole_year_'+outname_flabel+'_newData.png'
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()


def plotOMI(OMI_data,start_date,end_date,save=False,trend_type='standard',file_type='XR123',season='',minlat=30.):
    if(file_type=='NXAR'):
        title_flabel = 'No XTrack and All Rows'
        outname_flabel = 'noX_allRows'
    elif(file_type=='XAR'):
        title_flabel = 'XTrack and All Rows'
        outname_flabel = 'xtrack_allRows'
    elif(file_type=='XR123'):
        title_flabel = 'XTrack and Rows 1-23'
        outname_flabel = 'xtrack_rows1to23'
    trend_label=''
    if(trend_type=='theil-sen'):
        trend_label='_theilSen'
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False

    # If only summer months are being analyzed, remove all data except 
    # in summer
    if(season=='spring'):
        spring = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    OMI_data[lkey].pop(tkey)


    # Find the lowest lat in the file
    #lowest_lat = 50.

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    # Set up the polar stereographic map
    fig1 = plt.figure()
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,         \
                resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    

    cmap = plt.get_cmap('bwr')
    #bcmap = plt.cm.set_cmap(cmap)
    if(summer is True):
        v_min = -0.350 # Summer values
        mid_val = 0
        v_max = 0.900
    else:
        v_min = -0.350  # Whole-year values
        mid_val = 0
        v_max = 0.900

    if(minlat>30):
        v_max = 0.6
        mid_val = 0
        v_min = -0.6

    # Center the colorbar on zero
    norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
    #norm = Normalize(vmin=vmin,vmax=vmax)
    mapper = ScalarMappable(norm=norm,cmap=cmap)

    # Set up initial values for the analysis regions
    canada_avg = 0.
    canada_counts = 0
    siberia_avg = 0.
    siberia_counts = 0 
    eus_avg = 0.
    eus_counts = 0
    wus_avg = 0.
    wus_counts = 0
    ea_avg = 0.
    ea_counts = 0
    wa_avg = 0.
    wa_counts = 0
    europe_avg = 0.
    europe_counts = 0


    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    max_slope = -10.
    min_slope = 10.
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = (str(int(lat_ranges[i]))+'x'+                            \
                       str(int(lon_ranges[j])))

            # If no data are present for the curent lat/lon box, fill it with
            # black
            keylist = [ky for ky in OMI_data.keys()]
            if(dictkey not in OMI_data.keys()):
                color=(0,0,0,0)
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(OMI_data[dictkey].keys())==0):
                color=(0,0,0,0)
            else:
                avgs = [OMI_data[dictkey][date]['avg'] for date in \
                    sorted(OMI_data[dictkey].keys())]
                # Check the current max and min
                x_vals = np.arange(0,len(OMI_data[dictkey].keys()))
                temp_dates = np.array(sorted(OMI_data[dictkey].keys()))
                # Find the slope of the line of best fit for the time series of
                # average data
                if(trend_type=='standard'): 
                    slope, intercept, r_value, p_value, std_err = \
                        stats.linregress(x_vals,avgs)
                    slope *= len(x_vals)
                else:
                    #The slope
                    S=0
                    sm=0
                    nx = len(avgs)
                    num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
                    Sn=np.zeros(num_d)
                    for si in range(0,nx-1):
                        for sj in range(si+1,nx):
                            # Find the slope between the two points
                            Sn[sm] = (avgs[si]-avgs[sj])/(si-sj) 
                            sm=sm+1
                        # Endfor
                    # Endfor
                    Snsorted=sorted(Sn)
                    sm=int(num_d/2.)
                    print(dictkey,len(Snsorted))
                    if(len(Snsorted)==1):
                        color=(0,0,0,0) 
                    else:
                        if(2*sm    == num_d):
                            slope=0.5*(Snsorted[sm]+Snsorted[sm+1])
                        if((2*sm)+1 == num_d): 
                            slope=Snsorted[sm+1]
                        slope = slope*len(avgs)

                if(slope>max_slope):
                    max_slope=slope
                elif(slope<min_slope):
                    min_slope=slope
        
                color = mapper.to_rgba(slope)

                # Add value to analysis regions (if possible)
                if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                    canada_avg+=slope
                    canada_counts+=1
                elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                    siberia_avg+=slope
                    siberia_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                    eus_avg+=slope
                    eus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                    wus_avg+=slope
                    wus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                    ea_avg+=slope
                    ea_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                    wa_avg+=slope
                    wa_counts+=1
                elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                    europe_avg+=slope
                    europe_counts+=1
                


            # Find the y coordinates of the current LatxLon box
            y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
            # Find the x coordinates of the current LatxLon box
            x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
            # Convert x and y into map coordinates
            mx, my = m(x,y)
            mxy = zip(mx,my)
            pair1 = (mx[0],my[0])
            pair2 = (mx[1],my[1])
            pair3 = (mx[2],my[2])
            pair4 = (mx[3],my[3])
            color2 = rgb2hex(color)
            #color2 = rgb2hex(color[:2]+color[3:])
            # Plot the box on the map using color
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=color2,         \
                edgecolor=color2)
            plt.gca().add_patch(poly)
    
    #start_date = str(temp_dates[0].decode("utf-8"))
    #end_date = str(temp_dates[-1].decode("utf-8"))
    minmax_diff = (max_slope-min_slope) 
    minmax_range = np.arange(min_slope,max_slope,minmax_diff/8.0)
    print("min slope = ",min_slope)
    print("max slope = ",max_slope)
    print("Range = ",minmax_range)

    if(canada_counts>0):
        canada_avg = canada_avg/canada_counts
        print("Canada       avg = ",canada_avg,"  counts = ",canada_counts)
    if(siberia_counts>0):
        siberia_avg = siberia_avg/siberia_counts
        print("Siberia      avg = ",siberia_avg,"  counts = ",siberia_counts)
    if(eus_counts>0):
        eus_avg = eus_avg/eus_counts
        print("Eastern US   avg = ",eus_avg,"  counts = ",eus_counts)
    if(wus_counts>0):
        wus_avg = wus_avg/wus_counts
        print("Western US   avg = ",wus_avg,"  counts = ",wus_counts)
    if(ea_counts>0):
        ea_avg = ea_avg/ea_counts
        print("Eastern Asia avg = ",ea_avg,"  counts = ",ea_counts)
    if(wa_counts>0):
        wa_avg = wa_avg/wa_counts
        print("Western Asia avg = ",wa_avg,"  counts = ",wa_counts)
    if(europe_counts>0):
        europe_avg = europe_avg/europe_counts
        print("Europe       avg = ",europe_avg,"  counts = ",europe_counts)

    title_string = 'Change in OMI Ultraviolet Average Aerosol Index\n'+        \
        title_flabel+'\n'+str(start_date)+' to '+str(end_date)
        ##start_date+' to '+end_date+'\nXTrack and all rows'
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    if(summer is True):
        title_string = title_string+'\nJune, July, August'
    if(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    if(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)
    # Set up the colorbar
    cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    plt.xticks(rotation=45,fontsize=6)
    fig1.canvas.mpl_connect('button_press_event',onclick)
    if(save is True):
        if(spring is True):
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_spring_'+outname_flabel+'_newData.png'
        elif(summer is True):     
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_summer_'+outname_flabel+'_newData.png'
        elif(autumn is True):     
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_autumn_'+outname_flabel+'_newData.png'
        elif(winter is True):    
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_winter_'+outname_flabel+'_newData.png'
        else:                   
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_whole_year_'+outname_flabel+'_newData.png'
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 
# Plotting functions (NEW)
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# ptype is either 'trend', 'climo', or 'monthclimo'
# title is the plot title
def plotOMI_spatial(pax, plat, plon, pdata, ptype, ptitle = '', plabel = '', \
        vmin = None, vmax = None, colorbar_label_size = 14, minlat = 65.):

    if(vmin == None):
        vmin = np.nanmin(pdata)
    if(vmax == None):
        vmax = np.nanmax(pdata)

    if(ptype == 'trend'):
        colormap = plt.cm.bwr
    else:
        colormap = plt.cm.jet

    # Make copy of OMI_data array
    local_data  = np.copy(pdata)
    # Grid the data, fill in white space
    cyclic_data,cyclic_lons = add_cyclic_point(local_data,plon[0,:])
    plat,plon = np.meshgrid(plat[:,0],cyclic_lons)   
  
    # Mask any missing values
    mask_AI = np.ma.masked_where(cyclic_data < -998.9, cyclic_data)
    mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)

    #pax.gridlines()
    pax.coastlines(resolution='50m')
    mesh = pax.pcolormesh(plon, plat,\
            mask_AI.T,transform = datacrs,\
            cmap = colormap,vmin=vmin,vmax=vmax)
    pax.set_extent([-180,180,minlat,90],datacrs)
    pax.set_boundary(circle, transform=pax.transAxes)
    #ax.set_xlim(-3430748.535086173,3430748.438879491)
    #ax.set_ylim(-3413488.8763307533,3443353.899053069)
    #cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.2),\
    cbar = plt.colorbar(mesh,\
        ax = pax, orientation='vertical',shrink = 0.8, extend = 'both')
    #cbar10.set_label('UV Aerosol Index',weight='bold',fontsize=colorbar_label_size)
    #cbar.ax.tick_params(labelsize=14)
    cbar.set_label(plabel,fontsize=colorbar_label_size,weight='bold')
    pax.set_title(ptitle)


# Designed to work with the netCDF data
def plotOMI_MonthTrend(OMI_data,month_idx=None,save=False,\
        trend_type='standard',minlat=65.,return_trend=False, \
        title = '', label = '', pax = None):
    version = OMI_data['VERSION']
    trend_label=''
    if(trend_type=='theil-sen'):
        trend_label='_theilSen'

    if(month_idx == None):
        month_adder = ''
        do_month = False
        v_max = 0.5
        v_min = -0.5
    else:
        month_adder = '_month'
        do_month = True
        if(version == 'V003'):
            index_jumper = 12
        else:
            index_jumper = 6 
        v_max = 0.7
        v_min = -0.7
    # --------------------------------------------------------------
    #
    # Use calcOMI_grid_trend to calculate the trends in the AI data
    #
    # --------------------------------------------------------------
    ai_trends = calcOMI_grid_trend(OMI_data, month_idx, trend_type, \
        minlat)

    # --------------------------------------------------------------
    #
    # Plot the calculated trends on a figure
    #
    # --------------------------------------------------------------
    colormap = plt.cm.bwr

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(OMI_data['DATES'][month_idx::index_jumper][0],'%Y%m')
    end_date   = datetime.strptime(OMI_data['DATES'][month_idx::index_jumper][-1],'%Y%m')
    
    # Make figure title
    month_string = ''
    if(do_month == True):
        month_string = start_date.strftime('%B') + ' '

    if(title == ''):
        title = 'OMI AI ' + month_string + 'Trends ('+version+')\n'+\
            start_date.strftime('%b. %Y') + ' - ' + end_date.strftime('%b. %Y')
    if(label == ''):
        label = 'UV Aerosol Index Trend'

    # Call plotOMI_spatial to add the data to the figure
    # Make figure
    if(pax is None):
        # Make figure
        plt.close('all')
        fig1 = plt.figure(figsize = (8,8))
        ax = plt.axes(projection = mapcrs)

        plotOMI_spatial(ax, OMI_data['LAT'], OMI_data['LON'], ai_trends, 'trend', \
            ptitle = title, plabel = label, \
            vmin = v_min, vmax = v_max, minlat = minlat)

        fig1.tight_layout()

        if(save == True):
            month_adder = ''
            if(do_month == True):
                month_adder = '_' + start_date.strftime('%b') 
            out_name = 'omi_ai_trend' + month_adder + '_' + \
                start_date.strftime('%Y%m') + '_' + end_date.strftime('%Y%m') + \
                '_' + version + '.png'
            plt.savefig(out_name,dpi=300)
            print("Saved image",out_name)
        else:
            plt.show()
    else:
        plotOMI_spatial(pax, OMI_data['LAT'], OMI_data['LON'], ai_trends, 'trend', \
            ptitle = title, plabel = label, \
            vmin = v_min, vmax = v_max, minlat = minlat)

    #if(return_trend == True):
    #    return ai_trends


# Plots a single month of OMI climatology data (assumed to be from the 
# netCDF file).
# May work as a stand-alone function call or to plot on a premade
# axis by specifying the 'pax' argument
def plotOMI_NCDF_SingleMonth(OMI_data,time_idx,minlat=65, pax = None, \
        save=False):

    version = OMI_data['VERSION']
    if((version == 'VSJ2') | (version == 'VSJ4')):
        data_type = '(Perturbation)'
        label_adder = 'perturbation'
    else:
        data_type = '(Screened)'
        label_adder = ''

    # Make copy of OMI_data array
    local_data  = np.copy(OMI_data['AI'][time_idx,:,:])
    local_count = np.copy(OMI_data['OB_COUNT'][time_idx,:,:])

    # Mask any missing values
    mask_AI = np.ma.masked_where(local_count == 0, local_data)
    mask_AI = np.ma.masked_where(mask_AI == -999.9, mask_AI)

    # Mask any data below the threshold latitude
    mask_AI = np.ma.masked_where(OMI_data['LAT'] < minlat,mask_AI)

    # Make figure title
    first_date = OMI_data['DATES'][time_idx]
    title = 'OMI AI '+data_type+'\n'+label_dict[version]+'\n'+first_date

    if(pax is None):
        # Make figure
        plt.close('all')
        fig1 = plt.figure(figsize = (8,8))
        ax = plt.axes(projection = mapcrs)

        plotOMI_spatial(ax, OMI_data['LAT'], OMI_data['LON'], mask_AI, 'climo',\
            ptitle = title, plabel = 'UV Aerosol Index', \
            vmin = -1.0, vmax = 1.5, minlat = minlat)

        fig1.tight_layout()

        if(save == True):
            outname = 'omi_ai_single_month_' + first_date + '_'+version+'.png'
            plt.savefig(outname,dpi=300)
            print("Saved image",outname)
        else:
            plt.show()

    else:
        plotOMI_spatial(pax, OMI_data['LAT'], OMI_data['LON'], mask_AI, 'climo',\
            ptitle = title, plabel = 'UV Aerosol Index', \
            vmin = -1.0, vmax = 1.5, minlat = minlat)
    ##!#ax.gridlines()
    ##!#ax.coastlines(resolution='50m')
    ##!#mesh = ax.pcolormesh(OMI_data['LON'], OMI_data['LAT'],mask_AI,transform = datacrs,cmap = colormap,\
    ##!#        vmin = -1.0, vmax = 1.5)
    ##!#ax.set_boundary(circle, transform=ax.transAxes)
    ##!#ax.set_extent([-180,180,minlat,90],datacrs)
    ##!##ax.set_xlim(-3430748.535086173,3430748.438879491)
    ##!##ax.set_ylim(-3413488.8763307533,3443353.899053069)
    ##!#cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
    ##!#    aspect=50,shrink = 0.845)
    ##!#cbar.ax.tick_params(labelsize=14)
    ##!#cbar.set_label('UV Aerosol Index'+label_adder,fontsize=16,weight='bold')
    ##!#ax.set_title(title)


# Plots a total OMI climatology (assumed to be from the 
# netCDF file).
# May work as a stand-alone function call or to plot on a premade
# axis by specifying the 'pax' argument
def plotOMI_NCDF_Climo(OMI_data,start_idx=0,end_idx=None,season = '',minlat=60,\
                       title = '', pax = None, save=False):

    version = OMI_data['VERSION']

    season_dict = {
        'spring': '\nMAM',\
        'summer': '\nJJA',\
        'autumn': '\nSON',\
        'winter': '\nDJF',\
        '': ''
    }

    if(end_idx == None):
        end_idx = len(OMI_data['MONTH'])   
 
    # Make copy of OMI_data array
    local_data = np.copy(OMI_data['AI'])

    #start_date = datetime(year=2004,month=10,day=1)
    start_date = datetime.strptime(OMI_data['DATES'][0],"%Y%m")

    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False

    month_objects = []
    keepers = np.arange(1,13)
    if(season=='spring'):
        spring = True 
        keepers = [3,4,5]
    elif(season=='summer'):
        summer = True 
        keepers = [6,7,8]
    elif(season=='autumn'):
        autumn = True 
        keepers = [9,10,11]
    elif(season=='winter'):
        winter = True 
        keepers = [12,1,2]
   
    for m_idx in range(len(OMI_data['MONTH'])):
        new_date = datetime.strptime(OMI_data['DATES'][m_idx],'%Y%m')
        #new_date = start_date + relativedelta(months=m_idx)
        if(new_date.month not in keepers):  
            local_data[m_idx,:,:] = -999.9
        else:
            month_objects.append(new_date)

    # Mask any missing values
    mask_AI = np.ma.masked_where(((OMI_data['OB_COUNT'] == -99) | \
                (OMI_data['OB_COUNT'] == 0)), local_data)
    mask_AI = np.ma.masked_where(local_data == -999.9, mask_AI)

    # Calculate climatology between desired indices
    OMI_climo = np.nanmean(mask_AI[start_idx:end_idx,:,:],axis=0)

    print(np.min(OMI_climo),np.max(OMI_climo))

    # Make figure title
    first_date = month_objects[0].strftime("%Y%m")
    last_date = month_objects[-1].strftime("%Y%m")
    print(month_objects[0].strftime("%Y%m"),month_objects[-1].strftime("%Y%m"))
    if(title == ''):
        title = 'OMI AI Climatology\n'+first_date + ' - ' + last_date + \
            season_dict[season]

    # Make figure
    if(pax is None):
        plt.close('all')
        fig1 = plt.figure(figsize = (8,8))

        fig1.tight_layout()

        ax = plt.axes(projection = mapcrs)
        plotOMI_spatial(ax, OMI_data['LAT'], OMI_data['LON'], OMI_climo, 'climo',\
            ptitle = title, plabel = 'UV Aerosol Index', \
            vmin = -1.0, vmax = 1.0, minlat = minlat)

        if(save == True):
            season_adder = ''
            if(len(season.strip()) != 0):
                season_adder = '_' + season.strip()
            outname = 'omi_ai_climo_' + first_date + '_' + last_date + season_adder + \
                '_' + version + '.png'
            plt.savefig(outname,dpi=300)
            print("Saved image",outname)
        else:
            plt.show()
    else: 
        plotOMI_spatial(pax, OMI_data['LAT'], OMI_data['LON'], OMI_climo, 'climo',\
            ptitle = title, plabel = 'UV Aerosol Index', \
            vmin = -1.0, vmax = 1.0, minlat = minlat)
    ##!#ax.gridlines()
    ##!#ax.coastlines(resolution='50m')
    ##!#mesh = ax.pcolormesh(OMI_data['LON'], OMI_data['LAT'],OMI_climo,transform = datacrs,cmap = colormap,\
    ##!#        vmin = -1.0, vmax = 1.0)
    ##!#ax.set_extent([-180,180,minlat,90],datacrs)
    ##!#ax.set_boundary(circle, transform=ax.transAxes)
    ##!##ax.set_xlim(-3430748.535086173,3430748.438879491)
    ##!##ax.set_ylim(-3413488.8763307533,3443353.899053069)
    ##!##cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
    ##!##    aspect=50,shrink = 0.905,label='Aerosol Index')
    ##!#cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
    ##!#    aspect=50,shrink = 0.845)
    ##!#cbar.ax.tick_params(labelsize=14)
    ##!#cbar.set_label('UV Aerosol Index',fontsize=16,weight='bold')
    ##!#ax.set_title(title)

# NOTE: As of 11/20/2021, this function has not been updated to work with
# the "new" netCDF OMI data VJZ211 and VSJ4, which only contain data between
# April and September.
# The same as the above function, but plots all 4 seasons in a single 4-panel
# plot.
#def plotOMI_NCDF_Climo_FourPanel(OMI_data,start_idx=0,end_idx=169,minlat=60,\
def plotOMI_NCDF_Climo_SpringSummer(OMI_data,start_idx=0,end_idx=96,minlat=65,\
                       save=False):

    ##!## Make copy of OMI_data array
    ##!#local_spring = np.copy(OMI_data['AI'])
    ##!#local_summer = np.copy(OMI_data['AI'])
    ##!#local_autumn = np.copy(OMI_data['AI'])
    ##!#local_winter = np.copy(OMI_data['AI'])

    ##!#start_date = datetime.strptime(OMI_data['DATES'][0],"%Y%m")
    ##!##start_date = datetime(year=2004,month=10,day=1)

    ##!#month_objects = []
    ##!#spring_keepers = [3,4,5]
    ##!#summer_keepers = [6,7,8]
    ##!#autumn_keepers = [9,10,11]
    ##!#winter_keepers = [12,1,2]
   
    ##!#for m_idx in OMI_data['MONTH']:
    ##!#    new_date = start_date + relativedelta(months=m_idx)
    ##!#    if(new_date.month not in spring_keepers):  
    ##!#        local_spring[m_idx,:,:] = -999.9
    ##!#    if(new_date.month not in summer_keepers):  
    ##!#        local_summer[m_idx,:,:] = -999.9
    ##!#    if(new_date.month not in autumn_keepers):  
    ##!#        local_autumn[m_idx,:,:] = -999.9
    ##!#    if(new_date.month not in winter_keepers):  
    ##!#        local_winter[m_idx,:,:] = -999.9

    ##!## Mask any missing values
    ##!#mask_spring_AI = np.ma.masked_where(local_spring == -999.9, local_spring)
    ##!#mask_summer_AI = np.ma.masked_where(local_summer == -999.9, local_summer)
    ##!#mask_autumn_AI = np.ma.masked_where(local_autumn == -999.9, local_autumn)
    ##!#mask_winter_AI = np.ma.masked_where(local_winter == -999.9, local_winter)

    ##!## Calculate climatology between desired indices
    ##!#OMI_spring_climo = np.nanmean(mask_spring_AI[start_idx:end_idx,:,:],axis=0)
    ##!#OMI_summer_climo = np.nanmean(mask_summer_AI[start_idx:end_idx,:,:],axis=0)
    ##!#OMI_autumn_climo = np.nanmean(mask_autumn_AI[start_idx:end_idx,:,:],axis=0)
    ##!#OMI_winter_climo = np.nanmean(mask_winter_AI[start_idx:end_idx,:,:],axis=0)

    # Make figure title
    first_date = datetime.strptime(OMI_data['DATES'][0], "%Y%m")
    last_date  = datetime.strptime(OMI_data['DATES'][-1], "%Y%m")
    #print(month_objects[0].strftime("%Y%m"),month_objects[-1].strftime("%Y%m"))
    title = 'OMI AI Seasonal Climatology\n'+first_date.strftime("%b. %Y") + \
        ' - ' + last_date.strftime("%b. %Y")


    # Make figure
    plt.close('all')
    fig1 = plt.figure(1,figsize=(10,5))
    ax0 = fig1.add_subplot(1,2,1,projection = mapcrs)
    ax1 = fig1.add_subplot(1,2,2,projection = mapcrs)
    #gs = gridspec.GridSpec(nrows=1, ncols=2)
    #gs = gridspec.GridSpec(nrows=2, ncols=2, hspace = 0.10, wspace = 0.06)

    plt.suptitle(title)


    vmax_r = 1.5
    vmin_r = -1.0
    # Make spring plot
    # ----------------
    #ax0 = plt.subplot(gs[0,0],projection=mapcrs)
    plotOMI_NCDF_Climo(OMI_data,start_idx=0,end_idx=None,season = 'spring',\
        title = 'MAM', minlat=minlat,pax = ax0, save=False)

    ##!#ax0.set_extent([-180,180,60,90],ccrs.PlateCarree())
    ##!#ax0.gridlines()
    ##!#ax0.coastlines(resolution='50m')
    ##!#mesh0 = ax0.pcolormesh(OMI_data['LON'], OMI_data['LAT'],OMI_spring_climo,\
    ##!#        transform = datacrs,cmap = colormap,vmin = vmin_r, vmax = vmax_r)
    ##!#ax0.set_boundary(circle, transform=ax0.transAxes)
    ##!#cbar0 = plt.colorbar(mesh0,ticks = np.arange(-2.0,4.1,0.5), \
    ##!#    orientation='horizontal',pad=0,aspect=50,shrink = 0.905,\
    ##!#    label='Aerosol Index')
    ##!#ax0.text(0., 1.02, 'A', transform = ax0.transAxes, size=15, weight = 'bold')
    ##!#ax0.set_title('Spring (MAM)')

    # Make summer plot
    # ----------------
    #ax1 = plt.subplot(gs[0,1],projection=mapcrs)
    plotOMI_NCDF_Climo(OMI_data,start_idx=0,end_idx=None,season = 'summer',\
        title = 'JJA', minlat=minlat, pax = ax1, save=False)
    ##!#ax1 = plt.subplot(gs[0,1],projection=mapcrs)
    ##!#ax1.set_extent([-180,180,60,90],ccrs.PlateCarree())
    ##!#ax1.gridlines()
    ##!#ax1.coastlines(resolution='50m')
    ##!#mesh1 = ax1.pcolormesh(OMI_data['LON'], OMI_data['LAT'],OMI_summer_climo,\
    ##!#        transform = datacrs,cmap = colormap,vmin = vmin_r, vmax = vmax_r)
    ##!#ax1.set_boundary(circle, transform=ax1.transAxes)
    ##!#cbar1 = plt.colorbar(mesh1,ticks = np.arange(-2.0,4.1,0.5), \
    ##!#    orientation='horizontal',pad=0,aspect=50,shrink = 0.905,\
    ##!#    label='Aerosol Index')
    ##!#ax1.text(0., 1.02, 'B', transform = ax1.transAxes, size=15, weight = 'bold')
    ##!#ax1.set_title('Summer (JJA)')

    ##!## Make autumn plot
    ##!## ----------------
    ##!#ax2 = plt.subplot(gs[1,0],projection=mapcrs)
    ##!#ax2.set_extent([-180,180,60,90],ccrs.PlateCarree())
    ##!#ax2.gridlines()
    ##!#ax2.coastlines(resolution='50m')
    ##!#mesh2 = ax2.pcolormesh(OMI_data['LON'], OMI_data['LAT'],OMI_autumn_climo,\
    ##!#        transform = datacrs,cmap = colormap,vmin = vmin_r, vmax = vmax_r)
    ##!#ax2.set_boundary(circle, transform=ax2.transAxes)
    ##!#cbar2 = plt.colorbar(mesh2,ticks = np.arange(-2.0,4.1,0.5), \
    ##!#    orientation='horizontal',pad=0,aspect=50,shrink = 0.905,\
    ##!#    label='Aerosol Index')
    ##!#ax2.text(0., 1.02, 'C', transform = ax2.transAxes, size=15, weight = 'bold')
    ##!#ax2.set_title('Autumn (SON)')

    ##!## Make winter plot
    ##!## ----------------
    ##!#ax3 = plt.subplot(gs[1,1],projection=mapcrs)
    ##!#ax3.set_extent([-180,180,60,90],ccrs.PlateCarree())
    ##!#ax3.gridlines()
    ##!#ax3.coastlines(resolution='50m')
    ##!#mesh3 = ax3.pcolormesh(OMI_data['LON'], OMI_data['LAT'],OMI_winter_climo,\
    ##!#        transform = datacrs,cmap = colormap,vmin = vmin_r, vmax = vmax_r)
    ##!#ax3.set_boundary(circle, transform=ax3.transAxes)
    ##!#cbar3 = plt.colorbar(mesh3,ticks = np.arange(-2.0,4.1,0.5), \
    ##!#    orientation='horizontal',pad=0,aspect=50,shrink = 0.905,\
    ##!#    label='Aerosol Index')
    ##!#ax3.text(0., 1.02, 'D', transform = ax3.transAxes, size=15, weight = 'bold')
    ##!#ax3.set_title('Winter (DJF)')

    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')

    fig1.tight_layout()

    if(save == True):
        #outname = 'omi_ai_climo_fourpanel_' + first_date.strftime("%Y%m") \
        outname = 'omi_ai_climo_twopanel_' + first_date.strftime("%Y%m") \
            + '_' + last_date.strftime("%Y%m") + '_' + OMI_data['VERSION'] + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Plot a monthly climatology 
def plotOMI_MonthClimo(OMI_data,month_idx,minlat = 60, pax = None, \
        save=False):

    version = OMI_data['VERSION']

    if(version == 'VSJ2'):
        max_AI = 0.5
        min_AI = -0.5
        colormap = plt.cm.bwr
    else:
        max_AI = 1.5
        min_AI = -1.5
        colormap = plt.cm.jet

    min_AI = -1.0
    max_AI = 1.0

    # Make copy of OMI_data array
    if('MONTH_CLIMO' not in OMI_data.keys()):
        OMI_data = calcOMI_MonthClimo(OMI_data)

    local_data  = np.copy(OMI_data['MONTH_CLIMO'][month_idx,:,:])
    # Mask any missing values
    mask_AI = np.ma.masked_where(local_data == -999.9, local_data)
    mask_AI = np.ma.masked_where(OMI_data['LAT'] < minlat, mask_AI)

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(OMI_data['DATES'][month_idx],'%Y%m')
    end_date   = datetime.strptime(OMI_data['DATES'][-1],'%Y%m')

    # Make figure title
    #date_month = datetime(year = 1,month = month_idx+1, day = 1).strftime('%B')
    date_month = start_date.strftime('%B')
    title = 'OMI AI ' + date_month + ' Climatology ('+version+')\n'+\
        start_date.strftime('%b. %Y') + ' - ' + end_date.strftime('%b. %Y')

    if(pax is None):
        # Make figure
        plt.close('all')
        fig1 = plt.figure(figsize = (8,8))
        ax = plt.axes(projection = mapcrs)

        plotOMI_spatial(ax, OMI_data['LAT'], OMI_data['LON'], mask_AI, 'climo', \
            ptitle = title, plabel = 'UV Aerosol Index', \
            vmin = min_AI, vmax = max_AI, minlat = minlat)

        fig1.tight_layout()

        if(save == True):
            out_name = 'omi_ai_month_climo_' + date_month + '_' + version + '.png'
            plt.savefig(out_name,dpi=300)
            print("Saved image",out_name)
        else:
            plt.show()

    else:
        plotOMI_spatial(pax, OMI_data['LAT'], OMI_data['LON'], mask_AI, 'climo', \
            ptitle = title, plabel = 'UV Aerosol Index', \
            vmin = min_AI, vmax = max_AI, minlat = minlat)

    ##!#ax.gridlines()
    ##!#ax.coastlines(resolution='50m')
    ##!#mesh = ax.pcolormesh(OMI_data['LON'], OMI_data['LAT'],\
    ##!#        mask_AI,transform = datacrs,\
    ##!#        cmap = colormap, vmin = min_AI, vmax = max_AI)
    ##!#ax.set_extent([-180,180,minlat,90],datacrs)
    ##!#ax.set_boundary(circle, transform=ax.transAxes)
    ##!##ax.set_xlim(-3430748.535086173,3430748.438879491)
    ##!##ax.set_ylim(-3413488.8763307533,3443353.899053069)
    ##!#cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
    ##!#    aspect=50,shrink = 0.845)
    ##!#cbar.ax.tick_params(labelsize=14)
    ##!#cbar.set_label('UV Aerosol Index',fontsize=16,weight='bold')
    ##!#ax.set_title(title)


# Generate a six-panel figure comparing the climatology and trend between 3
# versions of the OMI data.
def plotOMI_Compare_ClimoTrend(OMI_data1,OMI_data2,OMI_data3,month_idx,\
        trend_type = 'standard', minlat=65,save=False):

    colormap = plt.cm.jet

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    if(OMI_data1['VERSION'] == 'V003'):
        index_jumper = 12
    else:
        index_jumper = 6 

    # Get the labels figured out
    new_label_dict = {
        'VBS1': 'Control',
        'VBS0': 'Control',
        'VJZ29': 'Screening Method',
        'VJZ211': 'Screening Method',
        'VSJ2': 'Perturbation Method',
        'VSJ4': 'Perturbation Method'
    } 
    colorbar_label_size = 13
    axis_title_size = 14.5
    row_label_size = 14.5

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(OMI_data1['DATES'][month_idx],'%Y%m')

    #fig = plt.figure()
    plt.close('all')
    fig = plt.figure(figsize=(16,10))
    plt.suptitle('OMI Comparisons: '+start_date.strftime("%B"),y=0.95,\
        fontsize=18,fontweight=4,weight='bold')
    gs = gridspec.GridSpec(nrows=2, ncols=3, hspace = 0.001, wspace = 0.15)

    # - - - - - - - - - - - - - - - - - - - - -
    # Plot the climatologies along the top row
    # - - - - - - - - - - - - - - - - - - - - -
       
    # Plot DATA1 climos
    # -----------------
    ##!## Make copy of OMI_data array
    local_data1  = np.copy(OMI_data1['MONTH_CLIMO'][month_idx,:,:])
    local_data2  = np.copy(OMI_data2['MONTH_CLIMO'][month_idx,:,:])
    local_data3  = np.copy(OMI_data3['MONTH_CLIMO'][month_idx,:,:])
    mask_AI1 = np.ma.masked_where(local_data1 == -999.9, local_data1)
    mask_AI2 = np.ma.masked_where(local_data2 == -999.9, local_data2)
    mask_AI3 = np.ma.masked_where(local_data3 == -999.9, local_data3)
    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,OMI_data1['LON'][0,:])
    ##!#plat,plon = np.meshgrid(OMI_data1['LAT'][:,0],cyclic_lons)   
  
    # Mask any missing values
    #mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)
    ax00 = plt.subplot(gs[0,0],projection=mapcrs)
    ax01 = plt.subplot(gs[0,1],projection=mapcrs)
    ax02 = plt.subplot(gs[0,2],projection=mapcrs)
    plotOMI_spatial(ax00, OMI_data1['LAT'], OMI_data1['LON'], mask_AI1, \
        'climo', ptitle = new_label_dict[OMI_data1['VERSION']] + '\n', \
        plabel = 'UV Aerosol Index', \
        vmin = -1.0, vmax = 1.0, minlat = minlat)
    print("here, ",new_label_dict[OMI_data2['VERSION']] + '\n')
    plotOMI_spatial(ax01, OMI_data2['LAT'], OMI_data2['LON'], mask_AI2, \
        'climo', ptitle = new_label_dict[OMI_data2['VERSION']] + '\n', \
        plabel = 'UV Aerosol Index', \
        vmin = -1.0, vmax = 1.0, minlat = minlat)
    plotOMI_spatial(ax02, OMI_data3['LAT'], OMI_data3['LON'], mask_AI3, \
        'climo', ptitle = new_label_dict[OMI_data3['VERSION']] + '\n', \
        plabel = 'UVAI Perturbation', \
        vmin = -1.0, vmax = 1.0, minlat = minlat)

    # Plot trends
    ax10 = plt.subplot(gs[1,0],projection=mapcrs)
    ax11 = plt.subplot(gs[1,1],projection=mapcrs)
    ax12 = plt.subplot(gs[1,2],projection=mapcrs)
    plotOMI_MonthTrend(OMI_data1,month_idx=month_idx,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = None, pax = ax10)
    ax01.set_title(OMI_data1['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data2,month_idx=month_idx,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = None, pax = ax11)
    plotOMI_MonthTrend(OMI_data3,month_idx=month_idx,\
        trend_type=trend_type,label = 'AI Pert. Trend (AI/Study Period)',\
        minlat=65.,title = None, pax = ax12)

    plot_subplot_label(ax00, '(a)')
    plot_subplot_label(ax01, '(b)')
    plot_subplot_label(ax02, '(c)')
    plot_subplot_label(ax10, '(d)')
    plot_subplot_label(ax11, '(e)')
    plot_subplot_label(ax12, '(f)')

    ##!#ax00.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    ##!#ax00.set_boundary(circle, transform=ax00.transAxes)
    ##!##ax00.gridlines()
    ##!#ax00.coastlines(resolution='50m')
    ##!#mesh00 = ax00.pcolormesh(plon,plat,mask_AI.T,\
    ##!#        transform=ccrs.PlateCarree(),vmin=-1.0,vmax=1.0,cmap=colormap)
    ##!#ax00.set_title(new_label_dict[OMI_data1['VERSION']]+'\n',weight='bold',\
    ##!#    fontsize=axis_title_size)
    ##!#cbar00 = plt.colorbar(mesh00,ticks = np.arange(-2.0,4.1,0.5),orientation='vertical',\
    ##!#    extend='both',shrink=0.8)
    ##!#cbar00.set_label('UV Aerosol Index',weight='bold',fontsize=colorbar_label_size)
    fig.text(0.10, 0.70, 'Climatology', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size)
    #cbar00 = plt.colorbar(mesh00,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
    #    aspect=50,shrink = 0.845)
    #cbar00.ax.tick_params(labelsize=14)
    #cbar00.set_label('UV Aerosol Index',fontsize=16,weight='bold')

    ##!## Plot DATA2 climos
    ##!## -----------------
    ##!## Make copy of OMI_data array
    ##!#local_data  = np.copy(OMI_data2['MONTH_CLIMO'][month_idx,:,:])
    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,OMI_data2['LON'][0,:])
    ##!#plat,plon = np.meshgrid(OMI_data2['LAT'][:,0],cyclic_lons)   
  
    ##!## Mask any missing values
    ##!#mask_AI = np.ma.masked_where(cyclic_data == -999.9, cyclic_data)
    ##!#mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)

    ##!#ax10 = plt.subplot(gs[0,1],projection=mapcrs)
    ##!#ax10.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    ##!#ax10.set_boundary(circle, transform=ax10.transAxes)
    ##!##ax10.gridlines()
    ##!#ax10.coastlines(resolution='50m')
    ##!#mesh10 = ax10.pcolormesh(plon,plat,mask_AI.T,\
    ##!#        transform=ccrs.PlateCarree(),vmin=-1.0,vmax=1.0,cmap=colormap)
    ##!#ax10.set_title(new_label_dict[OMI_data2['VERSION']]+'\n',weight='bold',\
    ##!#        fontsize=axis_title_size)
    ##!#cbar10 = plt.colorbar(mesh10,ticks = np.arange(-2.0,4.1,0.5),orientation='vertical',\
    ##!#    extend='both',shrink=0.8)
    ##!#cbar10.set_label('UV Aerosol Index',weight='bold',fontsize=colorbar_label_size)
    ##!##cbar10 = plt.colorbar(mesh10,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
    ##!##    aspect=50,shrink = 0.845)
    ##!##cbar10.ax.tick_params(labelsize=14)
    ##!##cbar10.set_label('UV Aerosol Index',fontsize=16,weight='bold')
    ##!# 
    ##!## Plot DATA3 climatology 
    ##!## -----------------
    ##!## Make copy of OMI_data array
    ##!#local_data  = np.copy(OMI_data3['MONTH_CLIMO'][month_idx,:,:])
    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,OMI_data3['LON'][0,:])
    ##!#plat,plon = np.meshgrid(OMI_data3['LAT'][:,0],cyclic_lons)   
  
    ##!## Mask any missing values
    ##!#mask_AI = np.ma.masked_where(cyclic_data == -999.9, cyclic_data)
    ##!#mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)

    ##!#ax20 = plt.subplot(gs[0,2],projection=mapcrs)
    ##!#ax20.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    ##!#ax20.set_boundary(circle, transform=ax20.transAxes)
    ##!##ax20.gridlines()
    ##!#ax20.coastlines(resolution='50m')
    ##!#mesh20 = ax20.pcolormesh(plon,plat,mask_AI.T,\
    ##!#        transform=ccrs.PlateCarree(),vmin=-1.0,vmax=1.0,cmap=colormap)
    ##!#ax20.set_title(new_label_dict[OMI_data3['VERSION']]+'\n',weight='bold',\
    ##!#    fontsize=axis_title_size)
    ##!#cbar20 = plt.colorbar(mesh20,ticks = np.arange(-2.0,4.1,0.5),orientation='vertical',\
    ##!#    extend='both',shrink=0.8)
    ##!##cbar20.ax.tick_params(labelsize=14)
    ##!#cbar20.set_label('AI Perturbation',weight='bold',fontsize=colorbar_label_size)

    ##!## Plot Trends
    ##!## -----------
    ##!#colormap = plt.cm.bwr

    ##!## Plot DATA1 trends
    ##!## Make copy of OMI_data array
    ##!#local_data   = np.copy(OMI_data1['AI'][month_idx::index_jumper,:,:])
    ##!#local_counts = np.copy(OMI_data1['OB_COUNT'][month_idx::index_jumper,:,:])
    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,OMI_data1['LON'][0,:])
    ##!#cyclic_counts,cyclic_lons = add_cyclic_point(local_counts,OMI_data1['LON'][0,:])
    ##!#plat,plon = np.meshgrid(OMI_data1['LAT'][:,0],cyclic_lons)   
  
    ##!#local_mask = np.ma.masked_where(cyclic_counts == 0, cyclic_data)
    ##!#local_mask = np.ma.masked_where(local_mask == -999.9, local_mask)
    ##!#ai_trends = np.zeros(cyclic_data.shape[1:])
    ##!## Loop over all the keys and print the regression slopes 
    ##!## Grab the averages for the key
    ##!#for i in range(0,len(lat_ranges)):
    ##!#    for j in range(0,len(lon_ranges)):
    ##!#        # Check the current max and min
    ##!#        x_vals = np.arange(0,len(local_mask[:,i,j]))
    ##!#        # Find the slope of the line of best fit for the time series of
    ##!#        # average data
    ##!#        if(trend_type=='standard'): 
    ##!#            slope, intercept, r_value, p_value, std_err = \
    ##!#                stats.linregress(x_vals,local_mask[:,i,j])
    ##!#            ai_trends[i,j] = slope * len(x_vals)
    ##!#        else:
    ##!#            print("Ignoring MK trend for now")

   
    ##!#ai_trends = np.ma.masked_where(plat.T < minlat, ai_trends)

    ##!#ax01 = plt.subplot(gs[1,0],projection=mapcrs)
    ##!#ax01.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    ##!#ax01.set_boundary(circle, transform=ax01.transAxes)
    fig.text(0.10, 0.30, 'Trend', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size)
    ##!##ax01.gridlines()
    ##!#ax01.coastlines(resolution='50m')
    ##!#mesh01 = ax01.pcolormesh(plon,plat,ai_trends.T,\
    ##!#        transform=ccrs.PlateCarree(),vmin=-0.7,vmax=0.7,cmap=colormap)
    ##!##ax01.set_title(OMI_data1['VERSION']+ '\n\n')
    ##!#cbar01 = plt.colorbar(mesh01,ticks = np.arange(-1.0,1.1,0.2),orientation='vertical',\
    ##!#    extend='both',shrink=0.8)
    ##!#cbar01.set_label('AI Trend (AI/Study Period)',weight='bold',\
    ##!#    fontsize=colorbar_label_size)
    ##!##cbar01 = plt.colorbar(mesh01,ticks = np.arange(-2.0,4.1,0.2),orientation='horizontal',pad=0,\
    ##!##    aspect=50,shrink = 0.845)
    ##!##cbar01.ax.tick_params(labelsize=14)
    ##!##cbar01.set_label('UV Aerosol Index',fontsize=16,weight='bold')
    ##!# 
    ##!## Plot DATA2 trends
    ##!## -----------------
    ##!## Make copy of OMI_data array
    ##!#local_data   = np.copy(OMI_data2['AI'][month_idx::index_jumper,:,:])
    ##!#local_counts = np.copy(OMI_data2['OB_COUNT'][month_idx::index_jumper,:,:])
    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,OMI_data2['LON'][0,:])
    ##!#cyclic_counts,cyclic_lons = add_cyclic_point(local_counts,OMI_data2['LON'][0,:])
    ##!#plat,plon = np.meshgrid(OMI_data2['LAT'][:,0],cyclic_lons)   

    ##!#local_mask = np.ma.masked_where(cyclic_counts == 0, cyclic_data)
    ##!#local_mask = np.ma.masked_where(local_mask == -999.9, local_mask)
    ##!#ai_trends = np.zeros(cyclic_data.shape[1:])
    ##!## Loop over all the keys and print the regression slopes 
    ##!## Grab the averages for the key
    ##!#for i in range(0,len(np.arange(np.min(OMI_data2['LAT']),90))):
    ##!#    for j in range(0,len(lon_ranges)):
    ##!#        # Check the current max and min
    ##!#        x_vals = np.arange(0,len(local_mask[:,i,j]))
    ##!#        # Find the slope of the line of best fit for the time series of
    ##!#        # average data
    ##!#        if(trend_type=='standard'): 
    ##!#            slope, intercept, r_value, p_value, std_err = \
    ##!#                stats.linregress(x_vals,local_mask[:,i,j])
    ##!#            ai_trends[i,j] = slope * len(x_vals)
    ##!#        else:
    ##!#            print("Ignoring MK trend for now")

   
    ##!#ai_trends = np.ma.masked_where(plat.T < minlat, ai_trends)

    ##!#ax11 = plt.subplot(gs[1,1],projection=mapcrs)
    ##!#ax11.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    ##!#ax11.set_boundary(circle, transform=ax11.transAxes)
    ##!##ax11.gridlines()
    ##!#ax11.coastlines(resolution='50m')
    ##!#mesh11 = ax11.pcolormesh(plon,plat,ai_trends.T,\
    ##!#        transform=ccrs.PlateCarree(),vmin=-0.7,vmax=0.7,cmap=colormap)
    ##!##ax11.set_title(OMI_data2['VERSION']+ ' Trend')
    ##!#cbar11 = plt.colorbar(mesh11,ticks = np.arange(-1.0,1.1,0.2),orientation='vertical',\
    ##!#    extend='both',shrink=0.8)
    ##!#cbar11.set_label('AI Trend (AI/Study Period)',weight='bold',\
    ##!#    fontsize=colorbar_label_size)
    ##!##$cbar11 = plt.colorbar(mesh11,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
    ##!##$    aspect=50,shrink = 0.845)
    ##!##$cbar11.ax.tick_params(labelsize=14)
    ##!##$cbar11.set_label('UV Aerosol Index',fontsize=16,weight='bold')
    ##!# 
    ##!## - - - - - - - - - - - - - - - - -
    ##!## Plot the DATA3 trends
    ##!## - - - - - - - - - - - - - - - - - 

    ##!## Make copy of OMI_data array
    ##!#local_data   = np.copy(OMI_data3['AI'][month_idx::index_jumper,:,:])
    ##!#local_counts = np.copy(OMI_data3['OB_COUNT'][month_idx::index_jumper,:,:])
    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,OMI_data3['LON'][0,:])
    ##!#cyclic_counts,cyclic_lons = add_cyclic_point(local_counts,OMI_data3['LON'][0,:])
    ##!#plat,plon = np.meshgrid(OMI_data3['LAT'][:,0],cyclic_lons)   

    ##!#local_mask = np.ma.masked_where(cyclic_counts == 0, cyclic_data)
    ##!#local_mask = np.ma.masked_where(local_mask == -999.9, local_mask)
    ##!#ai_trends = np.zeros(cyclic_data.shape[1:])
    ##!## Loop over all the keys and print the regression slopes 
    ##!## Grab the averages for the key
    ##!#for i in range(0,len(np.arange(np.min(OMI_data3['LAT']),90))):
    ##!#    for j in range(0,len(lon_ranges)):
    ##!#        # Check the current max and min
    ##!#        x_vals = np.arange(0,len(local_mask[:,i,j]))
    ##!#        # Find the slope of the line of best fit for the time series of
    ##!#        # average data
    ##!#        if(trend_type=='standard'): 
    ##!#            slope, intercept, r_value, p_value, std_err = \
    ##!#                stats.linregress(x_vals,local_mask[:,i,j])
    ##!#            ai_trends[i,j] = slope * len(x_vals)
    ##!#        else:
    ##!#            print("Ignoring MK trend for now")

   
    ##!#ai_trends = np.ma.masked_where(plat.T < minlat, ai_trends)
    ##!#   
    ##!#ax21 = plt.subplot(gs[1,2],projection=mapcrs)
    ##!#ax21.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    ##!#ax21.set_boundary(circle, transform=ax21.transAxes)
    ##!##ax21.gridlines()
    ##!#ax21.coastlines(resolution='50m')
    ##!#mesh21 = ax21.pcolormesh(plon,plat,ai_trends.T,\
    ##!#        transform=ccrs.PlateCarree(),vmin=-0.7,vmax=0.7,cmap=colormap)
    ##!##ax21.set_title(OMI_data3['VERSION']+ ' Trend')
    ##!#cbar21 = plt.colorbar(mesh21,ticks = np.arange(-1.0,1.1,0.2),orientation='vertical',\
    ##!#    extend='both',shrink=0.8)
    ##!#cbar21.set_label('AI Pert. Trend (AI/Study Period)',weight='bold',\
    ##!#    fontsize=colorbar_label_size)
    #cbar20 = plt.colorbar(mesh20,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
    #    aspect=50,shrink = 0.845)
    #cbar20.ax.tick_params(labelsize=14)

    #fig.tight_layout()

    outname = 'omi_TEST_ai_comps_'+start_date.strftime("%b")+'_'+\
        OMI_data1['VERSION']+'v'+OMI_data2['VERSION']+\
        'v'+OMI_data3['VERSION']+'.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_omi_da(OMI_da_nc,save=False):

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #
    # Plot the gridded data
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    mapcrs = ccrs.NorthPolarStereo()
    datacrs = ccrs.PlateCarree()
    colormap = plt.cm.jet
    
    # Set up the polar stereographic projection map
    fig1, ax = plt.subplots(figsize=(8,8))
    ax = plt.axes(projection = ccrs.NorthPolarStereo(central_longitude = 0.))
    ax.gridlines()
    ax.coastlines(resolution = '50m')
   
    UVAI   = OMI_da_nc.variables['AI'][:,:]
    counts = OMI_da_nc.variables['AI_COUNT'][:,:]
    lat    = OMI_da_nc.variables['lat'][:,:]
    lon    = OMI_da_nc.variables['lon'][:,:]
    plot_time = OMI_da_nc.variables['AI'].description[-21:-11]
 
    # Use meshgrid to convert the 1-d lat/lon arrays into 2-d, which is needed
    # for pcolormesh.
    #plot_lat, plot_lon = np.meshgrid(lat_ranges,lon_ranges)
    mask_UVAI = np.ma.masked_where(counts == 0, UVAI)
    
    plt.title('OMI DA AI ' + plot_time)
    mesh = ax.pcolormesh(lon, lat,mask_UVAI,transform = datacrs,cmap = colormap,\
            vmin = -2.0, vmax = 3.1)
            #vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])
    
    # Center the figure over the Arctic
    ax.set_extent([-180,180,60,90],ccrs.PlateCarree())
    
    # Depending on the desired variable, set the appropriate colorbar ticks
    tickvals = np.arange(-2.0,4.1,0.5)
    
    cbar = plt.colorbar(mesh,ticks = tickvals,orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label='AI')
    
    if(save == True):
        out_name = 'omi_da_ai_'+\
            plot_time+'.png'
        plt.savefig(out_name)
        print('Saved image '+out_name)
    else:
        plt.show()

###
###time_dict = {
###    '20200319': {
###        'start': 15,
###        'synop_max': 36
###    },
###    '20190411': {
###        'start': 12,
###        'synop_max': 18
###    }
###}
###
###first_back = time_dict['20200319']['start']-6
###secnd_back = time_dict['20200319']['start']-12
###
###if(first_back % 6 == 0):
###    end_time = time_dict['20200319']['synop_max'] + 1
###else:
###    end_time = 19
###
###for ftime in range(0,end_time):
###    file_name = ......'/hrrr.t' + str(first_back).zfill(2)+'z.wrfsfcf' + str(ftime).zfill(2) +' 
###    hrrr.t06z.wrfsfcf17z.20200319.grib2'

# NOTE: designed to be run with an OMI single-swath dictionary from
# readOMI_single_swath
def plotOMI_hrly(OMI_data_hrly, minlat = 60, pax = None, save=False):

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    local_data = np.copy(OMI_data_hrly['AI'])

    plot_lat, plot_lon = np.meshgrid(OMI_data_hrly['lat'],OMI_data_hrly['lon'])
    mask_AI = np.ma.masked_where(OMI_data_hrly['AI_count'] == 0, local_data)

    ##!## Determine the percentage of grid boxes that are actually filled
    ##!## with values.
    ##!#total_boxes = mask_AI.size
    ##!#total_good = mask_AI.count()
    ##!#pcnt_good = (total_good / total_boxes) * 100.
    ##!#print("Total_boxes = ",total_boxes,"Total good = ",total_good)
    ##!#print("Percent good = ",pcnt_good)

    if(pax is None):
        plt.close('all') 
        fig = plt.figure(figsize=(8,8))
        ax = plt.axes(projection = mapcrs)
        plotOMI_spatial(ax, OMI_data_hrly['LAT'], OMI_data_hrly['LON'], mask_AI, \
            'climo', ptitle = new_label_dict[OMI_data_hrly['VERSION']] + '\n', \
            plabel = 'UV Aerosol Index', \
            vmin = -2.0, vmax = 4.0, minlat = minlat)
   
    ##!#plt.close('all') 
    ##!#fig = plt.figure(figsize=(8,8))
    ##!#ax = plt.axes(projection = mapcrs)
    ##!#ax.gridlines()
    ##!#ax.coastlines(resolution = '50m')
    ##!#plt.title('OMI AI ' + OMI_data_hrly['date'])
    ##!##plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
    ##!#mesh = ax.pcolormesh(plot_lon, plot_lat,mask_AI,transform = datacrs,cmap = colormap,\
    ##!#        vmin = -2.0,vmax = 4.0)
    ##!#ax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    ##!#        #vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])
    ##!#cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
    ##!#    aspect=50,shrink = 0.845,label='UV Aerosol Index')
    
        if(save == True):
            out_name = 'omi_single_pass_uvai_' + OMI_data_hrly['date'] + \
                '_rows_0to' + str(OMI_data_hrly['row_max']) + '.png'
            plt.savefig(out_name,dpi=300)
            print('Saved image '+out_name)
        else:
            plt.show()
    else:
        plotOMI_spatial(pax, OMI_data_hrly['LAT'], OMI_data_hrly['LON'], mask_AI, \
            'climo', ptitle = new_label_dict[OMI_data_hrly['VERSION']] + '\n', \
            plabel = 'UV Aerosol Index', \
            vmin = -2.0, vmax = 4.0, minlat = minlat)

def plot_OMI_v_MODIS(MODIS_data,OMI_single_swath,save=False):

    # Pull out data from dictionaries
    local_modis = np.copy(MODIS_data['data'])[0,:,:].T
    local_omi   = np.ma.MaskedArray.copy(OMI_single_swath['AI_algae'])

    # shape = (1441, 121)
    plot_lat, plot_lon = np.meshgrid(MODIS_data['lat'],MODIS_data['lon'])

    # Mask any data outside the Finland region

    ##!#local_modis[(local_modis == -9) |\
    ##!#            (np.argwhere(np.isnan(local_modis))) |\
    ##!#            ((plot_lon < 17.) | (plot_lon > 54)) |\
    ##!#            ((plot_lat < 68.) | (plot_lat > 76)) |\
    ##!#            (local_modis > 1.2) |\
    ##!#            (local_omi == np.nan)] = np.nan
    ##!#local_omi[(local_modis == -9) |\
    ##!#            (np.argwhere(np.isnan(local_modis))) |\
    ##!#            ((plot_lon < 17.) | (plot_lon > 54)) |\
    ##!#            ((plot_lat < 68.) | (plot_lat > 76)) |\
    ##!#            (local_modis > 1.2) |\
    ##!#            (local_omi == np.nan)] = np.nan

    ##!#local_modis = np.ma.masked_invalid(local_modis)
    ##!#local_omi   = np.ma.masked_invalid(local_omi)

    # Ignore any data outside the Barents Sea box
    # -------------------------------------------
    local_modis = np.ma.masked_where(local_modis == -9, local_modis)
    local_modis = np.ma.masked_invalid(local_modis)
    local_modis = np.ma.masked_where(((plot_lon < 17.) | (plot_lon > 54)), local_modis)
    local_modis = np.ma.masked_where(((plot_lat < 68.) | (plot_lat > 76)), local_modis)
    local_modis = np.ma.masked_where((local_modis > 1.2), local_modis)
    local_omi   = np.ma.masked_where(((plot_lon < 17.) | (plot_lon > 54)), local_omi)
    local_omi   = np.ma.masked_where(((plot_lat < 68.) | (plot_lat > 76)), local_omi)
    #local_modis = np.ma.masked_where(local_omi == np.nan,local_modis)  

    # Get rid of any masked data
    local_modis = local_modis[~local_omi.mask]
    local_omi   = local_omi[~local_omi.mask]

    local_omi   = local_omi[~local_modis.mask]
    local_modis = local_modis[~local_modis.mask]


    # Make a datetime object for the title
    # ------------------------------------
    splitter = MODIS_data['titles'][0].split('/')[-1]
    plot_date = datetime.strptime(splitter[1:5],'%Y') + \
        relativedelta(days = int(splitter[5:8])-1)

    final_modis = local_modis.flatten()
    final_omi   = local_omi.flatten()
    print("Pearson:  ",pearsonr(final_modis,final_omi))
    print("Spearman: ",spearmanr(final_modis,final_omi))

    # Make the figure
    # --------------
    plt.close('all')
    fig1 = plt.figure(figsize=(7,6))
    plt.scatter(local_modis,local_omi,color='black')
    plt.plot(np.unique(final_modis), np.poly1d(np.polyfit(final_modis,final_omi,1))(np.unique(final_modis)),\
        color='red')
    plt.xlim(-0.001, 0.013)
    plt.title(plot_date.strftime('%Y%m%d'))
    plt.xlabel('MODIS ' + MODIS_data['ptitle'] + '\n['+MODIS_data['label']+']')
    plt.ylabel('OMI AI')

    if(save == True):
        outname = 'omi_algae_v_modis_'+MODIS_data['grabber'] + \
            '_' + plot_date.strftime('%Y%m%d') + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_time_diff(jz28,jz2,month):
    diffs = jz28 - jz2
    fig1 = plt.figure()
    plt.plot(diffs)
    plt.ylim(-0.15,0.35)
    plt.title('VJZ28 - VJZ2: ' + month)
    plt.savefig('omi_series_diff_'+month+'_jz28jz2_75N150.png')

# Plots a scatter of AI trends from one datatype to those from another  
# --------------------------------------------------------------------
def plot_compare_trends(OMI_data1,OMI_data2,month, pax = None, \
        trend_type = 'standard', minlat = 65., plot_trend = True, save=False):

    ##OMI_trend1 = plotOMI_MonthTrend(OMI_data1,month_idx=month,save=True,\
    ##            trend_type='standard',season='',minlat=65.,return_trend=True)
    ##OMI_trend2 = plotOMI_MonthTrend(OMI_data2,month_idx=month,save=True,\
    ##            trend_type='standard',season='',minlat=65.,return_trend=True)

    OMI_trend1 = calcOMI_grid_trend(OMI_data1, month, trend_type, \
        minlat)
    OMI_trend2 = calcOMI_grid_trend(OMI_data2, month, trend_type, \
        minlat)

    # Convert the index to a string using datetime
    if(month != None):
        dt_obj = datetime.strptime(OMI_data1['DATES'][month],"%Y%m")
        title = 'OMI AI ' + dt_obj.strftime("%b") + " Trend Comparison"
        outname = 'omi_ai_trend_comp_'+dt_obj.strftime("%b")+'_'+\
            OMI_data1['VERSION']+'v'+OMI_data2['VERSION']+'.png'
    else:
        title = 'OMI AI Trend Comparison'
        outname = 'omi_ai_trend_comp_'+\
            OMI_data1['VERSION']+'v'+OMI_data2['VERSION']+'.png'

    OMI_trend2[np.where(np.isnan(OMI_trend1) | np.isnan(OMI_trend2))] = -999.
    OMI_trend1[np.where(np.isnan(OMI_trend1) | np.isnan(OMI_trend2))] = -999.
    mask_trend1 = np.array(OMI_trend1[(OMI_trend1 != 0) & \
        (OMI_trend2 != 0) & (OMI_trend1 != -999.) & \
        (OMI_trend2 != -999.) & (OMI_data1['LAT'] >= minlat)])
    mask_trend2 = np.array(OMI_trend2[(OMI_trend1 != 0) & \
        (OMI_trend2 != 0) & (OMI_trend1 != -999.) & \
        (OMI_trend2 != -999.) & (OMI_data2['LAT'] >= minlat)])
    #mask_trend1 = np.ma.masked_where(OMI_trend1 == 0,OMI_trend1)
    #mask_trend2 = np.ma.masked_where(OMI_trend2 == 0,OMI_trend2)

    print(mask_trend1.shape,mask_trend2.shape)

    rval_p = pearsonr(mask_trend1,mask_trend2)[0]
    rval_s = spearmanr(mask_trend1,mask_trend2)[0]
    print("Pearson:  ",rval_p)
    print("Spearman: ",rval_s)

    xy = np.vstack([mask_trend1,mask_trend2])
    z = stats.gaussian_kde(xy)(xy)

    ##!## Plot a somewhat-robust best fit line using Huber Regression
    ##!## -----------------------------------------------------------
    ##!#x_scaler,y_scaler = StandardScaler(), StandardScaler()
    ##!#x_train = x_scaler.fit_transform(mask_trend1[...,None])
    ##!#y_train = y_scaler.fit_transform(mask_trend2[...,None])

    ##!#model = HuberRegressor(epsilon=1)
    ##!#model.fit(x_train,y_train)
    ##!#
    ##!#test_x = np.array([np.min(mask_trend1),np.max(mask_trend1)])
    ##!#print(test_x)
    ##!#predictions = y_scaler.inverse_transform(\
    ##!#    model.predict(x_scaler.transform(test_x[...,None])))

    # One to one line stuff
    xs = np.arange(np.min(mask_trend1),np.max(mask_trend1),0.1)

    no_ax = True
    if(pax is None):
        no_ax = False
        plt.close('all')
        fig1 = plt.figure()
        pax = fig1.add_subplot(1,1,1)
    pax.scatter(mask_trend1,mask_trend2,c=z,s=8)
    #pax.plot(test_x,predictions,color='tab:green',linestyle='--',label='Huber Fit')

    if(plot_trend):
        plot_trend_line(pax, mask_trend1, mask_trend2, color='tab:green', linestyle = '-', \
            slope = 'theil-sen')
    ### Plot an unrobust fit line using linear regression
    ### -------------------------------------------------
    ##pax.plot(np.unique(mask_trend1),np.poly1d(np.polyfit(mask_trend1,\
    ##    mask_trend2,1))(np.unique(mask_trend1)),color='tab:orange',\
    ##    linestyle='--',label='Polyfit Fit')
    # Plot a one-to-one line
    pax.plot(xs,xs,label='1-1',color='tab:red')
    if((month == 0) | (month == 1)):
        pax.set_xlim(-0.5,0.3)
        pax.set_ylim(-0.5,0.3)
    elif((month == 2)):
        pax.set_xlim(-0.6,0.5)
        pax.set_ylim(-0.6,0.5)
    elif((month == 3)):
        pax.set_xlim(-0.5,0.5)
        pax.set_ylim(-0.5,0.5)
    elif((month == 4)):
        pax.set_xlim(-0.5,0.7)
        pax.set_ylim(-0.5,0.7)
    elif((month == 5)):
        pax.set_xlim(-0.5,0.5)
        pax.set_ylim(-0.5,0.5)
    else:
        pax.set_xlim(-0.3,0.3)
        pax.set_ylim(-0.3,0.3)
    pax.legend(loc = 'lower right')
    pax.set_xlabel(OMI_data1['VERSION'])
    pax.set_ylabel(OMI_data2['VERSION'])
    # Add the correlations to the graph
    ##!#x_pos = (plt.gca().get_xlim()[1] - \
    ##!#    plt.gca().get_xlim()[0])*0.7 + \
    ##!#    plt.gca().get_xlim()[0]
    ##!#y_pos = (plt.gca().get_ylim()[1] - \
    ##!#    plt.gca().get_ylim()[0])*0.05 + \
    ##!#    plt.gca().get_ylim()[0]
    ##!#plt.text(x_pos,y_pos,\
    ##!#    'Pearson r     = '+str(np.round(rval_p,3)))
    ##!#x_pos = (plt.gca().get_xlim()[1] - \
    ##!#    plt.gca().get_xlim()[0])*0.7 + \
    ##!#    plt.gca().get_xlim()[0]
    ##!#y_pos = (plt.gca().get_ylim()[1] - \
    ##!#    plt.gca().get_ylim()[0])*0.1 + \
    ##!#    plt.gca().get_ylim()[0]
    ##!#plt.text(x_pos,y_pos,\
    ##!#    'Spearman r = '+str(np.round(rval_s,3)))
    pax.set_title(title+'\n Pearson Correlation = ' + str(np.round(rval_p,3)))
    if(not no_ax):
        if(save == True):
            plt.savefig(outname)
        else:
            plt.show()
    #return mask_trend1,mask_trend2

# Plot scatter comparisons of two trend types for each of the months
def plotOMI_trend_scatter_comp(OMI_data1, OMI_data2, minlat = 65.,\
        trend_type = 'standard', plot_trend = False, save = False):

    # ----------------------------------------------------
    #  Create a figure to hold all 6 scatter plots
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (14,10))
    ax0 = fig1.add_subplot(2,3,1)
    ax1 = fig1.add_subplot(2,3,2)
    ax2 = fig1.add_subplot(2,3,3)
    ax3 = fig1.add_subplot(2,3,4)
    ax4 = fig1.add_subplot(2,3,5)
    ax5 = fig1.add_subplot(2,3,6)

    # ----------------------------------------------------
    #
    # Plot the trend comparisons for each month
    #
    # ----------------------------------------------------
    plot_compare_trends(OMI_data1,OMI_data2,0, pax = ax0, \
        trend_type = trend_type, minlat = minlat, plot_trend = plot_trend, \
        save=False)
    plot_compare_trends(OMI_data1,OMI_data2,1, pax = ax1, \
        trend_type = trend_type, minlat = minlat, plot_trend = plot_trend, \
        save=False)
    plot_compare_trends(OMI_data1,OMI_data2,2, pax = ax2, \
        trend_type = trend_type, minlat = minlat, plot_trend = plot_trend, \
        save=False)
    plot_compare_trends(OMI_data1,OMI_data2,3, pax = ax3, \
        trend_type = trend_type, minlat = minlat, plot_trend = plot_trend, \
        save=False)
    plot_compare_trends(OMI_data1,OMI_data2,4, pax = ax4, \
        trend_type = trend_type, minlat = minlat, plot_trend = plot_trend, \
        save=False)
    plot_compare_trends(OMI_data1,OMI_data2,5, pax = ax5, \
        trend_type = trend_type, minlat = minlat, plot_trend = plot_trend, \
        save=False)

    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')
    plot_subplot_label(ax2, '(c)')
    plot_subplot_label(ax3, '(d)')
    plot_subplot_label(ax4, '(e)')
    plot_subplot_label(ax5, '(f)')

    #plt.suptitle(date_str)
    fig1.tight_layout()

    if(save):
        outname = 'omi_combined_scatter_compare_' + OMI_data1['VERSION'] + \
            'v' + OMI_data2['VERSION'] + '_min' + str(int(minlat)) + '.png'
        fig1.savefig(outname, dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# ----------------------------------------------------------------------------
#
# Compare 1x1 degree single-month data from OMI and CERES
#
# ----------------------------------------------------------------------------
def plot_compare_OMI_CERES_grid(OMI_data, CERES_data, midx, minlat=65, \
        max_AI = -200., only_sea_ice = False, save=False):
#def plot_compare_OMI_CERES_hrly(OMI_hrly,CERES_hrly,minlat=65,save=False):

    if('/home/bsorenson/Research/CERES/' not in sys.path):
        sys.path.append('/home/bsorenson/Research/CERES/')

    from gridCERESLib import plotCERES_Month

    # Step 1: Set up figure to have 3 panels
    # --------------------------------------
    plt.close('all')
    #fig, axs = plt.subplots(1,3,subplot_kw={'projection': mapcrs})
    fig1 = plt.figure(figsize=(17,5))
    ax0 = fig1.add_subplot(1,3,1,projection = mapcrs)
    ax1 = fig1.add_subplot(1,3,2,projection = mapcrs)
    ax2 = fig1.add_subplot(1,3,3)

    # Flip the CERES data to convert the longitudes from 0 - 360 to -180 - 180
    # ------------------------------------------------------------------------
    local_omi   = np.copy(OMI_data['AI'][midx])
    local_lon   = np.copy(CERES_data['lon'])
    local_lat   = np.copy(CERES_data['lat'])
    local_ceres = np.copy(CERES_data['data'][midx])
    over_180    = np.where(CERES_data['lon'][0,:] >= 179.9999)
    under_180   = np.where(CERES_data['lon'][0,:] < 179.9999)

    print(local_lat.shape, local_ceres.shape)

    for ii in range(local_lon.shape[0]):
        local_lon[ii,:] = np.concatenate([local_lon[ii,:][over_180] - 360.,\
            local_lon[ii,:][under_180]])
        local_lat[ii,:] = np.concatenate([local_lat[ii,:][over_180],\
            local_lat[ii,:][under_180]])
        local_ceres[ii,:] = np.concatenate([local_ceres[ii,:][over_180],\
            local_ceres[ii,:][under_180]])
   
    print("Before removal, ",local_ceres.shape)
 
    # First, mask any OMI data that are outside the bounds of the CERES data
    # Assume that the OMI data are larger than the CERES data
    # ----------------------------------------------------------------------
    local_omi[np.where(np.isnan(local_omi))] = -999.
    where_matching = np.where((OMI_data['LAT'][:,0] >= \
        (np.min(local_lat) - 0.5)) & (OMI_data['LAT'][:,0] <= \
        (np.max(local_lat) - 0.5)))[0]
    local_omi = local_omi[where_matching, :]

    mask_omi = np.array(local_omi[  (local_omi != 0) & (local_ceres != 0) \
        & (local_omi != -999.) & (local_ceres != -999.)])
    mask_ceres = np.array(local_ceres[(local_omi != 0) & (local_ceres != 0) \
        & (local_omi != -999.) & (local_ceres != -999.)])

    # Step 3: Plot OMI and CERES data in first two panels
    # ---------------------------------------------------
    plotOMI_NCDF_SingleMonth(OMI_data, midx ,minlat= minlat, pax = ax0)
    plotCERES_Month(CERES_data, midx, pax = ax1, minlat = 65)

    # Plot the scatter comparison
    # ---------------------------
    xy = np.vstack([mask_omi,mask_ceres])
    z = stats.gaussian_kde(xy)(xy)

    ax2.scatter(local_omi, local_ceres, c = z, s = 6)
    plot_trend_line(ax2, local_omi, local_ceres, color='tab:green', linestyle = '-', \
        slope = 'theil-sen')

    fig1.tight_layout()

    ##!## Step 4: Plot scatter OMI and CERES comparison in third panel
    ##!## ------------------------------------------------------------

    ##!## Make gridded lat and lon arrays to use for coloring the plot
    ##!#glat, glon = np.meshgrid(OMI_hrly['lat'],OMI_hrly['lon'])
    ##!#
    ##!## Mask any empty boxes
    ##!#mask_AI = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
    ##!#    (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),OMI_hrly['AI'])
    ##!#mask_flux = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
    ##!#    (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),CERES_hrly['data'])
    ##!#mask_sza  = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
    ##!#    (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),OMI_hrly['SZA'])
    ##!##mask_lat  = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
    ##!##    (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),glat)
    ##!###!## Convert the index to a string using datetime
    ##!###!#if(month != None):
    ##!###!#    dt_obj = datetime.strptime(OMI_data['DATES'][month],"%Y%m")
    ##!###!#    title = 'OMI AI / CERES ' + CERES_data['param'] + '\n'+ \
    ##!###!#        dt_obj.strftime("%b") + " Trend Comparison"
    ##!###!#    outname = 'omi_ceres_trend_comp_'+dt_obj.strftime("%b")+'_'+\
    ##!###!#        OMI_data['VERSION']+'vCERES.png'
    ##!###!#else:
    ##!###!#    title = 'OMI AI / CERES ' + CERES_data['param'] + ' Trend Comparison'
    ##!###!#    outname = 'omi_ceres_trend_comp_'+\
    ##!###!#        OMI_data1['VERSION']+'vCERES.png'


    ##!##print("Pearson:  ",pearsonr(mask_AI,mask_flux))
    ##!##print("Spearman: ",spearmanr(mask_AI,mask_flux))

    ##!###!#xy = np.vstack([mask_AI,mask_flux])
    ##!###!#z = stats.gaussian_kde(xy)(xy)


    ##!###!## One to one line stuff
    ##!###!#xs = np.arange(np.min(mask_AI),np.max(mask_AI),0.1)

    ##!#scat = ax2.scatter(mask_AI,mask_flux,c=mask_sza,s=8)
    ##!#cbar = plt.colorbar(scat,ax=ax2,orientation='vertical',\
    ##!#    label='Solar Zenith Angle',pad=0.02,aspect=50)
    ##!#ax2.set_title(OMI_hrly['date'])
    ##!#ax2.set_xlabel('OMI AI')
    ##!#ax2.set_ylabel('CERES ' + CERES_hrly['param'])
    ##!#plt.tight_layout()
    outname = 'omi_ceres_compare_'+OMI_data['DATES'][midx]+'.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()
    #return mask_AI,mask_flux


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Single-swath plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# NOTE: this is designed to be run with the brand new readOMI_swath_hdf
def plotOMI_single_swath(pax, OMI_hrly, pvar = 'UVAI', minlat = 65., \
        vmin = None, vmax = None, title = '', label = '', \
        circle_bound = True, gridlines = True, save=False):

    variable = 'UVAerosolIndex'

    if(vmin is None):
        vmin = var_dict[variable]['min']
    if(vmax is None):
        vmax = var_dict[variable]['max']

    if(gridlines):
        pax.gridlines()
    pax.coastlines(resolution = '50m')
    if(title == ''):
        title = 'OMI UVAI ' + OMI_hrly['dtype'] + ' ' + OMI_hrly['date']
    pax.set_title(title)
    if(pvar == 'GPQF'):
        # Convert the GPQF flag values to ice flags
        # -----------------------------------------
        local_GPQF = np.full(OMI_hrly[pvar].shape,np.nan)

        local_GPQF[(OMI_hrly[pvar] > 0) & (OMI_hrly[pvar] < 101)] = 9
        local_GPQF[OMI_hrly[pvar] == 0]                           = 8
        local_GPQF[OMI_hrly[pvar] == 101]                         = 10
        local_GPQF[OMI_hrly[pvar] == 103]                         = 11
        local_GPQF[OMI_hrly[pvar] == 104]                         = 12
        local_GPQF[OMI_hrly[pvar] == 124]                         = 13

        mask_GPQF = np.ma.masked_invalid(local_GPQF)
        cmap = plt.get_cmap('jet', np.max(mask_GPQF) - np.min(mask_GPQF)+1)

        cbar_labels = ['Shallow ocean','Land','Shallow inland water',\
                        'Ocean coastline / lake shoreline',\
                        'ephemeral (intermittent) water','Deep inland water',\
                        'Continental shelf ocean','Deep ocean',\
                        'Snow-free land','Sea ice',\
                        'Permanent ice\n(Greenland, Antarctica)','Dry snow',\
                        'Ocean','Mixed pixels\nat coastline',\
                        'Suspect ice value','Corners']

        mesh = pax.pcolormesh(OMI_hrly['LON'], OMI_hrly['LAT'], mask_GPQF,\
                transform = datacrs, cmap = cmap,\
                vmin = np.nanmin(mask_GPQF) - 0.5, \
                vmax = np.nanmax(mask_GPQF) + 0.5,\
                shading='auto')

        #cbar = plt.colorbar(mesh,ax = pax, orientation='horizontal',pad=0,\
        cbar = plt.colorbar(mesh,ax = pax, orientation='vertical',\
            shrink = 0.8, ticks = np.arange(int(np.nanmin(mask_GPQF)), \
            int(np.nanmax(mask_GPQF)) + 1))
            #shrink = 0.8, ticks = np.arange(np.nanmin(mask_GPQF), \
            #np.nanmax(mask_GPQF) + 1))
        #cbar.ax.set_xticks(np.arange(int(np.nanmin(mask_GPQF)),int(np.nanmax(mask_GPQF)) + 1))
        print(cbar_labels[int(np.nanmin(mask_GPQF)):int(np.nanmax(mask_GPQF))+1])
        cbar.ax.set_yticklabels(cbar_labels[int(np.nanmin(mask_GPQF)):int(np.nanmax(mask_GPQF))+1],\
            fontsize=10,weight = 'bold', rotation=0)
            #fontsize=8,rotation=35)
        
    else:
        mesh = pax.pcolormesh(OMI_hrly['LON'], OMI_hrly['LAT'], OMI_hrly[pvar],\
                transform = datacrs, cmap = colormap,\
                vmin = vmin, vmax = vmax,\
                shading='auto')

        if(label == ''):
            if(OMI_hrly['dtype'] == 'shawn'):
                label = 'UVAI Perturbation'
            else:
                label = 'UV Aerosol Index'

        #cbar.set_label(variable,fontsize=16,weight='bold')
        tickvals = np.arange(-2.0,4.1,0.5)
        cbar = plt.colorbar(mesh,ax = pax, ticks = tickvals,\
            orientation='vertical', shrink = 0.8, extend = 'both')
        cbar.set_label(label,fontsize = 14, weight='bold')
        cbar.ax.tick_params(labelsize=14)

    # Center the figure over the Arctic
    pax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    pax.add_feature(cfeature.BORDERS)
    pax.add_feature(cfeature.STATES)
    if(circle_bound):
        pax.set_boundary(circle, transform=pax.transAxes)
    # Depending on the desired variable, set the appropriate colorbar ticks
    ##!#if(variable == 'NormRadiance'):
    ##!#    tickvals = np.arange(0.0,1.010,0.10)
    ##!#elif(variable == 'Reflectivity'):
    ##!#    tickvals = np.arange(0.0,1.1,0.10)
    ##!#elif(variable == 'FinalAlgorithmFlags'):
    ##!#    tickvals = np.arange(0,9)
    ##!#elif(variable == 'SolarZenithAngle'):
    ##!#    tickvals = np.arange(0,90,10)
    ##!#else:
    
    #cbar10.set_label('UV Aerosol Index',weight='bold',fontsize=colorbar_label_size)
    #cbar.ax.tick_params(labelsize=14)

# Plot just a single swath on a 1-panel figure
# --------------------------------------------
def plotOMI_single_swath_figure(date_str, dtype = 'control',  \
        only_sea_ice = False, minlat = 65., save = False):

    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    if(dtype == 'shawn'):
        OMI_base  = readOMI_swath_shawn(date_str, latmin = minlat)
    else:
        OMI_base  = readOMI_swath_hdf(date_str, dtype, \
            only_sea_ice = only_sea_ice, latmin = minlat)

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (6,6))
    mapcrs = ccrs.Robinson()
    ax0 = fig1.add_subplot(1,1,1, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    plotOMI_single_swath(ax0, OMI_base, title = dtype.title(), \
        circle_bound = False)

    ax0.set_extent([-180., -140., -40., 0.], datacrs)

    plt.suptitle(date_str)

    #fig1.tight_layout()

    if(save):
        outname = 'omi_single_swath_figure_' + date_str + '_' + \
            dtype + '.png'
        fig1.savefig(outname, dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Plot four panels: two different single swaths, as well as
# data zoomed in over Greenland to show the high pixels
# ---------------------------------------------------------
def plotOMI_single_swath_multiple(date_str, dtype = 'control',  \
        only_sea_ice = False, minlat = 65., save = False):

    dates = ['201206141245','201206141920']
    date_str = dates[0]
    #dates = ['200806141226','200806141902']
    #dates = ['200807261124','200807261303']
    #dates = ['200804221027','200804221345']
    #dates = ['200804221027','200804221206']

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (10,10))
    #mapcrs = ccrs.Robinson()
    mapcrs2 = ccrs.NorthPolarStereo(central_longitude = -40)
    ax0 = fig1.add_subplot(2,2,1, projection = mapcrs)
    ax1 = fig1.add_subplot(2,2,2, projection = mapcrs)
    ax2 = fig1.add_subplot(2,2,3, projection = mapcrs2)
    ax3 = fig1.add_subplot(2,2,4, projection = mapcrs2)

    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    if(dtype == 'shawn'):
        OMI_base1 = readOMI_swath_shawn(dates[0], latmin = minlat - 5)
        OMI_base2 = readOMI_swath_shawn(dates[1], latmin = minlat - 5)
    else:
        OMI_base1 = readOMI_swath_hdf(dates[0], dtype, \
            only_sea_ice = only_sea_ice, latmin = minlat - 5)
        OMI_base2 = readOMI_swath_hdf(dates[1], dtype, \
            only_sea_ice = only_sea_ice, latmin = minlat - 5)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    plotOMI_single_swath(ax0, OMI_base1, title = dates[0], \
        circle_bound = True, gridlines = False)
    plotOMI_single_swath(ax1, OMI_base2, title = dates[1], \
        circle_bound = True, gridlines = False)
    plotOMI_single_swath(ax2, OMI_base1, title = dates[0] + '\nZoomed', \
        circle_bound = False, gridlines = False)
    plotOMI_single_swath(ax3, OMI_base2, title = dates[1] + '\nZoomed', \
        circle_bound = False, gridlines = False)

    ax2.set_extent([-70., -10., 65., 87.], datacrs)
    ax3.set_extent([-70., -10., 65., 87.], datacrs)

    #plt.suptitle(date_str)
    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')
    plot_subplot_label(ax2, '(c)', backgroundcolor = 'white')
    plot_subplot_label(ax3, '(d)', backgroundcolor = 'white')

    fig1.tight_layout()

    if(save):
        outname = 'omi_single_swath_figure_multiple_' + date_str + '_' + \
            dtype + '.png'
        fig1.savefig(outname, dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# 2 panel plot showing control UVAI and ground classification
# -----------------------------------------------------------
def plotOMI_single_ground(date_str, only_sea_ice = False, minlat = 65., \
        zoom = False, multi_panel = False, save = False):

    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    OMI_base  = readOMI_swath_hdf(date_str, 'control', \
        only_sea_ice = only_sea_ice, latmin = minlat)

    circle_bound = True
    zoom_add = ''
    if(zoom):
        zoom_add = '_zoom'
        circle_bound = False
        mapcrs = ccrs.NorthPolarStereo(central_longitude = -40)
        #mapcrs = ccrs.Robinson()
    #mapcrs = ccrs.NorthPolarStereo()
    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    if(multi_panel):
        fig1 = plt.figure(figsize = (10,10))
        ax0 = fig1.add_subplot(2,2,1, projection = mapcrs)
        ax1 = fig1.add_subplot(2,2,2, projection = mapcrs)
    else:
        fig1 = plt.figure(figsize = (11,5))
        ax0 = fig1.add_subplot(1,2,1, projection = mapcrs)
        ax1 = fig1.add_subplot(1,2,2, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    gridlines = False
    plotOMI_single_swath(ax0, OMI_base, title = 'Control', \
        circle_bound = circle_bound, gridlines = gridlines)
    plotOMI_single_swath(ax1, OMI_base, title = 'Screened', pvar = 'GPQF', \
        circle_bound = circle_bound, gridlines = gridlines)

    if(zoom):
        ax0.set_extent([-70., -10., 65., 87.], datacrs)
        ax1.set_extent([-70., -10., 65., 87.], datacrs)

    plt.suptitle(date_str)
    plot_subplot_label(ax0, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax1, '(b)', backgroundcolor = 'white')

    fig1.tight_layout()

    if(save):
        outname = 'omi_single_swath_ground_' + date_str + zoom_add + '.png'
        fig1.savefig(outname, dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Compare the control AI with the two methods
# -------------------------------------------
def plotOMI_single_3panel(date_str, only_sea_ice = False, minlat = 65., \
        save = False):

    # ----------------------------------------------------
    # Read in each of the 3 data types
    # ----------------------------------------------------
    OMI_base  = readOMI_swath_hdf(date_str, 'control', \
        only_sea_ice = only_sea_ice, latmin = minlat)
    OMI_JZ    = readOMI_swath_hdf(date_str, 'JZ', \
        only_sea_ice = only_sea_ice, latmin = minlat)
    OMI_shawn = readOMI_swath_shawn(date_str, latmin = minlat)

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (15,5))
    ax0 = fig1.add_subplot(1,3,1, projection = mapcrs)
    ax1 = fig1.add_subplot(1,3,2, projection = mapcrs)
    ax2 = fig1.add_subplot(1,3,3, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    plotOMI_single_swath(ax0, OMI_base, title = 'Control')
    plotOMI_single_swath(ax1, OMI_JZ, title = 'Screened')
    plotOMI_single_swath(ax2, OMI_shawn, pvar = 'UVAI_pert',\
        title = 'Perturbed')


    plt.suptitle(date_str)
    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')
    plot_subplot_label(ax2, '(c)')

    fig1.tight_layout()

    if(save):
        outname = 'omi_single_swath_compare_3panel_' + date_str + '.png'
        fig1.savefig(outname, dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

def plotOMI_shawn_3panel(date_str, only_sea_ice = False, minlat = 65.,\
        save = False):

    # ----------------------------------------------------
    # Read in the shawn data
    # ----------------------------------------------------
    OMI_shawn = readOMI_swath_shawn(date_str, latmin = minlat)

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (15,5))
    ax0 = fig1.add_subplot(1,3,1, projection = mapcrs)
    ax1 = fig1.add_subplot(1,3,2, projection = mapcrs)
    ax2 = fig1.add_subplot(1,3,3, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    plotOMI_single_swath(ax0, OMI_shawn, pvar = 'UVAI_raw', \
        label = 'UV Aerosol Index', title = 'Raw UVAI')
    plotOMI_single_swath(ax1, OMI_shawn, pvar = 'UVAI_climo', \
        label = 'UV Aerosol Index', title = 'UVAI Climatology')
    plotOMI_single_swath(ax2, OMI_shawn, pvar = 'UVAI_pert', \
        label = 'UVAI Perturbation', title = 'UVAI Perturbation')

    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')
    plot_subplot_label(ax2, '(c)')

    plt.suptitle(date_str)
    fig1.tight_layout()

    if(save):
        outname = 'omi_single_swath_shawn_3panel_' + date_str + '.png'
        fig1.savefig(outname, dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Plot a single swath of OMI data with total climatology subtracted
# mask_weakAI: removes AI values below 0.8 when plotting
def single_swath_anomaly_climo(OMI_data,swath_date,month_climo = True,\
        minlat = 60,row_max = 60,mask_weakAI = False, save=False): 
    # - - - - - - - - - - - - - - - -
    # Read in the single swath data
    # - - - - - - - - - - - - - - - -
    latmin = 60 
    
    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Set up values for gridding the AI data
    lat_gridder = latmin * 4.
    
    lat_ranges = np.arange(latmin,90.1,0.25)
    lon_ranges = np.arange(-180,180.1,0.25)
    base_path = '/home/bsorenson/data/OMI/H5_files/'
    year = swath_date[:4]
    date = swath_date[4:8]
    if(len(swath_date)==13):
        time = swath_date[9:]
    elif(len(swath_date)==12):
        time = swath_date[8:]
    else:
        time = ''
    total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+year+'m'+date+'t'+time+'*.he5',\
              shell=True).decode('utf-8').strip().split('\n')

    UVAI = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

    for fileI in range(len(total_list)):
        # read in data directly from HDF5 files
        print(total_list[fileI])
        data = h5py.File(total_list[fileI],'r')
        #ALBEDO = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo']
        #REFLECTANCE = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity']
        #CLD = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/CloudFraction']
        AI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
        #PIXEL= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/PixelQualityFlags']
        #MSMNT= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/MeasurementQualityFlags']
        #VZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/ViewingZenithAngle']
        #SZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SolarZenithAngle']
        #RAZ  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle']
        LAT = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
        LON = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
        #GRND = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags']
    
        #albedo = ALBEDO[:,:,0]   
        #reflectance = REFLECTANCE[:,:,0]   
        counter = 0
        #AI = AI[:,:,0]   
        # Loop over the values and rows 
        #for i in range(0,int(CBA2)):
        #for i in range(albedo.shape[0]):
        for i in range(AI.shape[0]):
            for j in range(0,row_max):
                #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
                if(AI[i,j]>-20):
                #if(plotAI[i,j]>-20):
                    # Only plot if XTrack flag is met
                    if((XTRACK[i,j] == 0) | ((XTRACK[i,j] & 4 == 4))):
                        # Print values to text file
                        if(LAT[i,j] > minlat):
                            counter+=1
    #                        fout.write("{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.6f} {9:.6f} {10:.6f} {11:.6f}\n".format(\
    #                            LAT[i,j],LON[i,j],AI[i,j],0.5,SZA[i,j],VZA[i,j],RAZ[i,j], \
    #                            ALBEDO[i,j,0],ALBEDO[i,j,1],REFLECTANCE[i,j,0],\
    #                            REFLECTANCE[i,j,1],CLD[i,j]))
    
    
                        #if((plotXTrack[i,j] == 0) | ((plotXTrack[i,j] & 4 == 4))):
                            index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                            index2 = int(np.floor(LON[i,j]*4 + 720.))
                            #index1 = int(np.floor(plotLAT[i,j]*4 + 360.))
                            #index2 = int(np.floor(plotLON[i,j]*4 + 720.))
                            
                            if(index1 < 0): index1 = 0
                            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                            if(index2 < 0): index2 = 0                                                                                            
                            if(index2 > 1439): index2 = 1439
                       
                            #diff = reflectance[i,j] - albedo[i,j]
                            #if(diff<min_diff):
                            #    min_diff = diff
                            #if(diff>max_diff):
                            #    max_diff = diff
                       
                            #if(diff<0.2): 
                            #    UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
                            UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
                            #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                                #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                            count[index2, index1] = count[index2,index1] + 1
    
    # Calculate the row-average AI for the secondary plot
    mask_avgs = np.nanmean(np.ma.masked_where(AI[:,:] < -20, AI[:,:]),axis=0)


    if(month_climo == True):
        # Determine which month to use based on the swath date
        monther = int(date[:2])-1
        local_climo = OMI_data['MONTH_CLIMO'][monther,:,:].T
        title_adder = 'Monthly Anomaly '
        file_adder = 'mnth_anom_'
    else:
        local_data = np.copy(OMI_data['AI'][:,:,:])
        local_mask = np.ma.masked_where(local_data == -999.9, local_data)
        local_climo = np.nanmean(local_mask, axis=0).T
        title_adder = 'Anomaly '
        file_adder = 'anom_'

    # Interpolate 2D data
    # Find where the high-res data fits in to the climo data
    begin_x = np.where(OMI_data['LAT'][:,0] == lat_ranges[0])[0][0]
    high_lats = np.arange(OMI_data['LAT'][begin_x,0], OMI_data['LAT'][-1,0]+0.25, 0.25)
    high_lons = np.arange(OMI_data['LON'][0,0], OMI_data['LON'][0,-1] + 0.25, 0.25)

    high_climo = np.zeros((high_lons.shape[0], high_lats.shape[0]))

    for yi in range(len(OMI_data['LON'][0,:])):
        for xj in range(len(OMI_data['LAT'][begin_x:,0])):
            high_climo[yi*4:yi*4+4, xj*4:xj*4+4] = local_climo[yi, xj+begin_x]

    # Find out where the high res climo data fits in to the
    # single swath data. The single swath data are assumed to be 
    # of a larger size than the climo data
    end_xj = np.where(lat_ranges == high_lats[-1])[0][0]
    end_yi = np.where(lon_ranges == high_lons[-1])[0][0]

    # Calculate anomalies
    UVAI_anomaly = UVAI[0:end_yi+1, 0:end_xj+1] - high_climo  

    # Mask grid boxes where the ob counts are zero
    plot_lat, plot_lon = np.meshgrid(high_lats,high_lons)
    mask_UVAI_anom = np.ma.masked_where(count[0:end_yi+1, 0:end_xj+1] == 0, UVAI_anomaly)
  
    ## Mask any data below AI values of 0.8
    #mask_UVAI_anom = np.ma.masked_where(mask_UVAI_anom < 0.8, mask_UVAI_anom)
 
    plt.close('all')
    fig1 = plt.figure(figsize=(8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
 
    plt.title('OMI AI '+title_adder + swath_date)
    #plt.title('OMI Reflectivity - Surface Albedo '+swath_date)
    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI_anom,transform = datacrs,cmap = colormap,\
            vmin = -2.0, vmax = 3.0)
    ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    #ax.set_xlim(-4170748.535086173,4167222.438879491)
    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label='Aerosol Index Anomaly')
    ##cax = fig.add_axes([0.16,0.075,0.7,0.025])
    ##cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    ##cb.ax.set_xlabel('Aerosol Index')
    #cb.ax.set_xlabel('Reflectivity - Surface Albedo')
    #out_name = 'omi_single_pass_ai_200804270052_to_0549_composite_rows_0to'+str(row_max)+'.png'       
    if(save == True):
        out_name = 'omi_single_pass_ai_'+file_adder+swath_date+'_rows_0to'+str(row_max)+'.png'       
        plt.savefig(out_name)
        print('Saved image '+out_name)
    else: 
        plt.show()
    #axs[1].plot(mask_avgs)
    #axs[1].set_xlabel('Sensor Row')
    #axs[1].set_ylabel('Row Average Aerosol Index')
    

# Plot a single swath of OMI data with single-time swath climatology subtracted
# single_swath = the swath that will be corrected (format = YYYYMMDDHHMM)
# climo_date   = the date on the swath climatology directory (format = MMDDHHMM)
def single_swath_anomaly_time(single_swath,climo_date,minlat = 60,row_max = 60): 

    # - - - - - - - - - - - - - - - -
    # Read in the single swath data
    # - - - - - - - - - - - - - - - -
    latmin = 60 
    
    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Set up values for gridding the AI data
    lat_gridder = latmin * 4.
    
    lat_ranges = np.arange(latmin,90.1,0.25)
    lon_ranges = np.arange(-180,180.1,0.25)
    climo_base_path = '/home/bsorenson/data/OMI/swath_anomaly_files/'
    single_base_path = '/home/bsorenson/data/OMI/H5_files/'

    # Grab the path date/time for finding the swath climatology files
    year = single_swath[0:4]
    date = single_swath[4:8]
    if(len(single_swath)==13):
        time = single_swath[9:]
    elif(len(single_swath)==12):
        time = single_swath[8:]
    else:
        time = ''

    total_climo_list = subprocess.check_output('ls '+climo_base_path+climo_date+'/*.he5',\
              shell=True).decode('utf-8').strip().split('\n')

    total_single_list = subprocess.check_output('ls '+single_base_path+'OMI-Aura_L2-OMAERUV_'+\
              year+'m'+date+'t'+time+'*.he5',shell=True).decode('utf-8').strip().split('\n')

    UVAI_climo = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_climo = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    UVAI_single = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_single = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

    print("Reading path climatology data")
    for fileI in range(len(total_climo_list)):
        # read in data directly from HDF5 files
        print(total_climo_list[fileI])
        data = h5py.File(total_climo_list[fileI],'r')
        #ALBEDO = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo']
        #REFLECTANCE = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity']
        #CLD = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/CloudFraction']
        AI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
        #PIXEL= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/PixelQualityFlags']
        #MSMNT= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/MeasurementQualityFlags']
        #VZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/ViewingZenithAngle']
        #SZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SolarZenithAngle']
        #RAZ  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle']
        LAT = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
        LON = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
        #GRND = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags']
    
        #albedo = ALBEDO[:,:,0]   
        #reflectance = REFLECTANCE[:,:,0]   
        counter = 0
        #AI = AI[:,:,0]   
        # Loop over the values and rows 
        #for i in range(0,int(CBA2)):
        #for i in range(albedo.shape[0]):
        for i in range(AI.shape[0]):
            for j in range(0,row_max):
                #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
                if(AI[i,j]>-20):
                #if(plotAI[i,j]>-20):
                    # Only plot if XTrack flag is met
                    if((XTRACK[i,j] == 0) | ((XTRACK[i,j] & 4 == 4))):
                        # Print values to text file
                        if(LAT[i,j] > minlat):
                            counter+=1
    #                        fout.write("{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.6f} {9:.6f} {10:.6f} {11:.6f}\n".format(\
    #                            LAT[i,j],LON[i,j],AI[i,j],0.5,SZA[i,j],VZA[i,j],RAZ[i,j], \
    #                            ALBEDO[i,j,0],ALBEDO[i,j,1],REFLECTANCE[i,j,0],\
    #                            REFLECTANCE[i,j,1],CLD[i,j]))
    
    
                        #if((plotXTrack[i,j] == 0) | ((plotXTrack[i,j] & 4 == 4))):
                            index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                            index2 = int(np.floor(LON[i,j]*4 + 720.))
                            #index1 = int(np.floor(plotLAT[i,j]*4 + 360.))
                            #index2 = int(np.floor(plotLON[i,j]*4 + 720.))
                            
                            if(index1 < 0): index1 = 0
                            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                            if(index2 < 0): index2 = 0                                                                                            
                            if(index2 > 1439): index2 = 1439
                       
                            #diff = reflectance[i,j] - albedo[i,j]
                            #if(diff<min_diff):
                            #    min_diff = diff
                            #if(diff>max_diff):
                            #    max_diff = diff
                       
                            #if(diff<0.2): 
                            #    UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
                            UVAI_climo[index2, index1] = (UVAI_climo[index2,index1]*count_climo[index2,index1] + \
                                AI[i,j])/(count_climo[index2,index1]+1)
                            #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                                #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                            count_climo[index2, index1] = count_climo[index2,index1] + 1
        data.close() 

    # Grab the single-pass data
    print("Reading single-swath data")
    for fileI in range(len(total_single_list)):
        # read in data directly from HDF5 files
        print(total_single_list[fileI])
        data2 = h5py.File(total_single_list[fileI],'r')
        AI  = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
        LAT = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
        LON = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
        XTRACK = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
    
        counter = 0
        # Loop over the values and rows 
        for i in range(AI.shape[0]):
            for j in range(0,row_max):
                if(AI[i,j]>-20):
                #if(plotAI[i,j]>-20):
                    # Only plot if XTrack flag is met
                    if((XTRACK[i,j] == 0) | ((XTRACK[i,j] & 4 == 4))):
                        # Print values to text file
                        if(LAT[i,j] > minlat):
                            counter+=1
    
                            index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                            index2 = int(np.floor(LON[i,j]*4 + 720.))
                            
                            if(index1 < 0): index1 = 0
                            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                            if(index2 < 0): index2 = 0                                                                                            
                            if(index2 > 1439): index2 = 1439
                       
                            UVAI_single[index2, index1] = (UVAI_single[index2,index1]*count_single[index2,index1] + \
                                AI[i,j])/(count_single[index2,index1]+1)
                            count_single[index2, index1] = count_single[index2,index1] + 1
        data2.close() 
    ##!## Calculate the row-average AI for the secondary plot
    ##!#mask_avgs = np.nanmean(np.ma.masked_where(AI[:,:] < -20, AI[:,:]),axis=0)

    # Calculate anomalies
    UVAI_anomaly = UVAI_single - UVAI_climo  

    plot_lat, plot_lon = np.meshgrid(lat_ranges,lon_ranges)
    mask_UVAI_climo = np.ma.masked_where(count_climo[:,:] == 0, UVAI_climo)
    mask_UVAI_anom = np.ma.masked_where(count_single[:,:] == 0, UVAI_anomaly)

    # Plot the climatology
    plt.close('all')
    fig1 = plt.figure(figsize=(8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
 
    plt.title('OMI Aerosol Index Climatology '+single_swath)
    #plt.title('OMI Reflectivity - Surface Albedo '+swath_date)
    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI_climo,transform = datacrs,cmap = colormap,\
            vmin = -2.0, vmax = 3.0)
    ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    #ax.set_xlim(-4170748.535086173,4167222.438879491)
    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label='Aerosol Index')
    ##cax = fig.add_axes([0.16,0.075,0.7,0.025])
    ##cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    ##cb.ax.set_xlabel('Aerosol Index')
    #cb.ax.set_xlabel('Reflectivity - Surface Albedo')
    #out_name = 'omi_single_pass_ai_200804270052_to_0549_composite_rows_0to'+str(row_max)+'.png'       
    out_name = 'omi_single_pass_ai_single_climo_'+single_swath+'_rows_0to'+str(row_max)+'.png'       
    #out_name = 'omi_single_pass_refl_albedo_diff_'+swath_date+'_rows_0to'+str(row_max)+'.png'       
    plt.savefig(out_name)
    print('Saved image '+out_name)
  
    # Plot the anomalies
    plt.close('all')
    fig1 = plt.figure(figsize=(8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
 
    plt.title('OMI Aerosol Index Anomaly '+single_swath)
    #plt.title('OMI Reflectivity - Surface Albedo '+swath_date)
    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI_anom,transform = datacrs,cmap = colormap,\
            vmin = -1.0, vmax = 1.5)
    ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    #ax.set_xlim(-4170748.535086173,4167222.438879491)
    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label='Aerosol Index Anomaly')
    ##cax = fig.add_axes([0.16,0.075,0.7,0.025])
    ##cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    ##cb.ax.set_xlabel('Aerosol Index')
    #cb.ax.set_xlabel('Reflectivity - Surface Albedo')
    #out_name = 'omi_single_pass_ai_200804270052_to_0549_composite_rows_0to'+str(row_max)+'.png'       
    out_name = 'omi_single_pass_ai_single_anom_'+single_swath+'_rows_0to'+str(row_max)+'.png'       
    #out_name = 'omi_single_pass_refl_albedo_diff_'+swath_date+'_rows_0to'+str(row_max)+'.png'       
    plt.savefig(out_name)
    print('Saved image '+out_name)
    
    #axs[1].plot(mask_avgs)
    #axs[1].set_xlabel('Sensor Row')
    #axs[1].set_ylabel('Row Average Aerosol Index')
    
    plt.show()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Comparison with CERES
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def plot_OMI_CERES_trend_compare(OMI_data, CERES_data,month,minlat=65,\
        trend_type = 'standard',save=False):

    if('/home/bsorenson/Research/CERES' not in sys.path):
        sys.path.append('/home/bsorenson/Reserach/CERES')
    #importlib.reload(gridCERESLib) 
    from gridCERESLib import calcCERES_grid_trend, plotCERES_spatial, \
        plotCERES_MonthTrend

    # Calculate trends for both the OMI and CERES data
    # ------------------------------------------------
    OMI_trend   = calcOMI_grid_trend(OMI_data, month, trend_type, \
        minlat)
    CERES_trend = calcCERES_grid_trend(CERES_data, month, trend_type, \
        minlat)

    print('before shapes', OMI_trend.shape, CERES_trend.shape)
    #for x, y in zip(CERES_data['lat'][:, 10], CERES_trend[:, 10]):
    #    print(x, y)

    # Convert the index to a string using datetime
    if(month != None):
        dt_obj = datetime.strptime(OMI_data['DATES'][month],"%Y%m")
        title = 'OMI AI / CERES ' + CERES_data['param'] + '\n'+ \
            dt_obj.strftime("%b") + " Trend Comparison"
        outname = 'omi_ceres_trend_comp_'+dt_obj.strftime("%b")+'_'+\
            OMI_data['VERSION']+'vCERES_min'+str(int(minlat))+'.png'
    else:
        title = 'OMI AI / CERES ' + CERES_data['param'] + ' Trend Comparison'
        outname = 'omi_ceres_trend_comp_'+\
            OMI_data1['VERSION']+'vCERES_min' + str(int(minlat)) + '.png'

    # Flip the CERES data to convert the longitudes from 0 - 360 to -180 - 180
    # ------------------------------------------------------------------------
    local_lon = np.copy(CERES_data['lon'])
    local_lat = np.copy(CERES_data['lat'])
    over_180  = np.where(CERES_data['lon'][0,:] >= 179.9999)
    under_180 = np.where(CERES_data['lon'][0,:] < 179.9999)

    for ii in range(local_lon.shape[0]):
        local_lon[ii,:] = np.concatenate([local_lon[ii,:][over_180] - 360.,\
            local_lon[ii,:][under_180]])
        local_lat[ii,:] = np.concatenate([local_lat[ii,:][over_180],\
            local_lat[ii,:][under_180]])
        CERES_trend[ii,:] = np.concatenate([CERES_trend[ii,:][over_180],\
            CERES_trend[ii,:][under_180]])
  
    print("Before removal, ",OMI_trend.shape)
 
    # First, mask any OMI data that are outside the bounds of the CERES data
    # Assume that the OMI data are larger than the CERES data
    # ----------------------------------------------------------------------
    OMI_trend[np.where(np.isnan(OMI_trend))] = -999.
    where_matching = np.where((OMI_data['LAT'][:,0] >= \
        np.min(local_lat) - 0.5 ) & (OMI_data['LAT'][:,0] <= \
        np.max(local_lat) - 0.5) )[0]
    #where_matching = np.where((OMI_data['LAT'][:,0] >= \
    #    (np.min(local_lat) - 0.5)) & (OMI_data['LAT'][:,0] <= \
    #    (np.max(local_lat) - 0.5)))[0]
    OMI_trend = OMI_trend[where_matching, :]

    mask_trend1 = np.array(OMI_trend[  (OMI_trend != 0) & (CERES_trend != 0) \
        & (OMI_trend != -999.) & (CERES_trend != -999.)])
    mask_trend2 = np.array(CERES_trend[(OMI_trend != 0) & (CERES_trend != 0) \
        & (OMI_trend != -999.) & (CERES_trend != -999.)])

    print('after shapes', mask_trend1.shape, mask_trend2.shape)
    print("average OMI trend", np.nanmean(mask_trend1))
    print("average CERES trend", np.nanmean(mask_trend2))

    print("Pearson:  ",pearsonr(mask_trend1,mask_trend2))
    print("Spearman: ",spearmanr(mask_trend1,mask_trend2))

    xy = np.vstack([mask_trend1,mask_trend2])
    z = stats.gaussian_kde(xy)(xy)

    # Set up the figure. 
    # 3 panels: OMI trend, CERES trends, scatter compare
    # ---------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (14,5))
    ax0 = fig1.add_subplot(1,3,1, projection = mapcrs)
    ax1 = fig1.add_subplot(1,3,2, projection = mapcrs)
    ax2 = fig1.add_subplot(1,3,3)

    # Plot the OMI and CERES trends
    # -----------------------------
    #plotOMI_spatial(ax0, OMI_data['LAT'], OMI_data['LON'], OMI_trend, 'trend', \
    #    ptitle = title, plabel = 'UV Aerosol Index Trend', \
    #    vmin = None, vmax = None, minlat = minlat)
    plotOMI_MonthTrend(OMI_data,month_idx=month,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=minlat,title = "OMI Trend", pax = ax0)
    plotCERES_MonthTrend(CERES_data,month_idx=month,\
        trend_type=trend_type,\
        minlat=minlat,pax = ax1)

    ##plotCERES_spatial(ax1, CERES_data['lat'], CERES_data['lon'], CERES_trend, 'trend', \
    ##    ptitle = title, plabel = 'CERES ' + CERES_data['param'] + ' Trend', \
    ##    vmin = None, vmax = None, minlat = minlat)

    

    ax2.scatter(mask_trend1,mask_trend2,c=z,s=8)
    plot_trend_line(ax2, mask_trend1, mask_trend2, color='tab:green', linestyle = '-', \
        slope = 'both')
    ##!#plt.plot(test_x,predictions,color='tab:green',linestyle='--',label='Huber Fit')
    ##!## Plot an unrobust fit line using linear regression
    ##!## -------------------------------------------------
    ##!#plt.plot(np.unique(mask_trend1),np.poly1d(np.polyfit(mask_trend1,\
    ##!#    mask_trend2,1))(np.unique(mask_trend1)),color='tab:orange',\
    ##!#    linestyle='--',label='Polyfit Fit')

    ##if((month == 0) | (month == 1)):
    ##    plt.xlim(-0.5,0.3)
    ##    plt.ylim(-0.5,0.3)
    ##elif((month == 2)):
    ##    plt.xlim(-0.6,0.5)
    ##    plt.ylim(-0.6,0.5)
    ##elif((month == 3)):
    ##    plt.xlim(-0.5,0.5)
    ##    plt.ylim(-0.5,0.5)
    ##elif((month == 4)):
    ##    plt.xlim(-0.5,0.7)
    ##    plt.ylim(-0.5,0.7)
    ##elif((month == 5)):
    ##    plt.xlim(-0.5,0.5)
    ##    plt.ylim(-0.5,0.5)
    ##else:
    ##    plt.xlim(-0.3,0.3)
    ##    plt.ylim(-0.3,0.3)
    ax2.legend()
    ax2.set_xlabel(OMI_data['VERSION'])
    ax2.set_ylabel('CERES')
    ax2.set_title(title + '\nMinlat = ' + str(minlat))

    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')
    plot_subplot_label(ax2, '(c)')

    fig1.tight_layout()

    if(save == True):
        fig1.savefig(outname)
        print("Saved image",outname)
    else:
        plt.show()
    #return mask_trend1,mask_trend2

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Comparison with CERES
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
