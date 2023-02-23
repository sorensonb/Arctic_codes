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

import matplotlib as mpl
import numpy as np
import sys
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
#import matplotlib.colors as cm
import matplotlib.colors as mc
import matplotlib.cm as cm
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.path as mpath
import matplotlib.patches as mpatches
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
from scipy.signal import argrelextrema, find_peaks
from glob import glob
import os
import importlib

home_dir = os.environ['HOME']

sys.path.append(home_dir)
from python_lib import circle, plot_subplot_label, plot_lat_circles
from python_lib import plot_trend_line, plot_subplot_label, plot_figure_text, \
    nearest_gridpoint, aerosol_event_dict, init_proj, \
    convert_radiance_to_temp, format_coord, circle, plot_point_on_map

# Bits 
#  0 - Missing
#  1 - Bad Pixel
#  2 - Processing Error
#  3 - Transient Pixel Warning
#  4 - RTS Pixel Warning
#  5 - Saturation Possibility Warning
#  6 - Noise Calculation Warning
#  7 - Dark Current Warning
#  8 - Offset Warning
#  9 - Exposure Smear Warning
#  10 - Stray Light Warning
#  11 - 13 Reserved 
#  14 - Dead Pixel Identifcation (OML1BCAL only)
#  15 - Dead Pixel Identifcation Error (OML1BCAL only)
def get_pxqf_flags(value, flag_idx = None):
    bit_vals = np.array([int(format(value,"016b")[-idx],2) \
        for idx in range(1, 17)])
    if(flag_idx == None):
        return bit_vals
    else:
        return bit_vals[flag_idx]

# Bits 0 to 3 together contain the land/water flags:
#      0 - shallow ocean
#      1 - land
#      2 - shallow inland water
#      3 - ocean coastline/lake shoreline
#      4 - ephemeral (intermittent) water
#      5 - deep inland water
#      6 - continental shelf ocean
#      7 - deep ocean
#      8-14 - not used
#      15 - error flag for land/water
#  Bits 4 to 6 are flags that are set to 0 for FALSE, or 1 for TRUE:
#      Bit 4 - sun glint possibility flag
#      Bit 5 - solar eclipse possibility flag
#      Bit 6 - geolocation error flag
#  Bit 7 is reserved for future use (currently set to 0)
#  Bits 8 to 14 together contain the snow/ice flags (based on NISE):
#      0 - snow-free land
#      1-100 - sea ice concentration (percent)
#      101 - permanent ice (Greenland, Antarctica)
#      102 - not used
#      103 - dry snow
#      104 - ocean (NISE-255)
#      105-123 - reserved for future use
#      124 - mixed pixels at coastline (NISE-252)
#      125 - suspect ice value (NISE-253)
#      126 - corners undefined (NISE-254)
#      127 - error
#  Bit 15 - NISE nearest neighbor filling flag
#  (See Section 6.2 of Reference 4 for more details.)
def get_gpqf_flags(value, flag_idx = None):
    bit_vals = np.full(7,np.nan)
    bit_vals[0] = int(format(value,"016b")[-4:],2) # land/water flags
    bit_vals[1] = int(format(value,"016b")[-5],2)  # sun glint possiblity
    bit_vals[2] = int(format(value,"016b")[-6],2)  # solar eclipse possiblity
    bit_vals[3] = int(format(value,"016b")[-7],2)  # geolocation error flag
    bit_vals[4] = int(format(value,"016b")[-8],2)  # future use
    bit_vals[5] = int(format(value,"016b")[-15:-8],2) # NISE snow/ice flags
    bit_vals[6] = int(format(value,"016b")[-16],2) # 
    if(flag_idx == None):
        return bit_vals
    else:
        return bit_vals[flag_idx]

def get_ice_flags(value):
    return get_gpqf_flags(value, flag_idx = 5)
    ##!#return int(format(value,"016b")[-15:-8],2)

# The cross track quality flags assigned to each pixel in
#  OMI L1B data. Flags indicate detection of the OMI row
#  anomaly and if the effect has been corrected.
#  Bits 0 to 2 together indicate row anomaly status:
#  0 - Not affected
#  1 - Affected, Not corrected, do not use
#  2 - Slightly affected, not corrected, use with caution
#  3 - Affected, corrected, use with caution
#  4 - Affected, corrected, use pixel
#  5 - Not used
#  6 - Not used
#  7 - Error during anomaly detection processing
#  Bit 3 - Reserved for future use.
#  Bit 4 - Possibly affected by wavelength shift
#  Bit 5 - Possibly affected by blockage
#  Bit 6 - Possibly affected by stray sunlight
#  Bit 7 - Possibly affected by stray earthshine
def get_xtrack_flags(value, flag_idx = None):
    bit_vals = np.full(6,np.nan)
    bit_vals[0] = int(format(value,"08b")[-3:],2) # row anomaly status
    bit_vals[1] = int(format(value,"08b")[-4],2)  # future use
    bit_vals[2] = int(format(value,"08b")[-5],2)  # possible wavelength shift
    bit_vals[3] = int(format(value,"08b")[-6],2)  # possible blockage
    bit_vals[4] = int(format(value,"08b")[-7],2)  # possible stray sunlight
    bit_vals[5] = int(format(value,"08b")[-8],2)  # possible stray earthshine
    if(flag_idx == None):
        return bit_vals
    else:
        return bit_vals[flag_idx]

##!## Compute a circle in axes coordinates, which we can use as a boundary
##!## for the map. We can pan/zoom as much as we like - the boundary will be
##!## permanently circular.
##!#theta = np.linspace(0, 2*np.pi, 100)
##!#center, radius = [0.5, 0.5], 0.5
##!#verts = np.vstack([np.sin(theta), np.cos(theta)]).T
##!#circle = mpath.Path(verts * radius + center)

# Set up mapping variables 
datacrs = ccrs.PlateCarree() 
colormap = plt.cm.jet
mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

label_dict = {
    'V003': 'Only XTrack == 0',
    'VJZ2': 'No Snow-free Land', # no snow-free land, all available rows
    'VJZ28': 'No Snow-free Land\nOnly 49, 50, 56 - 60', # 49, 50, 56 - 60, no snow-free land
    'VJZ282': 'No Snow-free Land, Only Pure Good Rows',
    'VJZ29': 'Include Snow-free Land\n49, 50, 56 - 60', #49, 50, 56 - 60, snow-free land
    'VJZ211': 'Include Snow-free Land\n56 - 60', # 55 - 60, snow-free land
    'VJZ212': 'Include Snow-free Land\n56 - 60, CLD < 0.2', # 55 - 60, snow-free land
    'VJZ4': 'XTrack == 0, not 4',
    'VJZ5': 'AI >= 0',
    'VBS0': 'No Bad Row Screening',
    'VBS1': 'Bad Row Screening Only',
    'VBS2': 'Only Rows 1-22',
    'VSJ2': 'Perturbation Version 2',
    'VSJ4': 'Perturbation Version 4',
    'VSJ42': 'Perturbation Version 4, CLD < 0.2'
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

sea_dict = {
    'esib': {
        'lat': np.array([90., 50, 50, 90.]), 
        'lon': np.array([176., 176., 145., 145.])
    },
    'laptev': {
        'lat': np.array([90., 50, 50, 90.]),
        'lon': np.array([145., 145., 100., 100.])
    },
    'barents': {
        'lat': np.array([90., 50, 50, 90.]),
        'lon': np.array([100., 100., 55., 55.])
    },
    'kara': {
        'lat': np.array([90., 50, 50, 90.]),
        'lon': np.array([55., 55., 12., 12.])
    },\
    'greenland': {
        'lat': np.array([90., 50, 50, 90.]),
        'lon': np.array([12., 12., -45., -45.])
    },\
    'baffin': {
        'lat': np.array([90., 50, 50, 90.]),
        'lon': np.array([-45., -45., -80., -80.])
    },\
    'canada': {
        'lat': np.array([90., 50, 50, 90.]),
        'lon': np.array([-80., -80., -125., -125.])
    },\
    'beaufort': {
        'lat': np.array([90., 50, 50, 90.]),
        'lon': np.array([-125., -125., -156., -156.])
    },
    'chukchi': {
        'lat': np.array([90., 50, 50, 90.]), 
        'lon': np.array([-184., -184., -156., -156.])
    }
}

sea_double_dict = {
    'esib_laptev': {
        'lat': np.array([90., 50, 50, 90.]), 
        'lon': np.array([176., 176., 100., 100.])
    },
    'kara_barents': {
        'lat': np.array([90., 50, 50, 90.]),
        'lon': np.array([100., 100., 12., 12.])
    },
    'greenland_baffin': {
        'lat': np.array([90., 50, 50, 90.]),
        'lon': np.array([12., 12., -80., -80.])
    },\
    'canada_beaufort_chukchi': {
        'lat': np.array([90., 50, 50, 90.]),
        'lon': np.array([-80., -80., 176., 176.])
    }
}


def init_sea_poly(in_sea, color = 'k', linewidth = 2):

    ##!#poly_corners = np.zeros((len(sea_dict[in_sea]['lat']), 2), np.float64)
    ##!#poly_corners[:,0] = sea_dict[in_sea]['lon']
    ##!#poly_corners[:,1] = sea_dict[in_sea]['lat']
    poly_corners = np.zeros((len(sea_double_dict[in_sea]['lat']), 2), np.float64)
    poly_corners[:,0] = sea_double_dict[in_sea]['lon']
    poly_corners[:,1] = sea_double_dict[in_sea]['lat']
    poly = mpatches.Polygon(poly_corners, closed = True, ec=color, \
        fill = False, lw = linewidth, fc = None, transform = ccrs.PlateCarree())
    return poly
   
def plot_arctic_regions(pax, linewidth = 2):
    ###pax.add_patch(init_sea_poly('chukchi', linewidth = linewidth))
    ###pax.add_patch(init_sea_poly('esib', linewidth = linewidth))
    ###pax.add_patch(init_sea_poly('laptev', linewidth = linewidth))
    ###pax.add_patch(init_sea_poly('barents', linewidth = linewidth))
    ###pax.add_patch(init_sea_poly('kara', linewidth = linewidth))
    ###pax.add_patch(init_sea_poly('greenland', linewidth = linewidth))
    ###pax.add_patch(init_sea_poly('baffin', linewidth = linewidth))
    ###pax.add_patch(init_sea_poly('canada', linewidth = linewidth))
    ###pax.add_patch(init_sea_poly('beaufort', linewidth = linewidth))
 
    pax.add_patch(init_sea_poly('esib_laptev', linewidth = linewidth))
    pax.add_patch(init_sea_poly('kara_barents', linewidth = linewidth))
    pax.add_patch(init_sea_poly('greenland_baffin', linewidth = linewidth))
    pax.add_patch(init_sea_poly('canada_beaufort_chukchi', linewidth = linewidth))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Miscellaneous functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

##!#def plot_subplot_label(ax, label, xval = None, yval = None, transform = None, \
##!#        color = 'black', backgroundcolor = None, fontsize = 14, \
##!#        location = 'upper_left'):
##!#
##!#    if(location == 'upper_left'):
##!#        y_lim = 0.90
##!#        x_lim = 0.05
##!#    elif(location == 'lower_left'):
##!#        y_lim = 0.05
##!#        x_lim = 0.05
##!#    elif(location == 'upper_right'):
##!#        y_lim = 0.90
##!#        x_lim = 0.90
##!#    elif(location == 'lower_right'):
##!#        y_lim = 0.05
##!#        x_lim = 0.90
##!#
##!#    if(xval is None):
##!#        xval = ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * x_lim
##!#    if(yval is None):
##!#        yval = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * y_lim
##!#    print('Xval = ',xval, 'Yval = ',yval)
##!#
##!#    if(transform is None):
##!#        if(backgroundcolor is None):
##!#            ax.text(xval,yval,label, \
##!#                color=color, weight='bold', \
##!#                fontsize=fontsize)
##!#        else:
##!#            ax.text(xval,yval,label, \
##!#                color=color, weight='bold', \
##!#                fontsize=fontsize, backgroundcolor = backgroundcolor)
##!#    else:
##!#        if(backgroundcolor is None):
##!#            ax.text(xval,yval,label, \
##!#                color=color, weight='bold', \
##!#                transform = transform, fontsize=fontsize)
##!#        else:
##!#            ax.text(xval,yval,label, \
##!#                color=color, weight='bold', \
##!#                transform = transform, fontsize=fontsize, \
##!#                backgroundcolor = backgroundcolor)

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

##!## Plots latitude circles on a given axis
##!## --------------------------------------
##!#def plot_lat_circles(pax, lat_circles):
##!#
##!#    if(len(lat_circles) > 5):
##!#        print("WARNING: more than 5 latitude circles selected")
##!#        print("    Only plotting the first 5")
##!#        lat_circles = lat_circles[:5]
##!#
##!#    colors = ['tab:blue','tab:red','tab:purple','tab:olive','tab:cyan']
##!#    lon = np.arange(-180, 181)
##!#    for ii, lat_val in enumerate(lat_circles):
##!#        lats = np.full(lon.shape, lat_val)
##!#
##!#        pax.plot(lon, lats, linewidth = 2.5, transform = datacrs, color = 'k')
##!#        pax.plot(lon, lats, transform = datacrs, color = colors[ii])
##!#        pax.text(-180, lat_val + 3, str(int(lat_val)) + ' $^{o}$N', \
##!#            horizontalalignment = 'center', weight = 'bold', \
##!#            color = colors[ii], transform = datacrs)



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
        if((len(templine)>1) & (templine[0] != 'Date')):
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
#     skiprows - array-like containing row indices (0-based) to not include
#                in the data
def readOMI_swath_hdf(plot_time, dtype, only_sea_ice = False, \
        only_ice = False, no_ice = False, latmin = 65, skiprows = None):

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

    base_path = home_dir + '/data/OMI/H5_files/'
    total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+\
        year+'m'+date+'t'+time+'*.he5',shell=True).decode('utf-8').strip().split('\n')

    print(total_list[0])
    data = h5py.File(total_list[0],'r')

    LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
    LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
    XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
    UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
    GPQF   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags'][:,:]
    PXQF   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/PixelQualityFlags'][:,:]
    AZM    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle'][:,:]
    LATcrnr = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/FoV75CornerLatitude'][:,:,:]
    LONcrnr = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/FoV75CornerLongitude'][:,:,:]
    TIME_O = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Time'][:]
    GPQF_decode = np.array([[get_gpqf_flags(val) for val in GPQFrow] \
        for GPQFrow in GPQF])
    XTRK_decode = np.array([[get_xtrack_flags(val) for val in XTRKrow] \
        for XTRKrow in XTRACK])
    PXQF_decode = np.array([[[get_pxqf_flags(val) for val in PXQFline] \
        for PXQFline in PXQFrow] for PXQFrow in PXQF])

    # Process the time
    # ----------------
    base_date_str = '199301010000'
    base_date = datetime.strptime(base_date_str, '%Y%m%d%H%M')
    TIME = np.array([base_date + timedelta(seconds = time) for time in TIME_O])

    base_date_str2 = '197001010000'
    base_date2 = datetime.strptime(base_date_str2, '%Y%m%d%H%M')
    TIME_2 = np.array([(time - base_date2).days + \
        ((time - base_date2).seconds / 86400.) for time in TIME])

    # If the user wants certain lines to not be included, set the AI values
    # to -9e5 to ensure they are removed
    # ---------------------------------------------------------------------
    if(skiprows is not None):
        UVAI[:,skiprows] = -9e5  

    # Mask the AI data where the values are either missing or the Xtrack flags
    # are not zero
    # ------------------------------------------------------------------------
    mask_UVAI = np.ma.masked_where((XTRACK != 0) | \
        (UVAI < -2e5) | (LAT < latmin), UVAI)

    if(only_sea_ice):
        mask_UVAI = np.ma.masked_where((GPQF_decode[:,:,5] < 1) | \
            (GPQF_decode[:,:,5] > 100), mask_UVAI)
    elif(only_ice):
        mask_UVAI = np.ma.masked_where((GPQF_decode[:,:,5] < 1) | \
            (GPQF_decode[:,:,5] > 101), mask_UVAI)
    elif(no_ice):
        mask_UVAI = np.ma.masked_where((GPQF_decode[:,:,5] >= 1) & \
            (GPQF_decode[:,:,5] < 103), mask_UVAI)

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

    OMI_swath['LAT']    = LAT
    OMI_swath['LON']    = LON
    OMI_swath['AZM']    = AZM
    OMI_swath['GPQF']   = GPQF_decode
    OMI_swath['PXQF']   = PXQF_decode
    OMI_swath['XTRACK'] = XTRK_decode
    OMI_swath['XTRACK_raw'] = XTRACK
    OMI_swath['LATcrnr']    = LATcrnr
    OMI_swath['LONcrnr']    = LONcrnr
    OMI_swath['TIME']   = np.repeat(TIME, 60).reshape(TIME.shape[0], 60)
    OMI_swath['TIME_2'] = np.repeat(TIME_2, 60).reshape(TIME_2.shape[0], 60)
    OMI_swath['base_date'] = base_date_str

    if(dtype == 'JZ'):

        mask_UVAI[~((((GPQF_decode[:,:,5] >= 0) & \
            (GPQF_decode[:,:,5] <= 101)) | (GPQF_decode[:,:,5] == 104)) & \
                      (AZM > 100))] = np.nan
        mask_UVAI = np.ma.masked_invalid(mask_UVAI)

    OMI_swath['UVAI'] = mask_UVAI

    data.close()

    return OMI_swath

def readOMI_swath_shawn(plot_time, latmin = 65., \
        shawn_path = home_dir + '/data/OMI/shawn_files/'):

    # Convert the year, month, and day of the base date to
    # datetime
    # ----------------------------------------------------
    base_date = datetime.strptime(plot_time[:8], '%Y%m%d')

    # Read in the OMI data for this file, for matching purposes
    # ---------------------------------------------------------
    OMI_data = readOMI_swath_hdf(plot_time, 'control', latmin = latmin)

    # Open the shawn file
    # -------------------
    print(shawn_path + plot_time)
    df_shawn = pd.read_csv(shawn_path + plot_time, delim_whitespace = True, \
        header = None)

    # Extract the needed variables
    # ----------------------------
    s_lat    = pd.to_numeric(df_shawn[0], errors = 'coerce').values
    s_lon    = pd.to_numeric(df_shawn[1], errors = 'coerce').values
    s_airaw  = pd.to_numeric(df_shawn[2], errors = 'coerce').values
    s_aiclim = pd.to_numeric(df_shawn[3], errors = 'coerce').values
    s_aipert = pd.to_numeric(df_shawn[4], errors = 'coerce').values
    #s_time   = pd.to_numeric(df_shawn[5], errors = 'coerce').values * 86400.
    s_sza    = pd.to_numeric(df_shawn[6], errors = 'coerce').values
    s_vza    = pd.to_numeric(df_shawn[7], errors = 'coerce').values
    s_raz    = pd.to_numeric(df_shawn[8], errors = 'coerce').values

    # Convert the SZA and VZA from the cosine value to the actual
    # value in degrees
    # -----------------------------------------------------------
    s_sza = np.degrees(np.arccos(s_sza))
    s_vza = np.degrees(np.arccos(s_vza))

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
    #sfile_TIME   = np.full(OMI_data['LAT'].shape, np.nan)
    sfile_SZA    = np.full(OMI_data['LAT'].shape, np.nan)
    sfile_VZA    = np.full(OMI_data['LAT'].shape, np.nan)
    sfile_RAZ    = np.full(OMI_data['LAT'].shape, np.nan)

    for slat, slon, sair, saic, saip, ssza, svza, sraz, in \
            zip(s_lat, s_lon, s_airaw, s_aiclim, s_aipert, \
                s_sza, s_vza, s_raz):
        match_loc = np.where((slat == OMI_data['LAT']) & \
            (slon == OMI_data['LON']))
        #print(slat, slon, sair, OMI_data['LAT'][match_loc], \
        #    OMI_data['LON'][match_loc], OMI_data['UVAI'][match_loc], match_loc)
        #sfile_LAT[match_loc] = slat
        #sfile_LON[match_loc] = slon
        sfile_AIraw[match_loc] = sair
        sfile_AIclim[match_loc] = saic
        sfile_AIpert[match_loc] = saip
        #print(base_date + timedelta(seconds = stime + 10), OMI_data['TIME'][match_loc])

        #sfile_TIME[match_loc]   = base_date + timedelta(stime)
        sfile_SZA[match_loc]    = ssza
        sfile_VZA[match_loc]    = svza
        sfile_RAZ[match_loc]    = sraz

    ## Convert the seconds per day to absolute time
    ## --------------------------------------------
    #sfile_TIME[~np.isnan(sfile_TIME)] = np.array([base_date + timedelta(seconds = stime2) \
    #    for stime2 in sfile_TIME[~np.isnan(sfile_TIME)]])

    #sfile_TIME = np.array([[base_date + timedelta(seconds = stime2) \
    #    for stime2 in stime1] for stime1 in sfile_TIME])

    OMI_swath = {}
    OMI_swath['date'] = plot_time
    OMI_swath['dtype'] = 'shawn'
    OMI_swath['LAT'] = OMI_data['LAT']
    OMI_swath['LON'] = OMI_data['LON']
    OMI_swath['LATcrnr'] = OMI_data['LATcrnr']
    OMI_swath['LONcrnr'] = OMI_data['LONcrnr']
    #OMI_swath['TIME'] = sfile_TIME
    #OMI_swath['TIME'] = abs_time
    OMI_swath['TIME']   = OMI_data['TIME']
    OMI_swath['TIME_2'] = OMI_data['TIME_2']
    OMI_swath['UVAI_raw']   = np.ma.masked_invalid(sfile_AIraw)
    OMI_swath['UVAI_climo'] = np.ma.masked_invalid(sfile_AIclim)
    OMI_swath['UVAI_pert']  = np.ma.masked_invalid(sfile_AIpert)
    OMI_swath['SZA'] = np.ma.masked_invalid(sfile_SZA)
    OMI_swath['VZA'] = np.ma.masked_invalid(sfile_VZA)
    OMI_swath['RAZ'] = np.ma.masked_invalid(sfile_RAZ)
   
    # Mask the data south of the desired latmin
    # -----------------------------------------
    OMI_swath['UVAI_raw']   = np.ma.masked_where(OMI_data['LAT'] < latmin, sfile_AIraw)
    OMI_swath['UVAI_climo'] = np.ma.masked_where(OMI_data['LAT'] < latmin, sfile_AIclim)
    OMI_swath['UVAI_pert']  = np.ma.masked_where(OMI_data['LAT'] < latmin, sfile_AIpert)

    return OMI_swath

# NOTE: this is the old way of reading Shawn data. Grids everything
# into the 0.25 x 0.25 degree grid
def readOMI_swath_shawn_old(plot_time, latmin = 65., resolution = 0.25):

    # This is the path that points to the HDF5 OMI files. This must be changed
    # if running on a new system.
    base_path = home_dir + '/data/OMI/shawn_files/'
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #
    # Select the needed data files
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
    # Extract date information to use when finding the data files
    year = plot_time[:4]
    date = plot_time[4:8]
    if(len(plot_time)==12):
        time = plot_time[8:12]
    else:
        time = ''
    
    total_list = subprocess.check_output('ls '+base_path+year+date+time+\
        '*',shell=True).decode('utf-8').strip().split('\n')
    
    # Set up values for gridding the AI data
    latmin = 65 
    lat_gridder = latmin * (1. / resolution)
    max_idx = int(360 / resolution)
    multiplier = int(1 / resolution)
    adder = int(180 / resolution)
    
    lat_ranges = np.arange(latmin,90, resolution)
    lon_ranges = np.arange(-180,180, resolution)
    
    # Set up blank grid arrays to hold the counts and the data
    UVAI = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #
    # Grid the data
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    for fileI in range(len(total_list)):
        # read in data directly from HDF5 files
        print(total_list[fileI])
        with open(total_list[fileI]) as infile:
            for line in infile:
                templine = line.strip().split()
                lat     = float(templine[0])
                lon     = float(templine[1])
                cleanAI = float(templine[4])
                if((cleanAI>-2e5)):
                    #if((j != 52) & (PDATA[i,j]>-2e5)):
                    # Only plot if XTrack flag is met
                    if(lat > latmin):
    
                        index1 = int(np.floor(lat*multiplier - lat_gridder))
                        index2 = int(np.floor(lon*multiplier + adder))
                          
                        if(index1 < 0): index1 = 0
                        if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                        if(index2 < 0): index2 = 0                                                                                            
                        if(index2 > (max_idx - 1)): index2 = max_idx - 1
                        #if(index2 > 1439): index2 = 1439
    
                        UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + cleanAI)/(count[index2,index1]+1)
                        count[index2, index1] = count[index2,index1] + 1
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    #plot_lat, plot_lon = np.meshgrid(lat_ranges,lon_ranges)
    mask_UVAI = np.ma.masked_where(count == 0, UVAI)

    OMI_data = {}
    OMI_data['date']      = plot_time
    OMI_data['dtype']     = 'shawn'
    OMI_data['LAT']       = lat_ranges
    OMI_data['LON']       = lon_ranges
    #OMI_data['LAT']       = plot_lat 
    #OMI_data['LON']       = plot_lon
    OMI_data['AI'] = mask_UVAI
    OMI_data['counts']    = count

    return OMI_data

# Reads a single swath of OMI data using the old method
# -----------------------------------------------------
def readOMI_single_swath(plot_time,row_max, row_min = 0, only_sea_ice = True,\
        latmin=65, resolution = 0.25, dtype = 'control', coccolith = False):


    if(dtype == 'JZ'):
        row_min = 56

    # Set up values for gridding the AI data
    lat_gridder = latmin * (1. / resolution)
    max_idx = int(360 / resolution)
    multiplier = int(1 / resolution)
    adder = int(180 / resolution)
    
    lat_ranges = np.arange(latmin,90., resolution)
    lon_ranges = np.arange(-180.,180., resolution)


    if(coccolith):
        g_NRAD_354     = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        g_NRAD_388     = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        g_NRAD_500     = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        g_REFL_354     = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        g_REFL_388     = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        g_SALB_354     = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_NRAD_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_NRAD_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_NRAD_500 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_REFL_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_REFL_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
        count_SALB_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    g_UVAI_354      = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_UVAI_354  = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    g_SZA           = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_SZA       = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

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
    base_path = home_dir + '/data/OMI/H5_files/'
    total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+year+'m'+date+'t'+time+'*.he5',\
              shell=True).decode('utf-8').strip().split('\n')

    for fileI in range(len(total_list)):
        # read in data directly from HDF5 files
        data = h5py.File(total_list[fileI],'r')
        print(total_list[fileI])
    
        LAT     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
        LON     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
        UVAI    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
        XTRACK  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
        GPQF    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags'][:,:]
        SZA     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SolarZenithAngle'][:,:]
        if(coccolith):
            NRAD    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/NormRadiance'][:,:]
            REFL    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity'][:,:]
            SALB    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo'][:,:]
    
        #albedo = ALBEDO[:,:,0]   
        #reflectance = REFLECTANCE[:,:,0]   
        counter = 0
        #AI = AI[:,:,0]   
        # Loop over the values and rows 
        #for i in range(0,int(CBA2)):
        #for i in range(albedo.shape[0]):
        for i in range(UVAI.shape[0]):
            for j in range(row_min,row_max):
                if(j == 52):
                    continue
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



                        index1 = int(np.floor(LAT[i,j]*multiplier - lat_gridder))
                        index2 = int(np.floor(LON[i,j]*multiplier + adder))
                          
                        if(index1 < 0): index1 = 0
                        if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                        if(index2 < 0): index2 = 0                                                                                            
                        if(index2 > (max_idx - 1)): index2 = max_idx - 1

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

        data.close()
  
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

    g_UVAI_354     = np.ma.masked_where(count_UVAI_354 == 0, g_UVAI_354)

    OMI_single_dict = {}
    OMI_single_dict['LAT'] = lat_ranges
    OMI_single_dict['LON'] = lon_ranges
    OMI_single_dict['AI'] = g_UVAI_354
    OMI_single_dict['AI_count'] = count_UVAI_354
    OMI_single_dict['SZA'] = g_SZA
    OMI_single_dict['SZA_count'] = count_SZA
    OMI_single_dict['date'] = plot_time
    OMI_single_dict['row_max'] = row_max
    OMI_single_dict['dtype'] = dtype
    if(coccolith):
        OMI_single_dict['NRAD500'] = g_NRAD_500
        OMI_single_dict['NRAD500_count'] = count_NRAD_500
        OMI_single_dict['AI_algae'] = mask_UVAI354
        OMI_single_dict['NRAD500_algae'] = mask_NRAD500

    return OMI_single_dict

def readOMI_NCDF(infile=home_dir + '/Research/OMI/omi_ai_V003_2005_2020.nc',\
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

def writeOMI_toNCDF(OMI_data,file_name, minlat = 65.):
    lat_ranges = np.arange(minlat,90.0,1.0)
    lon_ranges = np.arange(-180.,180.,1.0)
   
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
    MONTH = nc.createVariable('MONTH','i2',('nmth'))
    MONTH.description = 'Months since January 2005'
    LAT = nc.createVariable('Latitude','f4',('dlat','dlon'))
    LAT.description = 'Latitude (degrees North)'
    LAT.units = 'Degrees'
    LON = nc.createVariable('Longitude','f4',('dlat','dlon'))
    LON.description = 'Longitude (-180 - 180)'
    LON.units = 'Degrees'
   
    # Create a variable for the AI data, dimensioned using all three
    # dimensions.
    # --------------------------------------------------------------  
    AI = nc.createVariable('AI','f4',('nmth','dlat','dlon'))
    AI.description = 'Monthly Averaged Aerosol Index (using JZ method)'
    OB_COUNT = nc.createVariable('OB_COUNT','i2',('nmth','dlat','dlon'))
    OB_COUNT.description = '# of OMI AI measurements used in each monthly average'

    # Fill in dimension variables one-by-one.
    # NOTE: not sure if you can insert the entire data array into the
    # dimension (for example, doing MONTH = times), so this could be
    # something for you to try. Might make this faster if it works
    # ---------------------------------------------------------------
    for i in range(num_time):
        MONTH[i] = times[i]
    for i in range(num_lat):
        for j in range(num_lon):
            LAT[i,j]=lat_ranges[i] + 0.5
            LON[i,j]=lon_ranges[j] + 0.5

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
    outfile = home_dir + '/Research/OMI/omi_ai_da_'+da_time.strftime("%Y%m%d%H") + '.nc'
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

# Writes a single Shawn file to HDF5. 
def write_shawn_to_HDF5(OMI_base, save_path = './', minlat = 65., \
        shawn_path = home_dir + '/data/OMI/shawn_files/'):

    if(isinstance(OMI_base, str)):
        # Read the swath
        # --------------
        OMI_base  = readOMI_swath_shawn(OMI_base, latmin = minlat,\
            shawn_path = shawn_path)

    # Convert the filename object to datetime
    # ---------------------------------------
    file_date = OMI_base['date']
    dt_date_str = datetime.strptime(file_date, '%Y%m%d%H%M')
  
    # Create a new HDF5 dataset to write to the file
    # ------------------------------------------------
    outfile = save_path + 'omi_shawn_' + file_date + '.hdf5'
    dset = h5py.File(outfile,'w')
 
    dset.create_dataset('latitude',  data = OMI_base['LAT'][:,:])
    dset.create_dataset('longitude', data = OMI_base['LON'][:,:])
    dset.create_dataset('lat_crnr',  data = OMI_base['LATcrnr'][:,:,:])
    dset.create_dataset('lon_crnr',  data = OMI_base['LONcrnr'][:,:,:])
    dset.create_dataset('uvai_pert', data = OMI_base['UVAI_pert'][:,:])
    dset.create_dataset('uvai_raw',  data = OMI_base['UVAI_raw'][:,:])
    dset.create_dataset('time',      data = OMI_base['TIME_2'][:,:])
    dset.create_dataset('sza',       data = OMI_base['SZA'][:,:])
    dset.create_dataset('vza',       data = OMI_base['VZA'][:,:])
    dset.create_dataset('azm',       data = OMI_base['RAZ'][:,:])

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)  

# Writes a single Shawn file to netCDF. The file must be on the local
# computer under home_dir + /data/OMI/shawn_files/
def write_shawn_to_NCDF(filename,latmin, save_path = './', shawn_path = \
        home_dir + '/data/OMI/shawn_files/'):

    # Convert the filename object to datetime
    # ---------------------------------------
    file_date = filename.strip().split('/')[-1]
    dt_date_str = datetime.strptime(file_date, '%Y%m%d%H%M')

    # Read in the Shawn file
    # ----------------------
    data = readOMI_swath_shawn(filename, latmin = latmin, \
        shawn_path = shawn_path)

    # Don't use swaths that contain entirely masked data. This will
    # focus the netCDF file on only the chunk of data north of 65 N.
    # --------------------------------------------------------------
    data['UVAI_pert'] = np.ma.masked_invalid(data['UVAI_pert'])
    data['UVAI_raw']  = np.ma.masked_invalid(data['UVAI_raw'])
    data['SZA']  = np.ma.masked_invalid(data['SZA'])
    data['VZA']  = np.ma.masked_invalid(data['VZA'])
    data['RAZ']  = np.ma.masked_invalid(data['RAZ'])
    tmp = np.array([False in line for line in data['UVAI_pert'].mask])
    good_idx0 = np.where(tmp == True)

    final_pert = data['UVAI_pert'][good_idx0,:].squeeze()
    final_raw  = data['UVAI_raw'][good_idx0,:].squeeze()
    final_lat  = data['LAT'][good_idx0,:].squeeze() 
    final_lon  = data['LON'][good_idx0,:].squeeze() 
    final_sza  = data['SZA'][good_idx0,:].squeeze() 
    final_vza  = data['VZA'][good_idx0,:].squeeze() 
    final_raz  = data['RAZ'][good_idx0,:].squeeze() 

    # Create a new netCDF dataset to write to the file
    # ------------------------------------------------
    outfile = save_path + 'omi_uvai_perturbed_'+ file_date + '.nc'
    nc = Dataset(outfile,'w',format='NETCDF4')
  
    # Dimensions = lon, lat
    # Create the sizes of each dimension in the file. In this case,
    # the dimensions are "# of latitude" and "# of longitude"
    # -------------------------------------------------------------
    num_x = final_pert.shape[0]
    num_y = final_pert.shape[1]
    
    # Use the dimension size variables to actually create dimensions in 
    # the file.
    # ----------------------------------------------------------------- 
    n_x  = nc.createDimension('nx',num_x)
    n_y  = nc.createDimension('ny',num_y)

    # Create variables for the three dimensions. Note that since these
    # are variables, they are still given 'dimensions' using 'dlat', 
    # and 'dlon'. Latitude and longitude are each 2-d grids, so they are 
    # given 2 dimensions (dlat, dlon).
    # ------------------------------------------------------------------
    LON=nc.createVariable('Longitude','f4',('nx','ny'))
    LON.description='Longitude (-180 to 180)'
    LON.units='Degrees'
    LAT=nc.createVariable('Latitude','f4',('nx','ny'))
    LAT.description='Latitude (-90 to 90)'
    LAT.units='Degrees'

    # Create a variable for the AI data, dimensioned using both 
    # dimensions.
    # ---------------------------------------------------------  
    AI_raw = nc.createVariable('UVAI_raw','f4',('nx','ny'))
    AI_raw.description='Raw OMI OMAERUV UltraViolet Aerosol Index data ' + \
        'before the quality control methods. These ' + \
        'raw AI data were obtained from the online archive housed at the ' + \
        'GES DISC OMAERUV dataset page ' + \
        '(https://disc.gsfc.nasa.gov/datasets/OMAERUV_003/summary)'

    AI_pert = nc.createVariable('UVAI_perturbed','f4',('nx','ny'))
    AI_pert.description='Perturbed OMI OMAERUV UltraViolet Aerosol Index data. ' + \
        'Calculated by subtracting the ' + \
        'climatology value (determined by binning the raw AI data per month by' + \
        ' viewing geometry, surface albedo, and ground pixel classification) ' + \
        'from the raw AI value contained in the GES DISC HDF5 OMI files.'
    AI_pert.units = 'None'

    SZA = nc.createVariable('SolarZenithAngle','f4',('nx','ny'))
    SZA.description='Solar Zenith Angle'
    SZA.units = 'Degrees'

    VZA = nc.createVariable('ViewingZenithAngle','f4',('nx','ny'))
    VZA.description='Viewing Zenith Angle'
    VZA.units = 'degrees'

    RAZ = nc.createVariable('RelativeAzimuthAngle','f4',('nx','ny'))
    RAZ.description='Relative Azimuth Angle (East of North)'
    RAZ.units = 'Degrees'

    # Fill in dimension variables one-by-one.
    # NOTE: not sure if you can insert the entire data array into the
    # dimension (for example, doing :
    #      LAT, LON = np.meshgrid(lat_ranges,lon_ranges) ),
    # so this could be
    # something for you to try. Might make this faster if it works
    # ---------------------------------------------------------------
    LON[:,:]      = final_lon[:,:]
    LAT[:,:]      = final_lat[:,:]
    AI_raw[:,:]   = final_raw[:,:]
    AI_pert[:,:]  = final_pert[:,:]
    SZA[:,:]      = final_sza[:,:]
    VZA[:,:]      = final_vza[:,:]
    RAZ[:,:]      = final_raz[:,:]

    # Save, write, and close the netCDF file
    # --------------------------------------
    nc.close()

    print("Saved file ",outfile)  

# Wrapper to automate the shawn file conversion
# ---------------------------------------------
def auto_shawn_convert(start_year = 2005, end_year = 2020, latmin = 65., run = False):

    if(not run):
        print("NOT RUNNING")

    # This is the path that points to the HDF5 OMI files. This must be changed
    # if running on a new system.
    data_path = home_dir + '/data/OMI/'
    #data_path = '/data/OMI/'
    hdf_path = data_path + 'H5_files/'
    shawn_path = data_path + 'shawn_files/'
    ltc4_path = shawn_path + 'ltc4/'
    net_path = shawn_path + 'netCDF/'
    
    base_date = datetime(year=start_year,month=4,day=1)
    end_date  = datetime(year=end_year,month=10,day=1)
    
    get_data = False
    first_time = True
    
    # Keep looping until the end date is reached
    # ------------------------------------------
    while(base_date < end_date):
    
        work_time       = base_date.strftime("%Y%m")
        data_time       = base_date.strftime("%Ym")
        data_time_shawn = base_date.strftime("%Y")
    
        # Copy the working year's OMI data from Raindrop to the local
        # computer.
        # NOTE: If this code is being run on Raindrop, this is not
        # necessary.
        # -----------------------------------------------------------
        if(get_data == True):
            get_data = False
            last_time       = (base_date - relativedelta(years = 1)).strftime("%Ym")
            last_time_shawn = (base_date - relativedelta(years = 1)).strftime("%Y")
            print("getting data")
            for m_idx in range(base_date.month,10):

                # Remove this month's data from last year
                if(first_time == False):

                    # Remove the HDF5 files
                    # ---------------------
                    print('rm '+hdf_path+'OMI*OMAERUV_'+last_time+\
                        str(m_idx).zfill(2)+'*')
                    if(run):
                        os.system('rm '+ltc4_path+'OMI*OMAERUV_'+last_time+\
                            str(m_idx).zfill(2)+'*')

                    # Remove the shawn files
                    # ----------------------
                    print('rm '+ltc4_path + last_time_shawn + str(m_idx).zfill(2) + '*')
                    if(run):
                        os.system('rm '+ltc4_path + last_time_shawn + str(m_idx).zfill(2) + '*')

                # Copy over the HDF files
                # -----------------------
                print('scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/H5_files/'+\
                    'OMI*OMAERUV_'+data_time+str(m_idx).zfill(2)  + '* ' + hdf_path)
                if(run):
                    os.system('scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/H5_files/'+\
                        'OMI*OMAERUV_'+data_time+str(m_idx).zfill(2)  + '* ' + hdf_path)

                # Copy over the shawn files
                # -------------------------
                print('scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/out_files-ltc4/'+\
                    data_time_shawn + str(m_idx).zfill(2)  + '* ' + ltc4_path)
                if(run):
                    os.system('scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/out_files-ltc4/'+\
                        data_time_shawn + str(m_idx).zfill(2)  + '* ' + ltc4_path)
            first_time = False
    
        print(base_date.strftime("%Y%m%d"))
        
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        #
        # Select the needed data files
        #
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
        ##!#date = plot_time[4:8]
        ##!#if(len(plot_time)==13):
        ##!#    time = plot_time[9:]
        ##!#elif(len(plot_time)==12):
        ##!#    time = plot_time[8:]
        ##!#else:
        ##!#    time = ''
    
        try:
            total_list = subprocess.check_output('ls ' + ltc4_path + \
                work_time + '*', \
                shell=True).decode('utf-8').strip().split('\n')
        
            # Loop over each of these files and run the netCDF4 generator
            # -----------------------------------------------------------
            for sfile in total_list:
                print(sfile)
                if(run):
                    write_shawn_to_NCDF(sfile.strip().split('/')[-1],latmin, save_path = net_path,\
                        shawn_path = ltc4_path)

            ##!## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
            ##!##
            ##!## Grid the data
            ##!##
            ##!## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
            ##!#whole_rows = np.zeros((len(total_list),60))
            ##!#total_AI   = np.full((len(total_list),2000,60),np.nan)
            ##!#total_X   = np.full((len(total_list),2000,60),np.nan)
            ##!##total_AI   = np.full((len(total_list),2000,60),-99.)
            ##!##total_X   = np.full((len(total_list),2000,60),-99.)
            ##!#for fileI in range(len(total_list)):
            ##!#    # read in data directly from HDF5 files
            ##!#    #print(total_list[fileI])
            ##!#    data = h5py.File(total_list[fileI],'r')
    
            ##!#    AI1     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
            ##!#    XTRACK1 = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
            ##!#    XTRK_decode = np.array([[get_xtrack_flags(val) for val in XTRKrow] for XTRKrow in XTRACK1])
            ##!#    #print(AI1.shape)
        
            ##!#    # Mask values with missing AI or XTrack values that are not 0 or 4
            ##!#    # NEW: After decoding the XTRACK values, only see where the row anomaly
            ##!#    #      flag isequal to zero
            ##!#    AI1[(AI1 < -2e5) | ((XTRACK1 != 0) & (XTRACK1 != 4))] = np.nan 
    
            ##!#    total_AI[fileI,:AI1.shape[0],:] = AI1[:,:]
    
            ##!#    # NEW: After decoding the XTRACK values, only see where the row anomaly
            ##!#    #      flag is equal to zero. Ignoring the other bit values
            ##!#    total_X[fileI,:XTRACK1.shape[0],:] = XTRK_decode[:,:,0]
        
            ##!#    data.close()
    
            ##!## Calculate the averages along each row in the current swath over the 
            ##!## Arctic. NOTE: This section should go in the file loop, but am
            ##!## too lazy to do it now.
            ##!## --------------------------------------------------------------------
            ##!#total_avgs1 = np.array([[np.nanmean(total_AI[file_idx,1230:1500,idx]) \
            ##!#    for idx in range(60)] for file_idx in range(len(total_list))])
    
            ##!## Calculate the average xtrack value along each row in the current 
            ##!## swath over the Arctic. NOTE: This section should go in the file loop,
            ##!## but am too lazy to do it now. This is used for identifying known
            ##!## contaminted rows.
            ##!## --------------------------------------------------------------------
            ##!#total_avgs_X1 = np.array([[np.nanmean(total_X[file_idx,1230:1500,idx]) \
            ##!#    for idx in range(60)] for file_idx in range(len(total_list))])
        
            ##!#total_stds1 = np.array([[np.nanstd(total_AI[file_idx,1230:1500,idx]) \
            ##!#    for idx in range(60)] for file_idx in range(len(total_list))])
        
            ##!## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
            ##!##
            ##!## A row is identified as "bad" if the average of that row between the 
            ##!## 1230 and 1500 indices across each swath during the day is more than
            ##!## 2 standard deviations away from the average of all rows between 1230 and
            ##!## 1500 for the day..
            ##!##
            ##!## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        
            ##!#day_avgs = np.nanmean(total_avgs1,axis=0)
            ##!#avg_day_avg = np.nanmean(day_avgs)
            ##!#day_std  = np.nanstd(day_avgs)
    
            ##!## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
            ##!##
            ##!## A row is identified as "contaminated" if it has at least 1 pixel within
            ##!## indices 1230 and 1500 with an Xtrack QC value that is not zero, meaning
            ##!## that it is not perfectly clean as defined by the flag. This is determined
            ##!## by taking the average of all the Xtrack QC values along indices 1230 to
            ##!## 1500 of each row, and if that average is not equal to 0, it has a non-zero
            ##!## pixel inside it and is therefore contaminated.
            ##!##
            ##!## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
            ##!#day_avgs_X = np.nanmean(total_avgs_X1,axis=0)
    
            ##!## If all daily average Xtrack values are missing, then 
            ##!## don't try to count the xtrack rows
            ##!## ----------------------------------------------------
            ##!#count_xtrack = True
            ##!#if(np.count_nonzero(~np.isnan(day_avgs_X)) == 0):
            ##!#    count_xtrack = False
        
            ##!## Deal with xtrack rows
            ##!## If a row has any XTRACK flags set, add to xtrack row list
            ##!## ---------------------------------------------------------
            ##!#bad_rows = []
            ##!#xtrack_rows = []
            ##!#for rowI in range(len(day_avgs)):
            ##!#    if((day_avgs[rowI] - avg_day_avg) > (day_std * 3)):
            ##!#        bad_rows.append(rowI+1)
            ##!#    if(count_xtrack):
            ##!#        # If a row has any XTRACK flags set, add to xtrack row list
            ##!#        if(day_avgs_X[rowI] != 0):
            ##!#            xtrack_rows.append(rowI+1)
            ##!#    #        print("Bad row number ", rowI+1) 
            
            
           
        except subprocess.CalledProcessError:
            print("ERROR: no data found for time",work_time)

        # Print the date and bad row values to the output file 
        base_date = base_date + relativedelta(months=1)
        #base_date = base_date + timedelta(months=1)
        if(base_date.month == 10):
            get_data = True
            base_date = datetime(year = base_date.year+1,month=4,day=1)
    
    # Delete the ending data
    # NOTE: IF this is being run on Raindrop, this MUST be removed
    #       ensure that the OMI data on Raindrop are not deleted.
    # ------------------------------------------------------------
    last_time = (base_date - relativedelta(years = 1)).strftime("%Ym")
    last_time_shawn = (base_date - relativedelta(years = 1)).strftime("%Y")
    for m_idx in range(base_date.month,10):
        # Remove this month's data from last year
        if(first_time == False):
            ##!#print('rm ' + ltc4_path + 'OMI*OMAERUV_' + last_time + \
            ##!#      str(m_idx).zfill(2) + '*')
            ##!#if(run):
            ##!#    os.system('rm ' + ltc4_path + 'OMI*OMAERUV_' + last_time + \
            ##!#          str(m_idx).zfill(2) + '*')

            # Remove the HDF5 files
            # ---------------------
            print('rm '+hdf_path+'OMI*OMAERUV_'+last_time+\
                str(m_idx).zfill(2)+'*')
            if(run):
                os.system('rm '+ltc4_path+'OMI*OMAERUV_'+last_time+\
                    str(m_idx).zfill(2)+'*')

            # Remove the shawn files
            # ----------------------
            print('rm '+ltc4_path + last_time_shawn + str(m_idx).zfill(2) + '*')
            if(run):
                os.system('rm '+ltc4_path + last_time_shawn + str(m_idx).zfill(2) + '*')


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Calculation functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Calculate the amount of total sampled area reduced by the row anomaly
# ---------------------------------------------------------------------
def calcOMI_row_area_reduce(date_str):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # Read in OMI swath
    # -----------------
    base_path = home_dir + '/data/OMI/H5_files/'
    filename = glob(base_path + dt_date_str.strftime('OMI-*%Ym%m%dt%H%M*.he5'))
    ##!#total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+\
    ##!#    year+'m'+date+'t'+time+'*.he5',shell=True).decode('utf-8').strip().split('\n')
    data = h5py.File(filename[0], 'r')
    areas = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/FoV75Area'][:]

    data.close()
    OMI_base  = readOMI_swath_hdf(date_str, 'control', latmin = 55)

    # Find out exactly which rows are row anomaly-affected
    # ----------------------------------------------------
    flag_avgs = np.mean(OMI_base['XTRACK'][:,:,0], axis = 0)
    flagged   = flag_avgs != 0
    unflagged = flag_avgs == 0

    # Extract the pixel areas and the row anomaly flags
    # -------------------------------------------------
    rows = np.arange(1, len(flag_avgs) + 1)
    total_area = np.sum(areas)
    total_clear_area = np.sum(areas[unflagged])

    print("Total area", total_area)
    print("Total clear area", total_clear_area)
    print("% reduced", 100. - (total_clear_area / total_area) * 100.)

    # Make a plot showing this
    # ------------------------
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(rows, areas, label = 'Unaffected')
    #ax.plot(np.arange(len(areas))[unflagged], areas[unflagged], color = 'tab:blue')
    ax.plot(rows[flagged], areas[flagged], color = 'red', \
        label = 'Affected')
    ax.legend()
    ax.set_xlabel('Row number')
    ax.set_ylabel('FoV Area [km$^{2}$]')
    ax.set_title(date_str)
    #ax.fill_between(range(len(areas)), np.min(areas), np.max(areas), \
    #    where = flagged, alpha = 0.5)

    ax.grid()
    plt.show()

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
    print('HERE:',OMI_data['DATES'][month_idx::index_jumper])
    local_data   = np.copy(OMI_data['AI'][month_idx::index_jumper,:,:])
    local_counts = np.copy(OMI_data['OB_COUNT'][month_idx::index_jumper,:,:])
    local_mask = np.ma.masked_where(local_counts == 0, local_data)
    local_mask = np.ma.masked_where((local_mask == -999.9) & \
        (OMI_data['LAT'] < minlat), local_mask)

    ai_trends = np.full(local_data.shape[1:], np.nan)
    ai_pvals  = np.full(local_data.shape[1:], np.nan)
    ai_uncert = np.full(local_data.shape[1:], np.nan)
    # Add uncertainty here?

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
                #slope, intercept, r_value, p_value, std_err, intcpt_stderr = \
                result = stats.linregress(x_vals,local_mask[:,i,j])
                ai_trends[i,j] = result.slope * len(x_vals)
                ai_pvals[i,j]  = result.pvalue
                ai_uncert[i,j] = result.stderr * len(x_vals)
            else:
                res = stats.theilslopes(local_mask[:,i,j], x_vals, 0.90)
                ai_trends[i,j] = res[0]*len(x_vals)

    ai_trends = np.ma.masked_where(OMI_data['LAT'] < minlat, ai_trends)
    ai_pvals  = np.ma.masked_where(OMI_data['LAT'] < minlat, ai_pvals)
    ai_uncert = np.ma.masked_where(OMI_data['LAT'] < minlat, ai_uncert)

    print('in trend calc')
    for x, y in zip(OMI_data['LAT'][:,10], ai_trends[:,10]):
        print(x,y)

    return ai_trends, ai_pvals, ai_uncert

   
# This function assumes the data is being read from the netCDF file
# NOTE: Assume user is using new OMI climo file which starts in January
def calcOMI_MonthClimo(OMI_data):

    # Set up arrays to hold monthly climatologies
    month_climo = np.zeros((6,OMI_data['AI'].shape[1],OMI_data['AI'].shape[2]))

    # Mask the monthly averages
    local_data   = np.copy(OMI_data['AI'][:,:,:])
    local_counts = np.copy(OMI_data['OB_COUNT'][:,:,:])
    local_mask = np.ma.masked_where(local_counts == 0, local_data)
    local_mask = np.ma.masked_where(local_mask == -999.9, local_mask)
 
    # Calculate monthly climatologies
    for m_i in range(6):
        month_climo[m_i,:,:] = np.nanmean(local_mask[m_i::6,:,:],axis=0)
        print("Month: ",OMI_data['DATES'][m_i][4:]) 
        ##!#month_climo[m_i,:,:] = np.nanmean(local_mask[m_i::12,:,:],axis=0)
        ##!#print("Month: ",OMI_data['DATES'][m_i][4:], OMI_data['DATES'][::12]) 
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

    base_path = home_dir + '/data/OMI/H5_files/'
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
        vmin = None, vmax = None, colorbar_label_size = 14, colorbar = True, \
        minlat = 65., pvals = None):

    if(vmin == None):
        vmin = np.nanmin(pdata)
    if(vmax == None):
        vmax = np.nanmax(pdata)

    if(ptype == 'trend'):
        #colormap = plt.cm.bwr
        colormap = plt.cm.get_cmap('bwr', 5)
    elif(ptype == 'uncert'):
        #colormap = plt.cm.plasma
        colormap = plt.cm.get_cmap('jet', 6)
    else:
        colormap = plt.cm.jet

    # Make copy of OMI_data array
    local_data  = np.copy(pdata)
    # Grid the data, fill in white space
    cyclic_data,cyclic_lons = add_cyclic_point(local_data,plon[0,:])
    plat2,plon2 = np.meshgrid(plat[:,0],cyclic_lons)   
  
    # Mask any missing values
    mask_AI = np.ma.masked_where(cyclic_data < -998.9, cyclic_data)
    mask_AI = np.ma.masked_where(plat2.T < minlat, mask_AI)

    # Plot lat/lon lines
    gl = pax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, \
        linewidth = 1, color = 'gray', alpha = 0.5, linestyle = '-',\
        y_inline = True, xlocs = range(-180, 180, 30), ylocs = range(70, 90, 10))
    #gl.top_labels = False
    #gl.bottom_labels = False
    #gl.left_labels = False
    #gl.right_labels = False
    #gl.xlines = False
    #gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 30))
    #gl.ylocator = mticker.FixedLocator(np.arange(70, 90, 10))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    #gl.xlabel_style = {'size': 5, 'color': 'gray'}
    gl.xlabel_style = {'color': 'gray', 'weight': 'bold'}
    gl.ylabel_style = {'size': 15, 'color': 'gray'}
    gl.ylabel_style = {'color': 'black', 'weight': 'bold'}
    #pax.gridlines()

    pax.coastlines(resolution='50m')
    mesh = pax.pcolormesh(plon2, plat2,\
            mask_AI.T,transform = datacrs,\
            cmap = colormap,vmin=vmin,vmax=vmax, shading = 'auto')
    if(pvals is not None):
        cyclic_pvals, cyclic_lons = add_cyclic_point(pvals, plon[0,:])
        print(pvals.shape, plat.shape, plat2.shape)
        mask_pvals = np.ma.masked_where((plat < minlat) | \
            (pvals > 0.05), pvals)
        pax.pcolor(plon, plat, mask_pvals, hatch = '...', alpha = 0.0, \
            shading = 'auto', transform = datacrs)

    pax.set_extent([-180,180,minlat,90],datacrs)
    pax.set_boundary(circle, transform=pax.transAxes)
    #ax.set_xlim(-3430748.535086173,3430748.438879491)
    #ax.set_ylim(-3413488.8763307533,3443353.899053069)
    #cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.2),\
    if(colorbar):
        cbar = plt.colorbar(mesh,\
            ax = pax, orientation='vertical',shrink = 0.8, extend = 'both')
        #cbar10.set_label('UV Aerosol Index',weight='bold',fontsize=colorbar_label_size)
        #cbar.ax.tick_params(labelsize=14)
        cbar.set_label(plabel,fontsize=colorbar_label_size,weight='bold')
    pax.set_title(ptitle)


# Designed to work with the netCDF data
def plotOMI_MonthTrend(OMI_data,month_idx=None,save=False,\
        trend_type='standard',minlat=65.,return_trend=False, colorbar = True, \
        title = '', label = '', colorbar_label_size = 14, pax = None, \
        show_pval = False, uncert_ax = None):
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
    ai_trends, ai_pvals, ai_uncert = calcOMI_grid_trend(OMI_data, month_idx, trend_type, \
    #ai_trends, ai_pvals = calcOMI_grid_trend(OMI_data, month_idx, trend_type, \
        minlat)
    if(not show_pval):
        ai_pvals = None
    else:
        print('month_idx = ',month_idx,' PVAL nanmean = ', np.nanmean(ai_pvals))

    if(uncert_ax is None):
        ai_uncert = None
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
            ptitle = title, plabel = label, colorbar = colorbar, \
            colorbar_label_size = colorbar_label_size, \
            vmin = v_min, vmax = v_max, minlat = minlat, pvals = ai_pvals)

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
            ptitle = title, plabel = label, colorbar = colorbar, \
            colorbar_label_size = colorbar_label_size, \
            vmin = v_min, vmax = v_max, minlat = minlat, pvals = ai_pvals)

    if(uncert_ax is not None):
        plotOMI_spatial(uncert_ax, OMI_data['LAT'], OMI_data['LON'], ai_uncert, 'uncert', \
            ptitle = title, plabel = label, colorbar = colorbar, \
            colorbar_label_size = colorbar_label_size, \
            vmin = 0, vmax = 0.3, minlat = minlat)

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

    local_data[local_data == 0.0] = np.nan

    # Mask any missing values
    mask_AI = np.ma.masked_where(local_count == 0, local_data)
    #mask_AI = np.ma.masked_where(mask_AI == -999.9, mask_AI)

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
        title = 'April and May', minlat=minlat,pax = ax0, save=False)

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
        title = 'June, July, August', minlat=minlat, pax = ax1, save=False)
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
        vmin = -1., vmax = 1., ptitle = None, plabel = None, 
        save=False):

    version = OMI_data['VERSION']

    #if(version == 'VSJ2'):
    #    vmax = 0.5
    #    vmin = -0.5
    #    colormap = plt.cm.bwr
    #else:
    #    vmax = 1.5
    #    vmin = -1.5
    #    colormap = plt.cm.jet

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
    if(ptitle is None):
        ptitle = 'OMI AI ' + date_month + ' Climatology ('+version+')\n'+\
            start_date.strftime('%b. %Y') + ' - ' + end_date.strftime('%b. %Y')
    if(plabel is None):
        plabel = 'UV Aerosol Index'

    if(pax is None):
        # Make figure
        plt.close('all')
        fig1 = plt.figure(figsize = (8,8))
        ax = plt.axes(projection = mapcrs)

        plotOMI_spatial(ax, OMI_data['LAT'], OMI_data['LON'], mask_AI, 'climo', \
            ptitle = ptitle, plabel = plabel, \
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
            ptitle = ptitle, plabel = plabel, \
            vmin = vmin, vmax = vmax, minlat = minlat)

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

    fig.text(0.10, 0.70, 'Climatology', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size)
    fig.text(0.10, 0.30, 'Trend', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size)

    #fig.tight_layout()

    outname = 'omi_TEST_ai_comps_'+start_date.strftime("%b")+'_'+\
        OMI_data1['VERSION']+'v'+OMI_data2['VERSION']+\
        'v'+OMI_data3['VERSION']+'.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Generate a 15-panel figure comparing the climatology and trend between 3
# versions of the OMI data for all 3 summer months
def plotOMI_Compare_ClimoTrend_summer(OMI_data1,OMI_data2,OMI_data3,\
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

    #fig = plt.figure()
    plt.close('all')
    fig = plt.figure(figsize=(24,14))
    #plt.suptitle('OMI Comparisons: '+start_date.strftime("%B"),y=0.95,\
    #    fontsize=18,fontweight=4,weight='bold')
    gs = gridspec.GridSpec(nrows=3, ncols=5, hspace = 0.001, wspace = 0.15)

    # - - - - - - - - - - - - - - - - - - - - -
    # Plot the climatologies along the top row
    # - - - - - - - - - - - - - - - - - - - - -
       
    # Plot DATA1 climos
    # -----------------
    ##!## Make copy of OMI_data array
    local_data1_Jun  = np.copy(OMI_data1['MONTH_CLIMO'][2,:,:])
    local_data1_Jul  = np.copy(OMI_data1['MONTH_CLIMO'][3,:,:])
    local_data1_Aug  = np.copy(OMI_data1['MONTH_CLIMO'][4,:,:])
    local_data2_Jun  = np.copy(OMI_data2['MONTH_CLIMO'][2,:,:])
    local_data2_Jul  = np.copy(OMI_data2['MONTH_CLIMO'][3,:,:])
    local_data2_Aug  = np.copy(OMI_data2['MONTH_CLIMO'][4,:,:])
    local_data3_Jun  = np.copy(OMI_data3['MONTH_CLIMO'][2,:,:])
    local_data3_Jul  = np.copy(OMI_data3['MONTH_CLIMO'][3,:,:])
    local_data3_Aug  = np.copy(OMI_data3['MONTH_CLIMO'][4,:,:])

    mask_AI1_Jun = np.ma.masked_where(local_data1_Jun == -999.9, local_data1_Jun)
    mask_AI1_Jul = np.ma.masked_where(local_data1_Jul == -999.9, local_data1_Jul)
    mask_AI1_Aug = np.ma.masked_where(local_data1_Aug == -999.9, local_data1_Aug)
    mask_AI2_Jun = np.ma.masked_where(local_data2_Jun == -999.9, local_data2_Jun)
    mask_AI2_Jul = np.ma.masked_where(local_data2_Jul == -999.9, local_data2_Jul)
    mask_AI2_Aug = np.ma.masked_where(local_data2_Aug == -999.9, local_data2_Aug)
    mask_AI3_Jun = np.ma.masked_where(local_data3_Jun == -999.9, local_data3_Jun)
    mask_AI3_Jul = np.ma.masked_where(local_data3_Jul == -999.9, local_data3_Jul)
    mask_AI3_Aug = np.ma.masked_where(local_data3_Aug == -999.9, local_data3_Aug)

    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,OMI_data1['LON'][0,:])
    ##!#plat,plon = np.meshgrid(OMI_data1['LAT'][:,0],cyclic_lons)   
  
    # Mask any missing values
    #mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)
    ax00 = plt.subplot(gs[0,0], projection=mapcrs)   # June climo original
    ax01 = plt.subplot(gs[0,1], projection=mapcrs)   # June climo screened
    ax02 = plt.subplot(gs[0,2], projection=mapcrs)   # June trend original
    ax03 = plt.subplot(gs[0,3], projection=mapcrs)   # June trend screened
    ax04 = plt.subplot(gs[0,4], projection=mapcrs)   # June trend perturbed
    ax10 = plt.subplot(gs[1,0], projection=mapcrs)   # July climo original
    ax11 = plt.subplot(gs[1,1], projection=mapcrs)   # July climo screened
    ax12 = plt.subplot(gs[1,2], projection=mapcrs)   # July trend original
    ax13 = plt.subplot(gs[1,3], projection=mapcrs)   # July trend screened
    ax14 = plt.subplot(gs[1,4], projection=mapcrs)   # July trend perturbed
    ax20 = plt.subplot(gs[2,0], projection=mapcrs)   # August climo original
    ax21 = plt.subplot(gs[2,1], projection=mapcrs)   # August climo screened
    ax22 = plt.subplot(gs[2,2], projection=mapcrs)   # August trend original
    ax23 = plt.subplot(gs[2,3], projection=mapcrs)   # August trend screened
    ax24 = plt.subplot(gs[2,4], projection=mapcrs)   # August trend perturbed

    # Plot the figures in the first row: June
    # ---------------------------------------
    plotOMI_spatial(ax00, OMI_data1['LAT'], OMI_data1['LON'], mask_AI1_Jun, \
        'climo', ptitle = ' ', \
        plabel = '', \
        vmin = -1.0, vmax = 1.0, minlat = minlat)
    plotOMI_spatial(ax01, OMI_data2['LAT'], OMI_data2['LON'], mask_AI2_Jun, \
        'climo', ptitle = ' ', \
        plabel = 'UV Aerosol Index', \
        vmin = -1.0, vmax = 1.0, minlat = minlat)
    plotOMI_MonthTrend(OMI_data1,month_idx=2,\
        trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax02)
    #ax02.set_title(OMI_data2['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data2,month_idx=2,\
        trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax03)
    #ax03.set_title(OMI_data1['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data3,month_idx=2,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = ' ', pax = ax04)
    #ax04.set_title(OMI_data3['VERSION']+ '\n\n')

    # Plot the figures in the second row: July
    # ----------------------------------------
    plotOMI_spatial(ax10, OMI_data1['LAT'], OMI_data1['LON'], mask_AI1_Jul, \
        'climo', ptitle = ' ', \
        plabel = '', \
        vmin = -1.0, vmax = 1.0, minlat = minlat)
    plotOMI_spatial(ax11, OMI_data2['LAT'], OMI_data2['LON'], mask_AI2_Jul, \
        'climo', ptitle = ' ', \
        plabel = 'UV Aerosol Index', \
        vmin = -1.0, vmax = 1.0, minlat = minlat)
    plotOMI_MonthTrend(OMI_data1,month_idx=3,\
        trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax12)
    #ax12.set_title(OMI_data1['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data2,month_idx=3,\
        trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax13)
    #ax13.set_title(OMI_data2['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data3,month_idx=3,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = ' ', pax = ax14)
    #ax14.set_title(OMI_data3['VERSION']+ '\n\n')

    # Plot the figures in the third row: August
    # -----------------------------------------
    plotOMI_spatial(ax20, OMI_data1['LAT'], OMI_data1['LON'], mask_AI1_Aug, \
        'climo', ptitle = ' ', \
        plabel = '', \
        vmin = -1.0, vmax = 1.0, minlat = minlat)
    plotOMI_spatial(ax21, OMI_data2['LAT'], OMI_data2['LON'], mask_AI2_Aug, \
        'climo', ptitle = ' ', \
        plabel = 'UV Aerosol Index', \
        vmin = -1.0, vmax = 1.0, minlat = minlat)
    plotOMI_MonthTrend(OMI_data1,month_idx=4,\
        trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax22)
    #ax22.set_title(OMI_data2['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data2,month_idx=4,\
        trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax23)
    #ax23.set_title(OMI_data1['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data3,month_idx=4,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = ' ', pax = ax24)
    #ax24.set_title(OMI_data3['VERSION']+ '\n\n')

    #plot_subplot_label(ax00, '(a)')
    #plot_subplot_label(ax01, '(b)')
    #plot_subplot_label(ax02, '(c)')
    #plot_subplot_label(ax10, '(d)')
    #plot_subplot_label(ax11, '(e)')
    #plot_subplot_label(ax12, '(f)')

    fig.text(0.10, 0.75, 'June', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size)
    fig.text(0.10, 0.5, 'July', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size)
    fig.text(0.10, 0.25, 'August', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size)

    fig.text(0.18, 0.88, 'Original\nClimatology', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.34, 0.88, 'Screened\nClimatology', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.497, 0.88, 'Original\nTrend', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.655, 0.88, 'Screened\nTrend', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.815, 0.88, 'Perturbed\nTrend', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)

    #fig.tight_layout()

    outname = 'omi_ai_comps_summer_'+\
        OMI_data1['VERSION']+'v'+OMI_data2['VERSION']+\
        'v'+OMI_data3['VERSION']+'.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Generate a 15-panel figure comparing the climatology and trend between 3
# versions of the OMI data for all months
def plotOMI_Compare_ClimoTrend_all(OMI_data1,OMI_data2,OMI_data3,\
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

    colorbar_label_size = 7
    axis_title_size = 8
    row_label_size = 10 
    #colorbar_label_size = 13
    #axis_title_size = 14.5
    #row_label_size = 14.5

    #fig = plt.figure()
    plt.close('all')
    fig = plt.figure(figsize=(11,13))
    #plt.suptitle('OMI Comparisons: '+start_date.strftime("%B"),y=0.95,\
    #    fontsize=18,fontweight=4,weight='bold')
    gs = gridspec.GridSpec(nrows=6, ncols=5, hspace = 0.001, wspace = 0.15)

    # - - - - - - - - - - - - - - - - - - - - -
    # Plot the climatologies along the top row
    # - - - - - - - - - - - - - - - - - - - - -
       
    # Plot DATA1 climos
    # -----------------
    ##!## Make copy of OMI_data array
    local_data1_Apr  = np.copy(OMI_data1['MONTH_CLIMO'][0,:,:])
    local_data1_May  = np.copy(OMI_data1['MONTH_CLIMO'][1,:,:])
    local_data1_Jun  = np.copy(OMI_data1['MONTH_CLIMO'][2,:,:])
    local_data1_Jul  = np.copy(OMI_data1['MONTH_CLIMO'][3,:,:])
    local_data1_Aug  = np.copy(OMI_data1['MONTH_CLIMO'][4,:,:])
    local_data1_Sep  = np.copy(OMI_data1['MONTH_CLIMO'][5,:,:])
    local_data2_Apr  = np.copy(OMI_data2['MONTH_CLIMO'][0,:,:])
    local_data2_May  = np.copy(OMI_data2['MONTH_CLIMO'][1,:,:])
    local_data2_Jun  = np.copy(OMI_data2['MONTH_CLIMO'][2,:,:])
    local_data2_Jul  = np.copy(OMI_data2['MONTH_CLIMO'][3,:,:])
    local_data2_Aug  = np.copy(OMI_data2['MONTH_CLIMO'][4,:,:])
    local_data2_Sep  = np.copy(OMI_data2['MONTH_CLIMO'][5,:,:])
    local_data3_Apr  = np.copy(OMI_data3['MONTH_CLIMO'][0,:,:])
    local_data3_May  = np.copy(OMI_data3['MONTH_CLIMO'][1,:,:])
    local_data3_Jun  = np.copy(OMI_data3['MONTH_CLIMO'][2,:,:])
    local_data3_Jul  = np.copy(OMI_data3['MONTH_CLIMO'][3,:,:])
    local_data3_Aug  = np.copy(OMI_data3['MONTH_CLIMO'][4,:,:])
    local_data3_Sep  = np.copy(OMI_data3['MONTH_CLIMO'][5,:,:])

    mask_AI1_Apr = np.ma.masked_where(local_data1_Apr == -999.9, local_data1_Apr)
    mask_AI1_May = np.ma.masked_where(local_data1_May == -999.9, local_data1_May)
    mask_AI1_Jun = np.ma.masked_where(local_data1_Jun == -999.9, local_data1_Jun)
    mask_AI1_Jul = np.ma.masked_where(local_data1_Jul == -999.9, local_data1_Jul)
    mask_AI1_Aug = np.ma.masked_where(local_data1_Aug == -999.9, local_data1_Aug)
    mask_AI1_Sep = np.ma.masked_where(local_data1_Sep == -999.9, local_data1_Sep)
    mask_AI2_Apr = np.ma.masked_where(local_data2_Apr == -999.9, local_data2_Apr)
    mask_AI2_May = np.ma.masked_where(local_data2_May == -999.9, local_data2_May)
    mask_AI2_Jun = np.ma.masked_where(local_data2_Jun == -999.9, local_data2_Jun)
    mask_AI2_Jul = np.ma.masked_where(local_data2_Jul == -999.9, local_data2_Jul)
    mask_AI2_Aug = np.ma.masked_where(local_data2_Aug == -999.9, local_data2_Aug)
    mask_AI2_Sep = np.ma.masked_where(local_data2_Sep == -999.9, local_data2_Sep)
    mask_AI3_Apr = np.ma.masked_where(local_data3_Apr == -999.9, local_data3_Apr)
    mask_AI3_May = np.ma.masked_where(local_data3_May == -999.9, local_data3_May)
    mask_AI3_Jun = np.ma.masked_where(local_data3_Jun == -999.9, local_data3_Jun)
    mask_AI3_Jul = np.ma.masked_where(local_data3_Jul == -999.9, local_data3_Jul)
    mask_AI3_Aug = np.ma.masked_where(local_data3_Aug == -999.9, local_data3_Aug)
    mask_AI3_Sep = np.ma.masked_where(local_data3_Sep == -999.9, local_data3_Sep)

    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,OMI_data1['LON'][0,:])
    ##!#plat,plon = np.meshgrid(OMI_data1['LAT'][:,0],cyclic_lons)   
  
    # Mask any missing values
    #mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)
    ax00 = plt.subplot(gs[0,0], projection=mapcrs)   # April climo original
    ax01 = plt.subplot(gs[0,1], projection=mapcrs)   # April climo screened
    ax02 = plt.subplot(gs[0,2], projection=mapcrs)   # April trend original
    ax03 = plt.subplot(gs[0,3], projection=mapcrs)   # April trend screened
    ax04 = plt.subplot(gs[0,4], projection=mapcrs)   # April trend perturbed
    ax10 = plt.subplot(gs[1,0], projection=mapcrs)   # May climo original
    ax11 = plt.subplot(gs[1,1], projection=mapcrs)   # May climo screened
    ax12 = plt.subplot(gs[1,2], projection=mapcrs)   # May trend original
    ax13 = plt.subplot(gs[1,3], projection=mapcrs)   # May trend screened
    ax14 = plt.subplot(gs[1,4], projection=mapcrs)   # May trend perturbed
    ax20 = plt.subplot(gs[2,0], projection=mapcrs)   # June climo original
    ax21 = plt.subplot(gs[2,1], projection=mapcrs)   # June climo screened
    ax22 = plt.subplot(gs[2,2], projection=mapcrs)   # June trend original
    ax23 = plt.subplot(gs[2,3], projection=mapcrs)   # June trend screened
    ax24 = plt.subplot(gs[2,4], projection=mapcrs)   # June trend perturbed
    ax30 = plt.subplot(gs[3,0], projection=mapcrs)   # July climo original
    ax31 = plt.subplot(gs[3,1], projection=mapcrs)   # July climo screened
    ax32 = plt.subplot(gs[3,2], projection=mapcrs)   # July trend original
    ax33 = plt.subplot(gs[3,3], projection=mapcrs)   # July trend screened
    ax34 = plt.subplot(gs[3,4], projection=mapcrs)   # July trend perturbed
    ax40 = plt.subplot(gs[4,0], projection=mapcrs)   # August climo original
    ax41 = plt.subplot(gs[4,1], projection=mapcrs)   # August climo screened
    ax42 = plt.subplot(gs[4,2], projection=mapcrs)   # August trend original
    ax43 = plt.subplot(gs[4,3], projection=mapcrs)   # August trend screened
    ax44 = plt.subplot(gs[4,4], projection=mapcrs)   # August trend perturbed
    ax50 = plt.subplot(gs[5,0], projection=mapcrs)   # September climo original
    ax51 = plt.subplot(gs[5,1], projection=mapcrs)   # September climo screened
    ax52 = plt.subplot(gs[5,2], projection=mapcrs)   # September trend original
    ax53 = plt.subplot(gs[5,3], projection=mapcrs)   # September trend screened
    ax54 = plt.subplot(gs[5,4], projection=mapcrs)   # September trend perturbed

    # Plot the figures in the first row: April
    # ---------------------------------------
    plotOMI_spatial(ax00, OMI_data1['LAT'], OMI_data1['LON'], mask_AI1_Apr, \
        'climo', ptitle = ' ', plabel = '', vmin = -1.0, vmax = 1.0, \
        minlat = minlat, colorbar = False)
    plotOMI_spatial(ax01, OMI_data2['LAT'], OMI_data2['LON'], mask_AI2_Apr, \
        'climo', ptitle = ' ', plabel = '', vmin = -1.0, vmax = 1.0, \
        minlat = minlat, colorbar_label_size = colorbar_label_size, colorbar = False)
    plotOMI_MonthTrend(OMI_data1,month_idx=0,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax02, colorbar = False, \
        colorbar_label_size = colorbar_label_size)
    #ax02.set_title(OMI_data2['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data2,month_idx=0,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax03, colorbar = False, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    #ax03.set_title(OMI_data1['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data3,month_idx=0,trend_type=trend_type,\
        label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = ' ', pax = ax04, colorbar = False, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    # Plot the figures in the first row: May
    # ---------------------------------------
    plotOMI_spatial(ax10, OMI_data1['LAT'], OMI_data1['LON'], mask_AI1_May, \
        'climo', ptitle = ' ', \
        plabel = '', colorbar = False, \
        vmin = -1.0, vmax = 1.0, minlat = minlat, \
        colorbar_label_size = colorbar_label_size)
    plotOMI_spatial(ax11, OMI_data2['LAT'], OMI_data2['LON'], mask_AI2_May, \
        'climo', ptitle = ' ', plabel = 'UV Aerosol Index',vmin = -1.0, \
        vmax = 1.0, minlat = minlat, colorbar = False, \
        colorbar_label_size = colorbar_label_size)
    plotOMI_MonthTrend(OMI_data1,month_idx=1,\
        trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax12, colorbar = False, \
        colorbar_label_size = colorbar_label_size)
    #ax02.set_title(OMI_data2['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data2,month_idx=1,\
        trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax13, colorbar = False, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    #ax03.set_title(OMI_data1['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data3,month_idx=1,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = ' ', pax = ax14, colorbar = False, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    # Plot the figures in the first row: June
    # ---------------------------------------
    plotOMI_spatial(ax20, OMI_data1['LAT'], OMI_data1['LON'], mask_AI1_Jun, \
        'climo', ptitle = ' ', \
        plabel = '', colorbar = False, \
        vmin = -1.0, vmax = 1.0, minlat = minlat, \
        colorbar_label_size = colorbar_label_size)
    plotOMI_spatial(ax21, OMI_data2['LAT'], OMI_data2['LON'], mask_AI2_Jun, \
        'climo', ptitle = ' ', \
        plabel = 'UV Aerosol Index', colorbar = False, \
        vmin = -1.0, vmax = 1.0, minlat = minlat, \
        colorbar_label_size = colorbar_label_size)
    plotOMI_MonthTrend(OMI_data1,month_idx=2,\
        trend_type=trend_type,label = ' ', colorbar = False, \
        minlat=65.,title = ' ', pax = ax22, \
        colorbar_label_size = colorbar_label_size)
    #ax02.set_title(OMI_data2['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data2,month_idx=2,\
        trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax23, colorbar = False, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    #ax03.set_title(OMI_data1['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data3,month_idx=2,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = ' ', pax = ax24, colorbar = False, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    #ax04.set_title(OMI_data3['VERSION']+ '\n\n')

    # Plot the figures in the second row: July
    # ----------------------------------------
    plotOMI_spatial(ax30, OMI_data1['LAT'], OMI_data1['LON'], mask_AI1_Jul, \
        'climo', ptitle = ' ', \
        plabel = '', \
        vmin = -1.0, vmax = 1.0, minlat = minlat, colorbar = False, \
        colorbar_label_size = colorbar_label_size)
    plotOMI_spatial(ax31, OMI_data2['LAT'], OMI_data2['LON'], mask_AI2_Jul, \
        'climo', ptitle = ' ', colorbar = False, \
        plabel = 'UV Aerosol Index', \
        vmin = -1.0, vmax = 1.0, minlat = minlat, \
        colorbar_label_size = colorbar_label_size)
    plotOMI_MonthTrend(OMI_data1,month_idx=3,\
        trend_type=trend_type,label = ' ',colorbar = False, \
        minlat=65.,title = ' ', pax = ax32, \
        colorbar_label_size = colorbar_label_size)
    #ax12.set_title(OMI_data1['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data2,month_idx=3,\
        trend_type=trend_type,label = ' ',colorbar = False, \
        minlat=65.,title = ' ', pax = ax33, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    #ax13.set_title(OMI_data2['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data3,month_idx=3,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = ' ', pax = ax34, colorbar = False, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    #ax14.set_title(OMI_data3['VERSION']+ '\n\n')

    # Plot the figures in the third row: August
    # -----------------------------------------
    plotOMI_spatial(ax40, OMI_data1['LAT'], OMI_data1['LON'], mask_AI1_Aug, \
        'climo', ptitle = ' ', \
        plabel = '', colorbar = False, \
        vmin = -1.0, vmax = 1.0, minlat = minlat, \
        colorbar_label_size = colorbar_label_size)
    plotOMI_spatial(ax41, OMI_data2['LAT'], OMI_data2['LON'], mask_AI2_Aug, \
        'climo', ptitle = ' ', colorbar = False, \
        plabel = 'UV Aerosol Index', \
        vmin = -1.0, vmax = 1.0, minlat = minlat, \
        colorbar_label_size = colorbar_label_size)
    plotOMI_MonthTrend(OMI_data1,month_idx=4,\
        trend_type=trend_type,label = ' ',colorbar = False, \
        minlat=65.,title = ' ', pax = ax42, \
        colorbar_label_size = colorbar_label_size)
    #ax22.set_title(OMI_data2['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data2,month_idx=4,\
        trend_type=trend_type,label = ' ',colorbar = False, \
        minlat=65.,title = ' ', pax = ax43, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    #ax23.set_title(OMI_data1['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data3,month_idx=4,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = ' ', pax = ax44, colorbar = False, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    #ax24.set_title(OMI_data3['VERSION']+ '\n\n')

    # Plot the figures in the third row: August
    # -----------------------------------------
    plotOMI_spatial(ax50, OMI_data1['LAT'], OMI_data1['LON'], mask_AI1_Sep, \
        'climo', ptitle = ' ', \
        plabel = '', colorbar = False, \
        vmin = -1.0, vmax = 1.0, minlat = minlat, \
        colorbar_label_size = colorbar_label_size)
    plotOMI_spatial(ax51, OMI_data2['LAT'], OMI_data2['LON'], mask_AI2_Sep, \
        'climo', ptitle = ' ', colorbar = False, \
        plabel = 'UV Aerosol Index', \
        vmin = -1.0, vmax = 1.0, minlat = minlat, \
        colorbar_label_size = colorbar_label_size)
    plotOMI_MonthTrend(OMI_data1,month_idx=5,\
        trend_type=trend_type,label = ' ',colorbar = False, \
        minlat=65.,title = ' ', pax = ax52, \
        colorbar_label_size = colorbar_label_size)
    #ax22.set_title(OMI_data2['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data2,month_idx=5,\
        trend_type=trend_type,label = ' ',colorbar = False, \
        minlat=65.,title = ' ', pax = ax53, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    #ax23.set_title(OMI_data1['VERSION']+ '\n\n')
    plotOMI_MonthTrend(OMI_data3,month_idx=5,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        minlat=65.,title = ' ', pax = ax54, colorbar = False, \
        colorbar_label_size = colorbar_label_size, show_pval = True)
    #ax24.set_title(OMI_data3['VERSION']+ '\n\n')

    #plot_subplot_label(ax00, '(a)')
    #plot_subplot_label(ax01, '(b)')
    #plot_subplot_label(ax02, '(c)')
    #plot_subplot_label(ax10, '(d)')
    #plot_subplot_label(ax11, '(e)')
    #plot_subplot_label(ax12, '(f)')

    fig.text(0.10, 0.82, 'April', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.692, 'May', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.565, 'June', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.435, 'July', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.305, 'August', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.18, 'September', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)

    fig.text(0.195, 0.89, 'Original\nClimatology', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.355, 0.89, 'Screened\nClimatology', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.51, 0.89, 'Original\nTrend', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.67, 0.89, 'Screened\nTrend', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.83, 0.89, 'Perturbed\nTrend', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)

    cax = fig.add_axes([0.13, 0.09, 0.29, 0.01])
    norm = mpl.colors.Normalize(vmin = -2.0, vmax = 3.0)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap = plt.cm.jet, norm = norm, \
        orientation = 'horizontal', extend = 'both')
    cb1.set_label('UV Aerosol Index', weight = 'bold')

    cax2 = fig.add_axes([0.448, 0.09, 0.45, 0.01])
    norm2 = mpl.colors.Normalize(vmin = -0.5, vmax = 0.5)
    cb2 = mpl.colorbar.ColorbarBase(cax2, cmap = plt.cm.bwr, norm = norm2, \
        orientation = 'horizontal', extend = 'both')
    cb2.set_label('AI Trend (2005 - 2020)', weight = 'bold')
    #fig.colorbar(norm, cax = cax, orientation = 'horizontal')

    #fig.tight_layout()

    outname = 'omi_ai_comps_all6_'+\
        OMI_data1['VERSION']+'v'+OMI_data2['VERSION']+\
        'v'+OMI_data3['VERSION']+'_' + trend_type + '_2.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Compare trends and single monthly averages, and time series of
# monthly averages for two given lat/lon points.
def plotOMI_TrendUncertainty_SingleMonth(dtype1, dtype2, month_idx,\
         minlat = 65., save = False):

    #minlat = 65.
    OMI_data1 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
        'omi_ai_' + dtype1 + '_2005_2020.nc', minlat = minlat)
    OMI_data2   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
        'omi_ai_' + dtype2 + '_2005_2020.nc', minlat = minlat)
    
    # Determine the month being plotted 
    # ---------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (9.5, 6))
    ax1 = fig1.add_subplot(2,3,1, projection = mapcrs) # JZ211 trend plot
    ax2 = fig1.add_subplot(2,3,4, projection = mapcrs) # VSJ4 trend plot
    ax5 = fig1.add_subplot(2,3,2, projection = mapcrs) # JZ211 uncert plot
    ax6 = fig1.add_subplot(2,3,5, projection = mapcrs) # VSJ4 uncert plot
    ax3 = fig1.add_subplot(2,3,3)                      # JZ211 scatter
    ax4 = fig1.add_subplot(2,3,6)                      # VSJ4 scatter

    # Plot the spatial trends
    plotOMI_MonthTrend(OMI_data1,month_idx=month_idx,trend_type='standard',label = ' ',\
        minlat=65.,title = dtype1, pax = ax1, colorbar = True, \
        colorbar_label_size = 10, show_pval = True, uncert_ax = ax5)
    plotOMI_MonthTrend(OMI_data2,month_idx=month_idx,trend_type='standard',label = ' ',\
        minlat=65.,title = dtype2, pax = ax2, colorbar = True, \
        colorbar_label_size = 10, show_pval = True, uncert_ax = ax6)

    ai_trends, ai_pvals, ai_uncert = \
        calcOMI_grid_trend(OMI_data1, month_idx, 'standard', minlat)
    mask_trends = np.ma.masked_invalid(ai_trends).compressed()
    mask_uncert = np.ma.masked_invalid(ai_uncert).compressed()
    xy = np.vstack([mask_trends, mask_uncert])
    z = stats.gaussian_kde(xy)(xy)       
    ax3.scatter(mask_trends, mask_uncert, c = z, s = 6) 
    ax3.set_xlabel('AI Trend')
    ax3.set_ylabel('AI Trend Uncertainty')

    # Plot one to one line
    xlim = ax3.get_xlim()
    ylim = ax3.get_ylim()
    xvals = np.arange(xlim[0] - 1, xlim[-1] + 1)
    ax3.plot(xvals, xvals, '-', color = 'tab:cyan', linewidth = 1.5)
    ax3.set_xlim(xlim)   
    ax3.set_ylim([0, xlim[1]])   
    ax3.grid()

    ai_trends, ai_pvals, ai_uncert = \
        calcOMI_grid_trend(OMI_data2, month_idx, 'standard', minlat)
    mask_trends = np.ma.masked_invalid(ai_trends).compressed()
    mask_uncert = np.ma.masked_invalid(ai_uncert).compressed()
    xy = np.vstack([mask_trends, mask_uncert])
    z = stats.gaussian_kde(xy)(xy)       
    ax4.scatter(mask_trends, mask_uncert, c = z, s = 6) 
    ax4.set_xlabel('AI Trend')
    ax4.set_ylabel('AI Trend Uncertainty')

    # Plot one to one line
    xlim = ax4.get_xlim()
    ylim = ax4.get_ylim()
    xvals = np.arange(xlim[0] - 1, xlim[-1] + 1, 0.01)
    ax4.plot(xvals, xvals, '-', color = 'tab:cyan', linewidth = 1.5)
    ax4.set_xlim(xlim)   
    ax4.set_ylim([0, xlim[1]])   
    ax4.grid()

    month_dict = {
        0: 'April',
        1: 'May',
        2: 'June',
        3: 'July',
        4: 'August',
        5: 'September',

    }

    plt.suptitle(month_dict[month_idx])
    
    fig1.tight_layout()

    outname = 'omi_ai_trend_uncert_month'+str(month_idx) + '_' + \
        OMI_data1['VERSION']+'v'+OMI_data2['VERSION']+'.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Compare the trends and trend uncertanties
def plotOMI_Compare_TrendUncertainty_all(dtype1, dtype2,\
        trend_type = 'standard', minlat=65,save=False):

    colormap = plt.cm.jet

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    OMI_data1 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
        'omi_ai_' + dtype1 + '_2005_2020.nc', minlat = minlat)
    OMI_data2   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
        'omi_ai_' + dtype2 + '_2005_2020.nc', minlat = minlat)

    colorbar_label_size = 7
    axis_title_size = 8
    row_label_size = 10 
    #colorbar_label_size = 13
    #axis_title_size = 14.5
    #row_label_size = 14.5

    #fig = plt.figure()
    plt.close('all')
    fig = plt.figure(figsize=(11,13))
    #plt.suptitle('OMI Comparisons: '+start_date.strftime("%B"),y=0.95,\
    #    fontsize=18,fontweight=4,weight='bold')
    gs = gridspec.GridSpec(nrows=6, ncols=4, hspace = 0.001, wspace = 0.15)

    # - - - - - - - - - - - - - - - - - - - - -
    # Plot the climatologies along the top row
    # - - - - - - - - - - - - - - - - - - - - -
    # Mask any missing values
    #mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)
    ax00 = plt.subplot(gs[0,0], projection=mapcrs)   # April trend type1
    ax01 = plt.subplot(gs[0,1], projection=mapcrs)   # April uncert type1
    ax02 = plt.subplot(gs[0,2], projection=mapcrs)   # April trend type2
    ax03 = plt.subplot(gs[0,3], projection=mapcrs)   # April uncert type2
    ax10 = plt.subplot(gs[1,0], projection=mapcrs)   # May trend type1
    ax11 = plt.subplot(gs[1,1], projection=mapcrs)   # May uncert type1
    ax12 = plt.subplot(gs[1,2], projection=mapcrs)   # May trend type2
    ax13 = plt.subplot(gs[1,3], projection=mapcrs)   # May uncert type2
    ax20 = plt.subplot(gs[2,0], projection=mapcrs)   # June trend type1
    ax21 = plt.subplot(gs[2,1], projection=mapcrs)   # June uncert type1
    ax22 = plt.subplot(gs[2,2], projection=mapcrs)   # June trend type2
    ax23 = plt.subplot(gs[2,3], projection=mapcrs)   # June uncert type2
    ax30 = plt.subplot(gs[3,0], projection=mapcrs)   # July trend type1
    ax31 = plt.subplot(gs[3,1], projection=mapcrs)   # July uncert type1
    ax32 = plt.subplot(gs[3,2], projection=mapcrs)   # July trend type2
    ax33 = plt.subplot(gs[3,3], projection=mapcrs)   # July uncert type2
    ax40 = plt.subplot(gs[4,0], projection=mapcrs)   # August trend type1
    ax41 = plt.subplot(gs[4,1], projection=mapcrs)   # August uncert type1
    ax42 = plt.subplot(gs[4,2], projection=mapcrs)   # August trend type2
    ax43 = plt.subplot(gs[4,3], projection=mapcrs)   # August uncert type2
    ax50 = plt.subplot(gs[5,0], projection=mapcrs)   # September trend type1
    ax51 = plt.subplot(gs[5,1], projection=mapcrs)   # September uncert type1
    ax52 = plt.subplot(gs[5,2], projection=mapcrs)   # September trend type2
    ax53 = plt.subplot(gs[5,3], projection=mapcrs)   # September uncert type2

    # Plot the figures in the first row: April
    # ---------------------------------------
    plotOMI_MonthTrend(OMI_data1,month_idx=0,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax00, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax01)
    plotOMI_MonthTrend(OMI_data2,month_idx=0,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax02, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax03)

    # Plot the figures in the first row: May
    # ---------------------------------------
    plotOMI_MonthTrend(OMI_data1,month_idx=1,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax10, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax11)
    plotOMI_MonthTrend(OMI_data2,month_idx=1,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax12, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax13)

    # Plot the figures in the first row: June
    # ---------------------------------------
    plotOMI_MonthTrend(OMI_data1,month_idx=2,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax20, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax21)
    plotOMI_MonthTrend(OMI_data2,month_idx=2,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax22, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax23)

    # Plot the figures in the second row: July
    # ----------------------------------------
    plotOMI_MonthTrend(OMI_data1,month_idx=3,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax30, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax31)
    plotOMI_MonthTrend(OMI_data2,month_idx=3,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax32, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax33)

    # Plot the figures in the third row: August
    # -----------------------------------------
    plotOMI_MonthTrend(OMI_data1,month_idx=4,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax40, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax41)
    plotOMI_MonthTrend(OMI_data2,month_idx=4,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax42, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax43)

    # Plot the figures in the third row: September
    # -----------------------------------------
    plotOMI_MonthTrend(OMI_data1,month_idx=5,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax50, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax51)
    plotOMI_MonthTrend(OMI_data2,month_idx=5,trend_type=trend_type,label = ' ',\
        minlat=65.,title = ' ', pax = ax52, colorbar = False, \
        colorbar_label_size = colorbar_label_size, uncert_ax = ax53)

    fig.text(0.10, 0.82, 'April', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.692, 'May', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.565, 'June', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.435, 'July', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.305, 'August', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.18, 'September', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)

    fig.text(0.200, 0.90, dtype1 + '\nTrend', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.410, 0.90, dtype1 + '\nUncertainty', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.62, 0.90, dtype2 + 'Trend', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.83, 0.90, dtype2 + '\nUncertainty', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)

    cax = fig.add_axes([0.13, 0.09, 0.38, 0.01])
    norm = mpl.colors.Normalize(vmin = -0.5, vmax = 0.5)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap = plt.cm.bwr, norm = norm, \
        orientation = 'horizontal', extend = 'both')
    cb1.set_label('AI Trend', weight = 'bold')

    cax2 = fig.add_axes([0.522, 0.09, 0.38, 0.01])
    norm2 = mpl.colors.Normalize(vmin = 0.0, vmax = 0.3)
    cmap = plt.cm.get_cmap('jet', 6)
    cb2 = mpl.colorbar.ColorbarBase(cax2, cmap = cmap, norm = norm2, \
        orientation = 'horizontal', extend = 'both')
    cb2.set_label('AI Trend Uncertainty', weight = 'bold')
    #fig.colorbar(norm, cax = cax, orientation = 'horizontal')

    #fig.tight_layout()

    outname = 'omi_ai_trend_uncert_all6_'+\
        OMI_data1['VERSION']+'v'+OMI_data2['VERSION']+'.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Generate a 15-panel figure comparing the climatology and trend between 3
# versions of the OMI data for all months
def plotOMI_Type_MonthClimo_all(minlat=65, save = False, \
        save_dir = home_dir + '/Research/OMI/monthly_images/combined_climo/'):

    colormap = plt.cm.jet

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    colorbar_label_size = 7
    axis_title_size = 8
    row_label_size = 10 

    data_types = ['VBS1','VJZ2','VJZ28','VJZ29','VJZ211','VSJ2','VSJ4']    

    # Read in all the data
    OMI_VBS1 = readOMI_NCDF(infile = home_dir + '/Research/OMI/' + \
        'omi_ai_' + data_types[0] + '_2005_2019.nc', minlat = minlat)
    OMI_VJZ2 = readOMI_NCDF(infile = home_dir + '/Research/OMI/' + \
        'omi_ai_' + data_types[1] + '_2005_2019.nc', minlat = minlat)
    OMI_VJZ28= readOMI_NCDF(infile = home_dir + '/Research/OMI/' + \
        'omi_ai_' + data_types[2] + '_2005_2019.nc', minlat = minlat)
    OMI_VJZ29= readOMI_NCDF(infile = home_dir + '/Research/OMI/' + \
        'omi_ai_' + data_types[3] + '_2005_2019.nc', minlat = minlat)
    OMI_VJZ211 = readOMI_NCDF(infile = home_dir + '/Research/OMI/' + \
        'omi_ai_' + data_types[4] + '_2005_2019.nc', minlat = minlat)
    OMI_VSJ2 = readOMI_NCDF(infile = home_dir + '/Research/OMI/' + \
        'omi_ai_' + data_types[5] + '_2005_2019.nc', minlat = minlat)
    OMI_VSJ4 = readOMI_NCDF(infile = home_dir + '/Research/OMI/' + \
        'omi_ai_' + data_types[6] + '_2005_2020.nc', minlat = minlat)
   
    OMI_VBS1['AI']   = np.ma.masked_where((OMI_VBS1['OB_COUNT'] == 0)  | (OMI_VBS1['OB_COUNT'] == -99)  , OMI_VBS1['AI'])
    OMI_VJZ2['AI']   = np.ma.masked_where((OMI_VJZ2['OB_COUNT'] == 0)  | (OMI_VJZ2['OB_COUNT'] == -99)  , OMI_VJZ2['AI'])
    OMI_VJZ28['AI']  = np.ma.masked_where((OMI_VJZ28['OB_COUNT'] == 0) | (OMI_VJZ28['OB_COUNT'] == -99) , OMI_VJZ28['AI'])
    OMI_VJZ29['AI']  = np.ma.masked_where((OMI_VJZ29['OB_COUNT'] == 0) | (OMI_VJZ29['OB_COUNT'] == -99) , OMI_VJZ29['AI'])
    OMI_VJZ211['AI'] = np.ma.masked_where((OMI_VJZ211['OB_COUNT'] == 0)| (OMI_VJZ211['OB_COUNT'] == -99), OMI_VJZ211['AI'])
    OMI_VSJ2['AI']   = np.ma.masked_where((OMI_VSJ2['OB_COUNT'] == 0)  | (OMI_VSJ2['OB_COUNT'] == -99)  , OMI_VSJ2['AI'])
    OMI_VSJ4['AI']   = np.ma.masked_where((OMI_VSJ4['OB_COUNT'] == 0)  | (OMI_VSJ4['OB_COUNT'] == -99)  , OMI_VSJ4['AI'])
 
    #for ii in range(len(data_types)):
    #print(data_types[ii])
    #fig = plt.figure()

    for jj in range(len(OMI_VJZ2['DATES'])):

        plt.close('all')
        fig = plt.figure(figsize=(12,6))
        ax1 = fig.add_subplot(2,4,1, projection = mapcrs) # VBS1
        ax2 = fig.add_subplot(2,4,2, projection = mapcrs) # VJZ2
        ax3 = fig.add_subplot(2,4,3, projection = mapcrs) # VJZ28
        ax4 = fig.add_subplot(2,4,4, projection = mapcrs) # VJZ29
        ax5 = fig.add_subplot(2,4,5, projection = mapcrs) # VJZ211
        ax6 = fig.add_subplot(2,4,6, projection = mapcrs) # VSJ2
        ax7 = fig.add_subplot(2,4,7, projection = mapcrs) # VSJ4

        # Read the data
        
        # Plot trend
        plotOMI_NCDF_SingleMonth(OMI_VJZ2,jj,minlat = minlat, pax = ax1, \
            save=False)
        plotOMI_NCDF_SingleMonth(OMI_VJZ28,jj,minlat = minlat, pax = ax2, \
            save=False)
        plotOMI_NCDF_SingleMonth(OMI_VJZ29,jj,minlat = minlat, pax = ax3, \
            save=False)
        plotOMI_NCDF_SingleMonth(OMI_VJZ211,jj,minlat = minlat, pax = ax4, \
            save=False)
        plotOMI_NCDF_SingleMonth(OMI_VBS1,jj,minlat = minlat, pax = ax5, \
            save=False)
        plotOMI_NCDF_SingleMonth(OMI_VSJ2,jj,minlat = minlat, pax = ax6, \
            save=False)
        plotOMI_NCDF_SingleMonth(OMI_VSJ4,jj,minlat = minlat, pax = ax7, \
            save=False)
        
        fig.tight_layout()

        if(save == True):
            outname = save_dir + 'omi_ai_monthly_comps_' + OMI_VJZ2['DATES'][jj] + '.png'
            fig.savefig(outname,dpi=300)
            print("Saved image",outname)
        else:
            plt.show()

# Generate a 15-panel figure comparing the climatology and trend between 3
# versions of the OMI data for all months
def plotOMI_Type_Climo_all(trend_type = 'standard', \
        minlat=65,save=False):

    colormap = plt.cm.jet

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    colorbar_label_size = 7
    axis_title_size = 8
    row_label_size = 10 

    #fig = plt.figure()
    plt.close('all')
    fig = plt.figure(figsize=(11,13))
    gs = fig.add_gridspec(nrows=6, ncols=4, hspace = 0.001, wspace = 0.15)

    data_types = ['VJZ2','VJZ28','VJZ29','VJZ211']    
 
    for ii in range(len(data_types)):
        print(data_types[ii])

        axi1 = plt.subplot(gs[0,ii], projection = mapcrs) # April
        axi2 = plt.subplot(gs[1,ii], projection = mapcrs) # May
        axi3 = plt.subplot(gs[2,ii], projection = mapcrs) # June
        axi4 = plt.subplot(gs[3,ii], projection = mapcrs) # July
        axi5 = plt.subplot(gs[4,ii], projection = mapcrs) # August
        axi6 = plt.subplot(gs[5,ii], projection = mapcrs) # September

        # Read the data
        
        OMI_data = readOMI_NCDF(infile = home_dir + '/Research/OMI/' + \
            'omi_ai_' + data_types[ii] + '_2005_2019.nc', minlat = minlat)
        
        # Plot trend
        plotOMI_MonthClimo(OMI_data, month_idx = 0,\
            plabel = ' ', minlat=65., ptitle = data_types[ii], pax = axi1)
        plotOMI_MonthClimo(OMI_data, month_idx = 1,\
            plabel = ' ', minlat=65., ptitle = ' ', pax = axi2)
        plotOMI_MonthClimo(OMI_data, month_idx = 2,\
            plabel = ' ', minlat=65., ptitle = ' ', pax = axi3)
        plotOMI_MonthClimo(OMI_data, month_idx = 3,\
            plabel = ' ', minlat=65., ptitle = ' ', pax = axi4)
        plotOMI_MonthClimo(OMI_data, month_idx = 4,\
            plabel = ' ', minlat=65., ptitle = ' ', pax = axi5)
        plotOMI_MonthClimo(OMI_data, month_idx = 5,\
            plabel = ' ', minlat=65., ptitle = ' ', pax = axi6)

    fig.text(0.10, 0.82, 'April', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.692, 'May', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.565, 'June', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.435, 'July', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.305, 'August', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.18, 'September', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)

    
    #fig.text(0.195, 0.89, label_dict[data_types[0]], ha='center', va='center', \
    #    rotation='horizontal',weight='bold',fontsize=row_label_size)
    #fig.text(0.355, 0.89, label_dict[data_types[1]], ha='center', va='center', \
    #    rotation='horizontal',weight='bold',fontsize=row_label_size)
    #fig.text(0.51, 0.89, label_dict[data_types[2]], ha='center', va='center', \
    #    rotation='horizontal',weight='bold',fontsize=row_label_size)
    #fig.text(0.67, 0.89, label_dict[data_types[3]], ha='center', va='center', \
    #    rotation='horizontal',weight='bold',fontsize=row_label_size)

    cax = fig.add_axes([0.13, 0.09, 0.29, 0.01])
    norm = mpl.colors.Normalize(vmin = -2.0, vmax = 3.0)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap = plt.cm.jet, norm = norm, \
        orientation = 'horizontal', extend = 'both')
    cb1.set_label('UV Aerosol Index', weight = 'bold')

    cax2 = fig.add_axes([0.448, 0.09, 0.45, 0.01])
    norm2 = mpl.colors.Normalize(vmin = -0.5, vmax = 0.5)
    cb2 = mpl.colorbar.ColorbarBase(cax2, cmap = plt.cm.bwr, norm = norm2, \
        orientation = 'horizontal', extend = 'both')
    cb2.set_label('AI Trend (2005 - 2020)', weight = 'bold')
    #fig.colorbar(norm, cax = cax, orientation = 'horizontal')

    #fig.tight_layout()

    outname = 'omi_ai_climo_comps_jz.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Generate a 15-panel figure comparing the climatology and trend between 3
# versions of the OMI data for all months
def plotOMI_Type_Trend_all(trend_type = 'standard', \
        minlat=65,save=False):

    colormap = plt.cm.jet

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    colorbar_label_size = 7
    axis_title_size = 8
    row_label_size = 10 

    #fig = plt.figure()
    plt.close('all')
    fig = plt.figure(figsize=(11,13))
    gs = fig.add_gridspec(nrows = 6, ncols = 6, hspace = 0.001, wspace = 0.15)

    data_types = ['VJZ2','VJZ28','VJZ29','VJZ211','VSJ2','VSJ4']    
 
    for ii in range(len(data_types)):
        print(data_types[ii])

        axi1 = plt.subplot(gs[0,ii], projection = mapcrs) # April
        axi2 = plt.subplot(gs[1,ii], projection = mapcrs) # May
        axi3 = plt.subplot(gs[2,ii], projection = mapcrs) # June
        axi4 = plt.subplot(gs[3,ii], projection = mapcrs) # July
        axi5 = plt.subplot(gs[4,ii], projection = mapcrs) # August
        axi6 = plt.subplot(gs[5,ii], projection = mapcrs) # September

        # Read the data
       
        if(data_types[ii] == 'VSJ4'):
            OMI_data = readOMI_NCDF(infile = home_dir + '/Research/OMI/' + \
                'omi_ai_' + data_types[ii] + '_2005_2020.nc', minlat = minlat)
        else:
            OMI_data = readOMI_NCDF(infile = home_dir + '/Research/OMI/' + \
                'omi_ai_' + data_types[ii] + '_2005_2019.nc', minlat = minlat)

        # Plot trend
        plotOMI_MonthTrend(OMI_data, month_idx = 0, trend_type = trend_type,\
            label = ' ', minlat=65., title = data_types[ii], pax = axi1, colorbar = False, \
            colorbar_label_size = colorbar_label_size)
        plotOMI_MonthTrend(OMI_data, month_idx = 1, trend_type = trend_type,\
            label = ' ', minlat=65., title = ' ', pax = axi2, colorbar = False, \
            colorbar_label_size = colorbar_label_size)
        plotOMI_MonthTrend(OMI_data, month_idx = 2, trend_type = trend_type,\
            label = ' ', minlat=65., title = ' ', pax = axi3, colorbar = False, \
            colorbar_label_size = colorbar_label_size)
        plotOMI_MonthTrend(OMI_data, month_idx = 3, trend_type = trend_type,\
            label = ' ', minlat=65., title = ' ', pax = axi4, colorbar = False, \
            colorbar_label_size = colorbar_label_size)
        plotOMI_MonthTrend(OMI_data, month_idx = 4, trend_type = trend_type,\
            label = ' ', minlat=65., title = ' ', pax = axi5, colorbar = False, \
            colorbar_label_size = colorbar_label_size)
        plotOMI_MonthTrend(OMI_data, month_idx = 5, trend_type = trend_type,\
            label = ' ', minlat=65., title = ' ', pax = axi6, colorbar = False, \
            colorbar_label_size = colorbar_label_size)

    fig.text(0.10, 0.82, 'April', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.692, 'May', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.565, 'June', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.435, 'July', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.305, 'August', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.18, 'September', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)

    #fig.text(0.195, 0.89, 'Original\nClimatology', ha='center', va='center', \
    #    rotation='horizontal',weight='bold',fontsize=row_label_size)
    #fig.text(0.355, 0.89, 'Screened\nClimatology', ha='center', va='center', \
    #    rotation='horizontal',weight='bold',fontsize=row_label_size)
    #fig.text(0.51, 0.89, 'Original\nTrend', ha='center', va='center', \
    #    rotation='horizontal',weight='bold',fontsize=row_label_size)
    #fig.text(0.67, 0.89, 'Screened\nTrend', ha='center', va='center', \
    #    rotation='horizontal',weight='bold',fontsize=row_label_size)
    #fig.text(0.83, 0.89, 'Perturbed\nTrend', ha='center', va='center', \
    #    rotation='horizontal',weight='bold',fontsize=row_label_size)

    cax = fig.add_axes([0.13, 0.09, 0.29, 0.01])
    norm = mpl.colors.Normalize(vmin = -2.0, vmax = 3.0)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap = plt.cm.jet, norm = norm, \
        orientation = 'horizontal', extend = 'both')
    cb1.set_label('UV Aerosol Index', weight = 'bold')

    cax2 = fig.add_axes([0.448, 0.09, 0.45, 0.01])
    norm2 = mpl.colors.Normalize(vmin = -0.5, vmax = 0.5)
    cb2 = mpl.colorbar.ColorbarBase(cax2, cmap = plt.cm.bwr, norm = norm2, \
        orientation = 'horizontal', extend = 'both')
    cb2.set_label('AI Trend (2005 - 2020)', weight = 'bold')
    #fig.colorbar(norm, cax = cax, orientation = 'horizontal')

    #fig.tight_layout()

    outname = 'omi_ai_trend_comps_' + trend_type + '.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Compare trends and single monthly averages, and time series of
# monthly averages for two given lat/lon points.
def plotOMI_month_avg_comp(dtype1, dtype2, month_idx, lat1, lon1, \
        lat2, lon2, minlat = 65.):

    #minlat = 65.
    OMI_data1 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
        'omi_ai_' + dtype1 + '_2005_2020.nc', minlat = minlat)
    OMI_data2   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/' + \
        'omi_ai_' + dtype2 + '_2005_2020.nc', minlat = minlat)
    
    # Determine the month being plotted 
    # ---------------------------------
    date_str = OMI_data1['DATES'][month_idx]
    dt_date_str = datetime.strptime(date_str, '%Y%m')
    base_date = datetime.strptime(OMI_data1['DATES'][0], '%Y%m')
    month_offset = dt_date_str.month - base_date.month
    years = np.array([int(datetime.strptime(tval,"%Y%m").year) \
        for tval in OMI_data1['DATES'][month_offset::6]])
 
    plt.close('all')
    fig1 = plt.figure(figsize = (8, 11))
    ax1 = fig1.add_subplot(3,2,1, projection = mapcrs) # JZ211 plot
    ax2 = fig1.add_subplot(3,2,2, projection = mapcrs) # VSJ4 plot
    ax5 = fig1.add_subplot(3,2,3, projection = mapcrs) # JZ211 plot
    ax6 = fig1.add_subplot(3,2,4, projection = mapcrs) # VSJ4 plot
    ax3 = fig1.add_subplot(3,2,5)                      # JZ211 line plot
    ax4 = fig1.add_subplot(3,2,6)                      # VSJ4 line plot

    # Plot the spatial data
    plotOMI_NCDF_SingleMonth(OMI_data1, month_idx, minlat = 65,\
        pax = ax1)
    plotOMI_NCDF_SingleMonth(OMI_data2, month_idx, minlat = 65,\
        pax = ax2)

    # Plot the spatial trends
    plotOMI_MonthTrend(OMI_data1,month_idx=month_offset,trend_type='standard',label = ' ',\
        minlat=65.,title = dtype1, pax = ax5, colorbar = True, \
        colorbar_label_size = 10, show_pval = True)
    plotOMI_MonthTrend(OMI_data2,month_idx=month_offset,trend_type='standard',label = ' ',\
        minlat=65.,title = dtype2, pax = ax6, colorbar = True, \
        colorbar_label_size = 10, show_pval = True)
    

    # Plot points on the map at each lat/lon
    plot_point_on_map(ax1, lat1, lon1, markersize = 10)
    plot_point_on_map(ax1, lat2, lon2, markersize = 10)
    plot_point_on_map(ax2, lat1, lon1, markersize = 10)
    plot_point_on_map(ax2, lat2, lon2, markersize = 10)

    plot_point_on_map(ax5, lat1, lon1, markersize = 10)
    plot_point_on_map(ax5, lat2, lon2, markersize = 10)
    plot_point_on_map(ax6, lat1, lon1, markersize = 10)
    plot_point_on_map(ax6, lat2, lon2, markersize = 10)

    # Plot the time series at each point
    p1_idxs = nearest_gridpoint(lat1, lon1, OMI_data1['LAT'], OMI_data1['LON'])
    p2_idxs = nearest_gridpoint(lat2, lon2, OMI_data1['LAT'], OMI_data1['LON'])
   
    xvals = years
    ax3.plot(xvals, OMI_data1['AI'][month_offset::6,p1_idxs[0][0], p1_idxs[1][0]])
    ax3.plot(xvals, OMI_data2['AI'][month_offset::6,p1_idxs[0][0], p1_idxs[1][0]])
    ax4.plot(xvals, OMI_data1['AI'][month_offset::6,p2_idxs[0][0], p2_idxs[1][0]])
    ax4.plot(xvals, OMI_data2['AI'][month_offset::6,p2_idxs[0][0], p2_idxs[1][0]])

    ax3.set_title('Point blue\n'   + str(lat1) + 'N, ' + str(lon1) + 'E')
    ax4.set_title('Point orange\n' + str(lat2) + 'N, ' + str(lon2) + 'E')

    plot_trend_line(ax3, xvals, OMI_data1['AI'][month_offset::6,p1_idxs[0][0], p1_idxs[1][0]], \
        color='tab:blue', linestyle = '-', slope = 'linear')
    plot_trend_line(ax3, xvals, OMI_data2['AI'][month_offset::6,p1_idxs[0][0], p1_idxs[1][0]], \
        color='tab:orange', linestyle = '-', slope = 'linear')
    plot_trend_line(ax4, xvals, OMI_data1['AI'][month_offset::6,p2_idxs[0][0], p2_idxs[1][0]], \
        color='tab:blue', linestyle = '-', slope = 'linear')
    plot_trend_line(ax4, xvals, OMI_data2['AI'][month_offset::6,p2_idxs[0][0], p2_idxs[1][0]], \
        color='tab:orange', linestyle = '-', slope = 'linear')

    fig1.tight_layout()

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

def plot_OMI_shawn_old(OMI_shawn, ax = None, minlat = 65, labelsize = 12,\
        labelticksize = 11):

    mapcrs = ccrs.NorthPolarStereo()
    datacrs = ccrs.PlateCarree()
    colormap = plt.cm.jet
    
    # Set up the polar stereographic projection map
    in_ax = True
    if(ax is None):
        in_ax = False
        fig1 = plt.figure(figsize=(8,8))
        if(minlat<45):
            ax = plt.axes(projection = ccrs.Miller())
        else:
            ax = plt.axes(projection = ccrs.NorthPolarStereo(central_longitude = 300.))
    #ax.gridlines()
    ax.coastlines(resolution = '50m')
    
    
    # Use meshgrid to convert the 1-d lat/lon arrays into 2-d, which is needed
    # for pcolormesh.
    plot_lat = OMI_shawn['LAT']
    plot_lon = OMI_shawn['LON']
    mask_UVAI = OMI_shawn['AI']
    #plot_lat, plot_lon = np.meshgrid(OMI_shawn['LAT'], OMI_shawn['LON'])
    #mask_UVAI = np.ma.masked_where(OMI_shawn['counts'] == 0, OMI_shawn['AI'])
  
    print(OMI_shawn['date']) 
    dt_date_str = datetime.strptime(OMI_shawn['date'], '%Y%m%d') 
    ax.set_title(dt_date_str.strftime('%d %B %Y\nOMI Daily Average'))
    plot_lat, plot_lon = np.meshgrid(plot_lat, plot_lon)
    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI,transform = datacrs,cmap = colormap,\
            vmin = -1.0, vmax = 3.0, shading='auto')

    gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, \
        linewidth = 1, color = 'gray', alpha = 1.0, linestyle = '-',\
        y_inline = True, xlocs = range(-180, 180, 30), ylocs = range(70, 90, 10))
    
    # Center the figure over the Arctic
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_boundary(circle, transform=ax.transAxes)
    
    # Depending on the desired variable, set the appropriate colorbar ticks
    tickvals = np.arange(-1.0,4.1,0.5)
    
    cbar = plt.colorbar(mesh,\
        ax = ax, orientation='vertical', extend = 'both')
    cbar.ax.tick_params(labelsize=labelticksize)
    cbar.set_label('UV Aerosol Index Perturbation',fontsize=labelsize,weight='bold')


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

    if(home_dir + '/Research/CERES/' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES/')

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
        labelsize = 11, \
        circle_bound = True, gridlines = True, colorbar = True, save=False):

    variable = 'UVAerosolIndex'

    #if(vmin is None):
    #    vmin = var_dict[variable]['min']
    #if(vmax is None):
    #    vmax = var_dict[variable]['max']

    in_ax = True 
    if(pax is None): 
        in_ax = False
        fig0 = plt.figure()
        pax = fig0.add_subplot(1,1,1, projection = ccrs.NorthPolarStereo())

    if(gridlines):
        pax.gridlines()
    pax.coastlines(resolution = '50m')
    if(title == ''):
        title = 'OMI UVAI ' + OMI_hrly['dtype'] + ' ' + OMI_hrly['date']
    pax.set_title(title)
    if((pvar != 'UVAI') & (OMI_hrly['dtype'] != 'shawn')):
        if(pvar == 'ice'): #if(pvar == 'GPQF'):
            # Convert the GPQF flag values to ice flags
            # -----------------------------------------
            work_GPQF = OMI_hrly['GPQF'][:,:,5]
            local_GPQF = np.full(work_GPQF.shape,np.nan)

            local_GPQF[(work_GPQF > 0) & (work_GPQF < 101)] = 9
            local_GPQF[ work_GPQF == 0]                     = 8
            local_GPQF[ work_GPQF == 101]                   = 10
            local_GPQF[ work_GPQF == 103]                   = 11
            local_GPQF[ work_GPQF == 104]                   = 12
            local_GPQF[ work_GPQF == 124]                   = 13

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
        elif(pvar == 'ground'):
            # Convert the GPQF flag values to land/water flags
            # ------------------------------------------------
            work_GPQF = OMI_hrly['GPQF'][:,:,0]
            local_GPQF = np.copy(work_GPQF)
            #local_GPQF = np.full(work_GPQF.shape,np.nan)
            local_GPQF[ work_GPQF > 8]                     = 16

            mask_GPQF = np.ma.masked_invalid(local_GPQF)
            cmap = plt.get_cmap('jet', np.max(mask_GPQF) - np.min(mask_GPQF)+1)
            cbar_labels = ['Shallow ocean','Land','Shallow inland water',\
                            'Ocean coastline / lake shoreline',\
                            'ephemeral (intermittent) water','Deep inland water',\
                            'Continental shelf ocean','Deep ocean',\
                            'Snow-free land','Sea ice',\
                            'Permanent ice\n(Greenland, Antarctica)','Dry snow',\
                            'Ocean','Mixed pixels\nat coastline',\
                            'Suspect ice value','Corners','N/A']
        elif(pvar == 'glint'):
            # Convert the GPQF flag values to glint flags 
            # -------------------------------------------
            work_GPQF = OMI_hrly['GPQF'][:,:,1]
            local_GPQF = np.copy(work_GPQF)
            #local_GPQF = np.full(work_GPQF.shape,np.nan)

            mask_GPQF = np.ma.masked_invalid(local_GPQF)
            cmap = plt.get_cmap('jet', np.max(mask_GPQF) - np.min(mask_GPQF)+1)
            cbar_labels = ['Unaffected','Possible Glint']

        elif(pvar == 'stray'):
            # Convert the GPQF flag values to glint flags 
            # -------------------------------------------
            work_GPQF = OMI_hrly['PXQF'][:,:,0,10]
            local_GPQF = np.copy(work_GPQF)
            #local_GPQF = np.full(work_GPQF.shape,np.nan)

            mask_GPQF = np.ma.masked_invalid(local_GPQF)
            cmap = plt.get_cmap('jet', np.max(mask_GPQF) - np.min(mask_GPQF)+1)
            cbar_labels = ['Unaffected','Stray Light Warning']

        elif(pvar == 'sunshine'):
            # Convert the GPQF flag values to glint flags 
            # -------------------------------------------
            work_GPQF = OMI_hrly['XTRACK'][:,:,4]
            local_GPQF = np.copy(work_GPQF)
            #local_GPQF = np.full(work_GPQF.shape,np.nan)

            mask_GPQF = np.ma.masked_invalid(local_GPQF)
            cmap = plt.get_cmap('jet', np.max(mask_GPQF) - np.min(mask_GPQF)+1)
            cbar_labels = ['Unaffected','Possible Stray Sunshine']

        elif(pvar == 'earthshine'):
            # Convert the GPQF flag values to glint flags 
            # -------------------------------------------
            work_GPQF = OMI_hrly['XTRACK'][:,:,5]
            local_GPQF = np.copy(work_GPQF)
            #local_GPQF = np.full(work_GPQF.shape,np.nan)

            mask_GPQF = np.ma.masked_invalid(local_GPQF)
            cmap = plt.get_cmap('jet', np.max(mask_GPQF) - np.min(mask_GPQF)+1)
            cbar_labels = ['Unaffected','Possible Stray Earthshine']

        mesh = pax.pcolormesh(OMI_hrly['LON'], OMI_hrly['LAT'], mask_GPQF,\
                transform = datacrs, cmap = cmap,\
                vmin = np.nanmin(mask_GPQF) - 0.5, \
                vmax = np.nanmax(mask_GPQF) + 0.5,\
                shading='auto')
        gl = pax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, \
            linewidth = 1, color = 'gray', alpha = 0.5, linestyle = '-',\
            y_inline = True, xlocs = range(-180, 180, 30), ylocs = range(70, 90, 10))

        if(colorbar):
                #cbar = plt.colorbar(mesh,ax = pax, orientation='horizontal',pad=0,\
                cbar = plt.colorbar(mesh,ax = pax, orientation='vertical',\
                    ticks = np.arange(int(np.nanmin(mask_GPQF)), \
                    int(np.nanmax(mask_GPQF)) + 1), pad = 0.04, fraction = 0.040)
                    #shrink = 0.8, ticks = np.arange(np.nanmin(mask_GPQF), \
                    #np.nanmax(mask_GPQF) + 1))
                #cbar.ax.set_xticks(np.arange(int(np.nanmin(mask_GPQF)),int(np.nanmax(mask_GPQF)) + 1))
                print(cbar_labels[int(np.nanmin(mask_GPQF)):int(np.nanmax(mask_GPQF))+1])
                cbar.ax.set_yticklabels(cbar_labels[int(np.nanmin(mask_GPQF)):\
                    int(np.nanmax(mask_GPQF))+1],fontsize=10,weight = 'bold', \
                    rotation=0)
                    #fontsize=8,rotation=35)
        
    else:
        if((OMI_hrly['dtype'] == 'shawn') & (pvar == 'UVAI')):
            pvar = 'UVAI_pert'
        mesh = pax.pcolormesh(OMI_hrly['LON'], OMI_hrly['LAT'], OMI_hrly[pvar],\
                transform = datacrs, cmap = colormap,\
                vmin = vmin, vmax = vmax,\
                shading='auto')
        gl = pax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, \
            linewidth = 1, color = 'gray', alpha = 0.5, linestyle = '-',\
            y_inline = True, xlocs = range(-180, 180, 30), ylocs = range(70, 90, 10))

        if(label == ''):
            if(OMI_hrly['dtype'] == 'shawn'):
                label = 'UVAI Perturbation'
            else:
                label = 'UV Aerosol Index'

        if(colorbar):
            #cbar.set_label(variable,fontsize=16,weight='bold')
            #tickvals = np.arange(vmin,vmax,1.0)
            cbar = plt.colorbar(mesh,ax = pax, \
                orientation='vertical', extend = 'both', fraction = 0.046,\
                pad = 0.04)
            cbar.set_label(label,fontsize = labelsize, weight='bold')
            cbar.ax.tick_params(labelsize= labelsize)

    # Center the figure over the Arctic
    pax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    pax.add_feature(cfeature.BORDERS)
    pax.add_feature(cfeature.STATES)
    if(circle_bound):
        pax.set_boundary(circle, transform=pax.transAxes)

    if(not in_ax):
        plt.show()

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
#     skiprows - 
# --------------------------------------------
def plotOMI_single_swath_figure(date_str, dtype = 'control',  \
        only_sea_ice = False, minlat = 65., skiprows = None, \
        lat_circles = None, save = False, zoom = False, \
        vmin = None, vmax = None, circle_bound = True, ax = None, \
        shawn_path = home_dir + '/data/OMI/shawn_files/'):


    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    if(isinstance(date_str, str)):
        if(dtype == 'shawn'):
            OMI_base  = readOMI_swath_shawn(date_str, latmin = minlat,\
                shawn_path = shawn_path)
        else:
            OMI_base  = readOMI_swath_hdf(date_str, dtype, \
                only_sea_ice = only_sea_ice, latmin = minlat, \
                skiprows = skiprows)
    else:
        OMI_base = date_str
        date_str = OMI_base['date']

    #print(OMI_base['UVAI'][1300,:])

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    in_ax = True 
    if(ax is None): 
        in_ax = False
        plt.close('all')
        fig1 = plt.figure(figsize = (6,6))
        mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = ccrs.Robinson()
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    if(zoom):
        circle_bound = False
    plotOMI_single_swath(ax, OMI_base, title = dtype.title(), \
    #plotOMI_single_swath(ax0, OMI_base, title = 'No row 53', \
        #circle_bound = False, gridlines = False)
        vmin = vmin, vmax = vmax, 
        circle_bound = circle_bound, gridlines = False)

    #ax0.set_extent([-180., , -40., 0.], datacrs)
    #ax0.set_extent([-180,180,-90,90], datacrs)
    if(zoom):
        #ax.set_extent([-90., -20., 75., 87.], datacrs)
        dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
        ax.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lat'][1]],\
                        datacrs)
    else:
        ax.set_extent([-180,180,minlat,90], datacrs)

    plt.suptitle(date_str)

    # ----------------------------------------------------
    # If the user wants circles along latitude lines,    
    # plot them here      
    # ----------------------------------------------------
    if(lat_circles is not None):
        plot_lat_circles(ax, lat_circles) 


    if(not in_ax):
        fig1.tight_layout()

        if(save):
            row_adder = ''
            if(skiprows is not None):
                row_adder = '_no'
                for row in skiprows:
                    row_adder = row_adder + 'r' + str(row + 1)
                
            outname = 'omi_single_swath_figure_' + date_str + '_' + \
                dtype + row_adder + '.png'
            fig1.savefig(outname, dpi=300)
            print("Saved image",outname)
        else:
            plt.show()

# Plot four panels: two different single swaths, as well as
# data zoomed in over Greenland to show the high pixels
# ---------------------------------------------------------
def plotOMI_single_swath_multiple(dtype = 'control',  \
        only_sea_ice = False, minlat = 65., save = False):

    dates = ['200804221027','200804221920','200804221027']
    #dates = ['201206141245','201206141920','200804221027']
    #dates = ['200704020051','200804020057','200704221519','200804221524']
    date_str = dates[0]

    dt_dates = [datetime.strptime(ddate,'%Y%m%d%H%M') for ddate in dates]

    #dates = ['200806141226','200806141902']
    #dates = ['200807261124','200807261303']
    #dates = ['200804221027','200804221345']
    #dates = ['200804221027','200804221206']

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (9,3))
    #mapcrs = ccrs.Robinson()
    #mapcrs = ccrs.LambertConformal()
    mapcrs  = ccrs.NorthPolarStereo(central_longitude = -40)
    ax0 = fig1.add_subplot(1,2,1, projection = mapcrs)
    ax1 = fig1.add_subplot(1,2,2, projection = mapcrs)
    #ax2 = fig1.add_subplot(2,2,3, projection = mapcrs)
    #ax3 = fig1.add_subplot(2,2,4, projection = mapcrs)
    #ax2 = fig1.add_subplot(2,2,1, projection = mapcrs2)
    #ax3 = fig1.add_subplot(2,2,2, projection = mapcrs2)

    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    if(dtype == 'shawn'):
        OMI_base1 = readOMI_swath_shawn(dates[0], latmin = minlat - 5)
        #OMI_base2 = readOMI_swath_shawn(dates[1], latmin = minlat - 5)
        #OMI_base3 = readOMI_swath_shawn(dates[2], latmin = minlat - 5)
        #OMI_base4 = readOMI_swath_shawn(dates[3], latmin = minlat - 5)
    else:
        OMI_base1 = readOMI_swath_hdf(dates[0], dtype, \
            only_sea_ice = only_sea_ice, latmin = minlat - 5)
        #OMI_base2 = readOMI_swath_hdf(dates[1], dtype, \
        #    only_sea_ice = only_sea_ice, latmin = minlat - 5)
        #OMI_base3 = readOMI_swath_hdf(dates[2], dtype, \
        #    only_sea_ice = only_sea_ice, latmin = minlat - 5)
        #OMI_base4 = readOMI_swath_hdf(dates[3], dtype, \
        #    only_sea_ice = only_sea_ice, latmin = minlat - 5)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    plotOMI_single_swath(ax0, OMI_base1, title = None, \
        circle_bound = False, gridlines = False, vmax = 3.0)
    #plotOMI_single_swath(ax1, OMI_base2, title = dates[1], \
    plotOMI_single_swath(ax1, OMI_base1, title = None, \
        pvar = 'ice', circle_bound = False, gridlines = False)
    #plotOMI_single_swath(ax2, OMI_base3, title = dates[2], \
    #plotOMI_single_swath(ax2, OMI_base3, title = dt_dates[2].strftime('%H:%M UTC %m/%d/%Y'), \
    #    circle_bound = True, gridlines = False, vmax = 3.0)
    #plotOMI_single_swath(ax3, OMI_base4, title = dates[3], \
    #plotOMI_single_swath(ax3, OMI_base4, title = dt_dates[3].strftime('%H:%M UTC %m/%d/%Y'), \
    #    circle_bound = True, gridlines = False, vmax = 3.0)
    ##!#plotOMI_single_swath(ax3, OMI_base3, title = dates[2], \
    ##!#    pvar = 'GPQF', circle_bound = False, gridlines = False)

    #ax2.set_extent([-180., 180., minlat, 90.], datacrs)
    #ax3.set_extent([-180., 180., minlat, 90.], datacrs)
    ax0.set_extent([-70., -10., 65., 87.], datacrs)
    ax1.set_extent([-70., -10., 65., 87.], datacrs)
    #ax2.set_extent([-80., -28., 75., 87.], datacrs)
    #ax3.set_extent([-80., -28., 75., 87.], datacrs)

    ##!## ----------------------------------------------------
    ##!## Draw boxes for the comparison ranges
    ##!## ----------------------------------------------------
    ##!#lat_corners = np.array([61, 73, 73, 61]) 
    ##!#lon_corners = np.array([-140, -140, -160, -160]) 
    ##!#poly_corners = np.zeros((len(lat_corners), 2), np.float64)
    ##!#poly_corners[:,0] = lon_corners
    ##!#poly_corners[:,1] = lat_corners
    ##!#poly = mpatches.Polygon(poly_corners, closed = True, ec = 'r', \
    ##!#    fill = False, lw = 5, fc = None, transform = datacrs)
    ##!#ax0.add_patch(poly)
    ##!#ax1.add_patch(poly)
    ##!###!#ax0.add_patch(mpatches.Rectangle(xy = [-120, 67], width = 40, \
    ##!###!#    height = 12, alpha = 0.9, transform = datacrs, facecolor = 'none',\
    ##!###!#    edgecolor = 'red', linewidth = 5))
    ##!###!#ax1.add_patch(mpatches.Rectangle(xy = [-120, 67], width = 40, \
    ##!###!#    height = 12, alpha = 0.9, transform = datacrs, facecolor = 'none',\
    ##!###!#    edgecolor = 'red', linewidth = 5))


    plt.suptitle(dt_dates[0].strftime('%d %B %Y, %H:%M:%S UTC'))
    plot_subplot_label(ax0, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax1, '(b)', backgroundcolor = 'white')

    fig1.tight_layout()

    if(save):
        outname = 'omi_single_swath_figure_multiple_v4.png'
        fig1.savefig(outname, dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Plot four panels: two different single swaths, as well as
# data zoomed in over Greenland to show the high pixels
# ---------------------------------------------------------
def plotOMI_single_swath_multiple_v3(dtype = 'control',  \
        only_sea_ice = False, minlat = 65., save = False):

    dates = ['200704020051','200804020057','200704221519','200804221524']
    date_str = dates[0]

    dt_dates = [datetime.strptime(ddate,'%Y%m%d%H%M') for ddate in dates]

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (9.5,8))
    #mapcrs = ccrs.Robinson()
    #mapcrs = ccrs.LambertConformal()
    mapcrs  = ccrs.NorthPolarStereo(central_longitude = -40)
    ax0 = fig1.add_subplot(2,2,1, projection = mapcrs)
    ax1 = fig1.add_subplot(2,2,2, projection = mapcrs)
    ax2 = fig1.add_subplot(2,2,3, projection = mapcrs)
    ax3 = fig1.add_subplot(2,2,4, projection = mapcrs)

    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    if(dtype == 'shawn'):
        OMI_base1 = readOMI_swath_shawn(dates[0], latmin = minlat - 5)
        OMI_base2 = readOMI_swath_shawn(dates[1], latmin = minlat - 5)
        OMI_base3 = readOMI_swath_shawn(dates[2], latmin = minlat - 5)
        OMI_base4 = readOMI_swath_shawn(dates[3], latmin = minlat - 5)
    else:
        OMI_base1 = readOMI_swath_hdf(dates[0], dtype, \
            only_sea_ice = only_sea_ice, latmin = minlat - 5)
        OMI_base2 = readOMI_swath_hdf(dates[1], dtype, \
            only_sea_ice = only_sea_ice, latmin = minlat - 5)
        OMI_base3 = readOMI_swath_hdf(dates[2], dtype, \
            only_sea_ice = only_sea_ice, latmin = minlat - 5)
        OMI_base4 = readOMI_swath_hdf(dates[3], dtype, \
            only_sea_ice = only_sea_ice, latmin = minlat - 5)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    colorbar = True
    plotOMI_single_swath(ax0, OMI_base1, title = dt_dates[0].strftime('%d %B %Y, %H:%M:%S UTC'), \
        circle_bound = True, gridlines = False, vmax = 3.0, colorbar = colorbar)
    plotOMI_single_swath(ax1, OMI_base2, title = dt_dates[1].strftime('%d %B %Y, %H:%M:%S UTC'), \
        circle_bound = True, gridlines = False, vmax = 3.0, colorbar = colorbar)
    plotOMI_single_swath(ax2, OMI_base3, title = dt_dates[2].strftime('%d %B %Y, %H:%M:%S UTC'), \
        circle_bound = True, gridlines = False, vmax = 3.0, colorbar = colorbar)
    plotOMI_single_swath(ax3, OMI_base4, title = dt_dates[3].strftime('%d %B %Y, %H:%M:%S UTC'), \
        circle_bound = True, gridlines = False, vmax = 3.0, colorbar = colorbar)
    ##!#plotOMI_single_swath(ax3, OMI_base3, title = dates[2], \
    ##!#    pvar = 'GPQF', circle_bound = False, gridlines = False)

    #cmap = mpl.cm.jet
    #norm = mpl.colors.Normalize(vmin = -2, vmax = 3)
    #
    #cbar_ax = fig1.add_axes([0.90, 0.05, 0.05, 0.9])
    #fig1.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmap), \
    #    cax = cbar_ax, orientation = 'vertical', label = 'UV Aerosol Index')

    ax0.set_extent([-180., 180., minlat, 90.], datacrs)
    ax1.set_extent([-180., 180., minlat, 90.], datacrs)
    ax2.set_extent([-180., 180., minlat, 90.], datacrs)
    ax3.set_extent([-180., 180., minlat, 90.], datacrs)
    #ax0.set_extent([-70., -10., 65., 87.], datacrs)
    #ax1.set_extent([-70., -10., 65., 87.], datacrs)
    #ax2.set_extent([-80., -28., 75., 87.], datacrs)
    #ax3.set_extent([-80., -28., 75., 87.], datacrs)

    ##!## ----------------------------------------------------
    ##!## Draw boxes for the comparison ranges
    ##!## ----------------------------------------------------
    ##!#lat_corners = np.array([61, 73, 73, 61]) 
    ##!#lon_corners = np.array([-140, -140, -160, -160]) 
    ##!#poly_corners = np.zeros((len(lat_corners), 2), np.float64)
    ##!#poly_corners[:,0] = lon_corners
    ##!#poly_corners[:,1] = lat_corners
    ##!#poly = mpatches.Polygon(poly_corners, closed = True, ec = 'r', \
    ##!#    fill = False, lw = 5, fc = None, transform = datacrs)
    ##!#ax0.add_patch(poly)
    ##!#ax1.add_patch(poly)
    ##!###!#ax0.add_patch(mpatches.Rectangle(xy = [-120, 67], width = 40, \
    ##!###!#    height = 12, alpha = 0.9, transform = datacrs, facecolor = 'none',\
    ##!###!#    edgecolor = 'red', linewidth = 5))
    ##!###!#ax1.add_patch(mpatches.Rectangle(xy = [-120, 67], width = 40, \
    ##!###!#    height = 12, alpha = 0.9, transform = datacrs, facecolor = 'none',\
    ##!###!#    edgecolor = 'red', linewidth = 5))


    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')
    plot_subplot_label(ax2, '(c)')
    plot_subplot_label(ax3, '(d)')

    fig1.tight_layout()

    if(save):
        outname = 'omi_single_swath_figure_multiple_v3.png'
        fig1.savefig(outname, dpi=300)
        print("Saved image",outname)
    else:
        plt.show()
# Plot four panels: two different single swaths, as well as
# data zoomed in over Greenland to show the high pixels
# ---------------------------------------------------------
def plotOMI_single_swath_multiple_v4(dtype = 'control',  \
        only_sea_ice = False, minlat = 65., save = False):

    dates = ['200804221027','200804221920','200804221027']
    date_str = dates[0]

    dt_dates = [datetime.strptime(ddate,'%Y%m%d%H%M') for ddate in dates]

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (9,3))
    mapcrs  = ccrs.NorthPolarStereo(central_longitude = -40)
    ax0 = fig1.add_subplot(1,2,1, projection = mapcrs)
    ax1 = fig1.add_subplot(1,2,2, projection = mapcrs)

    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    if(dtype == 'shawn'):
        OMI_base1 = readOMI_swath_shawn(dates[0], latmin = minlat - 5)
    else:
        OMI_base1 = readOMI_swath_hdf(dates[0], dtype, \
            only_sea_ice = only_sea_ice, latmin = minlat - 5)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    plotOMI_single_swath(ax0, OMI_base1, title = None, \
        circle_bound = False, gridlines = False, vmax = 3.0)
    plotOMI_single_swath(ax1, OMI_base1, title = None, \
        pvar = 'ice', circle_bound = False, gridlines = False)

    ax0.set_extent([-70., -10., 65., 87.], datacrs)
    ax1.set_extent([-70., -10., 65., 87.], datacrs)

    plt.suptitle(dt_dates[0].strftime('%d %B %Y, %H:%M:%S UTC'))
    plot_subplot_label(ax0, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax1, '(b)', backgroundcolor = 'white')

    fig1.tight_layout()

    if(save):
        outname = 'omi_single_swath_figure_multiple_v4.png'
        fig1.savefig(outname, dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# 2 panel plot showing control UVAI and ice classification
# (formerly "ground")
# -----------------------------------------------------------
def plotOMI_single_flags(date_str, pvar = 'ice', \
        only_sea_ice = False, minlat = 65., \
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
    else:
        if(minlat > 50):
            mapcrs = ccrs.NorthPolarStereo()
            figsize = (11,5)
        else:
            mapcrs = ccrs.Robinson()
            figsize = (11,5)
    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    if(multi_panel):
        fig1 = plt.figure(figsize = (10,10))
        ax0 = fig1.add_subplot(2,2,1, projection = mapcrs)
        ax1 = fig1.add_subplot(2,2,2, projection = mapcrs)
    else:
        fig1 = plt.figure(figsize = figsize)
        ax0 = fig1.add_subplot(1,2,1, projection = mapcrs)
        ax1 = fig1.add_subplot(1,2,2, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    pvar_dict = {
        'stray': 'Possible Stray Light',\
        'earthshine': 'Possible Stray Earthshine', \
        'sunshine': 'Possible Stray Sunshine', \
        'glint': 'Possible Sun Glint',\
    }
    gridlines = False
    plotOMI_single_swath(ax0, OMI_base, title = 'Control', \
        circle_bound = circle_bound, gridlines = gridlines)
    plotOMI_single_swath(ax1, OMI_base, title = pvar_dict[pvar], \
        pvar = pvar, circle_bound = circle_bound, gridlines = gridlines)
    #plotOMI_single_swath(ax1, OMI_base, title = 'Screened', pvar = 'GPQF', \

    if(zoom):
        ax0.set_extent([-90., -20., 75., 87.], datacrs)
        ax1.set_extent([-90., -20., 75., 87.], datacrs)

    plt.suptitle(date_str)
    plot_subplot_label(ax0, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax1, '(b)', backgroundcolor = 'white')

    fig1.tight_layout()

    if(save):
        outname = 'omi_single_swath_' + pvar + '_' + date_str + \
            zoom_add + '.png'
        fig1.savefig(outname, dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Compare the control AI with the two methods
# -------------------------------------------
def plotOMI_single_multipanel(date_str, only_sea_ice = False, minlat = 65., \
        quad_panel = False, save = False):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # ----------------------------------------------------
    # Read in each of the 3 data types
    # ----------------------------------------------------
    OMI_base  = readOMI_swath_hdf(date_str, 'control', \
        only_sea_ice = only_sea_ice, latmin = minlat)
    OMI_JZ    = readOMI_swath_hdf(date_str, 'JZ', \
        only_sea_ice = only_sea_ice, latmin = minlat)
    OMI_shawn = readOMI_swath_shawn(date_str, latmin = minlat, \
        shawn_path = home_dir + '/data/OMI/shawn_files/ltc3/')

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    if(quad_panel):
        fig1 = plt.figure(figsize = (9.5,8))
        ax0 = fig1.add_subplot(2,2,1, projection = mapcrs)
        ax1 = fig1.add_subplot(2,2,2, projection = mapcrs)
        ax2 = fig1.add_subplot(2,2,4, projection = mapcrs)
        ax3 = fig1.add_subplot(2,2,3, projection = mapcrs)
    else:
        fig1 = plt.figure(figsize = (15,5))
        ax0 = fig1.add_subplot(1,3,1, projection = mapcrs)
        ax1 = fig1.add_subplot(1,3,2, projection = mapcrs)
        ax2 = fig1.add_subplot(1,3,3, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    plotOMI_single_swath(ax0, OMI_base, title = 'Original', gridlines = False)
    plotOMI_single_swath(ax1, OMI_JZ, title = 'Screened', gridlines = False)
    plotOMI_single_swath(ax2, OMI_shawn, pvar = 'UVAI_pert',\
        title = 'Perturbed', gridlines = False)
    if(quad_panel):
        plotOMI_single_swath(ax3, OMI_shawn, pvar = 'UVAI_climo',\
            title = 'Climatology', label = 'UVAI Climatology', gridlines = False)

    plt.suptitle(dt_date_str.strftime('%d %B %Y, %H:%M:%S UTC'))
    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')
    if(quad_panel):
        plot_subplot_label(ax2, '(d)')
        plot_subplot_label(ax3, '(c)')
    else:
        plot_subplot_label(ax2, '(c)')

    fig1.tight_layout()

    if(save):
        if(quad_panel):
            adder = '4panel'
        else:
            adder = '3panel'
        outname = 'omi_single_swath_compare_'+adder+'_' + date_str + '.png'
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

# Plot a daily average of OMI data
# --------------------------------
def plotOMI_daily(date_str, dtype = 'control',  \
        only_sea_ice = False, minlat = 65., skiprows = None, \
        lat_circles = None, ax = None, save = False, resolution = 0.25, \
        row_max = 60, row_min = 0, circle_bound = True, colorbar = True, \
        shawn_path = home_dir + '/data/OMI/shawn_files/'):

    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    if(dtype == 'shawn'):
        ##!#OMI_base  = readOMI_swath_shawn(date_str, latmin = minlat,\
        ##!#    shawn_path = shawn_path)
        OMI_base = readOMI_swath_shawn_old(date_str, latmin = minlat, 
            resolution = resolution)
    else:
        OMI_base = readOMI_single_swath(date_str, row_max, row_min = row_min,\
            only_sea_ice = only_sea_ice, dtype = dtype, latmin = minlat,\
            resolution = resolution)
        ##!#OMI_base  = readOMI_swath_hdf(date_str, dtype, \
        ##!#    only_sea_ice = only_sea_ice, latmin = minlat, \
        ##!#    skiprows = skiprows)

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    in_ax = True 
    if(ax is None): 
        in_ax = False
        plt.close('all')
        fig1 = plt.figure()
        #fig1 = plt.figure(figsize = (6,6))
        mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = ccrs.Robinson()
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    cmap = 'jet'
    y_val, x_val = np.meshgrid(OMI_base['LAT'], OMI_base['LON'])
    print(x_val.shape, y_val.shape, OMI_base['AI'].shape)
    mesh = ax.pcolormesh(x_val, y_val, OMI_base['AI'],\
            transform = datacrs, cmap = cmap,\
            vmin = -2.0, vmax = 4.0, \
            ##!#vmin = np.nanmin(mask_GPQF) - 0.5, \
            ##!#vmax = np.nanmax(mask_GPQF) + 0.5,\
            shading='auto')
    ax.coastlines()

    if(colorbar):
            #cbar = plt.colorbar(mesh,ax = pax, orientation='horizontal',pad=0,\
            cbar = plt.colorbar(mesh,ax = ax, orientation='vertical',\
                pad = 0.04, fraction = 0.040, label = 'UVAI')
            cbar.set_label('UVAI', fontsize = 12, weight='bold')
                #int(np.nanmax(mask_GPQF)) + 1), pad = 0.04, fraction = 0.040)
                #ticks = np.arange(int(np.nanmin(mask_GPQF)), \
                #shrink = 0.8, ticks = np.arange(np.nanmin(mask_GPQF), \
                #np.nanmax(mask_GPQF) + 1))
            #cbar.ax.set_xticks(np.arange(int(np.nanmin(mask_GPQF)),int(np.nanmax(mask_GPQF)) + 1))
            ##!#cbar.ax.set_yticklabels(cbar_labels[int(np.nanmin(mask_GPQF)):\
            ##!#    int(np.nanmax(mask_GPQF))+1],fontsize=10,weight = 'bold', \
            ##!#    rotation=0)
                #fontsize=8,rotation=35)
        
    if(circle_bound):
        ax.set_boundary(circle, transform=ax.transAxes)

    #ax0.set_extent([-180., , -40., 0.], datacrs)
    ax.set_extent([-180,180,minlat,90], datacrs)
    #ax0.set_extent([-180,180,-90,90], datacrs)

    plt.suptitle(date_str)

    # ----------------------------------------------------
    # If the user wants circles along latitude lines,    
    # plot them here      
    # ----------------------------------------------------
    if(lat_circles is not None):
        plot_lat_circles(ax0, lat_circles) 


    if(not in_ax):
        fig1.tight_layout()
        if(save):
            row_adder = ''
            if(skiprows is not None):
                row_adder = '_no'
                for row in skiprows:
                    row_adder = row_adder + 'r' + str(row + 1)
                
            outname = 'omi_single_swath_figure_' + date_str + '_' + \
                dtype + row_adder + '_latlines.png'
            fig1.savefig(outname, dpi=300)
            print("Saved image",outname)
        else:
            plt.show()

def plotOMI_daily_control_shawn(date_str, only_sea_ice = False, minlat = 65.,\
        lat_circles = None, ax = None, save = False, resolution = 0.25, \
        row_max = 60, row_min = 0, circle_bound = True, colorbar = True, \
        shawn_path = home_dir + '/data/OMI/shawn_files/', \
        save_dir = home_dir + '/Research/OMI/OMI_CERES_compares/'):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES')
    import gridCERESLib
    importlib.reload(gridCERESLib) 
    from gridCERESLib import plotCERES_daily_figure

    ##!## ----------------------------------------------------
    ##!## Read in the shawn and control data
    ##!## ----------------------------------------------------
    ##!#OMI_shawn = readOMI_swath_shawn_old(date_str, latmin = minlat, 
    ##!#        resolution = resolution)
    ##!#OMI_ctrl = readOMI_single_swath(date_str, row_max, row_min = row_min,\
    ##!#    only_sea_ice = only_sea_ice, dtype = 'control', latmin = minlat,\
    ##!#    resolution = resolution)

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (12,4))
    ax0 = fig1.add_subplot(1,3,1, projection = mapcrs)
    ax1 = fig1.add_subplot(1,3,2, projection = mapcrs)
    ax2 = fig1.add_subplot(1,3,3, projection = mapcrs)

    plotOMI_daily(date_str, dtype = 'control',  \
        only_sea_ice = False, minlat = minlat, skiprows = None, \
        lat_circles = None, ax = ax0, save = False, resolution = resolution, \
        row_max = 60, row_min = row_min, circle_bound = True, colorbar = False, \
        shawn_path = home_dir + '/data/OMI/shawn_files/')

    plotOMI_daily(date_str, dtype = 'shawn',  \
        only_sea_ice = False, minlat = minlat, skiprows = None, \
        lat_circles = None, ax = ax1, save = False, resolution = resolution, \
        row_max = 60, row_min = row_min, circle_bound = True, colorbar = True, \
        shawn_path = home_dir + '/data/OMI/shawn_files/')

    plotCERES_daily_figure(date_str, param = 'SWF',  \
        only_sea_ice = False, minlat = minlat, \
        lat_circles = None, ax = ax2, save = False, resolution = resolution, \
        circle_bound = circle_bound, colorbar = True)

    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')
    plot_subplot_label(ax2, '(c)')

    fig1.tight_layout()

    if(save):
        outname = save_dir + 'omi_ceres_daily_compare_' + date_str + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# Compare the quarter-degree coverage as a function of latitude for 
# single-swath daily averages (using the old read single swath function) 
# using all good rows and only the screened data
# ----------------------------------------------------------------------
def plot_Arctic_row_coverage_compare(date_str = '20180726', minlat = 65., \
        save = False):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d')
    
    # Read the data using all good rows
    # ---------------------------------
    OMI_data = readOMI_single_swath('20180726', 60, only_sea_ice = False, \
        latmin = minlat)

    # Read only the last 5 rows
    # -------------------------
    OMI_data2 = readOMI_single_swath('20180726', 60, only_sea_ice = False, \
        row_min = 55, latmin = minlat)

    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(OMI_data['LAT'], (np.count_nonzero(OMI_data['AI_count'], \
        axis = 0) / OMI_data['AI_count'][:,0].size) * 100., \
        label = 'Rows 1 - 22, 56 - 60')
   
    ax.plot(OMI_data2['LAT'], (np.count_nonzero(OMI_data2['AI_count'], \
        axis = 0) / OMI_data2['AI_count'][:,0].size) * 100., \
        label = 'Rows 56 - 60') 

    ax.grid()
    ax.set_xlabel('Latitude')
    ax.set_ylabel('%')
    ax.set_title('Percent coverage of quarter-degree grid boxes\n' + \
        dt_date_str.strftime('%Y-%m-%d') + ' daily average')
    ax.legend()
   
    if(save):
        outname = 'omi_coverage_compare_' + date_str + '_' + \
            str(int(minlat)) + 'to90.png' 
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

##!## Plot a single swath of OMI data with total climatology subtracted
##!## mask_weakAI: removes AI values below 0.8 when plotting
##!#def single_swath_anomaly_climo(OMI_data,swath_date,month_climo = True,\
##!#        minlat = 60,row_max = 60,mask_weakAI = False, save=False): 
##!#    # - - - - - - - - - - - - - - - -
##!#    # Read in the single swath data
##!#    # - - - - - - - - - - - - - - - -
##!#    latmin = 60 
##!#    
##!#    # Set up mapping variables 
##!#    datacrs = ccrs.PlateCarree() 
##!#    colormap = plt.cm.jet
##!#    if(minlat < 45):
##!#        mapcrs = ccrs.Miller()
##!#    else:
##!#        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)
##!#
##!#    # Set up values for gridding the AI data
##!#    lat_gridder = latmin * 4.
##!#    
##!#    lat_ranges = np.arange(latmin,90.1,0.25)
##!#    lon_ranges = np.arange(-180,180.1,0.25)
##!#    base_path = 'home_dir + /data/OMI/H5_files/'
##!#    year = swath_date[:4]
##!#    date = swath_date[4:8]
##!#    if(len(swath_date)==13):
##!#        time = swath_date[9:]
##!#    elif(len(swath_date)==12):
##!#        time = swath_date[8:]
##!#    else:
##!#        time = ''
##!#    total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+year+'m'+date+'t'+time+'*.he5',\
##!#              shell=True).decode('utf-8').strip().split('\n')
##!#
##!#    UVAI = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
##!#    count = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
##!#
##!#    for fileI in range(len(total_list)):
##!#        # read in data directly from HDF5 files
##!#        print(total_list[fileI])
##!#        data = h5py.File(total_list[fileI],'r')
##!#        #ALBEDO = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo']
##!#        #REFLECTANCE = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity']
##!#        #CLD = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/CloudFraction']
##!#        AI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
##!#        #PIXEL= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/PixelQualityFlags']
##!#        #MSMNT= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/MeasurementQualityFlags']
##!#        #VZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/ViewingZenithAngle']
##!#        #SZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SolarZenithAngle']
##!#        #RAZ  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle']
##!#        LAT = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
##!#        LON = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
##!#        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
##!#        #GRND = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags']
##!#    
##!#        #albedo = ALBEDO[:,:,0]   
##!#        #reflectance = REFLECTANCE[:,:,0]   
##!#        counter = 0
##!#        #AI = AI[:,:,0]   
##!#        # Loop over the values and rows 
##!#        #for i in range(0,int(CBA2)):
##!#        #for i in range(albedo.shape[0]):
##!#        for i in range(AI.shape[0]):
##!#            for j in range(0,row_max):
##!#                #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
##!#                if(AI[i,j]>-20):
##!#                #if(plotAI[i,j]>-20):
##!#                    # Only plot if XTrack flag is met
##!#                    if((XTRACK[i,j] == 0) | ((XTRACK[i,j] & 4 == 4))):
##!#                        # Print values to text file
##!#                        if(LAT[i,j] > minlat):
##!#                            counter+=1
##!#    #                        fout.write("{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.6f} {9:.6f} {10:.6f} {11:.6f}\n".format(\
##!#    #                            LAT[i,j],LON[i,j],AI[i,j],0.5,SZA[i,j],VZA[i,j],RAZ[i,j], \
##!#    #                            ALBEDO[i,j,0],ALBEDO[i,j,1],REFLECTANCE[i,j,0],\
##!#    #                            REFLECTANCE[i,j,1],CLD[i,j]))
##!#    
##!#    
##!#                        #if((plotXTrack[i,j] == 0) | ((plotXTrack[i,j] & 4 == 4))):
##!#                            index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
##!#                            index2 = int(np.floor(LON[i,j]*4 + 720.))
##!#                            #index1 = int(np.floor(plotLAT[i,j]*4 + 360.))
##!#                            #index2 = int(np.floor(plotLON[i,j]*4 + 720.))
##!#                            
##!#                            if(index1 < 0): index1 = 0
##!#                            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
##!#                            if(index2 < 0): index2 = 0                                                                                            
##!#                            if(index2 > 1439): index2 = 1439
##!#                       
##!#                            #diff = reflectance[i,j] - albedo[i,j]
##!#                            #if(diff<min_diff):
##!#                            #    min_diff = diff
##!#                            #if(diff>max_diff):
##!#                            #    max_diff = diff
##!#                       
##!#                            #if(diff<0.2): 
##!#                            #    UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
##!#                            UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
##!#                            #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
##!#                                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
##!#                                #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
##!#                            count[index2, index1] = count[index2,index1] + 1
##!#    
##!#    # Calculate the row-average AI for the secondary plot
##!#    mask_avgs = np.nanmean(np.ma.masked_where(AI[:,:] < -20, AI[:,:]),axis=0)
##!#
##!#
##!#    if(month_climo == True):
##!#        # Determine which month to use based on the swath date
##!#        monther = int(date[:2])-1
##!#        local_climo = OMI_data['MONTH_CLIMO'][monther,:,:].T
##!#        title_adder = 'Monthly Anomaly '
##!#        file_adder = 'mnth_anom_'
##!#    else:
##!#        local_data = np.copy(OMI_data['AI'][:,:,:])
##!#        local_mask = np.ma.masked_where(local_data == -999.9, local_data)
##!#        local_climo = np.nanmean(local_mask, axis=0).T
##!#        title_adder = 'Anomaly '
##!#        file_adder = 'anom_'
##!#
##!#    # Interpolate 2D data
##!#    # Find where the high-res data fits in to the climo data
##!#    begin_x = np.where(OMI_data['LAT'][:,0] == lat_ranges[0])[0][0]
##!#    high_lats = np.arange(OMI_data['LAT'][begin_x,0], OMI_data['LAT'][-1,0]+0.25, 0.25)
##!#    high_lons = np.arange(OMI_data['LON'][0,0], OMI_data['LON'][0,-1] + 0.25, 0.25)
##!#
##!#    high_climo = np.zeros((high_lons.shape[0], high_lats.shape[0]))
##!#
##!#    for yi in range(len(OMI_data['LON'][0,:])):
##!#        for xj in range(len(OMI_data['LAT'][begin_x:,0])):
##!#            high_climo[yi*4:yi*4+4, xj*4:xj*4+4] = local_climo[yi, xj+begin_x]
##!#
##!#    # Find out where the high res climo data fits in to the
##!#    # single swath data. The single swath data are assumed to be 
##!#    # of a larger size than the climo data
##!#    end_xj = np.where(lat_ranges == high_lats[-1])[0][0]
##!#    end_yi = np.where(lon_ranges == high_lons[-1])[0][0]
##!#
##!#    # Calculate anomalies
##!#    UVAI_anomaly = UVAI[0:end_yi+1, 0:end_xj+1] - high_climo  
##!#
##!#    # Mask grid boxes where the ob counts are zero
##!#    plot_lat, plot_lon = np.meshgrid(high_lats,high_lons)
##!#    mask_UVAI_anom = np.ma.masked_where(count[0:end_yi+1, 0:end_xj+1] == 0, UVAI_anomaly)
##!#  
##!#    ## Mask any data below AI values of 0.8
##!#    #mask_UVAI_anom = np.ma.masked_where(mask_UVAI_anom < 0.8, mask_UVAI_anom)
##!# 
##!#    plt.close('all')
##!#    fig1 = plt.figure(figsize=(8,8))
##!#    ax = plt.axes(projection = mapcrs)
##!#    ax.gridlines()
##!#    ax.coastlines(resolution='50m')
##!# 
##!#    plt.title('OMI AI '+title_adder + swath_date)
##!#    #plt.title('OMI Reflectivity - Surface Albedo '+swath_date)
##!#    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI_anom,transform = datacrs,cmap = colormap,\
##!#            vmin = -2.0, vmax = 3.0)
##!#    ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
##!#    ax.set_xlim(-3430748.535086173,3430748.438879491)
##!#    ax.set_ylim(-3413488.8763307533,3443353.899053069)
##!#    #ax.set_xlim(-4170748.535086173,4167222.438879491)
##!#    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
##!#    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
##!#        aspect=50,shrink = 0.905,label='Aerosol Index Anomaly')
##!#    ##cax = fig.add_axes([0.16,0.075,0.7,0.025])
##!#    ##cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
##!#    ##cb.ax.set_xlabel('Aerosol Index')
##!#    #cb.ax.set_xlabel('Reflectivity - Surface Albedo')
##!#    #out_name = 'omi_single_pass_ai_200804270052_to_0549_composite_rows_0to'+str(row_max)+'.png'       
##!#    if(save == True):
##!#        out_name = 'omi_single_pass_ai_'+file_adder+swath_date+'_rows_0to'+str(row_max)+'.png'       
##!#        plt.savefig(out_name)
##!#        print('Saved image '+out_name)
##!#    else: 
##!#        plt.show()
##!#    #axs[1].plot(mask_avgs)
##!#    #axs[1].set_xlabel('Sensor Row')
##!#    #axs[1].set_ylabel('Row Average Aerosol Index')
    

##!## Plot a single swath of OMI data with single-time swath climatology subtracted
##!## single_swath = the swath that will be corrected (format = YYYYMMDDHHMM)
##!## climo_date   = the date on the swath climatology directory (format = MMDDHHMM)
##!#def single_swath_anomaly_time(single_swath,climo_date,minlat = 60,row_max = 60): 
##!#
##!#    # - - - - - - - - - - - - - - - -
##!#    # Read in the single swath data
##!#    # - - - - - - - - - - - - - - - -
##!#    latmin = 60 
##!#    
##!#    # Set up mapping variables 
##!#    datacrs = ccrs.PlateCarree() 
##!#    colormap = plt.cm.jet
##!#    if(minlat < 45):
##!#        mapcrs = ccrs.Miller()
##!#    else:
##!#        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)
##!#
##!#    # Set up values for gridding the AI data
##!#    lat_gridder = latmin * 4.
##!#    
##!#    lat_ranges = np.arange(latmin,90.1,0.25)
##!#    lon_ranges = np.arange(-180,180.1,0.25)
##!#    climo_base_path = 'home_dir + /data/OMI/swath_anomaly_files/'
##!#    single_base_path = 'home_dir + /data/OMI/H5_files/'
##!#
##!#    # Grab the path date/time for finding the swath climatology files
##!#    year = single_swath[0:4]
##!#    date = single_swath[4:8]
##!#    if(len(single_swath)==13):
##!#        time = single_swath[9:]
##!#    elif(len(single_swath)==12):
##!#        time = single_swath[8:]
##!#    else:
##!#        time = ''
##!#
##!#    total_climo_list = subprocess.check_output('ls '+climo_base_path+climo_date+'/*.he5',\
##!#              shell=True).decode('utf-8').strip().split('\n')
##!#
##!#    total_single_list = subprocess.check_output('ls '+single_base_path+'OMI-Aura_L2-OMAERUV_'+\
##!#              year+'m'+date+'t'+time+'*.he5',shell=True).decode('utf-8').strip().split('\n')
##!#
##!#    UVAI_climo = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
##!#    count_climo = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
##!#    UVAI_single = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
##!#    count_single = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
##!#
##!#    print("Reading path climatology data")
##!#    for fileI in range(len(total_climo_list)):
##!#        # read in data directly from HDF5 files
##!#        print(total_climo_list[fileI])
##!#        data = h5py.File(total_climo_list[fileI],'r')
##!#        #ALBEDO = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo']
##!#        #REFLECTANCE = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity']
##!#        #CLD = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/CloudFraction']
##!#        AI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
##!#        #PIXEL= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/PixelQualityFlags']
##!#        #MSMNT= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/MeasurementQualityFlags']
##!#        #VZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/ViewingZenithAngle']
##!#        #SZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SolarZenithAngle']
##!#        #RAZ  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle']
##!#        LAT = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
##!#        LON = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
##!#        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
##!#        #GRND = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags']
##!#    
##!#        #albedo = ALBEDO[:,:,0]   
##!#        #reflectance = REFLECTANCE[:,:,0]   
##!#        counter = 0
##!#        #AI = AI[:,:,0]   
##!#        # Loop over the values and rows 
##!#        #for i in range(0,int(CBA2)):
##!#        #for i in range(albedo.shape[0]):
##!#        for i in range(AI.shape[0]):
##!#            for j in range(0,row_max):
##!#                #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
##!#                if(AI[i,j]>-20):
##!#                #if(plotAI[i,j]>-20):
##!#                    # Only plot if XTrack flag is met
##!#                    if((XTRACK[i,j] == 0) | ((XTRACK[i,j] & 4 == 4))):
##!#                        # Print values to text file
##!#                        if(LAT[i,j] > minlat):
##!#                            counter+=1
##!#    #                        fout.write("{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.6f} {9:.6f} {10:.6f} {11:.6f}\n".format(\
##!#    #                            LAT[i,j],LON[i,j],AI[i,j],0.5,SZA[i,j],VZA[i,j],RAZ[i,j], \
##!#    #                            ALBEDO[i,j,0],ALBEDO[i,j,1],REFLECTANCE[i,j,0],\
##!#    #                            REFLECTANCE[i,j,1],CLD[i,j]))
##!#    
##!#    
##!#                        #if((plotXTrack[i,j] == 0) | ((plotXTrack[i,j] & 4 == 4))):
##!#                            index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
##!#                            index2 = int(np.floor(LON[i,j]*4 + 720.))
##!#                            #index1 = int(np.floor(plotLAT[i,j]*4 + 360.))
##!#                            #index2 = int(np.floor(plotLON[i,j]*4 + 720.))
##!#                            
##!#                            if(index1 < 0): index1 = 0
##!#                            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
##!#                            if(index2 < 0): index2 = 0                                                                                            
##!#                            if(index2 > 1439): index2 = 1439
##!#                       
##!#                            #diff = reflectance[i,j] - albedo[i,j]
##!#                            #if(diff<min_diff):
##!#                            #    min_diff = diff
##!#                            #if(diff>max_diff):
##!#                            #    max_diff = diff
##!#                       
##!#                            #if(diff<0.2): 
##!#                            #    UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
##!#                            UVAI_climo[index2, index1] = (UVAI_climo[index2,index1]*count_climo[index2,index1] + \
##!#                                AI[i,j])/(count_climo[index2,index1]+1)
##!#                            #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
##!#                                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
##!#                                #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
##!#                            count_climo[index2, index1] = count_climo[index2,index1] + 1
##!#        data.close() 
##!#
##!#    # Grab the single-pass data
##!#    print("Reading single-swath data")
##!#    for fileI in range(len(total_single_list)):
##!#        # read in data directly from HDF5 files
##!#        print(total_single_list[fileI])
##!#        data2 = h5py.File(total_single_list[fileI],'r')
##!#        AI  = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
##!#        LAT = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
##!#        LON = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
##!#        XTRACK = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
##!#    
##!#        counter = 0
##!#        # Loop over the values and rows 
##!#        for i in range(AI.shape[0]):
##!#            for j in range(0,row_max):
##!#                if(AI[i,j]>-20):
##!#                #if(plotAI[i,j]>-20):
##!#                    # Only plot if XTrack flag is met
##!#                    if((XTRACK[i,j] == 0) | ((XTRACK[i,j] & 4 == 4))):
##!#                        # Print values to text file
##!#                        if(LAT[i,j] > minlat):
##!#                            counter+=1
##!#    
##!#                            index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
##!#                            index2 = int(np.floor(LON[i,j]*4 + 720.))
##!#                            
##!#                            if(index1 < 0): index1 = 0
##!#                            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
##!#                            if(index2 < 0): index2 = 0                                                                                            
##!#                            if(index2 > 1439): index2 = 1439
##!#                       
##!#                            UVAI_single[index2, index1] = (UVAI_single[index2,index1]*count_single[index2,index1] + \
##!#                                AI[i,j])/(count_single[index2,index1]+1)
##!#                            count_single[index2, index1] = count_single[index2,index1] + 1
##!#        data2.close() 
##!#    ##!## Calculate the row-average AI for the secondary plot
##!#    ##!#mask_avgs = np.nanmean(np.ma.masked_where(AI[:,:] < -20, AI[:,:]),axis=0)
##!#
##!#    # Calculate anomalies
##!#    UVAI_anomaly = UVAI_single - UVAI_climo  
##!#
##!#    plot_lat, plot_lon = np.meshgrid(lat_ranges,lon_ranges)
##!#    mask_UVAI_climo = np.ma.masked_where(count_climo[:,:] == 0, UVAI_climo)
##!#    mask_UVAI_anom = np.ma.masked_where(count_single[:,:] == 0, UVAI_anomaly)
##!#
##!#    # Plot the climatology
##!#    plt.close('all')
##!#    fig1 = plt.figure(figsize=(8,8))
##!#    ax = plt.axes(projection = mapcrs)
##!#    ax.gridlines()
##!#    ax.coastlines(resolution='50m')
##!# 
##!#    plt.title('OMI Aerosol Index Climatology '+single_swath)
##!#    #plt.title('OMI Reflectivity - Surface Albedo '+swath_date)
##!#    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI_climo,\
##!#        transform = datacrs,cmap = colormap,\
##!#        vmin = -2.0, vmax = 3.0)
##!#    ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
##!#    ax.set_xlim(-3430748.535086173,3430748.438879491)
##!#    ax.set_ylim(-3413488.8763307533,3443353.899053069)
##!#    #ax.set_xlim(-4170748.535086173,4167222.438879491)
##!#    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
##!#    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),\
##!#        orientation='horizontal',pad=0,\
##!#        aspect=50,shrink = 0.905,label='Aerosol Index')
##!#    ##cax = fig.add_axes([0.16,0.075,0.7,0.025])
##!#    ##cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
##!#    ##cb.ax.set_xlabel('Aerosol Index')
##!#    #cb.ax.set_xlabel('Reflectivity - Surface Albedo')
##!#    #out_name = 'omi_single_pass_ai_200804270052_to_0549_composite_rows_0to'+\
##!#    #    str(row_max)+'.png'       
##!#    out_name = 'omi_single_pass_ai_single_climo_'+single_swath+'_rows_0to'+\
##!#    #    str(row_max)+'.png'       
##!#    #out_name = 'omi_single_pass_refl_albedo_diff_'+swath_date+'_rows_0to'+
##!#    #    str(row_max)+'.png'       
##!#    plt.savefig(out_name)
##!#    print('Saved image '+out_name)
##!#  
##!#    # Plot the anomalies
##!#    plt.close('all')
##!#    fig1 = plt.figure(figsize=(8,8))
##!#    ax = plt.axes(projection = mapcrs)
##!#    ax.gridlines()
##!#    ax.coastlines(resolution='50m')
##!# 
##!#    plt.title('OMI Aerosol Index Anomaly '+single_swath)
##!#    #plt.title('OMI Reflectivity - Surface Albedo '+swath_date)
##!#    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI_anom,\
##!#            transform = datacrs,cmap = colormap,\
##!#            vmin = -1.0, vmax = 1.5)
##!#    ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
##!#    ax.set_xlim(-3430748.535086173,3430748.438879491)
##!#    ax.set_ylim(-3413488.8763307533,3443353.899053069)
##!#    #ax.set_xlim(-4170748.535086173,4167222.438879491)
##!#    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
##!#    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),\
##!#        orientation='horizontal',pad=0,\
##!#        aspect=50,shrink = 0.905,label='Aerosol Index Anomaly')
##!#    ##cax = fig.add_axes([0.16,0.075,0.7,0.025])
##!#    ##cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
##!#    ##cb.ax.set_xlabel('Aerosol Index')
##!#    #cb.ax.set_xlabel('Reflectivity - Surface Albedo')
##!#    #out_name = 'omi_single_pass_ai_200804270052_to_0549_composite_rows_0to'+\
##!#    #    str(row_max)+'.png'       
##!#    out_name = 'omi_single_pass_ai_single_anom_'+single_swath+'_rows_0to'+\
##!#        str(row_max)+'.png'       
##!#    #out_name = 'omi_single_pass_refl_albedo_diff_'+swath_date+'_rows_0to'+\
##!#    #    str(row_max)+'.png'       
##!#    plt.savefig(out_name)
##!#    print('Saved image '+out_name)
##!#    
##!#    #axs[1].plot(mask_avgs)
##!#    #axs[1].set_xlabel('Sensor Row')
##!#    #axs[1].set_ylabel('Row Average Aerosol Index')
##!#    
##!#    plt.show()

def plot_row_bias(minlat = 65., dataset = 'equator', save = False):

    # Open the 2006 CSCI training dataset
    # -----------------------------------
    net_data = Dataset(home_dir + '/CSCI/CSCI_543/final_project/train_dataset_2006.nc','r')
    if(dataset == 'equator'):
        net_data2 = Dataset(home_dir + '/Research/OMI/'+\
            'train_dataset_2006_equator.nc','r')
   
    # Read in OMI data from the two single swaths
    # ------------------------------------------- 
    dates = ['200604221050','200604221408']
    #dates = ['200804221027','200804221206']
    dt_dates = [datetime.strptime(ddate,'%Y%m%d%H%M') for ddate in dates]
    dtype = 'control'
    OMI_base1 = readOMI_swath_hdf(dates[0], dtype, \
        only_sea_ice = False, latmin = minlat - 5)
    OMI_base2 = readOMI_swath_hdf(dates[1], dtype, \
        only_sea_ice = False, latmin = minlat - 5)

    ##!## Open a sample OMI swath
    ##!## -----------------------
    ##!#data = h5py.File('home_dir + /data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2006m0629t0036-o10394_v003-2017m0720t195132.he5','r')

    # Calculate the average of the OMI AI data along the rows
    # for the climatology. The first "axis = 0" averages along
    # all the different swath images, so that the final dimension
    # is 60 x 60, with 60 lines and 60 rows. The second average
    # averages along the lines, so that the final average is 
    # 1 x 60, with one average for each row
    # -----------------------------------------------------------
    avg_AI_climo = np.nanmean(np.nanmean(net_data['AI'], axis = 0), axis = 0)
    if(dataset == 'equator'):
        avg_AI_climo2 = np.nanmean(np.nanmean(net_data2['AI'], axis = 0), axis = 0)
    #std_AI_climo = np.nanstd(np.nanstd(net_data['AI'], axis = 0), axis = 0)

    # Calculate the average azimuth angle value for each row 
    # over the approximate Arctic (1300:1600). Assume the viewing
    # geometry is about the same for all OMI swaths
    # -----------------------------------------------------------
    ##!#AZM = data[\
    ##!#    'HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle'][:,:]
    ##!#AI  = data[\
    ##!#    'HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]

    mask_AZM = np.ma.masked_where(OMI_base1['AZM'] < -9e4, OMI_base1['AZM'])
    #AZM = np.ma.masked_where(AI < -9e4, AZM)
    #AI  = np.ma.masked_where(AI < -9e4, AI)

    ##!#avg_AZM = np.nanmean(AZM[1300:1600,:],\
    ##!#    axis = 0)   
    avg_AZM = np.nanmean(mask_AZM[1300:1600,:],\
        axis = 0)   
    #avg_AI  = np.nanmean(AI[1300:1600,:],\
    #    axis = 0)   

    # Generate figure
    # ---------------
    mapcrs2 = ccrs.NorthPolarStereo()
    if(dataset == 'equator'):
        fig = plt.figure(figsize = (8,4.5))
        ax0 = fig.add_subplot(1,1,1)
    else:
        fig = plt.figure(figsize = (8,8.5))
        gs  = fig.add_gridspec(nrows = 2, ncols = 2, hspace = 0.3)
        ax0 = plt.subplot(gs[0,:])
        ax1 = plt.subplot(gs[1,0], projection = mapcrs2)
        ax2 = plt.subplot(gs[1,1], projection = mapcrs2)

    if(dataset == 'equator'):
        print('Arctic: ', np.min(net_data['LAT']), np.max(net_data['LAT']))
        print('Equator: ', np.min(net_data2['LAT']), np.max(net_data2['LAT']))
        ln1 = ax0.plot(np.arange(1,len(avg_AI_climo)+1), avg_AI_climo, \
            label =   'AI:  2$^{o}$ S - 90$^{o}$ N')
        ln2 = ax0.plot(np.arange(1,len(avg_AI_climo2)+1), avg_AI_climo2, \
            label = 'AI: 65$^{o}$ S - 67$^{o}$ N', color = 'tab:green')
    else:
        ln1 = ax0.plot(np.arange(1,len(avg_AI_climo)+1), avg_AI_climo, label = 'AI')
    #ax.plot(avg_AI_climo + std_AI_climo, label = 'AI', color = 'tab:blue', linestyle = '--')
    #ax.plot(std_AI_climo, label = 'AI', color = 'tab:blue', linestyle = '-.')
    #ax.plot(avg_AI,       label = 'AI', color = 'tab:blue', linestyle = '--')
    ax02 = ax0.twinx()
    ln3 = ax02.plot(np.arange(1, len(avg_AI_climo) + 1), avg_AZM, \
        color = 'tab:orange', label = 'AZM')
    ax0.set_ylabel('Raw OMI AI')
    ax02.set_ylabel('Relative Azimuth Angle [$^{o}$]')
    ax0.set_xlabel('Row number')
    
    ax0.tick_params(axis='y', colors='tab:blue')
    ax02.tick_params(axis='y', colors='tab:orange')

    ax0.yaxis.label.set_color('tab:blue')
    ax02.yaxis.label.set_color('tab:orange')
    ax0.grid()
    ax0.set_title('Average AI and Relative Azimuth Angle\n01 April 2006 - 30 September 2006')

    if(dataset != 'equator'):
        # Plot single-swath data
        # ----------------------
        plotOMI_single_swath(ax1, OMI_base1, title = 'Raw OMI UVAI\n' + \
            dt_dates[0].strftime('%d %B %Y, %H:%M:%S UTC'), circle_bound = True, \
            gridlines = False, vmax = 3.0, colorbar = False)
        #plotOMI_single_swath(ax1, OMI_base2, title = dates[1], \
        plotOMI_single_swath(ax2, OMI_base2, title = 'Raw OMI UVAI\n' + \
            dt_dates[1].strftime('%d %B %Y, %H:%M:%S UTC'), circle_bound = True, \
            gridlines = False, vmax = 3.0)

        plot_subplot_label(ax0, '(a)')
        plot_subplot_label(ax1, '(b)')
        plot_subplot_label(ax2, '(c)')

        lns = ln1 + ln3
    else:
        lns = ln1 + ln2 + ln3

    labs = [l.get_label() for l in lns]
    if(dataset == 'equator'):
        ax0.legend(lns, labs, loc='upper center', bbox_to_anchor=(0.5, -0.20),
              fancybox=True, shadow=False, ncol=3)

    net_data.close()
    if(dataset == 'equator'):
        net_data2.close()
        fig.tight_layout()
    ##!#data.close()

    #fig.subplots_adjust(hspace = 0.2)

    if(save):
        if(dataset == 'equator'):
            outname = 'omi_row_bias_v3.png'
        else:
            outname = 'omi_row_bias_v2.png'
        fig.savefig(outname, dpi = 300)
        print("Saved iamge", outname)
    else:
        plt.show()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Comparison with CERES
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def plot_OMI_CERES_trend_compare(OMI_data, CERES_data,month,ax0 = None, \
        ax1 = None, ax2 = None, minlat=65,\
        trend_type = 'standard', titles = True, region_avgs = True, \
        save=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Reserach/CERES')
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
    ceres_param = CERES_data['param'].split('_')[1]
    if(month != None):
        dt_obj = datetime.strptime(OMI_data['DATES'][month],"%Y%m")
        title = 'OMI AI / CERES ' + CERES_data['param'] + '\n'+ \
            dt_obj.strftime("%b") + " Trend Comparison"
        outname = 'omi_ceres_trend_comp_'+dt_obj.strftime("%b")+'_'+\
            OMI_data['VERSION']+'vCERES_'+ceres_param+'_min'+str(int(minlat))+'.png'
    else:
        title = 'OMI AI / CERES ' + CERES_data['param'] + ' Trend Comparison'
        outname = 'omi_ceres_trend_comp_'+\
            OMI_data1['VERSION']+'vCERES_'+ceres_param+'_min' + str(int(minlat)) + '.png'
 
    label = 'AI Trend (AI/Study Period)'
    omi_title = 'OMI Trend'
    if(not titles):
        omi_title = ' ' 
        title = ' '
        label = ' ' 
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


    if(region_avgs):
        # Grab the matching OMI lon
        # -------------------------
        matching_OMI_lon = OMI_data['LON'][where_matching]

        # Mask the data while retaining the shape
        # ---------------------------------------
        mask_OMI_trend = np.ma.masked_where((OMI_trend == -999.) | \
            (OMI_trend == 0) | (CERES_trend == 0) | (CERES_trend == -999.), \
            OMI_trend)
        mask_CERES_trend = np.ma.masked_where((OMI_trend == -999.) | \
            (OMI_trend == 0) | (CERES_trend == 0) | (CERES_trend == -999.), \
            CERES_trend)

        # Calculate the averages
        # ----------------------
        for region in sea_dict.keys():
            print(region)
            if(region == 'chukchi'):
                omi_region_avg = np.nanmean(mask_OMI_trend[np.where( \
                    (matching_OMI_lon < sea_dict[region]['lon'][0]) | \
                    (matching_OMI_lon >= sea_dict[region]['lon'][2]))])
                ceres_region_avg = np.nanmean(mask_CERES_trend[np.where( \
                    (local_lon < sea_dict[region]['lon'][0]) | \
                    (local_lon >= sea_dict[region]['lon'][2]))])

                print('omi avg   - ',omi_region_avg)
                print('ceres avg - ',ceres_region_avg)
                
            if(region != 'chukchi'):
                omi_region_avg = np.nanmean(mask_OMI_trend[np.where( \
                    (matching_OMI_lon >= sea_dict[region]['lon'][2]) & \
                    (matching_OMI_lon <= sea_dict[region]['lon'][0]))])
                ceres_region_avg = np.nanmean(mask_CERES_trend[np.where( \
                    (local_lon >= sea_dict[region]['lon'][2]) & \
                    (local_lon <= sea_dict[region]['lon'][0]))])

                print('omi avg   - ',omi_region_avg)
                print('ceres avg - ',ceres_region_avg)

        for region in sea_double_dict.keys():
            print(region)
            if(region == 'canada_beaufort_chukchi'):
                omi_region_avg = np.nanmean(mask_OMI_trend[np.where( \
                    (matching_OMI_lon <  sea_double_dict[region]['lon'][0]) | \
                    (matching_OMI_lon >= sea_double_dict[region]['lon'][2]))])
                ceres_region_avg = np.nanmean(mask_CERES_trend[np.where( \
                    (local_lon <  sea_double_dict[region]['lon'][0]) | \
                    (local_lon >= sea_double_dict[region]['lon'][2]))])

                print('omi avg   - ',omi_region_avg)
                print('ceres avg - ',ceres_region_avg)
                
            else:
                omi_region_avg = np.nanmean(mask_OMI_trend[np.where( \
                    (matching_OMI_lon >= sea_double_dict[region]['lon'][2]) & \
                    (matching_OMI_lon <= sea_double_dict[region]['lon'][0]))])
                ceres_region_avg = np.nanmean(mask_CERES_trend[np.where( \
                    (local_lon >= sea_double_dict[region]['lon'][2]) & \
                    (local_lon <= sea_double_dict[region]['lon'][0]))])

                print('omi avg   - ',omi_region_avg)
                print('ceres avg - ',ceres_region_avg)

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
    local_figure = False
    if((ax0 is None) & (ax1 is None) & (ax2 is None)):
        local_figure = True
        plt.close('all')
        fig1 = plt.figure(figsize = (14,5))
        ax0 = fig1.add_subplot(1,3,1, projection = mapcrs)
        ax1 = fig1.add_subplot(1,3,2, projection = mapcrs)
        ax2 = fig1.add_subplot(1,3,3)
        if(region_avgs):
            plot_arctic_regions(ax0, linewidth = 2)
            plot_arctic_regions(ax1, linewidth = 2)

    # Plot the OMI and CERES trends
    # -----------------------------
    plotOMI_MonthTrend(OMI_data,month_idx=month,\
        trend_type=trend_type,label = label,\
        minlat=minlat,title = omi_title, pax = ax0)
    plotCERES_MonthTrend(CERES_data,month_idx=month,\
        trend_type=trend_type,\
        minlat=minlat,pax = ax1)


    ax2.scatter(mask_trend1,mask_trend2,c=z,s=8)
    plot_trend_line(ax2, mask_trend1, mask_trend2, color='tab:green', linestyle = '-', \
        slope = 'theil-sen')
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
    #ax2.legend()
    ax2.set_xlabel(OMI_data['VERSION'])
    ax2.set_ylabel('CERES')
    ax2.set_title(title + '\nMinlat = ' + str(minlat))

    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')
    plot_subplot_label(ax2, '(c)')

    if(local_figure):
        fig1.tight_layout()
        if(save == True):
            fig1.savefig(outname)
            print("Saved image",outname)
        else:
            plt.show()
    #return mask_trend1,mask_trend2

def plot_OMI_CERES_trend_compare_summer(minlat=72,\
        ceres_type = 'sw', trend_type = 'standard', titles = False, \
        region_avgs = True, save=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        print("Appending CERES path")
        sys.path.append(home_dir + '/Research/CERES')

    #importlib.reload(gridCERESLib) 
    from gridCERESLib import readgridCERES

    # Read in the OMI and CERES data
    # ------------------------------
    OMI_data = readOMI_NCDF(infile = \
        home_dir + '/Research/OMI/omi_ai_VJZ211_2005_2020.nc',\
        start_date = 200504, end_date = 202010, minlat = minlat)
    CERES_data =  readgridCERES(200504,202010,'toa_'+ceres_type+\
        '_all_mon', minlat = minlat + 0.5, season = 'sunlight')

    return OMI_data, CERES_data

    # Set up total figure layout
    # --------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (12,10))
    ax1 = fig1.add_subplot(3,3,1, projection = mapcrs)
    ax2 = fig1.add_subplot(3,3,4, projection = mapcrs)
    ax3 = fig1.add_subplot(3,3,7, projection = mapcrs)
    ax4 = fig1.add_subplot(3,3,2, projection = mapcrs)
    ax5 = fig1.add_subplot(3,3,5, projection = mapcrs)
    ax6 = fig1.add_subplot(3,3,8, projection = mapcrs)
    ax7 = fig1.add_subplot(3,3,3)
    ax8 = fig1.add_subplot(3,3,6)
    ax9 = fig1.add_subplot(3,3,9)

    if(region_avgs):
        plot_arctic_regions(ax1, linewidth = 2)
        plot_arctic_regions(ax2, linewidth = 2)
        plot_arctic_regions(ax3, linewidth = 2)
        plot_arctic_regions(ax4, linewidth = 2)
        plot_arctic_regions(ax5, linewidth = 2)
        plot_arctic_regions(ax6, linewidth = 2)

    # Process data for each month
    # ---------------------------
    plot_OMI_CERES_trend_compare(OMI_data, CERES_data,2,ax0 = ax1, \
        ax1 = ax4, ax2 = ax7, minlat=minlat,\
        trend_type = trend_type, titles = titles, save=False)
    plot_OMI_CERES_trend_compare(OMI_data, CERES_data,3,ax0 = ax2, \
        ax1 = ax5, ax2 = ax8, minlat=minlat,\
        trend_type = trend_type, titles = titles, save=False)
    plot_OMI_CERES_trend_compare(OMI_data, CERES_data,4,ax0 = ax3, \
        ax1 = ax6, ax2 = ax9, minlat=minlat,\
        trend_type = trend_type, titles = titles, save=False)
    
    fig1.tight_layout()
    plt.show()

# Compare the OMI, CERES, and MODIS data over the Arctic
# ------------------------------------------------------
def plot_compare_OMI_CERES_MODIS_NSIDC(modis_date_str, ch1, \
        omi_dtype = 'shawn', minlat = 65., zoom = False, save = False):

    if(home_dir + '/Research/CERES/' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES/')
    if(home_dir + '/Research/MODIS/obs_smoke_forcing/' not in sys.path):
        sys.path.append(home_dir + '/Research/MODIS/obs_smoke_forcing/')
    if(home_dir + '/Research/NSIDC/' not in sys.path):
        sys.path.append(home_dir + '/Research/NSIDC/')
    from gridCERESLib import plotCERES_hrly_figure
    from MODISLib import plot_MODIS_channel
    from NSIDCLib import plotNSIDC_daily_figure

    file_date_dict = {
        '201807052125': {
            'OMI': '201807052034', 
            'CERES': '2018070521', 
        },
        '201808241435': {
            'OMI': '201808241343', 
            'CERES': '2018082414', 
        },
        '201908110125': {
            'OMI': '201908110033', 
            'CERES': '2019081101', 
        },
        '201908110440': {
            'OMI': '201908110351', 
            'CERES': '2019081104', 
        },
    }

    # Make the overall figure
    # -----------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (14,5))
    ax1 = fig1.add_subplot(2,3,1, projection = mapcrs)
    ax2 = fig1.add_subplot(2,3,2, projection = mapcrs)
    ax3 = fig1.add_subplot(2,3,3, projection = mapcrs)
    ax4 = fig1.add_subplot(2,3,4, projection = mapcrs)
    ax5 = fig1.add_subplot(2,3,5, projection = mapcrs)
    ax6 = fig1.add_subplot(2,3,6, projection = mapcrs)
   
    # Plot the MODIS true-color and channel data
    # ------------------------------------------
    plot_MODIS_channel(modis_date_str, 'true_color', swath = True, \
        zoom = zoom, ax = ax1)
    plot_MODIS_channel(modis_date_str, ch1, swath = True, \
        zoom = zoom, ax = ax2, vmax = 0.4)
    #plot_MODIS_channel(modis_date_str, ch2, swath = True, \
    #    zoom = zoom, ax = ax3)

    # Plot the NSIDC data
    # -------------------
    plotNSIDC_daily_figure(modis_date_str[:8], minlat = minlat, \
        zoom = zoom, ax = ax3, gridlines = False, save = False)

    # Plot the OMI data
    # -----------------
    plotOMI_single_swath_figure(file_date_dict[modis_date_str]['OMI'], \
            dtype = omi_dtype, only_sea_ice = False, minlat = minlat, \
            ax = ax4, skiprows = [52], lat_circles = None, save = False, \
            zoom = zoom)
    
    # Plot the CERES data
    # -------------------
    plotCERES_hrly_figure(file_date_dict[modis_date_str]['CERES'], 'SWF',  \
        minlat = minlat, lat_circles = None, ax = ax5, title = 'SWF',\
        grid_data = True, zoom = zoom, vmax = 450, vmin = None, save = False)
    plotCERES_hrly_figure(file_date_dict[modis_date_str]['CERES'], 'LWF',  \
        minlat = minlat, lat_circles = None, ax = ax6, title = 'LWF',\
        grid_data = True, zoom = zoom, vmax = 275, vmin = 150, save = False)

    if(save):
        outname = 'omi_ceres_modis_nsidc_compare_' + modis_date_str + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# OMI fort out plotting stuff 
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def read_OMI_fort_out(file_name, min_lat, vtype = 'areas', max_lat = None, \
        all_data = False):
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
    # Extract date information from the file name
    name_split = file_name.strip().split('/')[-1].split('_')
    dtype      = name_split[1]
    vtype      = name_split[2]  # counts or areas?
    start_year = name_split[3]
    end_year   = name_split[4]
    thresh     = name_split[5].split('.')[0]
    ai_thresh  = float(int(thresh)/100)
        
    if(vtype == 'areas'):
        data_str = 'Area'
        #divider = 1e5
        #axis_label = 'Area of high AI [10$^{5}$ km$^{2}$]'
    elif(vtype == 'counts'):
        data_str = 'Cnt'
        divider = 1
        axis_label = 'Count of quarter-degree lat/lon grids of high AI'
    else:
        print("INVALID VARIABLE TYPE. Must be areas or counts")
        return
    
    if(int(min_lat) == 75):
        interval = 1
        h_val = 1
    elif(int(min_lat) > 75):
        interval = 0.2
        h_val = 1
        #h_val = 0.1
        divider = 1e4
        axis_label = 'Area of high AI [10$^{4}$ km$^{2}$]'
    else:
        interval = 2
        h_val = 1
        divider = 1e5
        axis_label = 'Area of high AI [10$^{5}$ km$^{2}$]'
    
    print(ai_thresh)
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Read area data from the input file
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    in_data = pd.read_csv(file_name, delim_whitespace=True)
    dates  = in_data['Date'].values
    if(len(str(dates[0])) == 10):
        time_fmt = "%Y%m%d%H"
        day_avgs = False
    elif(len(str(dates[0])) == 8):
        time_fmt = "%Y%m%d"
        day_avgs = True
    else:
        print("INVALID DATE FORMAT")
        return

    dt_dates = [datetime.strptime(str(tmpdate),time_fmt) \
        for tmpdate in dates]
    
    # Convert the areas from square kilometers to 1e5 square kilometers
    count65  = in_data[data_str + str(int(min_lat))].values / divider
    x_range = np.arange(len(dates))
   
    # If the user only wants the areas between two lat bins, do so here
    if(max_lat is not None):
        print("Reading max lat data now")
        count_max  = in_data[data_str + str(int(max_lat))].values / divider
        count65 = count65 - count_max
 
    # Calculate daily totals
    #daily_counts_65 = [np.average(tmparr) for tmparr in \
    if(not day_avgs):
        daily_counts_65 = [np.sum(tmparr) for tmparr in \
            np.array_split(count65,len(count65)/4)]
        daily_dt_dates = dt_dates[::4]
        daily_dates = dates[::4]/100
        daily_xrange = x_range[::4]
    else:
        # Just use the daily totals from the file
        daily_counts_65 = count65
        daily_dates = dates
        daily_dt_dates = dt_dates
    
    # Add up the areas for each year
    daily_dt_dates = np.array(daily_dt_dates)
    daily_counts_65 = np.array(daily_counts_65)
    
    years = np.arange(int(start_year),int(end_year) + 1)
    ##!#yearly_totals = np.zeros(len(years))
    
    ##!## Create a 2d array to hold all the area data from each day
    ##!#day_values = np.full((int(datetime(year=2008,month=9,day=30).strftime('%j')) \
    ##!#    - int(datetime(year=2008,month=4,day=1).strftime('%j'))+1, len(years)+1),-9.)
    ##!#
    ##!## Insert each day's area data into the array
    ##!#for xx in range(len(daily_dt_dates)):
    ##!#    indices = int(daily_dt_dates[xx].strftime('%j')) - \
    ##!#        int(datetime(year=daily_dt_dates[xx].year,month=4,\
    ##!#        day=1).strftime('%j')), daily_dt_dates[xx].year - \
    ##!#        daily_dt_dates[0].year 
    ##!#    day_values[indices] = daily_counts_65[xx]
    
    ##!## Mask any -9s, which are caused by days without data
    ##!#mask_values = np.ma.masked_where(day_values == -9,day_values)
    # Calculate the mean and standard deviation of the areas for each
    # day of the year.
    ##!#mask_day_avgs = np.nanmean(mask_values,axis=1)
    ##!#mask_day_stds = np.nanstd(mask_values,axis=1)
    ##!#
    ##!#mean_total_std = np.nanstd(mask_values)
    ##!##mean_std = np.nanmean(mask_day_stds)
    ##!#
    ##!## Calculate the
    ##!## This adds and subtracts the standard deviation for each day
    ##!## ---------------------------------------------------------------------
    ##!#lower_range = mask_day_avgs - mask_day_stds * 0.75
    ##!#upper_range = mask_day_avgs + mask_day_stds * 0.75
    ##!#
    ##!## This adds and subtracts the total standard deviation of all the areas
    ##!## ---------------------------------------------------------------------
    ##!##lower_range = mask_day_avgs - mean_total_std * 0.75
    ##!##upper_range = mask_day_avgs + mean_total_std * 0.75
    ##!#
    ##!## Use a 1-sigma check to mask any daily_counts_65 values that are
    ##!## outside of the average
    ##!#new_daily_counts = np.full((len(daily_dt_dates)),-9.)
    ##!#for xx in range(len(daily_dt_dates)):
    ##!#    indices = int(daily_dt_dates[xx].strftime('%j')) - \
    ##!#        int(datetime(year=daily_dt_dates[xx].year,month=4,\
    ##!#        day=1).strftime('%j')), daily_dt_dates[xx].year - \
    ##!#        daily_dt_dates[0].year 
    ##!#    if((daily_counts_65[xx] - mask_day_avgs[indices[0]]) > 1. * \
    ##!#    #if((daily_counts_65[xx] - mask_day_avgs[indices[0]]) < 0.5 * \
    ##!#        mask_day_stds[indices[0]]):
    ##!#        #daily_counts_65[xx] = -9
    ##!#        new_daily_counts[xx] = daily_counts_65[xx]
    
    daily_counts_65 = np.ma.masked_where(daily_counts_65 == -9,daily_counts_65)
    
    # Create an output dictionary to hold the data
    # --------------------------------------------
    omi_fort_dict = {}
    omi_fort_dict['daily_data'] = daily_counts_65
    omi_fort_dict['daily_dt_dates'] = daily_dt_dates
    omi_fort_dict['years'] = years
    omi_fort_dict['interval'] = interval
    omi_fort_dict['h_val'] = h_val
    omi_fort_dict['axis_label'] = axis_label
    omi_fort_dict['ai_thresh'] = ai_thresh
    omi_fort_dict['dtype']      = dtype
    omi_fort_dict['vtype']      = vtype  # counts or areas?
    omi_fort_dict['start_year'] = start_year
    omi_fort_dict['end_year']   = end_year

    return omi_fort_dict   

def calc_aerosol_event_dates(infile, minlat = 70., vtype = 'areas', \
        max_lat = None):

    # Read the daily data
    # -------------------
    omi_fort_dict = read_OMI_fort_out(infile, minlat, vtype = vtype, \
        max_lat = max_lat)

    outname = 'omi_area_event_dates_minlat' + str(int(minlat)) + '.txt'
   
    with open(outname, 'w') as fout: 
        # Loop over the dates and figure out which ones have area over the
        # threshold
        # ----------------------------------------------------------------
        for ii in range(len(omi_fort_dict['daily_data'])):
            local_date = omi_fort_dict['daily_dt_dates'][ii]
            local_area = omi_fort_dict['daily_data'][ii]
    
            if(local_area > omi_fort_dict['h_val']):
                fout.write("{0} {1:.3f}\n".format(\
                    local_date.strftime('%Y%m%d'), local_area))

                print(local_date.strftime('%Y%m%d'), omi_fort_dict['daily_data'][ii])


def plot_OMI_fort_out_func(infile, ax = None, min_lat = 70., vtype = 'areas', \
        max_lat = None, save = False):

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Read the data from the file
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    omi_fort_dict = read_OMI_fort_out(infile, min_lat, vtype = vtype)

    if(max_lat is not None):
        print("Reading max lat data now")
        omi_fort_dict_max = read_OMI_fort_out(infile, max_lat, vtype = vtype)

        if(max_lat > 75):
            # Handle the case where the very small data north of 80 has been
            # multiplied by 10.
            omi_fort_dict_max['daily_data'] = omi_fort_dict_max['daily_data'] / 10.      
 
        omi_fort_dict['daily_data'] = omi_fort_dict['daily_data'] - \
            omi_fort_dict_max['daily_data'] 

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Plot the individual omi_fort_dict['years'] of area time series
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
  
    in_ax = True 
    if(ax is None): 
        in_ax = False
        fig0 = plt.figure(figsize = (10,4))
        ax = fig0.add_subplot(1,1,1)
    plot_c = cm.turbo((omi_fort_dict['years']-np.min(omi_fort_dict['years']))/\
        (np.max(omi_fort_dict['years'])-np.min(omi_fort_dict['years'])))
    
    # Set up an array to hold the event counts
    # ----------------------------------------
    d_val = 4
    #omi_fort_dict['omi_fort_dict['h_val']'] = 1
    if(omi_fort_dict['interval'] == 0.2):
        max_val = np.round(np.max(omi_fort_dict['daily_data']), 1)+omi_fort_dict['interval']
    else:
        if(int(min_lat) == 75):
            max_val = int(np.max(omi_fort_dict['daily_data']))+2*omi_fort_dict['interval']
        else:
            max_val = int(np.max(omi_fort_dict['daily_data']))+omi_fort_dict['interval']
    
    event_sizes = np.arange(omi_fort_dict['h_val'], max_val , omi_fort_dict['interval'])
    events_yearly = np.zeros((omi_fort_dict['years'].shape[0], event_sizes.shape[0]))
    
    for ii, year in enumerate(omi_fort_dict['years']):
        ##!#yearly_totals[ii] = np.sum(omi_fort_dict['daily_data'][np.where( \
        ##!#    (omi_fort_dict['daily_dt_dates'] >= datetime(year,4,1)) & \
        ##!#    (omi_fort_dict['daily_dt_dates'] <= datetime(year,9,30)))])
    
        ##!## Test the local maxima
        ##!## ---------------------
        ##!#for jj in range(len(event_sizes) - 1):
        ##!#    # Find the peaks for this size
        ##!#    # ----------------------------
        ##!#    peaks, _ = find_peaks(omi_fort_dict['daily_data'][np.where( \
        ##!#        (omi_fort_dict['daily_dt_dates'] >= datetime(year,4,1)) & \
        ##!#        (omi_fort_dict['daily_dt_dates'] <= datetime(year,9,30)))], \
        ##!#        height = [event_sizes[jj], event_sizes[jj+1]], distance = d_val)
        ##!#    events_yearly[ii,jj] = len(peaks)        
     
        # Determine the counts of peaks greater than each height range
        # ------------------------------------------------------------
    
        # Calculate the days since April 1
        # --------------------------------
        loop_dates = omi_fort_dict['daily_dt_dates'][np.where( \
            (omi_fort_dict['daily_dt_dates'] >= datetime(year,4,1)) & \
            (omi_fort_dict['daily_dt_dates'] <= datetime(year,9,30)))]
    
        plot_dates = datetime(year = 2000, month = 4, day = 1) + \
            (loop_dates - datetime(year = year, month = 4, day = 1))
   
        ax.plot(plot_dates,\
            omi_fort_dict['daily_data'][np.where( \
            (omi_fort_dict['daily_dt_dates'] >= datetime(year,4,1)) & \
            (omi_fort_dict['daily_dt_dates'] <= datetime(year,9,30)))], c = plot_c[ii])
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    norm = mc.Normalize(vmin=np.min(omi_fort_dict['years']), \
        vmax = np.max(omi_fort_dict['years']))
   
    bounds = np.arange(np.min(omi_fort_dict['years']) - 0.5, \
        np.max(omi_fort_dict['years']) + 1.5) 
    ticks = np.arange(np.min(omi_fort_dict['years']), \
        np.max(omi_fort_dict['years']) + 1) 

    cmap = cm.turbo
    cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                 ax=ax, orientation='vertical', ticks = ticks[::3], \
                 boundaries = bounds, label='Year')
    
    #ax.plot(plot_dates, mask_day_avgs,label='avg',color='black')
    #ax.plot(plot_dates, upper_range,label='+avg σ',linestyle='--',color='black')
    ax.set_ylabel(omi_fort_dict['axis_label'])
    if(max_lat is not None):
        ax.set_title('High AI ' + vtype + '\nPerturbed AI threshold of '+\
            str(omi_fort_dict['ai_thresh']) + '\n' + str(int(min_lat)) + \
            '$^{o}$N - ' + str(int(max_lat))+'$^{o}$N')
    else: 
        ax.set_title('High AI ' + vtype + '\nPerturbed AI threshold of ' + \
            str(omi_fort_dict['ai_thresh']) + '\nNorth of ' + \
            str(int(min_lat)) + '$^{o}$N')
    ax.tick_params(axis='both', labelsize = 12)
    ax.grid()
    
    if(not in_ax):
        fig0.tight_layout()
        if(save == True):
            outname = "ai_" + vtype + "_indiv_yearly_" + dtype + "_thresh"+thresh+"_minlat"+min_lat+".png"
            fig0.savefig(outname,dpi=300)
            print("Saved image",outname) 
        else:
            plt.show() 
    #plt.legend()
    
    ##!## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    ##!##
    ##!## Plot the total yearly areas on a graph
    ##!##
    ##!## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    ##!#fig2 = plt.figure()
    ##!#ax4 = fig2.add_subplot(1,1,1)
    ##!#ax4.plot(omi_fort_dict['years'],yearly_totals)
    ##!#
    ##!#if(save == True):
    ##!#    outname = "ai_" + vtype + "_total_" + dtype + "_thresh"+thresh+"_minlat"+min_lat+".png"
    ##!#    fig2.savefig(outname,dpi=300)
    ##!#    print("Saved image",outname) 


# Plot the time series of OMI peaks with 'x's at peaks for the OMI fort out
# dictionary returned from read_OMI_fort_out
# ----------------------------------------------------------------------------
def plot_OMI_fort_out_time_series(omi_fort_dict, pax, minlat = 70., d_val = 4):
   
    if(np.max(omi_fort_dict['daily_data']) < 1.):
        multiplier = 10
        exponent = 4
    else:
        multiplier = 1
        exponent = 5

    label_str = 'Area of high AI [10$^{' + \
        str(exponent) + '}$ km$^{2}$]'

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Plot the total time series of daily areas with peaks
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    pax.plot(omi_fort_dict['daily_dt_dates'], \
        omi_fort_dict['daily_data'] * multiplier, \
        label='daily '+omi_fort_dict['dtype'])
 
    # Test peaks here
    # ---------------
    d_val = 4
    peaks, _ = find_peaks(omi_fort_dict['daily_data'] * multiplier, \
        height = omi_fort_dict['h_val'] * multiplier, distance = d_val)
    pax.plot(omi_fort_dict['daily_dt_dates'][peaks], \
        omi_fort_dict['daily_data'][peaks] * multiplier, 'x', color='black')
    
    #ax1.legend()
    #pax.set_ylabel(omi_fort_dict['axis_label'])
    pax.set_ylabel(label_str, weight = 'bold')
    #pax.set_title('Latitude > '+str(int(minlat))+'$^{o}$')
    #pax.set_title('AI ' + omi_fort_dict['vtype'].title() + ': Threshold of '+str(omi_fort_dict['ai_thresh'])+\
    #    '\nNorth of '+str(int(minlat))+'$^{o}$')
    pax.grid()


# Plot the bar chart of yearly peaks in each height range for the OMI fort out
# dictionary returned from read_OMI_fort_out
# ----------------------------------------------------------------------------
def plot_OMI_fort_out_peak_bar(omi_fort_dict, pax, minlat = 70., d_val = 4):

    if(omi_fort_dict['interval'] == 0.2):
        max_val = np.round(np.max(omi_fort_dict['daily_data']), 1)+omi_fort_dict['interval']
    else:
        if(int(minlat) == 75):
            max_val = int(np.max(omi_fort_dict['daily_data']))+2*omi_fort_dict['interval']
        else:
            max_val = int(np.max(omi_fort_dict['daily_data']))+omi_fort_dict['interval']
    event_sizes = np.arange(omi_fort_dict['h_val'], max_val , omi_fort_dict['interval'])
    events_yearly = np.zeros((omi_fort_dict['years'].shape[0], event_sizes.shape[0]))
    
    for ii, year in enumerate(omi_fort_dict['years']):
    
        # Test the local maxima
        # ---------------------
        for jj in range(len(event_sizes) - 1):
            # Find the peaks for this size
            # ----------------------------
            peaks, _ = find_peaks(omi_fort_dict['daily_data'][np.where( \
                (omi_fort_dict['daily_dt_dates'] >= datetime(year,4,1)) & \
                (omi_fort_dict['daily_dt_dates'] <= datetime(year,9,30)))], \
                height = [event_sizes[jj], event_sizes[jj+1]], distance = d_val)
            events_yearly[ii,jj] = len(peaks)        

    for ii in range(len(event_sizes)-1):
        if(event_sizes[ii] < 1):
            multiplier = 10
            exponent = 4
        else:
            multiplier = 1
            exponent = 5

        label_str = str(int(np.round(event_sizes[ii] * multiplier,1))) + ' - ' + \
            str(int(np.round(event_sizes[ii+1] * multiplier,1))) + ' * 10$^{' + \
            str(exponent) + '}$ km$^{2}$'
        pax.bar(omi_fort_dict['years'], events_yearly[:,ii], \
            bottom = np.sum(events_yearly[:,0:ii],axis=1), \
            label = label_str)
    pax.legend(fontsize = 10)
    pax.set_ylabel('# of AI peaks in each size range', weight = 'bold')
    #pax.set_title('Latitude > '+str(int(minlat))+'$^{o}$')
    #pax.set_title('AI ' + omi_fort_dict['vtype'].title() + ' : Threshold of '+str(omi_fort_dict['ai_thresh'])+\
    #    '\nNorth of '+str(int(minlat))+'$^{o}$')
    

# Plot a combination of both the total time series and the bar chart
# ------------------------------------------------------------------
def plot_OMI_fort_out_peaks(infile, minlat = 70., vtype = 'areas', save = False):

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Read the data from the file
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    omi_fort_dict = read_OMI_fort_out(infile, minlat, vtype = vtype)
    
    fig1 = plt.figure(figsize = (13,5)) 
    ax1 = fig1.add_subplot(1,2,1)
    ax2 = fig1.add_subplot(1,2,2)

    # Plot the time series with the peaks as 'x's
    # -------------------------------------------
    plot_OMI_fort_out_time_series(omi_fort_dict, ax1, minlat = minlat)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Plot the yearly event size box graph
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    plot_OMI_fort_out_peak_bar(omi_fort_dict, ax2, minlat = minlat)

    if(save):
        outname = 'ai_daily_areas_peaks_thresh'+thresh+"_minlat"+min_lat+'.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image",outname)
    
    else:    
        plt.show()
    
# Plot a combination of both the total time series and the bar chart
# ------------------------------------------------------------------
def plot_OMI_fort_out_two_lats(infile, minlat = 70., vtype = 'areas', save = False):

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Read the data from the file
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    omi_fort_dict_70 = read_OMI_fort_out(infile, 70., vtype = vtype, max_lat = 80.)
    omi_fort_dict_80 = read_OMI_fort_out(infile, 80., vtype = vtype)
    
    fig1 = plt.figure(figsize = (11,10)) 
    ax1 = fig1.add_subplot(2,2,1)
    ax2 = fig1.add_subplot(2,2,2)
    ax3 = fig1.add_subplot(2,2,3)
    ax4 = fig1.add_subplot(2,2,4)

    # Plot the time series with the peaks as 'x's
    # -------------------------------------------
    plot_OMI_fort_out_time_series(omi_fort_dict_70, ax1, minlat = 70.) 
    plot_OMI_fort_out_peak_bar(omi_fort_dict_70, ax2, minlat = 70.)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Plot the yearly event size box graph
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    plot_OMI_fort_out_time_series(omi_fort_dict_80, ax3, minlat = 80.)
    plot_OMI_fort_out_peak_bar(omi_fort_dict_80, ax4, minlat = 80.)

    # Add subplot labels
    # ------------------
    plot_subplot_label(ax1, '(a)', location = 'upper_left')
    plot_subplot_label(ax2, '(b)', location = 'upper_left')
    plot_subplot_label(ax3, '(c)', location = 'upper_left')
    plot_subplot_label(ax4, '(d)', location = 'upper_left')

    row_label_size = 14
    fig1.text(0.05, 0.705, '---------- 70$^{o}$ N - 80$^{o}$ N ----------', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size)
    fig1.text(0.05, 0.285, '--------------- > 80$^{o}$ N ----------------', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size)

    plt.suptitle('AI ' + vtype + ': Threshold of '+str(omi_fort_dict_70['ai_thresh']), \
        fontsize = 14, weight = 'bold')

    if(save):
        thresh = str(int(omi_fort_dict_70['ai_thresh']*100))
        outname = 'ai_fort_out_thresh'+thresh+'_minlat7080_split.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image",outname)
    
    else:    
        plt.show()

# Plot a three-panel figure. Panel 1 is the yearly time series of the OMI AI 
# areas. Panel 2 is the WorldView image. Panel 3 is the single-day average
# ----------------------------------------------------------------------------
def plot_combined_fort_out(plot_time, min_lat = 70., vtype = 'areas', \
        max_lat = None, save = False):

    if(len(plot_time) == 8):
        fmt_str = '%Y%m%d'
    elif(len(plot_time) == 10):
        fmt_str = '%Y%m%d%H'
    else:
        print("ERROR: invalid date/time format")
        return

    dt_date_str = datetime.strptime(plot_time, fmt_str)

    if(plot_time == '20190811'):
        fort_minlat = 70.
    elif(plot_time == '20170818'):
        fort_minlat = 80.

    # Set up the figure
    # -----------------
    fig = plt.figure(figsize = (9,8))
    ##!#mosaic = \
    ##!#    """
    ##!#    AA
    ##!#    BC
    ##!#    """
    ##!#axs = fig.subplot_mosaic(mosaic)

    #gs = fig.add_gridspec(nrows = 2, ncols = 2)
    gs = gridspec.GridSpec(nrows=2, ncols=2)
    mapcrs = ccrs.NorthPolarStereo(central_longitude = 320.)
    ax0 = plt.subplot(gs[0,:])
    ax1 = plt.subplot(gs[1,0])
    ax2 = plt.subplot(gs[1,1],projection=mapcrs)

    # Read the Shawn data the old way
    # ------------------------------
    OMI_shawn = readOMI_swath_shawn_old(plot_time, latmin = 65., resolution = 0.25)

    print("OMI_shawn lat shape", OMI_shawn['LAT'].shape)

    plot_OMI_shawn_old(OMI_shawn, ax = ax2, minlat = 65., labelsize = 11, \
        labelticksize = 10)

    # Read in the matching Worldview image
    # ------------------------------------
    base_path = home_dir + '/Research/OMI/'
    img_name = dt_date_str.strftime(base_path + 'snapshot-%Y-%m-%dT00_00_00Z.jpg') 
    print("Looking for ",img_name)

    pic = plt.imread(img_name)
    ax1.imshow(pic)
    ax1.axis('off')
    #ax1.set_title('MODIS Aqua True-Color\n' + dt_date_str.strftime('%d %B %Y'))
    plot_figure_text(ax1, 'MODIS Aqua True Color\n' + dt_date_str.strftime('%d %B %Y'), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 11, backgroundcolor = 'white', halign = 'right')
 
    # Read (and plot?) the yearly area time series
    # --------------------------------------------
    infile = home_dir + '/Research/OMI/shawn_analysis/count_analysis/' + \
        'omi_vsj4_areas_2005_2020_100.txt'
    plot_OMI_fort_out_func(infile, ax = ax0, min_lat = fort_minlat, vtype = 'areas', \
        max_lat = max_lat, save = False)

    # Add subplot labels
    # ------------------
    plot_subplot_label(ax0, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax1, '(b)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(c)')

    #fig.tight_layout()

    if(save):
        if(max_lat is not None):
            outname = 'omi_fort_out_combined_' + plot_time + '_max_lat.png'
        else:
            outname = 'omi_fort_out_combined_' + plot_time + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# OMI bad row anomaly table stuff
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def plot_bad_row_table(bad_row_file, xtrack_file = None, ax = None, \
        save = False):

    do_xtrack = False
    bad_color = 'jet'
    xtrack_add = ''
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
    # Extract date information from the file name
    name_split = bad_row_file.strip().split('/')[-1].split('_')
    start_year = name_split[-2]
    end_year   = name_split[-1].split('.')[0]
    
    #bad_dict = {}
    
    str_dates = []
    
    # Read bad row file
    # -----------------
    with open(bad_row_file,'r') as f_in:
        flines = f_in.readlines()
        total_array = np.zeros((len(flines),61))
        #int_dates = np.zeros(len(flines))
        for ii, line in enumerate(flines):
        #for line in f_in:
            templine = line.strip().split()
            #int_dates[ii] = int(templine[0])
            str_dates.append(templine[0][:4])
            #int_dates[ii] = int(templine[0])
            #bad_dict[templine[0]] = np.zeros(int(templine[1]))
            if(int(templine[1]) > 0):
                for jj in templine[2:]:
                    total_array[ii,int(jj)] = 1
                #bad_dict[templine[0]][ii] = int(templine[2+ii])
   
    do_xtrack = False 
    if(xtrack_file is not None):
        bad_color = 'bwr_r'
        do_xtrack = True
        str_dates_X = []
        
        # Read xtrack row file
        # -----------------
        with open(xtrack_file,'r') as f_in:
            flines = f_in.readlines()
            total_array_X = np.zeros((len(flines),61))
            #int_dates = np.zeros(len(flines))
            for ii, line in enumerate(flines):
            #for line in f_in:
                templine = line.strip().split()
                #int_dates[ii] = int(templine[0])
                str_dates_X.append(templine[0][:4])
                #int_dates[ii] = int(templine[0])
                #bad_dict[templine[0]] = np.zeros(int(templine[1]))
                if((int(templine[1]) > 0) & (int(templine[1]) < 60)):
                    for jj in templine[2:]:
                        total_array_X[ii,int(jj)] = 1
                    #bad_dict[templine[0]][ii] = int(templine[2+ii])
        
        x_vals_X = np.arange(total_array_X.shape[0])
        y_vals_X = np.arange(total_array_X.shape[1])
        mask_array_X = np.ma.masked_where(total_array_X == 0,total_array_X)
    
    
    #index = 0
    #for key in bad_dict:
    #    total_array[index,:(len(bad_dict[key]))] = bad_dict[key]
    #    index += 1
    
    x_vals = np.arange(total_array.shape[0])
    y_vals = np.arange(total_array.shape[1])
    mask_array = np.ma.masked_where(total_array == 0,total_array)

    in_ax = True
    if(ax is None):
        in_ax = False
        fig1, ax = plt.subplots(figsize=(10,4))
   
    ax.pcolormesh(x_vals,y_vals + 0.5,mask_array.T,shading='auto',cmap=bad_color)
    if(do_xtrack):
        xtrack_color = 'jet'
        ax.pcolormesh(x_vals_X,y_vals_X + 0.5,mask_array_X.T,shading='auto',cmap=xtrack_color)
    ax.set_yticks(y_vals,minor=True)
    ax.set_xticks(x_vals[::183] + 0.5)
    ax.set_xticklabels(str_dates[::183])
    ax.grid(which='minor',axis='y')
    ax.grid(which='major',axis='y',linewidth=1.0,color='black')
    ax.set_xlabel('Year')
    ax.set_ylabel('Row Number')
    ax.set_ylim(1, 60)

    # Add legend
    # ----------
    handles, labels = ax.get_legend_handles_labels()
    f_patch  = mpatches.Patch(color = 'darkblue', label = 'Flagged')
    uf_patch = mpatches.Patch(color = 'red', label = 'Unflagged')
    handles.extend([f_patch, uf_patch])
    ax.legend(handles = handles, loc = 'lower right', framealpha = 1)
   
    if(in_ax is False): 
        if(save):
            outname = 'bad_rows_'+start_year + '_'+end_year+xtrack_add + '_2.png'
            plt.savefig(outname,dpi=300)
            print("Saved image",outname)
        plt.show()

def plot_row_anomaly_combined(dtype = 'control', minlat = 65., save = False):
   
    date_str  = '201204102151'
    date_str2 = '201204102330'
    #date_str2 = '201807260244'
 
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    dt_date_str2 = datetime.strptime(date_str2, '%Y%m%d%H%M')

    # Set up the figure
    # -----------------
    fig = plt.figure(figsize = (9,11))
    #fig = plt.figure(figsize = (16,4))
    ##!#ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
    ##!#ax2 = fig.add_subplot(1,2,2)
    gs = fig.add_gridspec(nrows = 3, ncols = 2)
    ax1  = fig.add_subplot(gs[0,0], projection = mapcrs)   # true color    
    ax2  = fig.add_subplot(gs[0,1], projection = mapcrs)   # true color    
    ax4  = fig.add_subplot(gs[1,0]) # MODIS Ch 31
    ax5  = fig.add_subplot(gs[1,1], projection = mapcrs) # MODIS Ch 31
    ax3  = fig.add_subplot(gs[2,:]) # MODIS Ch 31

    # Read the single-swath OMI data
    # ------------------------------
    if(dtype == 'shawn'):
        OMI_base  = readOMI_swath_shawn(date_str, latmin = minlat,\
            skiprows = skiprows)
        OMI_base2  = readOMI_swath_shawn(date_str2, latmin = minlat,\
            skiprows = skiprows)
        OMI_base2_skip  = readOMI_swath_shawn(date_str2, latmin = minlat,\
            skiprows = [42, 43])
    else:
        OMI_base  = readOMI_swath_hdf(date_str, dtype, \
            only_sea_ice = False, latmin = minlat, \
            skiprows = None)
        OMI_base2  = readOMI_swath_hdf(date_str2, dtype, \
            only_sea_ice = False, latmin = minlat, \
            skiprows = None)
        OMI_base2_skip  = readOMI_swath_hdf(date_str2, dtype, \
            only_sea_ice = False, latmin = minlat, \
            skiprows = [42,43])
    
    # Plot the single-swath OMI data
    # ------------------------------
    plotOMI_single_swath(ax1, OMI_base, title = dtype.title(), \
        circle_bound = True, gridlines = False, colorbar = False)
    plotOMI_single_swath(ax2, OMI_base2, title = dtype.title(), \
        circle_bound = True, gridlines = False)
    plotOMI_single_swath(ax5, OMI_base2_skip, title = dtype.title(), \
        circle_bound = True, gridlines = False)

    ax1.set_extent([-180,180,minlat,90], datacrs)
    ax1.set_title(dt_date_str.strftime('Raw OMI UVAI\n%d %b %Y, %H:%M:%S UTC'))
    ax2.set_extent([-180,180,minlat,90], datacrs)
    ax2.set_title(dt_date_str2.strftime('Raw OMI UVAI\n%d %b %Y, %H:%M:%S UTC'))
    ax5.set_extent([-180,180,minlat,90], datacrs)
    ax5.set_title(dt_date_str2.strftime('Raw OMI UVAI\nUnflagged Bad Rows Removed\n%d %b %Y, %H:%M:%S UTC'))

    # Plot the row averages
    # ---------------------
    ax4.bar(np.arange(OMI_base2['UVAI'].shape[1]), np.average(OMI_base2['UVAI'][1230:1500,:], axis = 0))
    #ax4.bar(np.arange(OMI_base['UVAI'].shape[1]), np.average(OMI_base['UVAI'][1230:1500,:], axis = 0))
    ax4.axhline(0.0, color = 'black')
    ax4.set_xlim(0, 60)
    ax4.set_xlabel('Row Number', fontsize = 12)
    ax4.set_ylabel('Average AI over Arctic', fontsize = 12)
    ax4.set_title(dt_date_str2.strftime('Raw OMI UVAI\n%d %b %Y, %H:%M:%S UTC'))

    ##!#ax5.bar(np.arange(OMI_base2['UVAI'].shape[1]), np.average(OMI_base2['UVAI'][1230:1500,:], axis = 0))
    ##!#ax5.axhline(0.0, color = 'black')
    ##!#ax5.set_xlim(0, 60)
    ##!#ax5.set_xlabel('Row Number')
    ##!#ax5.set_ylabel('Average AI over Arctic')
    ##!#ax5.set_title(dt_date_str2.strftime('%H:%M UTC %d %b %Y'))

    # Plot the row anomaly table
    # -------------------------- 
    bad_row_file = 'row_anomaly_dates_20050401_20201001.txt'
    xtrack_file = 'row_anomaly_xtrack_dates_20050401_20201001.txt'
    plot_bad_row_table(bad_row_file, xtrack_file = xtrack_file, ax = ax3)
    ax3.set_title('Temporal Row Anomaly Flagging')

    plot_subplot_label(ax1, '(a)')
    plot_subplot_label(ax2, '(b)')
    plot_subplot_label(ax4, '(c)')
    plot_subplot_label(ax5, '(d)')
    plot_subplot_label(ax3, '(e)', backgroundcolor = 'white')
    
    fig.tight_layout()

    if(save):
        outname = 'combined_bad_row_'+date_str+'_stacked_v3.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_OMI_row_avg(date_str, plot_swath = False, minlat = 65., \
        save = True):
    
    if(len(date_str) == 4):
        str_fmt = '%Y'
        out_fmt = '%Ym'
    elif(len(date_str) == 6):
        str_fmt = '%Y%m'
        out_fmt = '%Ym%m'
    elif(len(date_str) == 8):
        str_fmt = '%Y%m%d'
        out_fmt = '%Ym%m%d'
    else:
        print("INVALID DATE STRING. RETURNING")
        sys.exit()
    
    dt_date_str = datetime.strptime(date_str, str_fmt)
    
    files = glob(dt_date_str.strftime(\
        '/home/bsorenson/data/OMI/H5_files/OMI-A*_' + out_fmt + '*.he5'))
    
    num_files = len(files)
    
    swath_ais = np.full((num_files, 2000, 60), np.nan)
    swath_lat = np.full((num_files, 2000, 60), np.nan)
    swath_xtk = np.full((num_files, 2000, 60, 6), np.nan)
    swath_xrw = np.full((num_files, 2000, 60), np.nan)
    
    latmin = 65
    for ii, fname in enumerate(files):
        dt_name = datetime.strptime(fname.strip().split('/')[-1][20:34],\
            '%Ym%m%dt%H%M')

        # Read and resample each file
        # ---------------------------
        OMI_data = readOMI_swath_hdf(dt_name.strftime('%Y%m%d%H%M'), 'control', \
            latmin = latmin, \
            skiprows = [52])
   
        if(plot_swath): 
            # Plot this swath
            # ---------------
            plotOMI_single_swath_figure(OMI_data, dtype = 'control', skiprows = [52], \
                vmin = -2, vmax = 4, minlat = minlat, save = True)
    
        swath_ais[ii, :OMI_data['UVAI'].shape[0], :]      = OMI_data['UVAI'][:,:]
        swath_lat[ii, :OMI_data['LAT'].shape[0], :]       = OMI_data['LAT'][:,:]
        swath_xtk[ii, :OMI_data['XTRACK'].shape[0], :,:]  = OMI_data['XTRACK'][:,:,:]
        swath_xrw[ii, :OMI_data['XTRACK'].shape[0], :]    = OMI_data['XTRACK_raw'][:,:]
    
    # Clean the data
    swath_ais = np.ma.masked_where(swath_ais < -50, swath_ais)
    swath_ais = np.ma.masked_invalid(swath_ais)
    swath_ais = np.ma.masked_where(swath_xrw != 0, swath_ais)
    swath_ais = np.ma.masked_where(swath_lat < latmin, swath_ais)
    
    # Average all the data for each day along each row, so there
    # is one set of row averages per day
    day_avgs = np.nanmean(swath_ais, axis = 1)
    # Average all the daily row averages
    row_avgs = np.nanmean(day_avgs, axis = 0)
   
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(row_avgs)
    ax.set_xlabel('Row #')
    ax.set_ylabel('Daily average AI')
    ax.set_title('OMI daily average UVAI by row\nNorth of 65\n' + date_str)
    ax.grid()
    fig.tight_layout()
    
    outname = 'omi_row_avgs_' + date_str + '.png'
    fig.savefig(outname, dpi = 300)
    print("Saved image", outname)
    
