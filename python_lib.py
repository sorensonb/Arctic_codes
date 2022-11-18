"""
  NAME:
    python_lib.py

  PURPOSE:
    Contain general useful functions across different projects

"""

import numpy as np
import matplotlib.path as mpath
from scipy import stats
from datetime import datetime
import cartopy.crs as ccrs
from glob import glob
from PIL import Image
import h5py
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
from math import sin, cos, sqrt, atan2, radians
from bs4 import BeautifulSoup
import requests
import json

# LAADS DAAC key
# --------------
laads_daac_key = 'Ymxha2Uuc29yZW5zb246WW14aGEyVXVjMjl5Wlc1emIyNUFkVzVr'+\
    'TG1Wa2RRPT06MTY2MzYwNjAyMDo1NWU4MjMwNjI2ZTVkZDE2MmU3MmQ3M2Q3OGQxN'+\
    'Dg3Y2NjZTU4ZTAx'

# Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

datacrs = ccrs.PlateCarree()
#mapcrs = ccrs.LambertConformal(central_longitude = -110, central_latitude = 40)

def format_coord(x, y):
    xarr = X[0,:]
    yarr = Y[:,0]
    if ((x > xarr.min()) & (x <= xarr.max()) & 
        (y > yarr.min()) & (y <= yarr.max())):
        col = np.searchsorted(xarr, x)-1
        row = np.searchsorted(yarr, y)-1
        z = Z[row, col]
        return f'x={x:1.4f}, y={y:1.4f}, z={z:1.4f}   [{row},{col}]'
    else:
        return f'x={x:1.4f}, y={y:1.4f}'

def find_distance_between_points(lat1, lon1, lat2, lon2):
    R = 6371. # km
    
    lat1 = radians(lat1) 
    lon1 = radians(lon1) 
    lat2 = radians(lat2) 
    lon2 = radians(lon2) 

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = sin(dlat / 2)**2. + cos(lat1) * cos(lat2) * sin(dlon / 2)**2.
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    distance = R * c

    return distance

# Find the gridpoint in the gridded lat/lon data that 
# corresponds to the station at slat and slon
# ---------------------------------------------------- 
def nearest_gridpoint(slat, slon, grid_lat, grid_lon):
    fun_c = ((slat - grid_lat)**2. + (slon - grid_lon)**2.)**0.5
    #fun_c = np.maximum(np.abs(grid_lat - slat), \
    #    np.abs(grid_lon - slon))
    m_idx = np.where(fun_c == np.min(fun_c))

    # NOTE: If broken when running with code other than VIIRS data,
    # comment this out
    # -------------------------------------------------------------
    if(len(m_idx[0]) > 1):
        print("ERROR:", m_idx, m_idx[0][0])
        m_idx = m_idx[0][0]
    
    return m_idx

def covariance(x,y):
    avg_x = np.average(x)
    avg_y = np.average(y)
    N = len(x)
    if(len(x)!=len(y)):
        print("ERROR: Arrays are not the same size.\nArray x has len=",len(x),\
              "\nArray y has len=",len(y))
    else:
        cov = (np.sum((x-avg_x)*(y-avg_y)))/(N-1)
        return cov

def correlation(x,y):
    avg_x = np.average(x)
    avg_y = np.average(y)
    N = len(x)
    if(len(x)!=len(y)):
        print("ERROR: Arrays are not the same size.\nArray x has len=",len(x),\
              "\nArray y has len=",len(y))
    else:
        cov = (np.sum((x-avg_x)*(y-avg_y)))/(N-1)
        std_x = np.std(x)
        std_y = np.std(y)
        corr = cov/(std_x*std_y)
        return(corr)

# Wavenumber assumed to be in cm-1
# Wavelength returned in μm 
def convert_wavenumber_to_wavelength(wavenum):
    return ( (1. / (wavenum * 1e2)) * 1e6)

# Wavelength assumed to be in μm 
# Wavenumber returned in cm-1
def convert_wavelength_to_wavenumber(wavelength):
    return (1. / (wavelength * 1e-6 )) * 1e-2

def convert_temp_to_radiance(lmbda, temp):

    # Define constants for converting radiances to temperatures
    lmbda = (1e-6) * lmbda # convert to meters
    c_const = 3e8
    h_const = 6.626e-34 # J*s
    k_const = 1.381e-23 # J/K

    data = (2. * h_const * c_const**2.) / \
        (lmbda** 5. * (np.exp( (h_const * c_const ) / \
        (lmbda * k_const * temp )) - 1.) * 1e6)

    #data = (h_const * c_const) / \
    #    (lmbda * k_const * np.log( ((2.0 * h_const * (c_const**2.0) ) / \
    #    ((lmbda**4.) * (lmbda / 1e-6) * radiance ) ) + 1.0 ) )

    return data

def convert_radiance_to_temp(lmbda, radiance):

    # Define constants for converting radiances to temperatures
    lmbda = (1e-6) * lmbda
    c_const = 3e8
    h_const = 6.626e-34 # J*s
    k_const = 1.381e-23 # J/K

    data = (h_const * c_const) / \
        (lmbda * k_const * np.log( ((2.0 * h_const * (c_const**2.0) ) / \
        ((lmbda**4.) * (lmbda / 1e-6) * radiance ) ) + 1.0 ) )

    return data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# identify_axes automatically ads plot labels on subplots generated using
# plt.subplot_mosaic. For example, here is how the figure is to be set up
#
# >>> mosaic = \
# ...   """
# ...   AAB
# ...   AAC
# ...   """
# ...
# >>> fig = plt.figure()
# >>> axs = fig.subplot_mosaic(mosaic)
#
# Then, to add labels to all the figures...
# 
# >>> identify_axes(ax)
#

def identify_axes(ax_dict, fontsize = 14, color = 'k', \
        backgroundcolor = None):
    if(backgroundcolor is None):
        kw = dict(ha='left', va = 'top', fontsize = fontsize, color = color)
    else:
        kw = dict(ha='left', va = 'top', fontsize = fontsize, color = color, \
            backgroundcolor = backgroundcolor)

    # Check if ax_dict is actually a dict. If not, then try to
    # extract the dictionary from the tuple
    # --------------------------------------------------------
    if(not isinstance(ax_dict, dict)):
        ax_dict = ax_dict[1]
    
    for k, ax in ax_dict.items():
        ax.text(0.08, 0.92, k, transform = ax.transAxes, **kw)
    
def init_proj(date_str):
    #mapcrs = Miller()
    if(date_str == None):
        mapcrs = ccrs.LambertConformal()
    else:
        dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

        mapcrs = ccrs.LambertConformal(central_longitude = \
            np.mean(aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']),\
            central_latitude = \
            np.mean(aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']))

    return mapcrs

def plot_theil_sen_trend(pax, xdata, ydata, color, linestyle, linewidth):
    res = stats.theilslopes(ydata, xdata, 0.90)
    print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))

    # Then, plot the trend line on the figure
    pax.plot(xdata, res[1] + res[0] * xdata, \
        color='k', linewidth = linewidth + 2.0, linestyle = linestyle)
    # Then, plot the trend line on the figure
    pax.plot(xdata, res[1] + res[0] * xdata, \
        color=color, linestyle = linestyle)
   
    return res[0], res[1]
 
def plot_lin_regress_trend(pax, xdata, ydata, color, linestyle, linewidth):
    # First, calculate the trend
    zdata = np.polyfit(xdata, ydata, 1)

    print("Lin Regress: {0}x + {1}".format(*zdata))

    # Then, plot the trend line on the figure
    pax.plot(np.unique(xdata), np.poly1d(zdata)(np.unique(xdata)), \
        color='k', linewidth = linewidth + 2.0, linestyle = linestyle)
    pax.plot(np.unique(xdata), np.poly1d(zdata)(np.unique(xdata)), \
        color=color, linestyle = linestyle)

    return zdata

def plot_trend_line(pax, xdata, ydata, color='black', linestyle = '-', \
        linewidth = 1.5, slope = 'theil-sen'):

    if(slope == 'theil-sen'):
        out = plot_theil_sen_trend(pax, xdata, ydata, color, linestyle, linewidth)
    elif(slope == 'both'):
        out = plot_theil_sen_trend(pax, xdata, ydata, color, linestyle, linewidth)
        out = plot_lin_regress_trend(pax, xdata, ydata, color, linestyle, linewidth)
    else:
        out = plot_lin_regress_trend(pax, xdata, ydata, color, linestyle, linewidth)

    return out

##!#def plot_trend_line(pax, xdata, ydata, color='black', linestyle = '-', \
##!#        slope = 'thiel-sen'):
##!#
##!#    if(slope == 'thiel-sen'):
##!#        res = stats.theilslopes(ydata, xdata, 0.95)
##!#        print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))
##!#
##!#        # Then, plot the trend line on the figure
##!#        pax.plot(xdata, res[1] + res[0] * xdata, \
##!#            color='k', linewidth = 2.5, linestyle = linestyle)
##!#        # Then, plot the trend line on the figure
##!#        pax.plot(xdata, res[1] + res[0] * xdata, \
##!#            color=color, linestyle = linestyle)
##!#    else:
##!#        # First, calculate the trend
##!#        zdata = np.polyfit(xdata, ydata, 1)
##!#
##!#        print("{0}x + {1}".format(*zdata))
##!#
##!#        # Then, plot the trend line on the figure
##!#        pax.plot(np.unique(xdata), np.poly1d(zdata)(np.unique(xdata)), \
##!#            color=color, linestyle = linestyle)

def plot_subplot_label(ax, label, xval = None, yval = None, transform = None, \
        color = 'black', backgroundcolor = None, fontsize = 14, \
        location = 'upper_left'):

    if(location == 'upper_left'):
        y_lim = 0.90
        #y_lim = 1.03
        x_lim = 0.05
    if(location == 'upper_upper_left'):
        y_lim = 1.03
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

def plot_figure_text(ax, text, xval = None, yval = None, transform = None, \
        color = 'black', fontsize = 12, backgroundcolor = 'white',\
        halign = 'left', location = 'lower_right'):

    if(location == 'lower_right'):
        if(xval is None):
            xval = ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.95
        if(yval is None):
            yval = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.05
    elif(location == 'upper_right'):
        if(xval is None):
            xval = ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.95
        if(yval is None):
            yval = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.90

    print('Xval = ',xval, 'Yval = ',yval)

    if(transform is None):
        ax.text(xval,yval,text, \
            color=color, weight='bold', \
            fontsize=fontsize, backgroundcolor = backgroundcolor, \
            horizontalalignment = halign)
    else:
        ax.text(xval,yval,text, \
            color=color, weight='bold', \
            transform = transform, fontsize=fontsize, \
            backgroundcolor = backgroundcolor, \
            horizontalalignment = halign)

# Plots latitude circles on a given axis
# --------------------------------------
def plot_lat_circles(pax, lat_circles):

    if(len(lat_circles) > 5):
        print("WARNING: more than 5 latitude circles selected")
        print("    Only plotting the first 5")
        lat_circles = lat_circles[:5]

    colors = ['tab:blue','tab:red','tab:purple','tab:olive','tab:cyan']
    lon = np.arange(-180, 181)
    for ii, lat_val in enumerate(lat_circles):
        lats = np.full(lon.shape, lat_val)

        pax.plot(lon, lats, linewidth = 2.5, transform = datacrs, color = 'k')
        pax.plot(lon, lats, transform = datacrs, color = colors[ii])
        pax.text(-180, lat_val + 3, str(int(lat_val)) + ' $^{o}$N', \
            horizontalalignment = 'center', weight = 'bold', \
            color = colors[ii], transform = datacrs)

def make_gif(frame_folder, gif_name):
    frames = [Image.open(image) for image in glob(f"{frame_folder}/*.png")]
    frame_one = frames[0]
    frame_one.save(gif_name, format = 'GIF', append_images = frames,\
        save_all = True, duration = 200, loop = 0)
    print("Saved gif",gif_name)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Downloading functions 
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# listFD grabs all files and directories that end with the provided ext at the
# provided url.
# ----------------------------------------------------------------------------
def listFD(url, ext=''):
    page = requests.get(url).text
    #print page
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') \
        if node.get('href').endswith(ext)]

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Event dictionary  
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
with open('/home/bsorenson/Research/Arctic_compares/json_comp_files.txt') as fin:
    aerosol_event_dict = json.load(fin)
##!#aerosol_event_dict = {
##!#    "2008-04-22": {
##!#        'Lat': [65., 80.0],
##!#        'Lon': [147., 228.],
##!#        '1935': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2008113.1935.061.2018033064331.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2008m0422t1841-o20060_v003-2017m0721t115822.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc',
##!#            'ceres_time': '19',
##!#            'swath': ['200804221935','200804221940', '200804221945'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [65., 90.0],
##!#            #'Lon': [-180., 180.],
##!#            'Lat': [65., 80.0],
##!#            'Lon': [147., 228.],
##!#            'modis_Lat': [65., 80.0],
##!#            'modis_Lon': [147., 228.],
##!#        },
##!#        '1940': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2008113.1940.061.2018033064402.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2008m0422t1841-o20060_v003-2017m0721t115822.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc',
##!#            'ceres_time': '19',
##!#            'swath': ['200804221935','200804221940', '200804221945'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 80.0],
##!#            'Lon': [147., 228.],
##!#            'modis_Lat': [65., 80.0],
##!#            'modis_Lon': [147., 228.],
##!#        },
##!#        '1945': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2008113.1945.061.2018033064241.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2008m0422t1841-o20060_v003-2017m0721t115822.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc',
##!#            'ceres_time': '19',
##!#            'swath': ['200804221935','200804221940', '200804221945'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 80.0],
##!#            'Lon': [147., 228.],
##!#            'modis_Lat': [65., 80.0],
##!#            'modis_Lon': [147., 228.],
##!#        },
##!#        '2110': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2008113.2110.061.2018033064303.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2008m0422t2020-o20061_v003-2017m0721t115818.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc',
##!#            'ceres_time': '21',
##!#            'swath': ['200804222110','200804222115', '200804222120'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 80.0],
##!#            'Lon': [147., 228.],
##!#            'modis_Lat': [65., 80.0],
##!#            'modis_Lon': [147., 228.],
##!#        },
##!#        '2115': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2008113.2115.061.2018033064300.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2008m0422t2020-o20061_v003-2017m0721t115818.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc',
##!#            'ceres_time': '21',
##!#            'swath': ['200804222110','200804222115', '200804222120'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 80.0],
##!#            'Lon': [147., 228.],
##!#            'modis_Lat': [65., 80.0],
##!#            'modis_Lon': [147., 228.],
##!#        },
##!#        '2120': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2008113.2120.061.2018033064306.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2008m0422t2020-o20061_v003-2017m0721t115818.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc',
##!#            'ceres_time': '21',
##!#            'swath': ['200804222110','200804222115', '200804222120'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 80.0],
##!#            'Lon': [147., 228.],
##!#            'modis_Lat': [65., 80.0],
##!#            'modis_Lon': [147., 228.],
##!#        },
##!#        '2250': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2008113.2250.061.2018033064306.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2008m0422t2159-o20062_v003-2017m0721t115805.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc',
##!#            'ceres_time': '22',
##!#            'swath': ['200804222250','200804222255', '200804222300'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 80.0],
##!#            'Lon': [147., 228.],
##!#            'modis_Lat': [65., 80.0],
##!#            'modis_Lon': [147., 228.],
##!#        },
##!#        '2255': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2008113.2255.061.2018033064354.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2008m0422t2159-o20062_v003-2017m0721t115805.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc',
##!#            'ceres_time': '22',
##!#            'swath': ['200804222250','200804222255', '200804222300'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 80.0],
##!#            'Lon': [147., 228.],
##!#            'modis_Lat': [65., 80.0],
##!#            'modis_Lon': [147., 228.],
##!#        },
##!#        '2300': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2008113.2300.061.2018033064357.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2008m0422t2159-o20062_v003-2017m0721t115805.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc',
##!#            'ceres_time': '22',
##!#            'swath': ['200804222250','200804222255', '200804222300'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 80.0],
##!#            'Lon': [147., 228.],
##!#            'modis_Lat': [65., 80.0],
##!#            'modis_Lon': [147., 228.],
##!#        },
##!#    },
##!#    "2016-05-15": {
##!#        'Lat': [65., 90.0],
##!#        'Lon': [-180., 180.],
##!#        '1435': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018236.1435.061.2018237153916.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2016m0515t2104-o62948_v003-2017m0724t143803.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2016051500-2016051523.nc',
##!#            'ceres_time': '22',
##!#            'swath': ['201808241435','201808241440', '201808241445'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [65., 90.0],
##!#            'Lon': [-180., 180.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '1440': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018236.1440.061.2018237153927.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2016m0515t2104-o62948_v003-2017m0724t143803.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2016051500-2016051523.nc',
##!#            'ceres_time': '22',
##!#            'swath': ['201808241435','201808241440', '201808241445'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 90.0],
##!#            'Lon': [-180., 180.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '1445': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018236.1445.061.2018237153918.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2016m0515t2104-o62948_v003-2017m0724t143803.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2016051500-2016051523.nc',
##!#            'ceres_time': '22',
##!#            'swath': ['201808241435','201808241440', '201808241445'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 90.0],
##!#            'Lon': [-180., 180.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#    },
##!#    "2017-08-17": {
##!#        '0000': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2017231.1810.061.2018034125202.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2017m0817t2043-o69632_v003-2017m0818t223219.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2017081700-2017081723.nc',
##!#            'ceres_time': '21',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [80.0, 90.0],
##!#            'Lon': [-125.0, -75.0],
##!#            'modis_Lat': [80.0, 90.0],
##!#            'modis_Lon': [-125.0, -75.0],
##!#        },
##!#    },
##!#    "2017-08-18": {
##!#        '0000': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2017231.1810.061.2018034125202.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2017m0818t1630-o69644_v003-2017m0819t230655.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2017081800-2017081823.nc',
##!#            'ceres_time': '17',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [80.0, 90.0],
##!#            'Lon': [-95.0, -15.0],
##!#            'modis_Lat': [80.0, 90.0],
##!#            'modis_Lon': [-95.0, -15.0],
##!#        },
##!#    },
##!#    "2017-08-19": {
##!#        '1810': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2017231.1810.061.2018034125202.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2017m0819t1713-o69659_v003-2017m0820t234933.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/minlat_55/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2017081900-2017081923.nc',
##!#            'ceres_time': '18',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [80.0, 88.0],
##!#            'Lon': [-100.0, -20.0],
##!#            'modis_Lat': [82.0, 88.0],
##!#            'modis_Lon': [-105.0, -25.0],
##!#        },
##!#    },
##!#    "2018-07-05": {
##!#        'Lat': [64., 80.0],
##!#        'Lon': [145., 218.],
##!#        '1950': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018186.1950.061.2018187153425.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0705t1856-o74320_v003-2019m0802t151352.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018070500-2018070523.nc',
##!#            'ceres_time': '19',
##!#            'swath': ['201807051950','201807051955', '201807052000'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [64., 80.0],
##!#            'Lon': [145., 218.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '1955': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018186.1955.061.2018187153053.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0705t1856-o74320_v003-2019m0802t151352.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018070500-2018070523.nc',
##!#            'ceres_time': '19',
##!#            'swath': ['201807051950','201807051955', '201807052000'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [64., 80.0],
##!#            'Lon': [145., 218.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '2000': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018186.2000.061.2018187153028.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0705t1856-o74320_v003-2019m0802t151352.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018070500-2018070523.nc',
##!#            'ceres_time': '19',
##!#            'swath': ['201807051950','201807051955', '201807052000'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [64., 80.0],
##!#            'Lon': [145., 218.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '2125': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018186.2125.061.2018187153719.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0705t2034-o74321_v003-2019m0802t151255.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018070500-2018070523.nc',
##!#            'ceres_time': '21',
##!#            'swath': ['201807052125','201807052130', '201807052135'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [64., 80.0],
##!#            'Lon': [145., 218.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '2130': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018186.2130.061.2018187153814.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0705t2034-o74321_v003-2019m0802t151255.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018070500-2018070523.nc',
##!#            'ceres_time': '21',
##!#            'swath': ['201807052125','201807052130', '201807052135'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [64., 80.0],
##!#            'Lon': [145., 218.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '2135': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018186.2135.061.2018187153628.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0705t2034-o74321_v003-2019m0802t151255.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018070500-2018070523.nc',
##!#            'ceres_time': '21',
##!#            'swath': ['201807052125','201807052130', '201807052135'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [64., 80.0],
##!#            'Lon': [145., 218.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '2305': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018186.2305.061.2018187153333.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0705t2213-o74322_v003-2019m0802t151345.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018070500-2018070523.nc',
##!#            'ceres_time': '23',
##!#            'swath': ['201807052305','201807052310', '201807052315'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [64., 80.0],
##!#            'Lon': [145., 218.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '2310': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018186.2310.061.2018187153449.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0705t2213-o74322_v003-2019m0802t151345.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018070500-2018070523.nc',
##!#            'ceres_time': '23',
##!#            'swath': ['201807052305','201807052310', '201807052315'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [64., 80.0],
##!#            'Lon': [145., 218.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '2315': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018186.2315.061.2018187153325.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0705t2213-o74322_v003-2019m0802t151345.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018070500-2018070523.nc',
##!#            'ceres_time': '23',
##!#            'swath': ['201807052305','201807052310', '201807052315'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [64., 80.0],
##!#            'Lon': [145., 218.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#    },
##!#    "2018-08-24": {
##!#        'Lat': [65., 90.0],
##!#        'Lon': [-180., 180.],
##!#        '1435': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018236.1435.061.2018237153916.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0824t1343-o75045_v003-2019m0802t155420.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018082400-2018082423.nc',
##!#            'ceres_time': '14',
##!#            'swath': ['201808241435','201808241440', '201808241445'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [65., 90.0],
##!#            'Lon': [-180., 180.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '1440': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018236.1440.061.2018237153927.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0824t1343-o75045_v003-2019m0802t155420.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018082400-2018082423.nc',
##!#            'ceres_time': '14',
##!#            'swath': ['201808241435','201808241440', '201808241445'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 90.0],
##!#            'Lon': [-180., 180.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#        '1445': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2018236.1445.061.2018237153918.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2018m0824t1343-o75045_v003-2019m0802t155420.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2018082400-2018082423.nc',
##!#            'ceres_time': '14',
##!#            'swath': ['201808241435','201808241440', '201808241445'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [60.0, 90.0],
##!#            #'Lon': [-180.0, 180.0],
##!#            'Lat': [65., 90.0],
##!#            'Lon': [-180., 180.],
##!#            'modis_Lat': [60.0, 90.0],
##!#            'modis_Lon': [-180.0, 180.0],
##!#        },
##!#    },
##!#    "2019-08-11": {
##!#        'Lat': [70.0, 85.0],
##!#        'Lon': [105.0, 151.0],
##!#        '0125': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0125.061.2019223152835.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2019m0811t0033-o80163_v003-2019m0812t234149.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2019081100-2019081123.nc',
##!#            'ceres_time': '01',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [68.0, 82.0],
##!#            #'Lon': [95.0, 161.0],
##!#            'swath': ['201908110125','201908110130', '201908110135'],
##!#            'Lat': [70.0, 85.0],
##!#            'Lon': [105.0, 151.0],
##!#            'modis_Lat': [68.0, 82.0],
##!#            'modis_Lon': [95.0, 161.0],
##!#        },
##!#        '0130': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0130.061.2019223152859.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2019m0811t0033-o80163_v003-2019m0812t234149.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2019081100-2019081123.nc',
##!#            'ceres_time': '01',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [68.0, 82.0],
##!#            #'Lon': [95.0, 161.0],
##!#            'swath': ['201908110125','201908110130', '201908110135'],
##!#            'Lat': [70.0, 85.0],
##!#            'Lon': [105.0, 151.0],
##!#            'modis_Lat': [68.0, 82.0],
##!#            'modis_Lon': [95.0, 161.0],
##!#        },
##!#        '0135': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0135.061.2019223152904.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2019m0811t0033-o80163_v003-2019m0812t234149.he5',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2019081100-2019081123.nc',
##!#            'ceres_time': '01',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [68.0, 82.0],
##!#            #'Lon': [95.0, 161.0],
##!#            'swath': ['201908110125','201908110130', '201908110135'],
##!#            'Lat': [70.0, 85.0],
##!#            'Lon': [105.0, 151.0],
##!#            'modis_Lat': [68.0, 82.0],
##!#            'modis_Lon': [95.0, 161.0],
##!#        },
##!#        '0310': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': ['/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0305.061.2019223152906.hdf',\
##!#                      '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0310.061.2019223152911.hdf'], 
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2019m0811t0212-o80164_v003-2019m0812t234246.he5',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2019081100-2019081123.nc',
##!#            'ceres_time': '03',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [68.0, 82.0],
##!#            #'Lon': [95.0, 161.0],
##!#            'Lat': [70.0, 85.0],
##!#            'Lon': [105.0, 151.0],
##!#            'modis_Lat': [68.0, 82.0],
##!#            'modis_Lon': [95.0, 161.0],
##!#        },
##!#        '0440': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0440.061.2019223152242.hdf', 
##!#            ##!#'modis': ['/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0130.061.2019223152859.hdf', \
##!##           ##!#           '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0305.061.2019223152906.hdf', \
##!##           ##!#           '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0310.061.2019223152911.hdf', \
##!#            ##!#          '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0445.061.2019223152919.hdf'], 
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2019m0811t0351-o80165_v003-2019m0812t234332.he5',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2019081100-2019081123.nc',
##!#            'ceres_time': '04',
##!#            'swath': ['201908110440','201908110445', '201908110450'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [74.0, 80.0],
##!#            #'Lon': [120.0, 151.0],
##!#            'Lat': [70.0, 85.0],
##!#            'Lon': [105.0, 151.0],
##!#            'modis_Lat': [74.0, 80.0],
##!#            'modis_Lon': [110.0, 161.0],
##!#        },
##!#        '0445': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0445.061.2019223152919.hdf', 
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2019m0811t0351-o80165_v003-2019m0812t234332.he5',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2019081100-2019081123.nc',
##!#            'ceres_time': '04',
##!#            'swath': ['201908110440','201908110445', '201908110450'],
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            #'Lat': [74.0, 80.0],
##!#            #'Lon': [120.0, 151.0],
##!#            'Lat': [70.0, 85.0],
##!#            'Lon': [105.0, 151.0],
##!#            'modis_Lat': [74.0, 80.0],
##!#            'modis_Lon': [110.0, 161.0],
##!#        },
##!#        '0450': {
##!#            #'asos': 'asos_data_20210713.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0450.061.2019223152915.hdf', 
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2019m0811t0351-o80165_v003-2019m0812t234332.he5',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2019081100-2019081123.nc',
##!#            'ceres_time': '04',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'swath': ['201908110440','201908110445', '201908110450'],
##!#            'Lat': [74.0, 80.0],
##!#            'Lon': [120.0, 151.0],
##!#            'modis_Lat': [74.0, 80.0],
##!#            'modis_Lon': [110.0, 161.0],
##!#        },
##!#    },
##!#    "2021-07-13": {
##!#        '2110': {
##!#            'asos': 'asos_data_20210713.csv',
##!#            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [39.5, 42.0],
##!#            'Lon': [-122.0, -119.5],
##!#            'modis_Lat': [39.0, 42.5],
##!#            'modis_Lon': [-123., -119.]
##!#        }
##!#    },
##!#    "2021-07-20": {
##!#        '2125': {
##!#            'asos': 'asos_data_20210720.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021201.2125.061.2021202154814.hdf',
##!#            'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021201.2125.061.2021202225717.hdf',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/modis_comp/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072010-2021072021.nc',
##!#            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.20.214.L2.SUBS2RET.v6.0.32.0.G21202153435.hdf'],
##!#            'Lat': [39.5, 42.0],
##!#            'Lon': [-122.0, -119.5],
##!#            'data_lim': {
##!#                1:  [None, None],
##!#                2:  [None, None],
##!#                3:  [None, None],
##!#                4:  [None, None],
##!#                5:  [None, None],
##!#                7:  [None, None],
##!#                31: [270., 330.],
##!#                32: [270., 330.],
##!#                'wv_ir': [0.2, 2.5],
##!#            },
##!#            'modis_Lat': [39.5, 42.0],
##!#            'modis_Lon': [-122.0, -119.5]
##!#        }
##!#    },
##!#    "2021-07-21": {
##!#        '2030': {
##!#            'asos': 'asos_data_20210722_4.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021202.2030.061.2021203174050.hdf',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072109-2021072122.nc',
##!#            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.21.205.L2.SUBS2RET.v6.0.32.0.G21203185004.hdf'],
##!#            #'Lat': [39.5, 42.0],
##!#            #'Lon': [-122.0, -119.5],
##!#            'Lat': [39.5, 42.0],
##!#            'Lon': [-122.0, -119.5],
##!#            'data_lim': {
##!#                1:  [None, None],
##!#                2:  [None, None],
##!#                3:  [None, None],
##!#                4:  [None, None],
##!#                5:  [None, None],
##!#                6:  [None, None],
##!#                7:  [None, None],
##!#                31: [270., 330.],
##!#                32: [270., 330.],
##!#                'wv_ir': [0.2, 1.5],
##!#            },
##!#            'modis_Lat': [39.5, 42.0],
##!#            'modis_Lon': [-122.0, -119.5],
##!#            'goes_Lat': [39.5, 42.0],
##!#            'goes_Lon': [-122.0, -119.5],
##!#        }
##!#    },
##!#    "2021-07-22": {
##!#        '2110': {
##!#            'asos': 'asos_data_20210722_4.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
##!#            'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
##!#            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [39.5, 42.0],
##!#            'Lon': [-122.3, -119.5],
##!#            #'Lat': [39.5, 42.0],
##!#            #'Lon': [-122.0, -119.5],
##!#            'data_lim': {
##!#                1:  [5., 50.],
##!#                2:  [None, None],
##!#                3:  [None, None],
##!#                4:  [None, None],
##!#                5:  [None, None],
##!#                6:  [None, None],
##!#                7:  [None, None],
##!#                31: [270., 330.],
##!#                32: [270., 330.],
##!#                'wv_ir': [0.2, 1.5],
##!#                'SWF': [130., 250.],
##!#                'LWF': [300., 370.],
##!#            },
##!#            #'modis_Lat': [39.5, 42.0],
##!#            #'modis_Lon': [-122.0, -119.5]
##!#            'modis_Lat': [39.5, 42.0],
##!#            'modis_Lon': [-122.0, -119.5], 
##!#            #'goes_Lat': [39.5, 42.0],
##!#            #'goes_Lon': [-122.0, -119.5]
##!#            'goes_Lat': [39.5, 42.0],
##!#            'goes_Lon': [-122.0, -119.5],
##!#        }
##!#    },
##!#    "2021-07-23": {
##!#        'Lat': [39.5, 42.0],
##!#        'Lon': [-122.0, -119.5],
##!#        '2155': {
##!#            'asos': 'asos_data_20210722_2.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021204.2155.061.2021205153516.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
##!#            #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
##!#            'Lat': [39.5, 42.0],
##!#            'Lon': [-122.0, -119.5],
##!#            #'Lat': [39.5, 42.0],
##!#            #'Lon': [-122.0, -119.5],
##!#            'data_lim': {
##!#                1:  [0.05, 0.5],
##!#                2:  [None, None],
##!#                3:  [None, None],
##!#                4:  [None, None],
##!#                5:  [None, None],
##!#                7:  [None, None],
##!#                31: [270., 330.],
##!#                32: [270., 330.],
##!#                'wv_ir': [0.2, 1.5],
##!#            },
##!#            #'modis_Lat': [39.5, 42.0],
##!#            #'modis_Lon': [-122.0, -119.5]
##!#            'modis_Lat': [39.5, 42.0],
##!#            'modis_Lon': [-122.0, -119.5]
##!#        }
##!#    },
##!#    "2021-08-04": {
##!#        '2110': {
##!#            'asos': 'asos_data_20210806.csv',
##!#            #'asos': 'asos_nevada_20210806.csv',
##!#            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
##!#            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.04.206.L2.SUBS2RET.v6.0.32.0.G21217152448.hdf',\
##!#                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.04.207.L2.SUBS2RET.v6.0.32.0.G21217152904.hdf'],\
##!#            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
##!#            'Lat': [36.0, 39.0],
##!#            'Lon': [-118.0, -114.0],
##!#            'modis_Lat': [35.0, 40.0],
##!#            'modis_Lon': [-119., -113.]
##!#        }
##!#    },
##!#    "2021-08-05": {
##!#        '2120': {
##!#            'asos': 'asos_data_20210805.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021217.2120.061.2021218164201.hdf',
##!#            'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021217.2120.061.2021218165546.hdf',
##!#            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
##!#            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0805t2038-o90733_v003-2021m0807t014855.he5',
##!#            'Lat': [36.5, 39.0],
##!#            'Lon': [-118.0, -114.0],
##!#            'modis_Lat': [36.0, 39.0],
##!#            'modis_Lon': [-118., -114.]
##!#        },
##!#        '2125': {
##!#            'asos': 'asos_california_20210805.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021217.2125.061.2021218161010.hdf',
##!#            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
##!#            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
##!#            'Lat': [39.5, 42.0],
##!#            'Lon': [-122.0, -119.0],
##!#            'data_lim': {
##!#                1:  [None, None],
##!#                2:  [None, 1.0],
##!#                3:  [None, None],
##!#                4:  [None, None],
##!#                5:  [None, None],
##!#                7:  [None, 1.0],
##!#                31: [270., 330.],
##!#                32: [270., 330.],
##!#                'wv_ir': [0.2, 1.5],
##!#            },
##!#            'modis_Lat': [39.5, 42.0],
##!#            'modis_Lon': [-122., -119.]
##!#        }
##!#    },
##!#    "2021-08-06": {
##!#        '2025': {
##!#            'asos': 'asos_data_20210806.csv',
##!#            #'asos': 'asos_nevada_20210806.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
##!#            'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
##!#            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080609-2021080620.nc',
##!#            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.06.204.L2.SUBS2RET.v6.0.32.0.G21219130523.hdf',\
##!#                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.06.205.L2.SUBS2RET.v6.0.32.0.G21219130455.hdf'],
##!#            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
##!#            'Lat': [36.0, 39.0],
##!#            'Lon': [-118.0, -114.0],
##!#            'data_lim': {
##!#                1:  [None, None],
##!#                2:  [None, None],
##!#                3:  [None, None],
##!#                4:  [None, None],
##!#                5:  [None, None],
##!#                7:  [None, None],
##!#                31: [295., 330.],
##!#                32: [270., 330.],
##!#                'wv_ir': [0.2, 1.5],
##!#                'SWF': [250., 300.],
##!#                'LWF': [300., 370.],
##!#            },
##!#            'modis_Lat': [36.0, 39.0],
##!#            'modis_Lon': [-118., -114.]
##!#        }
##!#    },
##!#    "2021-08-07": {
##!#        '2110': {
##!#            'asos': 'asos_data_20210806.csv',
##!#            #'asos': 'asos_nevada_20210806.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021219.2110.061.2021220151612.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
##!#            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.07.212.L2.SUBS2RET.v6.0.32.0.G21220123225.hdf'],\
##!#            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
##!#            'Lat': [36.0, 39.0],
##!#            'Lon': [-118.0, -114.0],
##!#            'modis_Lat': [35.0, 40.0],
##!#            'modis_Lon': [-119., -113.]
##!#        }
##!#    },
##!#    "2021-08-08": {
##!#        '2110': {
##!#            'asos': 'asos_data_20210806.csv',
##!#            #'asos': 'asos_nevada_20210806.csv',
##!#            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
##!#            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
##!#            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.08.202.L2.SUBS2RET.v6.0.32.0.G21221124420.hdf',\
##!#                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.08.203.L2.SUBS2RET.v6.0.32.0.G21221185932.hdf'],\
##!#            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
##!#            'Lat': [36.0, 39.0],
##!#            'Lon': [-118.0, -114.0],
##!#            'modis_Lat': [35.0, 40.0],
##!#            'modis_Lon': [-119., -113.]
##!#        }
##!#    },
##!#    "2021-08-17": {
##!#        '2145': {
##!#            'asos': 'asos_data_20210817.csv',
##!#            'Lat': [38.0, 42.0],
##!#            'Lon': [-122.0, -117.0]
##!#        }
##!#    },
##!#    "2021-08-30": {
##!#        '2115': {
##!#            'asos': 'asos_data_20210830.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021242.2115.061.2021243183953.hdf',
##!#            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
##!#            'Lat': [38.0, 40.0],
##!#            'Lon': [-121.0, -118.5],
##!#            'data_lim': {
##!#                1:  [None, None],
##!#                2:  [None, None],
##!#                3:  [None, None],
##!#                4:  [None, None],
##!#                5:  [None, None],
##!#                7:  [None, None],
##!#                31: [270., 330.],
##!#                32: [270., 330.],
##!#                'wv_ir': [0.2, 1.5],
##!#            },
##!#        }
##!#    },
##!#    "2021-09-01": {
##!#        '2105': {
##!#            'asos': 'asos_data_20210830.csv',
##!#            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021244.2105.061.2021245152256.hdf',
##!#            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
##!#            'Lat': [38.0, 42.0],
##!#            'Lon': [-121.5, -118.0],
##!#            'data_lim': {
##!#                1:  [None, None],
##!#                2:  [None, None],
##!#                3:  [None, None],
##!#                4:  [None, None],
##!#                5:  [None, None],
##!#                7:  [None, None],
##!#                31: [270., 330.],
##!#                32: [270., 330.],
##!#                'wv_ir': [0.2, 1.5],
##!#            },
##!#            'modis_Lat': [38.0, 42.0],
##!#            'modis_Lon': [-121.5, -118.]
##!#        }
##!#    } 
##!#}
##!#
