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
import glob
from PIL import Image
import h5py
import matplotlib.pyplot as plt
import cartopy.feature as cfeature

# Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

datacrs = ccrs.PlateCarree()
#mapcrs = ccrs.LambertConformal(central_longitude = -110, central_latitude = 40)

# Find the gridpoint in the gridded lat/lon data that 
# corresponds to the station at slat and slon
# ---------------------------------------------------- 
def nearest_gridpoint(slat, slon, grid_lat, grid_lon):
    fun_c = ((slat - grid_lat)**2. + (slon - grid_lon)**2.)**0.5
    #fun_c = np.maximum(np.abs(grid_lat - slat), \
    #    np.abs(grid_lon - slon))
    m_idx = np.where(fun_c == np.min(fun_c))
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

def plot_VIIRS_granule(filename, ax):
    
    # Read the filename to determine the file type
    # --------------------------------------------
    total_split = filename.split('/')
    name_split = total_split[-1].split('.')
    dataset_name = name_split[0]
    if(dataset_name == 'VNP46A1'):
        viirs_ds = h5py.File(filename, 'r')
        lat_max = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/'].attrs['NorthBoundingCoord']
        lat_min = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/'].attrs['SouthBoundingCoord']
        lon_max = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/'].attrs['EastBoundingCoord']
        lon_min = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/'].attrs['WestBoundingCoord']
        dnb = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/DNB_At_Sensor_Radiance_500m']

        lats = np.linspace(lat_max[0], lat_min[0], dnb.shape[1])
        lons = np.linspace(lon_min[0], lon_max[0], dnb.shape[0])
        x, y = np.meshgrid(lons, lats)

    elif(dataset_name == 'VNP02DNB'):
        # Determine if geolocation data is present in same path
        # -----------------------------------------------------
        geoloc_list = glob.glob('/'.join(total_split[:-1]) + '/VNP03DNB.' + \
            '.'.join(name_split[1:4]) + '*') 
        if(len(geoloc_list) == 0):
            print("ERROR: no geolocation found for file",filename)
            return
        else:
            viirs_ds   = h5py.File(filename, 'r')
            viirs_ds_g = h5py.File(geoloc_list[0], 'r')

            # Extract the DNB values
            # ----------------------
            dnb = np.log10(viirs_ds['observation_data/DNB_observations'][:])

            # Extract the lat and lons
            # ------------------------
            x = viirs_ds_g['geolocation_data']['longitude'][:]
            y = viirs_ds_g['geolocation_data']['latitude'][:]

            viirs_ds.close()
            viirs_ds_g.close() 

    ax.pcolormesh(x,y, dnb, cmap = 'cividis', vmin = 0, vmax = 100, transform = ccrs.PlateCarree())
    viirs_ds.close()

def plot_VIIRS(filename):
    plt.close('all')
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(1,1,1, projection = mapcrs)
    if(isinstance(filename, list)):

        for ff in filename:
            plot_VIIRS_granule(ff,ax)

    elif(isinstance(filename, str)):
        plot_VIIRS_granule(filename,ax)

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES)
    plt.show()

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
        halign = 'left'):

    if(xval is None):
        print(len(text))
        xval = ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.95
    if(yval is None):
        yval = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.05
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
    frames = [Image.open(image) for image in glob.glob(f"{frame_folder}/*.png")]
    frame_one = frames[0]
    frame_one.save(gif_name, format = 'GIF', append_images = frames,\
        save_all = True, duration = 200, loop = 0)
    print("Saved gif",gif_name)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Event dictionary  
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
aerosol_event_dict = {
    "2017-08-17": {
        '0000': {
            #'asos': 'asos_data_20210713.csv',
            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2017231.1810.061.2018034125202.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2017m0817t2043-o69632_v003-2017m0818t223219.he5',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2017081700-2017081723.nc',
            'ceres_time': '21',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [80.0, 90.0],
            'Lon': [-125.0, -75.0],
            'modis_Lat': [80.0, 90.0],
            'modis_Lon': [-125.0, -75.0],
        },
    },
    "2017-08-18": {
        '0000': {
            #'asos': 'asos_data_20210713.csv',
            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2017231.1810.061.2018034125202.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2017m0818t1630-o69644_v003-2017m0819t230655.he5',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2017081800-2017081823.nc',
            'ceres_time': '17',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [80.0, 90.0],
            'Lon': [-95.0, -15.0],
            'modis_Lat': [80.0, 90.0],
            'modis_Lon': [-95.0, -15.0],
        },
    },
    "2017-08-19": {
        '1810': {
            #'asos': 'asos_data_20210713.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2017231.1810.061.2018034125202.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2017m0819t1713-o69659_v003-2017m0820t234933.he5',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/minlat_55/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2017081900-2017081923.nc',
            'ceres_time': '18',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [82.0, 88.0],
            'Lon': [-105.0, -25.0],
            'modis_Lat': [82.0, 88.0],
            'modis_Lon': [-105.0, -25.0],
        },
    },
    "2019-08-11": {
        '0130': {
            #'asos': 'asos_data_20210713.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0130.061.2019223152859.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2019m0811t0033-o80163_v003-2019m0812t234149.he5',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2019081100-2019081123.nc',
            'ceres_time': '01',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [68.0, 82.0],
            'Lon': [95.0, 161.0],
            'modis_Lat': [68.0, 82.0],
            'modis_Lon': [95.0, 161.0],
        },
        '0310': {
            #'asos': 'asos_data_20210713.csv',
            'modis': ['/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0305.061.2019223152906.hdf',\
                      '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0310.061.2019223152911.hdf'], 
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2019m0811t0212-o80164_v003-2019m0812t234246.he5',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2019081100-2019081123.nc',
            'ceres_time': '03',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [68.0, 82.0],
            'Lon': [95.0, 161.0],
            'modis_Lat': [68.0, 82.0],
            'modis_Lon': [95.0, 161.0],
        },
        '0445': {
            #'asos': 'asos_data_20210713.csv',
            'modis': ['/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2019223.0445.061.2019223152919.hdf'], 
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2019m0811t0351-o80165_v003-2019m0812t234332.he5',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2019081100-2019081123.nc',
            'ceres_time': '04',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [74.0, 80.0],
            'Lon': [110.0, 161.0],
            'modis_Lat': [74.0, 80.0],
            'modis_Lon': [110.0, 161.0],
        }
    },
    "2021-07-13": {
        '2110': {
            'asos': 'asos_data_20210713.csv',
            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            'modis_Lat': [39.0, 42.5],
            'modis_Lon': [-123., -119.]
        }
    },
    "2021-07-20": {
        '2125': {
            'asos': 'asos_data_20210720.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021201.2125.061.2021202154814.hdf',
            'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021201.2125.061.2021202225717.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/modis_comp/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072010-2021072021.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.20.214.L2.SUBS2RET.v6.0.32.0.G21202153435.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            'data_lim': {
                1:  [None, None],
                2:  [None, None],
                3:  [None, None],
                4:  [None, None],
                5:  [None, None],
                7:  [None, None],
                31: [270., 330.],
                32: [270., 330.],
                'wv_ir': [0.2, 2.5],
            },
            'modis_Lat': [39.5, 42.0],
            'modis_Lon': [-122.0, -119.5]
        }
    },
    "2021-07-21": {
        '2030': {
            'asos': 'asos_data_20210722_4.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021202.2030.061.2021203174050.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072109-2021072122.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.21.205.L2.SUBS2RET.v6.0.32.0.G21203185004.hdf'],
            #'Lat': [39.5, 42.0],
            #'Lon': [-122.0, -119.5],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            'data_lim': {
                1:  [None, None],
                2:  [None, None],
                3:  [None, None],
                4:  [None, None],
                5:  [None, None],
                7:  [None, None],
                31: [270., 330.],
                32: [270., 330.],
                'wv_ir': [0.2, 1.5],
            },
            'modis_Lat': [39.5, 42.0],
            'modis_Lon': [-122.0, -119.5]
        }
    },
    "2021-07-22": {
        '2110': {
            'asos': 'asos_data_20210722_4.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
            'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            #'Lat': [39.5, 42.0],
            #'Lon': [-122.0, -119.5],
            'data_lim': {
                1:  [0.05, 0.5],
                2:  [None, None],
                3:  [None, None],
                4:  [None, None],
                5:  [None, None],
                7:  [None, None],
                31: [270., 330.],
                32: [270., 330.],
                'wv_ir': [0.2, 1.5],
            },
            #'modis_Lat': [39.5, 42.0],
            #'modis_Lon': [-122.0, -119.5]
            'modis_Lat': [39.5, 42.0],
            'modis_Lon': [-122.0, -119.5]
        }
    },
    "2021-07-23": {
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        '2155': {
            'asos': 'asos_data_20210722_2.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021204.2155.061.2021205153516.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            #'Lat': [39.5, 42.0],
            #'Lon': [-122.0, -119.5],
            'data_lim': {
                1:  [0.05, 0.5],
                2:  [None, None],
                3:  [None, None],
                4:  [None, None],
                5:  [None, None],
                7:  [None, None],
                31: [270., 330.],
                32: [270., 330.],
                'wv_ir': [0.2, 1.5],
            },
            #'modis_Lat': [39.5, 42.0],
            #'modis_Lon': [-122.0, -119.5]
            'modis_Lat': [39.5, 42.0],
            'modis_Lon': [-122.0, -119.5]
        }
    },
    "2021-08-04": {
        '2110': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.04.206.L2.SUBS2RET.v6.0.32.0.G21217152448.hdf',\
                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.04.207.L2.SUBS2RET.v6.0.32.0.G21217152904.hdf'],\
            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'modis_Lat': [35.0, 40.0],
            'modis_Lon': [-119., -113.]
        }
    },
    "2021-08-05": {
        '2120': {
            'asos': 'asos_data_20210805.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021217.2120.061.2021218164201.hdf',
            'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021217.2120.061.2021218165546.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0805t2038-o90733_v003-2021m0807t014855.he5',
            'Lat': [36.5, 39.0],
            'Lon': [-118.0, -114.0],
            'modis_Lat': [36.0, 39.0],
            'modis_Lon': [-118., -114.]
        },
        '2125': {
            'asos': 'asos_california_20210805.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021217.2125.061.2021218161010.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.0],
            'data_lim': {
                1:  [None, None],
                2:  [None, 1.0],
                3:  [None, None],
                4:  [None, None],
                5:  [None, None],
                7:  [None, 1.0],
                31: [270., 330.],
                32: [270., 330.],
                'wv_ir': [0.2, 1.5],
            },
            'modis_Lat': [39.5, 42.0],
            'modis_Lon': [-122., -119.]
        }
    },
    "2021-08-06": {
        '2025': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080609-2021080620.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.06.204.L2.SUBS2RET.v6.0.32.0.G21219130523.hdf',\
                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.06.205.L2.SUBS2RET.v6.0.32.0.G21219130455.hdf'],
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'data_lim': {
                1:  [None, None],
                2:  [None, None],
                3:  [None, None],
                4:  [None, None],
                5:  [None, None],
                7:  [None, None],
                31: [270., 330.],
                32: [270., 330.],
                'wv_ir': [0.2, 1.5],
            },
            'modis_Lat': [36.0, 39.0],
            'modis_Lon': [-118., -114.]
        }
    },
    "2021-08-07": {
        '2110': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021219.2110.061.2021220151612.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.07.212.L2.SUBS2RET.v6.0.32.0.G21220123225.hdf'],\
            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'modis_Lat': [35.0, 40.0],
            'modis_Lon': [-119., -113.]
        }
    },
    "2021-08-08": {
        '2110': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.08.202.L2.SUBS2RET.v6.0.32.0.G21221124420.hdf',\
                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.08.203.L2.SUBS2RET.v6.0.32.0.G21221185932.hdf'],\
            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'modis_Lat': [35.0, 40.0],
            'modis_Lon': [-119., -113.]
        }
    },
    "2021-08-17": {
        '2145': {
            'asos': 'asos_data_20210817.csv',
            'Lat': [38.0, 42.0],
            'Lon': [-122.0, -117.0]
        }
    },
    "2021-08-30": {
        '2115': {
            'asos': 'asos_data_20210830.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021242.2115.061.2021243183953.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
            'Lat': [38.0, 40.0],
            'Lon': [-121.0, -118.5],
            'data_lim': {
                1:  [None, None],
                2:  [None, None],
                3:  [None, None],
                4:  [None, None],
                5:  [None, None],
                7:  [None, None],
                31: [270., 330.],
                32: [270., 330.],
                'wv_ir': [0.2, 1.5],
            },
        }
    },
    "2021-09-01": {
        '2105': {
            'asos': 'asos_data_20210830.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021244.2105.061.2021245152256.hdf',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
            'Lat': [38.0, 42.0],
            'Lon': [-121.5, -118.0],
            'data_lim': {
                1:  [None, None],
                2:  [None, None],
                3:  [None, None],
                4:  [None, None],
                5:  [None, None],
                7:  [None, None],
                31: [270., 330.],
                32: [270., 330.],
                'wv_ir': [0.2, 1.5],
            },
            'modis_Lat': [38.0, 42.0],
            'modis_Lon': [-121.5, -118.]
        }
    } 
}

