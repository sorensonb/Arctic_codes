"""
  NAME:
    GOESLib_lite
    
  PURPOSE:
    Read and plot GOES data using satpy. Also contains functions for 
        automatically downloading GOES data from the Amazon bucket. 

    NOTE: must be running a Conda environment with satpy installed. 

  SYNTAX:
    # = = = = = = = = = = =
    # Notes for downloading
    # = = = = = = = = = = =
    Before running the downloading scripts, the main of which is 
        auto_GOES_download(), you MUST have the Amazon AWS CLI 
        interface installed on your Linux terminal. To do this, 
        you must have sudo permissions on your computer. If you
        have sudo permissions, run the following in the command
        line:
    
    $ sudo apt install awscli

    Once this is installed, and if you have a Conda environment
        installed that has satpy and datetime (and the other essential
        packages listed in the "import" section), you may download
        GOES data automatically using auto_GOES_download. To run
        the function, use the following syntax:

    >>> from GOESLib_lite import *
    >>> begin_date = '202107201200' # begin_date in YYYYMMDDHHMM format
    >>> end_date   = '202107210000' # end_date in YYYYMMDDHHMM format
    >>> interval   = 30  # this is the number of minutes between each file
    >>> auto_GOES_download(begin_date, end_date, interval)
  
    # = = = = = = = = = = 
    # Notes for plotting
    # = = = = = = = = = =
    When running in a separate Python script, use the following syntax:
    >>> from GOESLib_lite import *
    >>> date_str = '202107202126'  # date_str in YYYYMMDDHHMM format, for whatever the time of the GOES file 
    >>> plot_GOES_satpy(date_str, 2)   # Plots for channel 2 data
    
    If you want to save the image, do:
    >>> plot_GOES_satpy(date_str, 2, save = True)

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>         - 2022/04/06:
        Written
 
"""
import numpy as np
import numpy.ma as ma
import sys
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy import find_files_and_readers
from satpy.scene import Scene
from satpy.writers import get_enhanced_image
import subprocess

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Set up global variables
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
datacrs = ccrs.PlateCarree()
data_dir = '/home/bsorenson/data/GOES/goes17_abi/'
debug = False

# - 0.64  (band 2, Red)
# - 2.25  (band 6, Cloud particle size)
# - 3.90  (band 7, Shortwave Window)
# - 6.18  (band 8, Upper-Level Water Vapor)
# - 6.95  (band 9, Mid-Level Water Vapor)
# - 7.34  (band 10, Lower-Level Water Vapor)
# - 10.35 (band 13, Clean IR Longwave Window)
channel_dict = {
    '1': {
        'name': 'Blue',
        'wavelength': 0.47
    },\
    '2': {
        'name': 'Red',
        'wavelength': 0.64
    },\
    '3': {
        'name': 'Veggie',
        'wavelength': 0.86
    },\
    '4': {
        'name': 'Cirrus',
        'wavelength': 1.38
    },\
    '5': {
        'name': 'Snow/Ice',
        'wavelength': 1.61
    },\
    '6': {
        'name': 'Cloud Particle Size',
        'wavelength': 2.25
    },\
    '7': {
        'name': 'Shortwave Window',
        'wavelength': 3.90
    },\
    '8': {
        'name': 'Upper-Level Water Vapor',
        'wavelength': 6.18
    },\
    '9': {
        'name': 'Mid-Level Water Vapor',
        'wavelength': 6.95
    },\
    '10': {
        'name': 'Lower-Level Water Vapor',
        'wavelength': 7.34
    },\
    '11': {
        'name': 'Cloud-Top Phase',
        'wavelength': 8.50
    },\
    '12': {
        'name': 'Ozone',
        'wavelength': 9.61
    },\
    '13': {
        'name': 'Clean IR Longwave Window',
        'wavelength': 10.35
    },\
    '14': {
        'name': 'IR Longwave Window',
        'wavelength': 11.20
    },\
    '15': {
        'name': 'Dirty IR Longwave Window',
        'wavelength': 12.30
    },\
    '16': {
        'name': 'CO2 Longwave Infrared',
        'wavelength': 13.30
    }
}

for key in channel_dict.keys():
    if(channel_dict[key]['wavelength'] is not None):
        channel_dict[key]['wavelength_label'] = \
            str(channel_dict[key]['wavelength']) + ' μm'
    else:
        channel_dict[key]['wavelength_label'] = ''

##!#plot_limits_dict = {
##!#    "2021-07-20": {
##!#        'asos': 'asos_data_20210720.csv',
##!#        'Lat': [39.5, 42.0],
##!#        'Lon': [-122.0, -119.5],
##!#        'data_lim': {
##!#            1:  [0.05, 0.5],
##!#            31: [270., 330.],
##!#        },
##!#        'goes_Lat': [39.5, 42.0],
##!#        'goes_Lon': [-122.0, -119.5]
##!#    }
##!#}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Download data
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# date_str: YYYYMMDDHHMM
# channels is a list containing the numbers of the desired channels
def download_GOES_bucket(date_str, sat = 'goes17', \
        channels = [2,6,8,9,10,13], dest_dir = \
        '/home/bsorenson/data/GOES/goes17_abi/'):

    # Convert the input date_str to datetime
    # --------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    #if((dt_date_str.minute == 0)):
    #    dt_date_str = dt_date_str - timedelta(minutes = 5)

    # Pull the entire file list for this date and time
    # ------------------------------------------------
    request_add = dt_date_str.strftime('s3://noaa-'+sat+'/ABI-L1b-RadC/2021/%j/%H/')

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = subprocess.run(['aws','s3','ls','--no-sign-request',\
            request_add], check=True, \
            capture_output=True).stdout.decode('utf-8').strip().split('\n') 
    except subprocess.CalledProcessError:
        print("ERROR: No ",sat," files for the input DTG",date_str)
        return

    # Remove the timestamps from the file strings
    # -------------------------------------------
    files_only = [tfile.strip().split()[3] for tfile in files]

    # Use the actual file times to get timestamps
    # -------------------------------------------
    file_dates = [datetime.strptime(tfile[27:40],'%Y%j%H%M%S') for tfile in files_only] 

    # Extract the channels
    # --------------------
    file_channels = np.array([int(tfile[19:21]) for tfile in files_only])

    # Figure out where the closest files to the desired time are
    # ----------------------------------------------------------
    time_diffs = np.array([abs((dt_date_str - ddate).total_seconds()) \
        for ddate in file_dates])

    # Select those files. Should result in all 16 channels for a single
    # time
    # -----------------------------------------------------------------
    files_found = np.array(files_only)[np.where(\
        time_diffs == np.min(time_diffs))]

    # Use the input channels to extract only the desired channels
    # -----------------------------------------------------------
    np_channels = np.array(channels) - 1
    good_files = files_found[np_channels]

    # Loop over the selected files and download them
    # ----------------------------------------------
    for gf in good_files:
        request_add = dt_date_str.strftime('s3://noaa-'+sat+\
            '/ABI-L1b-RadC/2021/%j/%H/' + gf)
        #print(request_add)
        cmnd_list = ['aws','s3','cp','--no-sign-request',\
            request_add, dest_dir]
        subprocess.run(cmnd_list, check=True, \
            capture_output=False)

# begin_date : YYYYMMDDHHMM
# end_date   : YYYYMMDDHHMM
# interval   : minutes between desired images
# sat        : default is 'goes17', but can set to 'goes16' or 
#              'goes18' to download other data.
# channels   : a list containing the channel numbers you want
#              to download. For example, [2] would download
#              just data for the 0.64 micron channel (channel 2)
#              [2,6,13] downloads data for the 0.64 micron channel,
#              the 2.25 micron channel, and the 10.35 micron channel.
# dest_dir   : the location where files are going to be downloaded to.
#              NOTE: This MUST be changed when running on a computer
#              that does not belong to bsorenson.
def auto_GOES_download(begin_date, end_date, interval, sat = 'goes17', \
        channels = [2,6,8,9,10,13], dest_dir = \
        '/home/bsorenson/data/GOES/goes17_abi/'):

    # Convert the input date_str to datetime
    # --------------------------------------
    begin_dt_date = datetime.strptime(begin_date,"%Y%m%d%H%M")
    end_dt_date   = datetime.strptime(end_date,"%Y%m%d%H%M")

    # Using the interval, get the desired file times for each
    # GOES image
    # -------------------------------------------------------
    num_seconds_in_minute = 60
    num_minutes = (end_dt_date - begin_dt_date).total_seconds() / \
        num_seconds_in_minute
    num_imgs = num_minutes / interval


    image_times = [begin_dt_date + timedelta(minutes = interval * ii) \
        for ii in range(int(num_imgs))]

    for ttime in image_times:
        print(ttime)
        download_GOES_bucket(ttime.strftime('%Y%m%d%H%M'), \
            sat = sat, channels = channels, dest_dir = dest_dir)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Miscellaneous plot-related functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# phi  = latitude of the point
# lbda = longitude of the point
def calc_zenith_angle(phi, lbda):

    # GOES-17 values
    sat_height = 35786.0234375 # km
    sat_lon    = -137.2
    r_e        = 6371.0 # km

    r_dist = r_e + sat_height
    
    # gamma_rad is the angle between the radius vectors of
    # point E (the satellite subpoint along the equator) and point
    # P (the point location)
    gamma_rad = np.arccos(np.cos(np.radians(phi)) * \
        np.cos(np.radians(sat_lon - lbda)))

    # d_dist is the distance from the point on the earth's surface 
    # to the satellite.
    d_dist = r_dist * (1. + (r_e / r_dist)**2. - 2. * \
        (r_e / r_dist) * np.cos(gamma_rad))**0.5 

    # Calculate the viewing zenith angle in degrees
    vza = np.degrees(np.arcsin((r_dist / d_dist) * np.sin(gamma_rad)))

    return vza

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

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

calib_dict = {
    0.64: 'reflectance',
    2.25: 'reflectance',
    3.90: 'radiance',
    6.18: 'brightness_temperature',
    6.95: 'brightness_temperature',
    7.34: 'brightness_temperature',
    10.35: 'brightness_temperature'
}   
label_dict = {
    0.64: '%',
    2.25: '%',
    3.90: 'mW m$^{-2}$ Sr$^{-1}$ (cm$^{-1}$)$^{-1}$',
    #6.18: 'mW m$^{-2}$ Sr$^{-1}$ (cm$^{-1}$)$^{-1}$',
    6.18: 'K',
    6.95: 'K',
    7.34: 'K',
    10.35: 'K'
}   
cmap_dict = {
    0.64:  'Greys_r',
    2.25:  'Greys_r',
    3.90:  'plasma',
    6.18:  'plasma',
    6.95:  'plasma',
    7.34:  'plasma',
    10.35: 'plasma'
}   
# Channel must be:
# - 0.64  (band 2, Red)
# - 2.25  (band 6, Cloud particle size)
# - 3.90  (band 7, Shortwave Window)
# - 6.18  (band 8, Upper-Level Water Vapor)
# - 6.95  (band 9, Mid-Level Water Vapor)
# - 7.34  (band 10, Lower-Level Water Vapor)
# - 10.35 (band 13, Clean IR Longwave Window)
def read_GOES_satpy(date_str, channel, zoom = True):


    # Extract the channel wavelength using the input string
    # -----------------------------------------------------
    channel = channel_dict[str(channel)]['wavelength']

    # Determine the correct GOES files associated with the date
    # ---------------------------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    dt_date_str_end = dt_date_str + timedelta(minutes = 10)

    # Use the Satpy find_files_and_readers to grab the files
    # ------------------------------------------------------
    files = find_files_and_readers(start_time = dt_date_str, \
        end_time = dt_date_str_end, base_dir = data_dir, reader = 'abi_l1b')

    print(files)

    # Extract the goes true-color plot limits
    # ----------------------------------------
    lat_lims = [29.5, 48.0]
    lon_lims = [-122.0, -100.5]

    # Use satpy (Scene) to open the file
    # ----------------------------------
    scn = Scene(reader = 'abi_l1b', filenames = files)

    # Load the desired channel data
    # -----------------------------
    scn.load([channel], calibration = [calib_dict[channel]])

    ## Set the map projection and center the data
    ## ------------------------------------------
    #my_area = scn[channel].attrs['area'].compute_optimal_bb_area({\
    #    'proj':'lcc', 'lon_0': lon_lims[0], 'lat_0': lat_lims[0], \
    #    'lat_1': lat_lims[0], 'lat_2': lat_lims[0]})
    #new_scn = scn.resample(my_area)

    ##!## Enhance the image for plotting
    ##!## ------------------------------
    ##!#var = get_enhanced_image(scn[channel]).data
    ##!#var = var.transpose('y','x','bands')

    # Zoom the image on the desired area
    # ----------------------------------
    if(zoom):
        scn = scn.crop(ll_bbox = (lon_lims[0] + 0.65, lat_lims[0], \
            lon_lims[1] - 0.65, lat_lims[1]))


    # Extract the lats, lons, and data
    # -----------------------------------------------------
    lons, lats = scn[channel].attrs['area'].get_lonlats()
    var = scn[channel].data

    # Extract the map projection from the data for plotting
    # -----------------------------------------------------
    crs = scn[channel].attrs['area'].to_cartopy_crs()

    # Extract the appropriate units
    # -----------------------------
    units = label_dict[channel]
    #units = scn[channel].units
    plabel = calib_dict[channel].title() + ' [' + units + ']'

    return var, crs, lons, lats, lat_lims, lon_lims, plabel
    #return var, crs, lat_lims, lon_lims

# channel must be an integer between 1 and 16
def plot_GOES_satpy(date_str, channel, ax = None, var = None, crs = None, \
        lons = None, lats = None, lat_lims = None, lon_lims = None, \
        vmin = None, vmax = None, ptitle = None, plabel = None, \
        labelsize = 10, colorbar = True, zoom=True,save=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    if(var is None): 
        var, crs, lons, lats, lat_lims, lon_lims, plabel = read_GOES_satpy(date_str, channel)

    # Plot the GOES data
    # ------------------
    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection=crs)

    ##!#ax.imshow(var.data, transform = crs, extent=(var.x[0], var.x[-1], \
    ##!#    var.y[-1], var.y[0]), vmin = vmin, vmax = vmax, origin='upper', \
    ##!#    cmap = 'Greys_r')
    #im1 = ax.imshow(var, transform = crs, vmin = vmin, vmax = vmax, \
    im1 = ax.pcolormesh(lons, lats, var, transform = datacrs, \
        vmin = vmin, vmax = vmax, \
        cmap = cmap_dict[channel_dict[str(channel)]['wavelength']], \
        shading = 'auto')
    ax.add_feature(cfeature.STATES)
    if(colorbar):
        cbar = plt.colorbar(im1, ax = ax, pad = 0.03, fraction = 0.052, \
            extend = 'both')
        cbar.set_label(plabel.replace('_',' '), size = labelsize, weight = 'bold')

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
    ##!#    ax.set_extent([lon_lims[0]+0.55,lon_lims[1]-0.6,lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
        zoom_add = '_zoom'
    else:
        zoom_add = ''

    # NOTE: commented out after removing the 'enhanced_image' code because
    #       it doesn't work now .
    ##!#ax.coastlines(resolution = '50m')
    ##!#ax.add_feature(cfeature.STATES)
    ##!#ax.add_feature(cfeature.BORDERS)
    if(ptitle is None):
        ax.set_title('GOES-17\n'+dt_date_str.strftime('%Y-%m-%d %H:%M'))
    else:
        ax.set_title(ptitle)

    if(not in_ax): 
        fig.tight_layout()
        if(save):
            outname = 'goes_' + date_str + zoom_add + '.png'
            plt.savefig(outname,dpi=300)
            print("Saved image",outname)
        else:
            plt.show()

def plot_GOES_satpy_6panel(date_str, ch1, ch2, ch3, ch4, ch5, ch6, \
        zoom = True, save = False):
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    plt.close('all')
    fig1 = plt.figure(figsize = (10,6))
    var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel0 = read_GOES_satpy(date_str, ch1)
    var1, crs1, lons1, lats1, lat_lims, lon_lims, plabel1 = read_GOES_satpy(date_str, ch2)
    var2, crs2, lons2, lats2, lat_lims, lon_lims, plabel2 = read_GOES_satpy(date_str, ch3)
    var3, crs3, lons3, lats3, lat_lims, lon_lims, plabel3 = read_GOES_satpy(date_str, ch4)
    var4, crs4, lons4, lats4, lat_lims, lon_lims, plabel4 = read_GOES_satpy(date_str, ch5)
    var5, crs5, lons5, lats5, lat_lims, lon_lims, plabel5 = read_GOES_satpy(date_str, ch6)

    ax0 = fig1.add_subplot(2,3,1, projection = crs0)
    ax1 = fig1.add_subplot(2,3,2, projection = crs1)
    ax2 = fig1.add_subplot(2,3,3, projection = crs2)
    ax3 = fig1.add_subplot(2,3,4, projection = crs3)
    ax4 = fig1.add_subplot(2,3,5, projection = crs4)
    ax5 = fig1.add_subplot(2,3,6, projection = crs5)

    ##!#ax1.set_title('GOES-17 Band ' + str(ch2) + '\n' + \
    ##!#    channel_dict[str(ch2)]['name'] + '\n' + \
    labelsize = 11
    font_size = 10
    plot_GOES_satpy(date_str, ch1, ax = ax0, var = var0, crs = crs0, \
        lons = lons0, lats = lats0, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = None, vmax = None, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch2, ax = ax1, var = var1, crs = crs0, \
        lons = lons1, lats = lats1, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = None, vmax = 30., ptitle = '', plabel = plabel1, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch3, ax = ax2, var = var2, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = 270, vmax = 330, ptitle = '', plabel = plabel2, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch4, ax = ax3, var = var3, crs = crs0, \
        lons = lons3, lats = lats3, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = None, vmax = None, ptitle = '', plabel = plabel3, \
        colorbar = True, labelsize = labelsize, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch5, ax = ax4, var = var4, crs = crs0, \
        lons = lons4, lats = lats4, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = None, vmax = None, ptitle = '', plabel = plabel4, \
        colorbar = True, labelsize = labelsize, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch6, ax = ax5, var = var5, crs = crs0, \
        lons = lons5, lats = lats5, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = None, vmax = None, ptitle = '', plabel = plabel5, \
        colorbar = True, labelsize = labelsize, zoom=True,save=False)

    plot_figure_text(ax0, 'GOES-17 ' + \
        str(channel_dict[str(ch1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax1, 'GOES-17 ' + \
        str(channel_dict[str(ch2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax2, 'GOES-17 ' + \
        str(channel_dict[str(ch3)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax3, 'GOES-17 ' + \
        str(channel_dict[str(ch4)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax4, 'GOES-17 ' + \
        str(channel_dict[str(ch5)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax5, 'GOES-17 ' + \
        str(channel_dict[str(ch6)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')

    plot_subplot_label(ax0,  '(a)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax1,  '(b)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax2,  '(c)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax3,  '(d)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax4,  '(e)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax5,  '(f)', backgroundcolor = 'white', fontsize = font_size)

    lon_stn = -120.7605
    lat_stn = 41.2098

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
    ##!#    ax0.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
    ##!#    ax1.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
    ##!#    ax2.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
    ##!#    ax3.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
    ##!#    ax4.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
    ##!#    ax5.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
        zoom_add = '_zoom'
    else:
        zoom_add = ''

    fig1.tight_layout()

    if(save):
        outname = 'goes17_'+date_str+'_6panel.png'
        fig1.savefig(outname, dpi = 300)
        print('Saved image', outname)
    else:
        plt.show()
    
