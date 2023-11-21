#!/usr/bin/env python

"""
  NAME:
    plot_OMI_fourpanel_falsecolor.py

  PURPOSE:
    Generate a four-panel OMI image, consisting of (a) the 500 nm normalized
        radiance, (b) the 388 nm normalized radiance, (c) the 354 nm
        normalized radiance, and (d) a false-color image using 
        these three channels. By adding "equalize" as a command line argument,
        a histogram equalization procedure is applied to each of the three
        channels in the false color image.

  SYNTAX:
    $ ./plot_OMI_fourpanel_falsecolor.py [equalize]

        equalize: optional argument. Adding "equalize" as an argument
            applies the histogram equalization procedure.

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2023/11/20:
      Written

"""

import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# BEGIN FUNCTION DEFINITIONS
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Perform nonlinear histogram equalizing stuff
def hist_equal(data):
    
    # Prep the data for the procedure by neatly handling missing pixels.    
    # Also, convert the data to integer values.
    # -----------------------------------------------------------------
    local_data = np.ma.masked_where(data < 0, data)
    local_data[(local_data.mask == True) | (local_data.data < 0)] = -9.
    local_data = local_data.data.astype(int)
    local_data = np.ma.masked_where(local_data < 0, local_data)
    work_data = local_data.compressed()
 
    # Determine the number of pixel values to work with
    # -------------------------------------------------
    pixel_vals = int(np.nanmax(work_data)) + 1 
    
    # Set up an array containing all unique values from 0 to
    # the maximum value in the input channel
    # ------------------------------------------------------ 
    bins = np.arange(pixel_vals+1) 
    values = np.full((pixel_vals), -9.)

    # Calculate the histogram of the input channel
    # --------------------------------------------
    for ii in range(len(values)):
        values[ii] = len(np.where(work_data == ii)[0])
    
    # Calculate the cumulative histogram
    # ----------------------------------
    cum_hist = np.cumsum(values)
 
    # Apply the normalization methods
    # ------------------------------- 
    new_hist = ((pixel_vals - 1.)/cum_hist[-1]) * cum_hist
    new_bright = np.round(new_hist,0)

    # Set up an output array
    # ----------------------
    new_data = np.copy(local_data)
    new_data = np.ma.masked_where(new_data < 0, new_data)

    # Loop over the new data and adjust using the new 
    # brightness values
    # -----------------------------------------------
    for ii in range(len(new_bright)):
        new_data[local_data == bins[ii]] = new_bright[ii]

    return new_data

# plot_subplot_label allows the user to easily add subplot labels to a figure
# ---------------------------------------------------------------------------
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

# plot_swath_on_map automates the plotting of the OMI data from a single
# swath. The function takes the following required arguments:
#
# - infile: string name consisting of the entire filename (path and file) 
#       to a single OMI HDF5 file.
# - ax?:    each individual figure subplot axis object
#
# Optional arguments include:
# - minlon/lat: the minimum latitude or longitude value to use in the plot
# - maxlon/lat: the maximum latitude or longitude value to use in the plot
# - equalize: By setting to True, applies histogram equalization to each
#               channel in the false color image.
# -------------------------------------------------------------------------
def plot_swath_on_map(infile, ax1, ax2, ax3, ax4, minlon = -60, maxlon = 50, \
        minlat = -20, maxlat = 40, equalize = True):

    # Open the file and extract the data
    # ----------------------------------
    data = h5py.File(infile)
    allrad =  data[\
        'HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/NormRadiance'][:,:,:]
    allrad = np.ma.masked_where(allrad <= 0, allrad)
    LAT    = data[\
        'HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
    LON    = data[\
        'HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
    data.close()

    # Mask the data that is outside of the desired lat/lon bounds
    # -----------------------------------------------------------
    rad354 = np.ma.masked_where( (LON < minlon) | (LON > maxlon)  | \
        (LAT < minlat) | (LAT > maxlat), allrad[:,:,0])
    rad388 = np.ma.masked_where( (LON < minlon) | (LON > maxlon)  | \
        (LAT < minlat) | (LAT > maxlat), allrad[:,:,1])
    rad500 = np.ma.masked_where( (LON < minlon) | (LON > maxlon)  | \
        (LAT < minlat) | (LAT > maxlat), allrad[:,:,2])


    # Scale the RGB values to 0 - 255.
    # --------------------------------
    allmaxs = np.array([np.max(rad354), np.max(rad388), np.max(rad500)])
    red   = rad500*(255./np.max(allmaxs)) 
    green = rad388*(255./np.max(allmaxs)) 
    blue  = rad354*(255./np.max(allmaxs)) 

    # If desired, perform the histogram equalization on each channel.
    # Doing so also converts the values to integers
    # ---------------------------------------------------------------
    if(equalize):
        red   = hist_equal(red)
        green = hist_equal(green)
        blue  = hist_equal(blue)
    
    # Create color scales for each RGB channel
    # ----------------------------------------
    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,110,160,210,240,255])
    plotred     = np.interp(red,old_val,ehn_val) / 255.
    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,110,160,200,230,240])
    plotgreen   = np.interp(green,old_val,ehn_val) / 255.
    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,100,150,210,240,255])
    plotblue    = np.interp(blue,old_val,ehn_val) / 255.
    
    # Combine the three RGB channels into 1 3-d array
    # -----------------------------------------------
    image = np.zeros((plotred.shape[0],plotred.shape[1],3))
    image[:,:,0] = plotred
    image[:,:,1] = plotgreen
    image[:,:,2] = plotblue
    
    # Convert the color values into a format usable for plotting
    # ----------------------------------------------------------
    colortuple = tuple(np.array([image[:,:,0].flatten(), \
        image[:,:,1].flatten(), image[:,:,2].flatten()]).transpose().tolist())
       
    # Plot the "Red" channel
    # ---------------------- 
    cmap = 'viridis'
    datacrs = ccrs.PlateCarree()
    ax1.pcolormesh(LON, LAT, red, transform = datacrs, shading = 'auto', cmap = cmap)
    ax1.coastlines()
    ax1.set_extent([minlon, maxlon, minlat, maxlat], datacrs)

    # Plot the "Green" channel
    # ------------------------ 
    ax2.pcolormesh(LON, LAT, green, transform = datacrs, shading = 'auto', cmap = cmap)
    ax2.coastlines()
    ax2.set_extent([minlon, maxlon, minlat, maxlat], datacrs)
    
    # Plot the "Blue" channel
    # ------------------------ 
    ax3.pcolormesh(LON, LAT, blue, transform = datacrs, shading = 'auto', cmap = cmap)
    ax3.coastlines()
    ax3.set_extent([minlon, maxlon, minlat, maxlat], datacrs)

    # Plot the false color image
    # --------------------------
    ax4.pcolormesh(LON,LAT,\
        image[:,:,0], color= colortuple, \
        shading='auto', transform = datacrs) 
    ax4.coastlines()
    ax4.set_extent([minlon, maxlon, minlat, maxlat], datacrs)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# END FUNCTION DEFINITIONS
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Check arguments
# ---------------
if('equalize' in sys.argv):
    equalize = True
else:
    equalize = False

# NOTE: The paths to each of these three files must be changed.
# Also, additional files can be added here to plot more OMI swaths
# ----------------------------------------------------------------
infiles = [\
    '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2007m0728t0927-o16137_v003-2017m0721t040315.he5',
    '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2007m0728t1106-o16138_v003-2017m0721t040329.he5',
    '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2007m0728t1245-o16139_v003-2017m0721t040251.he5'
    ]

# Set up the figure
# -----------------
fig = plt.figure(figsize = (8, 6.5))
ax1 = fig.add_subplot(2,2,1, projection = ccrs.PlateCarree())
ax2 = fig.add_subplot(2,2,2, projection = ccrs.PlateCarree())
ax3 = fig.add_subplot(2,2,3, projection = ccrs.PlateCarree())
ax4 = fig.add_subplot(2,2,4, projection = ccrs.PlateCarree())

# Loop over each desired OMI swath and plot the data on the figure
# ----------------------------------------------------------------
for tfile in infiles:
    plot_swath_on_map(tfile, ax1, ax2, ax3, ax4, equalize = equalize)

# Add subplot titles
# ------------------
if(equalize):
    title_adder = '\nHistogram Equalized'
else:
    title_adder = ''
ax1.set_title('500 nm Norm. Radiance' + title_adder)
ax2.set_title('388 nm Norm. Radiance' + title_adder)
ax3.set_title('354 nm Norm. Radiance' + title_adder)
ax4.set_title('False Color Image\n(R = 500 nm; G = 388 nm; B = 354 nm)' + title_adder)

# Add subplot labels
# ------------------
plot_subplot_label(ax1, '(a)', location = 'upper_left')
plot_subplot_label(ax2, '(b)', location = 'upper_left')
plot_subplot_label(ax3, '(c)', location = 'upper_left')
plot_subplot_label(ax4, '(d)', location = 'upper_left')

# Add an overall plot title
# -------------------------
plt.suptitle('OMI Normalized Radiances\n2007-07-28')

fig.tight_layout()

# Save the image
# --------------
if(equalize):
    outname = 'omi_normrad_falsecolor_equalize.png'
else:
    outname = 'omi_normrad_falsecolor.png'
fig.savefig(outname, dpi = 300)
print("Saved image", outname)
#plt.show()
