#!/usr/bin/env python

"""
  NAME:

  PURPOSE:


  SYNTAX:
    $ python write_aerosol_profiles.py

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2024/08/07
      Cleaned up the code to make it more readable
      Modified Logan Lee's code to just print the profiles without
          averaging.

    Logan Lee - 2020/05
      Modified

    Logan Lee - 2018/05
      Written

"""

#Program to average CATS vertical profiles
#Written by Logan Lee May 2018 Modified May 2020
#Modified by J.Zhang to output aerosol profiles

import numpy as np
import os
import h5py
import fnmatch
from datetime import datetime, timedelta
from glob import glob
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import matplotlib.colors as cm
import sys
import os

home_dir = os.environ['HOME']
cats_dir = home_dir + '/Research/CATS/'

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# BEGIN LOGICAL SWITCHES
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# When set to True, this logical switch makes the script only print the
# CATS layers above 10 km. If Falso, all above-ground layers are printed
# ----------------------------------------------------------------------
l_ONLY_STRAT = True

# When True, this switch prevents the script from printing aerosol layers
# that contain zero-value extinctions
# -----------------------------------------------------------------------
l_ONLY_NONZERO_LAYERS = True

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# END LOGICAL SWITCHES
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

filename  = cats_dir + 'CATS-ISS_L2O_D-M7.2-V3-01_05kmPro.2017-07-25T14-08-20T14-55-56UTC.hdf5'
#filename  = cats_dir + 'CATS-ISS_L2O_N-M7.2-V3-01_05kmPro.2017-07-25T19-34-08T19-48-01UTC.hdf5'

print(filename)

infile=h5py.File(filename,'r')

# Read variables from the HDF5 file
# ---------------------------------
catlat    = infile['/geolocation/CATS_Fore_FOV_Latitude'][:]
catlon    = infile['/geolocation/CATS_Fore_FOV_Longitude'][:]
time      = infile['/profile/Profile_UTC_Time'][:]
date      = infile['/profile/Profile_UTC_Date'][:]
aod       = infile['/profile/Aerosol_Optical_Depth_1064_Fore_FOV'][:]
aer_type  = infile['/profile/Aerosol_Type_Fore_FOV'][:]
cld_phase = infile['/profile/Cloud_Phase_Fore_FOV'][:]
feattype  = infile['/profile/Feature_Type_Fore_FOV'][:]
sky       = infile['/profile/Sky_Condition_Fore_FOV'][:]
qc        = infile['/profile/Extinction_QC_Flag_1064_Fore_FOV'][:]
extuncat  = infile['/profile/Extinction_Coefficient_Uncertainty_1064_Fore_FOV'][:]
featscore = infile['/profile/Feature_Type_Score_Fore_FOV'][:]
dnfcat    = infile['/profile/Day_Night_Flag'][:]
opac      = infile['/profile/Percent_Opacity_Fore_FOV'][:]
demsfc    = infile['/profile/DEM_Surface_Altitude_Fore_FOV'][:]
ext       = infile['/profile/Extinction_Coefficient_1064_Fore_FOV'][:]
bin_alt   = infile['/metadata_parameters/Bin_Altitude_Array'][:]

# Close the file
infile.close()

# Extract the mean lat, lon, and time values 
# ------------------------------------------
catlon = catlon[:,1]
catlat = catlat[:,1]
time = time[:,1]

# Construct datetime objects to hold the full date/time information
# -----------------------------------------------------------------
juldayc2 = np.array([datetime.strptime(str(tdate), '%Y%m%d') + \
    timedelta(days = ttime) for tdate, ttime in zip(date, time)])


mask_ext = np.ma.masked_where( (ext < 0) | (ext > 100), ext)
mask_ext = np.log10(mask_ext)


work_time = juldayc2
#work_time = time
work_alt = bin_alt

#fig = plt.figure(figsize = (12, 9))
fig = plt.figure(figsize = (12, 11))
axs = fig.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
ax1 = axs[0]
ax2 = axs[3]
ax3 = axs[1]
ax4 = axs[2]
#ax1 = fig.add_subplot(3,1,1)
#ax2 = fig.add_subplot(3,1,2)
#ax3 = fig.add_subplot(3,1,3)

mesh = ax1.pcolormesh(work_time, work_alt, mask_ext, shading = 'auto')
cbar = fig.colorbar(mesh, ax = ax1, pad = 0.04, fraction = 0.040)
ax1.axhline(10, linestyle = ':', color = 'k')
ax1.set_xlabel('Time')
ax1.set_ylabel('CATS bin height [km]')
ax1.set_title('Log Extinction')

cmap = cm.ListedColormap(['blue','cyan','grey','orange','black'])
cbar_labels = ['N/A','Cloud','Undet.','Aerosol','Atten.']
bounds = [0,1,2,3,4,5]
tickvals = (np.array(bounds[1:]) + np.array(bounds[:-1])) / 2
norm = cm.BoundaryNorm(bounds, cmap.N)
mesh = ax2.pcolormesh(work_time, work_alt, feattype, shading = 'auto', cmap = cmap)
cbar = fig.colorbar(ScalarMappable(cmap = cmap, norm = norm), ax = ax2, \
                    ticks = tickvals, pad = 0.04, fraction = 0.040)
cbar.ax.set_yticklabels(cbar_labels,fontsize=10, rotation=0)
#cbar = fig.colorbar(mesh, ax = ax2)
ax2.axhline(10, linestyle = ':', color = 'yellow')
ax2.set_xlabel('Time')
ax2.set_ylabel('CATS bin height [km]')
ax2.set_title('Feature Type')

cmap = cm.ListedColormap(['gray','blue','cyan','yellow','orange','green','red','black','brown'])
cbar_labels = ['N/A','Marine','Marine Mix.','Dust','Dust Mix','Clean/Bkgd','Pol. Cont.','Smoke','Volcanic']
bounds = [0,1,2,3,4,5,6,7,8,9]
tickvals = (np.array(bounds[1:]) + np.array(bounds[:-1])) / 2
norm = cm.BoundaryNorm(bounds, cmap.N)
mesh = ax3.pcolormesh(work_time, work_alt, aer_type, shading = 'auto', cmap = cmap)
cbar = fig.colorbar(ScalarMappable(cmap = cmap, norm = norm), ax = ax3, ticks = tickvals, pad = 0.04, fraction = 0.040)
cbar.ax.set_yticklabels(cbar_labels,fontsize=10, rotation=0)
#cbar = fig.colorbar(mesh, ax = ax3)
ax3.axhline(10, linestyle = ':', color = 'yellow')
ax3.set_xlabel('Time')
ax3.set_ylabel('CATS bin height [km]')
ax3.set_title('Aerosol Type')

cmap = cm.ListedColormap(['palegreen','blue','gray','white'])
cbar_labels = ['N/A','Water cloud','Unknown','Ice cloud']
bounds = [0,1,2,3,4]
tickvals = (np.array(bounds[1:]) + np.array(bounds[:-1])) / 2
norm = cm.BoundaryNorm(bounds, cmap.N)
mesh = ax4.pcolormesh(work_time, work_alt, cld_phase, shading = 'auto', cmap = cmap)
cbar = fig.colorbar(ScalarMappable(cmap = cmap, norm = norm), ax = ax4, \
    ticks = tickvals, pad = 0.04, fraction = 0.040)
cbar.ax.set_yticklabels(cbar_labels,fontsize=10, rotation=0)
#cbar = fig.colorbar(mesh, ax = ax4)
ax4.axhline(10, linestyle = ':', color = 'yellow')
ax4.set_xlabel('Time')
ax4.set_ylabel('CATS bin height [km]')
ax4.set_title('Cloud phase')

plt.suptitle(filename)

#
#sys.exit()

## Extract the mean lat, lon, and time values 
## ------------------------------------------
#catlon = catlon[:,1]
#catlat = catlat[:,1]
#time = time[:,1]
#
## Construct datetime objects to hold the full date/time information
## -----------------------------------------------------------------
#juldayc2 = np.array([datetime.strptime(str(tdate), '%Y%m%d') + \
#    timedelta(days = ttime) for tdate, ttime in zip(date, time)])

# Variables for QA procedure
# Valid range is from -10 to 10.
# -10 is confidently aerosol, 10 is confidently cloud, 0 could be either.
featscore2=np.copy(featscore).astype(float)  
featscore2[featscore2 < -10]  = np.nan   #Don't want undefined values

catcadmax    = np.nanmax(featscore2,axis=0) #Max featscore
catcadmin    = np.nanmin(featscore2,axis=0) #Min featscore
catextunmax  = np.nanmax(extuncat,axis=0)   #Max Extinction Uncertainty
catextunmin  = np.nanmin(extuncat,axis=0)   #Min Extintion Uncertainty
catextmax    = np.nanmax(ext,axis=0)        #Max Extinction
featscoremax = np.nanmax(featscore,axis=0)


#CATS QC
#Feature Type Values 0=, 1=, 2=, 3=, 4=
clouds=np.unique(np.where(feattype == 1)[1])  #Get the x-values of cloudy pixels (not the heights)
cloud_param=np.zeros(aod.shape)
cloud_param[clouds]=1   #clouds present = 1, no clouds = 0
#Excluding profiles with underfined feature types...    
unds=np.unique(np.where(feattype == 2)[1])
unds_param=np.zeros(aod.shape)
unds_param[unds]=1  #undefined types present = 1, not present = 0
qcmax=np.nanmax(qc,axis=0)  #QC Flag values, 

#Valid Values based on CATS QC
valid5 = (catextmax < 1.25) & (qcmax == 0) & (sky == 1) & \
         (cloud_param < 1) & (unds_param < 1) & (catextunmax < 10) & \
         (catcadmax < -1) & (catcadmin > -11) & (dnfcat != 1)
aod[~valid5]=-999.0

#catlat    = infile['/geolocation/CATS_Fore_FOV_Latitude'][:]
#catlon    = infile['/geolocation/CATS_Fore_FOV_Longitude'][:]
#time      = infile['/profile/Profile_UTC_Time'][:]
#date      = infile['/profile/Profile_UTC_Date'][:]
#aod       = infile['/profile/Aerosol_Optical_Depth_1064_Fore_FOV'][:]
#aer_type  = infile['/profile/Aerosol_Type_Fore_FOV'][:]
#cld_phase = infile['/profile/Cloud_Phase_Fore_FOV'][:]
#feattype  = infile['/profile/Feature_Type_Fore_FOV'][:]
#sky       = infile['/profile/Sky_Condition_Fore_FOV'][:]
#qc        = infile['/profile/Extinction_QC_Flag_1064_Fore_FOV'][:]
#extuncat  = infile['/profile/Extinction_Coefficient_Uncertainty_1064_Fore_FOV'][:]
#featscore = infile['/profile/Feature_Type_Score_Fore_FOV'][:]
#dnfcat    = infile['/profile/Day_Night_Flag'][:]
#opac      = infile['/profile/Percent_Opacity_Fore_FOV'][:]
#demsfc    = infile['/profile/DEM_Surface_Altitude_Fore_FOV'][:]
#ext       = infile['/profile/Extinction_Coefficient_1064_Fore_FOV'][:]
#bin_alt   = infile['/metadata_parameters/Bin_Altitude_Array'][:]



for ii in range(aod.shape[0]):
    plot_line = False
    if(valid5[ii] == True):
        print(ii, juldayc2[ii].strftime('%D %H:%M:%S'), catlat[ii], catlon[ii])

        for jj in range(bin_alt.shape[0]):
            if(  \
                (ext[jj,ii] < 1.25) & \
                (ext[jj,ii] > 0) & \
                (qc[jj,ii] <= 0) & \
                (sky[ii] == 1) & \
                (feattype[jj,ii] != 1) & \
                (feattype[jj,ii] != 2) & \
                (extuncat[jj,ii] < 10) & \
                (featscore[jj,ii] < -1) & \
                (featscore[jj,ii] > -11)): 

                # Now, check stratosphere stuff
                # -----------------------------
                #bin_alt_agl = bin_alt[jj]
                bin_alt_agl = bin_alt[jj] - demsfc[ii]
                if(bin_alt_agl >= 10):
                    plot_line = True
                    print("    LAYER", np.round(bin_alt_agl, 2), np.round(ext[jj,ii], 4), aer_type[jj,ii])
                    #print("    LAYER", np.round(bin_alt_agl, 2), np.round(ext[jj,ii], 4))

    if(plot_line):
        ax1.axvline(juldayc2[ii], linestyle = ':', color = 'black')
        ax2.axvline(juldayc2[ii], linestyle = ':', color = 'yellow')
        ax3.axvline(juldayc2[ii], linestyle = ':', color = 'yellow')
        ax4.axvline(juldayc2[ii], linestyle = ':', color = 'yellow')

fig.tight_layout()
titler = filename.strip().split('/')[-1].split('.')[2]
outname = 'cats_combined_data_' + titler + '.png'
fig.savefig(outname, dpi = 200)
print("Saved image", outname)
plt.show()

sys.exit()

#AGL CORRECT VERTICAL EXTINCTION PROFILES
if (aod[valid5].shape[0] > 0): 

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Now, print the profiles with aerosol
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

    # Figure out which indices contain valid AODs
    # -------------------------------------------
    aod_idxs = np.where(aod != -999.)
    num_good_profs = len(aod_idxs[0])

    # Subtract the DEM surface elevations from the CATS lidar altitudes
    # -----------------------------------------------------------------
    bin_alt_agl = \
        np.array([[talt - dsfc for dsfc in demsfc] for talt in bin_alt])

    # Mask any AGL altitudes that are below 0 or below 10
    # ---------------------------------------------------
    if(l_ONLY_STRAT):
        mask_alt_agl = np.ma.masked_where(bin_alt_agl < 10, bin_alt_agl)
    else:
        mask_alt_agl = np.ma.masked_where(bin_alt_agl < 0, bin_alt_agl)

    # Figure out the profiles that contain non-zero extinction value
    # in the desired layers
    # --------------------------------------------------------------
    mask_ext = np.ma.masked_where(ext <= 0, ext)

    test_ext = np.ma.masked_where(  (mask_ext.mask == True) | \
        (mask_alt_agl.mask == True), mask_ext)

    liner = np.array([False in test_ext[:,nn].mask \
        for nn in range(test_ext.shape[1])])

    time_idxs = np.where(liner == True)

    # See if there are any valid profiles that contain aerosols in the 
    # desired layers.
    # ----------------------------------------------------------------
    if(len(time_idxs[0]) > 0):

        # Determine the day for this file. If it's a new day, save the 
        # working file and open another one
        # ------------------------------------------------------------
        new_day = int(cats_list[ii].strip().split('/')[-1][42:44])
        if( (new_day != working_day)):

            # Only close the file object if it is not the first iteration
            # -----------------------------------------------------------
            if(working_day != -999):
                print("\nClosing current outfile")
                outfile.close()
           
            # Figure out the timestamp for this day
            # ------------------------------------- 
            working_day = new_day
            date_str = datetime.strptime(\
                cats_list[ii].strip().split('/')[-1][34:44], \
                '%Y-%m-%d').strftime('%Y%m%d')

            # Build an output filename
            # ------------------------
            if(l_ONLY_STRAT):
                strat_add = '_strat'
            else:
                strat_add = ''

            if(l_ONLY_NONZERO_LAYERS):
                aer_add = '_aer'
            else:
                aer_add = ''

            new_filename  = 'cats_aer_profiles_' + date_str + strat_add + \
                aer_add + '.txt'

            print("NEW FILENAME:", new_filename)

            # Open the file
            # -------------
            outfile = open(new_filename, 'w')


        for jj in time_idxs[0]:
            # Grab only the indices above ground, or above 10 km
            # --------------------------------------------------
            if(l_ONLY_STRAT):
                end_idx = np.where(mask_alt_agl >= 10)[0][-1] + 1
            else:
                end_idx = np.where(mask_alt_agl > 0)[0][-1] + 1

            out_str = '{0:20s} {1:8.7f} {2:9.4f} {3:9.4f} \n'.format(\
                juldayc2[jj].strftime('%Y%m%d-%H:%M:%S'), \
                time[jj], \
                float(catlat[jj]), float(catlon[jj]))
            print(out_str.strip())
            outfile.write(out_str)

            for kk in range(end_idx):
                if( (l_ONLY_NONZERO_LAYERS and ext[kk,jj] != 0.0) or (not l_ONLY_NONZERO_LAYERS)):
                    out_str = '{0:7.2f} {1:7.5f}\n'.format(float(bin_alt_agl[kk,jj]), float(ext[kk,jj]))
                    outfile.write(out_str)
