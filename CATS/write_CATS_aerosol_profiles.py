#!/usr/bin/env python

"""
  NAME:
    write_CATS_aerosol_profiles.py

  PURPOSE:
    Prints all CATS aerosol profiles in a day to a single text file.
    Generates files formatted following:

        cats_aer_profiles_YYYYMMDD.txt:
            - All layers (troposphere and stratosphere) both with and 
              without aerosol are printed. NOTE: these files will be 
              very, very large. 
            - Generated when write_aerosol_profiles.py is run with:
                - l_ONLY_STRAT = False
                - l_ONLY_NONZERO_LAYERS = False
        
        cats_aer_profiles_YYYYMMDD_strat.txt
            - All stratosphere profile layers are printed if there are 
              aerosols in ANY CATS profile layer (even if the aerosols 
              are in the troposphere). In this format, many of the printed 
              profiles will contain only 0-value extinctions in the 
              stratosphere because the only aerosols in the proflie
              were in the troposphere.
            - Generated when write_aerosol_profiles.py is run with:
                - l_ONLY_STRAT = True
                - l_ONLY_NONZERO_LAYERS = False
        
        cats_aer_profiles_YYYYMMDD_aer.txt:
            - Only aerosol layers anywhere in the column (troposphere or 
              stratosphere) are printed.
            - Generated when write_aerosol_profiles.py is run with:
                - l_ONLY_STRAT = False
                - l_ONLY_NONZERO_LAYERS = True
        
        cats_aer_profiles_YYYYMMDD_strat_aer.txt:
            - Only aerosol layers within the stratosphere are printed. 
              This type results in the smallest output files. 
            - Generated when write_aerosol_profiles.py is run with:
                - l_ONLY_STRAT = True
                - l_ONLY_NONZERO_LAYERS = True


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
import sys

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

#Read in data files
#cats_init = glob.glob('/data/CATS_V3/L2_Profiles/*/*.hdf5')
cats_init = glob('/home/WRF/CATS/Profile/2017/*.hdf5')

# Extract the date information from each file name, storing them inside
# datetime objects
# ----------------------------------------------------------------------
file_dates = np.array([datetime.strptime(tfile.strip().split('/')[-1][34:44], \
                       '%Y-%m-%d') for tfile in cats_init])

# Extract the months and years from each datetime object
# ------------------------------------------------------
file_months = np.array([fdate.month for fdate in file_dates])
file_years  = np.array([fdate.year for fdate in file_dates])

# Determine which files are from May - October of 2017
# ----------------------------------------------------
good_idxs = np.where( (file_years == 2017) & ( (file_months >= 5) & \
        (file_months <= 10)))

# Extract only those files from the overall CATS file list
# --------------------------------------------------------
cats_list = sorted(np.array(cats_init)[good_idxs])

working_day = -999
for ii in range(0,len(cats_list)):   #Loop through list of CATS files

    #Get variable data
    filename2=cats_list[ii]
    print(filename2)

    all_sum=np.zeros((200,4))
    all_count=np.zeros((200,4))


    infile=h5py.File(filename2,'r')

    # Read variables from the HDF5 file
    # ---------------------------------
    catlat    = infile['/geolocation/CATS_Fore_FOV_Latitude'][:]
    catlon    = infile['/geolocation/CATS_Fore_FOV_Longitude'][:]
    time      = infile['/profile/Profile_UTC_Time'][:]
    date      = infile['/profile/Profile_UTC_Date'][:]
    aod       = infile['/profile/Aerosol_Optical_Depth_1064_Fore_FOV'][:]
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
