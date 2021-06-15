#!/usr/bin/env python
"""
  NAME:

  PURPOSE:

  SYNTAX:

  EXAMPLE:

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2021/04/01:
      Written (modification of AerosolIndexPlotting_Procedure_August112014t0046.pro by rcontreras) 

"""

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color
import cartopy.crs as ccrs
import subprocess
import h5py

if(len(sys.argv)<2):
    print("SYNTAX: python plot_single_OMI.py date")
    #print("SYNTAX: python plot_single_OMI.py date variable row_max")
    print("\n        Accepted variable names:")
    sys.exit()

plot_time = sys.argv[1]

# This is the path that points to the HDF5 OMI files. This must be changed
# if running on a new system.
base_path = '/home/bsorenson/data/OMI/H5_files/'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Select the needed data files
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Extract date information to use when finding the data files
year = plot_time[:4]
date = plot_time[4:8]
if(len(plot_time)==13):
    time = plot_time[9:]
elif(len(plot_time)==12):
    time = plot_time[8:]
else:
    time = ''

total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+\
    year+'m'+date+'t'+time+'*.he5',shell=True).decode('utf-8').strip().split('\n')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Grid the data
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
whole_rows = np.zeros((len(total_list),60))
total_AI   = np.full((len(total_list),2000,60),-99.)
total_X   = np.full((len(total_list),2000,60),-99.)
for fileI in range(len(total_list)):
    # read in data directly from HDF5 files
    #print(total_list[fileI])
    data = h5py.File(total_list[fileI],'r')

    AI1     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
    XTRACK1 = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]

    #print(AI1.shape)

    # Mask values with missing AI or XTrack values that are not 0 or 4
    AI1[(AI1 < -2e5) | ((XTRACK1 != 0) & (XTRACK1 != 4))] = np.nan 

    total_AI[fileI,:AI1.shape[0],:] = AI1[:,:]

    #print(np.nanmin(total_AI[fileI,:,:]),np.min(total_AI[fileI,:,:]))
    total_X[fileI,:XTRACK1.shape[0],:] = XTRACK1

    # Calculate the averages along each row in the current swath over the 
    # Arctic.
    #avgs1 = np.array([np.nanmean(AI1[1230:1500,idx]) for idx in range(60)])

    #bad_rows = np.where(avgs1 > 1)[0]

    #print(total_list[fileI],bad_rows)

    data.close()

## Mask values with missing AI or XTrack values that are not 0 or 4
#total_AI[(total_AI[:,:,:] < -2e5) | ((total_X[:,:,:] != 0) & \
#    (total_X[:,:,:] != 4))] = np.nan 
total_AI[total_AI[:,:,:] <= -99.] = np.nan 

# Calculate the averages along each row in the current swath over the 
# Arctic. NOTE: This section should go in the file loop, but am
# too lazy to do it now.
total_avgs1 = np.array([[np.nanmean(total_AI[file_idx,1230:1500,idx]) \
    for idx in range(60)] for file_idx in range(len(total_list))])

total_stds1 = np.array([[np.nanstd(total_AI[file_idx,1230:1500,idx]) \
    for idx in range(60)] for file_idx in range(len(total_list))])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Plot the gridded data
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

day_avgs = np.nanmean(total_avgs1,axis=0)
avg_day_avg = np.nanmean(day_avgs)
day_std  = np.nanstd(day_avgs)

for fileI in range(len(day_avgs)):
    if((day_avgs[fileI] - avg_day_avg) > (day_std * 2)):
        print("Bad row at index ", fileI) 
