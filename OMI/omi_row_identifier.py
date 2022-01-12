#!/usr/bin/env python
"""
  NAME:
    omi_row_identifier.py

  PURPOSE:
    Identify row anomaly-affected rows that are not flagged by the XTRACK QC    
    flag. The identified bad row numbers (1-based) for each day are printed to
    a text file. Also identifies flagged bad rows and writes them to an output
    file. 

    PLEASE READ:

    Note that this is currently set up to automatically download the needed
    OMI files at each year in the loop. This was done to avoid overloading my
    disk space with all the OMI data at once. If this is being run on 
    Raindrop, where the data are already located, that step is unnecessary.
    An associated section at the bottom of the code which deletes the 
    downloaded data is not needed if being run on Raindrop. These two sections
    are commented out for now.

  SYNTAX:
    $ python omi_row_identifier.py

  EXAMPLE:

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2021/04/01:
      Written (modification of AerosolIndexPlotting_Procedure_August112014t0046.pro by rcontreras) 

"""

import sys
import numpy as np
import subprocess
import h5py
import os
import warnings
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta

# Returns the actual flag value from the overall Xtrack integer value
def get_xtrack_flags(value, flag_idx = None):
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

warnings.simplefilter(action = "ignore", category = RuntimeWarning)

if(len(sys.argv)<2):
    print("SYNTAX: python omi_row_identifier.py date")
    sys.exit()

start_year = 2005
end_year   = 2020

# This is the path that points to the HDF5 OMI files. This must be changed
# if running on a new system.
base_path = '/home/bsorenson/data/OMI/H5_files/'

base_date = datetime(year=start_year,month=4,day=1)
end_date  = datetime(year=end_year,month=10,day=1)

# Open additional bad row output file
# -----------------------------------
outname = "row_anomaly_dates_" + base_date.strftime("%Y%m%d") + '_' + \
    end_date.strftime("%Y%m%d") + '.txt'

# Open XTRACK row output file
# ---------------------------
outname_xtrack = "row_anomaly_xtrack_dates_" + base_date.strftime("%Y%m%d") + '_' + \
    end_date.strftime("%Y%m%d") + '.txt'

with open(outname,'w') as fout:
    with open(outname_xtrack, 'w') as fout2:

        get_data = True
        first_time = True

        # Keep looping until the end date is reached
        # ------------------------------------------
        while(base_date < end_date):
        
            plot_time = base_date.strftime("%Y%m%d")
            data_time = base_date.strftime("%Ym")

            bad_rows = []
            xtrack_rows = []
        
            ##!## Copy the working year's OMI data from Raindrop to the local
            ##!## computer.
            ##!## NOTE: If this code is being run on Raindrop, this is not
            ##!## necessary.
            ##!## -----------------------------------------------------------
            ##!#if(get_data == True):
            ##!#    get_data = False
            ##!#    last_time = (base_date - relativedelta(years = 1)).strftime("%Ym")
            ##!#    print("getting data")
            ##!#    for m_idx in range(base_date.month,10):
            ##!#        # Remove this month's data from last year
            ##!#        if(first_time == False):
            ##!#            print('rm '+base_path+'OMI*OMAERUV_'+last_time+\
            ##!#                str(m_idx).zfill(2)+'*')
            ##!#            os.system('rm '+base_path+'OMI*OMAERUV_'+last_time+\
            ##!#                str(m_idx).zfill(2)+'*')
            ##!#        print('scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/H5_files/'+\
            ##!#            'OMI*OMAERUV_'+data_time+str(m_idx).zfill(2)  + '* ' + base_path)
            ##!#        os.system('scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/H5_files/'+\
            ##!#            'OMI*OMAERUV_'+data_time+str(m_idx).zfill(2)  + '* ' + base_path)
            ##!#    first_time = False
        
            print(base_date.strftime("%Y%m%d"))
            
            
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
        
            try:
                total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+\
                    year+'m'+date+'t'+time+'*.he5',shell=True).decode('utf-8').strip().split('\n')
                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
                #
                # Grid the data
                #
                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
                whole_rows = np.zeros((len(total_list),60))
                total_AI   = np.full((len(total_list),2000,60),np.nan)
                total_X   = np.full((len(total_list),2000,60),np.nan)
                #total_AI   = np.full((len(total_list),2000,60),-99.)
                #total_X   = np.full((len(total_list),2000,60),-99.)
                for fileI in range(len(total_list)):
                    # read in data directly from HDF5 files
                    #print(total_list[fileI])
                    data = h5py.File(total_list[fileI],'r')
        
                    AI1     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
                    XTRACK1 = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
                    XTRK_decode = np.array([[get_xtrack_flags(val) for val in XTRKrow] for XTRKrow in XTRACK1])
                    #print(AI1.shape)
            
                    # Mask values with missing AI or XTrack values that are not 0 or 4
                    # NEW: After decoding the XTRACK values, only see where the row anomaly
                    #      flag isequal to zero
                    AI1[(AI1 < -2e5) | ((XTRACK1 != 0) & (XTRACK1 != 4))] = np.nan 
        
                    total_AI[fileI,:AI1.shape[0],:] = AI1[:,:]
        
                    # NEW: After decoding the XTRACK values, only see where the row anomaly
                    #      flag is equal to zero. Ignoring the other bit values
                    total_X[fileI,:XTRACK1.shape[0],:] = XTRK_decode[:,:,0]
            
                    data.close()
      
                # Calculate the averages along each row in the current swath over the 
                # Arctic. NOTE: This section should go in the file loop, but am
                # too lazy to do it now.
                # --------------------------------------------------------------------
                total_avgs1 = np.array([[np.nanmean(total_AI[file_idx,1230:1500,idx]) \
                    for idx in range(60)] for file_idx in range(len(total_list))])

                # Calculate the average xtrack value along each row in the current 
                # swath over the Arctic. NOTE: This section should go in the file loop,
                # but am too lazy to do it now. This is used for identifying known
                # contaminted rows.
                # --------------------------------------------------------------------
                total_avgs_X1 = np.array([[np.nanmean(total_X[file_idx,1230:1500,idx]) \
                    for idx in range(60)] for file_idx in range(len(total_list))])
            
                total_stds1 = np.array([[np.nanstd(total_AI[file_idx,1230:1500,idx]) \
                    for idx in range(60)] for file_idx in range(len(total_list))])
            
                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
                #
                # A row is identified as "bad" if the average of that row between the 
                # 1230 and 1500 indices across each swath during the day is more than
                # 2 standard deviations away from the average of all rows between 1230 and
                # 1500 for the day..
                #
                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
            
                day_avgs = np.nanmean(total_avgs1,axis=0)
                avg_day_avg = np.nanmean(day_avgs)
                day_std  = np.nanstd(day_avgs)

                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
                #
                # A row is identified as "contaminated" if it has at least 1 pixel within
                # indices 1230 and 1500 with an Xtrack QC value that is not zero, meaning
                # that it is not perfectly clean as defined by the flag. This is determined
                # by taking the average of all the Xtrack QC values along indices 1230 to
                # 1500 of each row, and if that averag eis not equal to 0, it has a non-zero
                # pixel inside it and is therefore contaminated.
                #
                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
                day_avgs_X = np.nanmean(total_avgs_X1,axis=0)
            
                # Deal with xtrack rows
                # If a row has any XTRACK flags set, add to xtrack row list
                # ---------------------------------------------------------
                bad_rows = []
                xtrack_rows = []
                for rowI in range(len(day_avgs)):
                    if((day_avgs[rowI] - avg_day_avg) > (day_std * 2)):
                        bad_rows.append(rowI+1)
                    # If a row has any XTRACK flags set, add to xtrack row list
                    if(day_avgs_X[rowI] != 0):
                        xtrack_rows.append(rowI+1)
                #        print("Bad row number ", rowI+1) 
                
                
               
            except subprocess.CalledProcessError:
                print("ERROR: no data found for time",plot_time)
  
            # Prep bad row strings
            # --------------------
            outstring = plot_time + ' ' + str(len(bad_rows))
            for row in bad_rows:
                outstring = outstring + ' ' + str(row)
            outstring = outstring + '\n'

            # Prep xtrack row strings
            # -----------------------
            outstring_X = plot_time + ' ' + str(len(xtrack_rows))
            for row in xtrack_rows:
                outstring_X = outstring_X + ' ' + str(row)
            outstring_X = outstring_X + '\n'
 
            # Print the date and bad row values to the output file 
            fout2.write(outstring_X)

            # Print the date and bad row values to the output file 
            fout.write(outstring)
            base_date = base_date + timedelta(days=1)
            if(base_date.month == 10):
                get_data = True
                base_date = datetime(year = base_date.year+1,month=4,day=1)
        
        ##!## Delete the ending data
        ##!## NOTE: IF this is being run on Raindrop, this MUST be removed
        ##!##       ensure that the OMI data on Raindrop are not deleted.
        ##!## ------------------------------------------------------------
        ##!#last_time = (base_date - relativedelta(years = 1)).strftime("%Ym")
        ##!#for m_idx in range(base_date.month,10):
        ##!#    # Remove this month's data from last year
        ##!#    if(first_time == False):
        ##!#        print('rm '+base_path+'OMI*OMAERUV_'+last_time+\
        ##!#              str(m_idx).zfill(2)+'*')
        ##!#        os.system('rm '+base_path+'OMI*OMAERUV_'+last_time+\
        ##!#              str(m_idx).zfill(2)+'*')

