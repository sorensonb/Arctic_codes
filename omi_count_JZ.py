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
import numpy as np
import subprocess
import h5py
import os
import warnings
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta

warnings.simplefilter(action = "ignore", category = RuntimeWarning)

if(len(sys.argv)<2):
    print("SYNTAX: python plot_single_OMI.py date")
    #print("SYNTAX: python plot_single_OMI.py date variable row_max")
    print("\n        Accepted variable names:")
    sys.exit()

plot_time = sys.argv[1]

start_year = 2018
end_year   = 2018

# This is the path that points to the HDF5 OMI files. This must be changed
# if running on a new system.
base_path = '/home/bsorenson/data/OMI/H5_files/'


base_date = datetime(year=start_year,month=4,day=1)
end_date  = datetime(year=end_year,month=10,day=1)

# Open output file
outname = "row_anomaly_dates_" + base_date.strftime("%Y%m%d") + '_' + \
    end_date.strftime("%Y%m%d") + '.txt'

with open(outname,'w') as fout:

    get_data = True
    first_time = True
    while(base_date < end_date):
    
        plot_time = base_date.strftime("%Y%m%d")
        data_time = base_date.strftime("%Ym")
        bad_rows = []
    
        # Copy data for this year 
        if(get_data == True):
            get_data = False
            last_time = (base_date - relativedelta(years = 1)).strftime("%Ym")
            print("getting data")
            for m_idx in range(base_date.month,10):
                # Remove this month's data from last year
                if(first_time == False):
                    print('rm '+base_path+'OMI*OMAERUV_'+last_time+\
                        str(m_idx).zfill(2)+'*')
                    os.system('rm '+base_path+'OMI*OMAERUV_'+last_time+\
                        str(m_idx).zfill(2)+'*')
                print('scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/H5_files/'+\
                    'OMI*OMAERUV_'+data_time+str(m_idx).zfill(2)  + '* ' + base_path)
                os.system('scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/H5_files/'+\
                    'OMI*OMAERUV_'+data_time+str(m_idx).zfill(2)  + '* ' + base_path)
            first_time = False
    
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
            # A row is identified as "bad" if the average of that row between the 
            # 1230 and 1500 indices across each swath during the day is more than
            # 2 standard deviations away from the average of all rows between 1230 and
            # 1500 for the day..
            #
            # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        
            day_avgs = np.nanmean(total_avgs1,axis=0)
            avg_day_avg = np.nanmean(day_avgs)
            day_std  = np.nanstd(day_avgs)
        
            bad_rows = []
            for rowI in range(len(day_avgs)):
                if((day_avgs[rowI] - avg_day_avg) > (day_std * 2)):
                    bad_rows.append(rowI+1)
            #        print("Bad row number ", rowI+1) 
           
        except subprocess.CalledProcessError:
            print("ERROR: no data found for time",plot_time)
  
        outstring = plot_time + ' ' + str(len(bad_rows))
        for row in bad_rows:
            outstring = outstring + ' ' + str(row)
        outstring = outstring + '\n'
 
        # Print the date and bad row values to the output file 
        fout.write(outstring)
        base_date = base_date + timedelta(days=1)
        if(base_date.month == 10):
            get_data = True
            base_date = datetime(year = base_date.year+1,month=4,day=1)
    
    # Delete the ending data
    
    last_time = (base_date - relativedelta(years = 1)).strftime("%Ym")
    for m_idx in range(base_date.month,10):
        # Remove this month's data from last year
        if(first_time == False):
            print('rm '+base_path+'OMI*OMAERUV_'+last_time+\
                  str(m_idx).zfill(2)+'*')
            os.system('rm '+base_path+'OMI*OMAERUV_'+last_time+\
                  str(m_idx).zfill(2)+'*')

