#!/usr/bin/env python
"""
  NAME:
    automate_OMI_download.py

  PURPOSE:
    Act as a script to automatically download archived OMI data from
    "https://aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMAERUV.003/"

  SYNTAX:
    python automate_OMI_download.py

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2018/11/09:
      Written

"""
import os
import sys
import numpy as np
import datetime

test = "https://ceres-tool.larc.nasa.gov/ord-tool/data1//CERES_2020-01-31:14283/dir1/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_200709-200709.nc"

base_cmnd = "wget https://ceres-tool.larc.nasa.gov/ord-tool/data1//CERES_2020-01-31:14283/dir1/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_"

years = np.arange(2000,2020)
months = np.arange(1,13)

for year in years:
    for month in months:
        date = str(year)+str(month).zfill(2)
        cmnd = base_cmnd+date+'-'+date+'.nc'
        print(cmnd)
        os.system(cmnd)

#### Find current date information
###now = datetime.datetime.now()
###c_jul_day = now.timetuple().tm_yday
###
#### DOESN'T WORK. Each call to 'wget' downloads all the data for that entire year
###
###cmnd = base_path+"2018/001"
###os.system(base_cmnd+' '+cmnd)
###
###sys.exit()

#for year in years:
#    if(year == now.year):
#        # The archive may not contain data all the way up 
#        date_max = c_jul_day-1
#        for i in range(1,date_max):
#            cmnd = base_path+str(year)+"/"+str(i).zfill(3)
#            print base_cmnd+' '+cmnd
#            os.system(base_cmnd+' '+cmnd)
#            # Execute command
#        
#            # Move newly-downloaded files to the main storage directory
#            print 'mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path
#            os.system('mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path)
#    else:
#        if(year in leap_years):
#            for i in range(1,3):
#            #for i in range(1,367):
#                cmnd = base_path+str(year)+"/"+str(i).zfill(3)
#                print base_cmnd+' '+cmnd
#                os.system(base_cmnd+' '+cmnd)
#                #print 'mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path
#                print 'mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path
#                os.system('mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path)
#        else:
#            for i in range(1,3):
#            #for i in range(1,366):
#                cmnd = base_path+str(year)+"/"+str(i).zfill(3)
#                print base_cmnd+' '+cmnd
#                os.system(base_cmnd+' '+cmnd)
#                print 'mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path
#                os.system('mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path)
#            
#    
