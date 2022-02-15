#!/usr/bin/env python
"""
  NAME:
    automate_GOES17_download.py

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
from datetime import datetime

if(len(sys.argv) != 2):
    print("SYNTAX: ./automate_GOES17_download.py YYYYMMDDHHMM")
    sys.exit()

date_str = sys.argv[1]
dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

#years = [2008]
#years = [2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]
#leap_years = [2004,2008,2012,2016,2020,2024]
base_path = "https://noaa-goes17.s3.amazonaws.com/index.html#ABI-L1b-RadC/"
base_file = 'OR_ABI-L1b-RadC-M6C'
#base_file = 'OR_ABI-L1b-RadC-M6C10_G17_s20212020601177_e20212020603561_c20212020603590.nc'
base_cmnd = "wget "
#base_cmnd = "wget --no-clobber --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies -np -r -A .he5 "
main_storage_path = "/home/bsorenson/data/GOES/goes17_abi/"
#local_storage_path = "aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMAERUV.003/"

###for year in years:
###    print(base_cmnd+' '+base_path+str(year)+'/')
###    os.system(base_cmnd+' '+base_path+str(year)+'/')

bands = np.array([2,6,7,8,9,10,13])

for band in bands:
    total_cmnd = base_cmnd + base_path + \
        dt_date_str.strftime('%Y/%j/%H/') + base_file + \
        str(band).zfill(2) + '_G17_s' + dt_date_str.strftime('%Y%m%d%H%M') + '*'
    print(total_cmnd)

sys.exit()
    

# Find current date information
#now = datetime.datetime.now()
#c_jul_day = now.timetuple().tm_yday

#### DOESN'T WORK. Each call to 'wget' downloads all the data for that entire year
###
###cmnd = base_path+"2018/001"
###os.system(base_cmnd+' '+cmnd)
###
###sys.exit()

for year in years:
    if(year == now.year):
        # The archive may not contain data all the way up 
        date_max = c_jul_day-1
        for i in range(1,date_max):
            cmnd = base_path+str(year)+"/"+str(i).zfill(3)
            print(base_cmnd+' '+cmnd)
            #os.system(base_cmnd+' '+cmnd)
            # Execute command
        
            # Move newly-downloaded files to the main storage directory
            print('mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path)
            #os.system('mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path)
    else:
        if(year in leap_years):
            #for i in range(1,3):
            for i in range(1,367):
                cmnd = base_path+str(year)+"/"+str(i).zfill(3)
                print(base_cmnd+' '+cmnd)
                #os.system(base_cmnd+' '+cmnd)
                #print 'mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path
                print('mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path)
                #os.system('mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path)
        else:
            #for i in range(1,3):
            for i in range(182,184):
            #for i in range(182,244):
            #for i in range(1,366):
                cmnd = base_path+str(year)+"/"+str(i).zfill(3)
                print(base_cmnd+' '+cmnd)
                os.system(base_cmnd+' '+cmnd)
                print('mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path)
                os.system('mv '+local_storage_path+str(year)+"/"+str(i).zfill(3)+"/*.he5 "+main_storage_path)
            
    
