#!/usr/bin/env python
"""
  NAME:

  PURPOSE:

  SYNTAX:

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2018/11/09:
      Written

"""
import os
import sys
import numpy as np
import datetime

if(len(sys.argv)!=2):
    print("Syntax: ./automate_CryoSat2_download.py file_list_file")
    print("         file_list_file = /home/bsorenson/Research/CryoSat2/nsidc-download_e19eb67b.txt")
    sys.exit()

base_cmnd = "curl -b ~/.urs_cookies -c ~/.urs_cookies -L -n -O "

##!## For downloading directly from the website
##!#base_url = "https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0051.001/"
##!#years = np.arange(1980,2000)
##!##months = np.arange(6,9)
##!#months = np.array([4,9])
##!#for year in years:
##!#    for month in months:
##!#        total_cmnd = base_cmnd + base_url + str(year)+'.'+str(month).zfill(2)+'.01/nt_'+str(year)+str(month).zfill(2)+'_f13_v1.1_n.bin'
##!#        print(total_cmnd)
##!#        os.system(total_cmnd)

# Open up  and read in file file
infile = sys.argv[1]
with open(infile,'r') as f:
    for line in f:
        templine = line.strip()
        if(templine[-2:]=='nc'):
            total_cmnd = base_cmnd+templine
            print(total_cmnd)
            os.system(total_cmnd)

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
