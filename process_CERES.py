#!/usr/bin/env python
"""
  NAME:

  PURPOSE:

  PYTHON VERSION:

  MODULES:
    
"""

import numpy as np
import sys
import os
import subprocess

if(len(sys.argv)!=2):
    print("SYNTAX: ./process_CERES.py year_block")
    print("        0 : 2001-2005")
    print("        1 : 2006-2010")
    print("        2 : 2011-2015")
    print("        3 : custom")
    sys.exit()

base_path = '/home/bsorenson/HighLatitudeStudy/CERES/'
exec_name = 'autoceres_exec'
#exec_name = 'misr_exec'
#avg_path = '30to90/'
#avg_base_name = base_path+avg_path+'misr_avg_aod_newData'
#directories = np.array(subprocess.check_output('ls '+data_path,shell=True).decode('utf-8').strip().split('\n'))
#temp_date = directories[0]

# Compile program
#print('gcc -o misr_exec /home/bsorenson/HighLatitudeStudy/MISR/Updated_Data/netcdf_helper_newest.c -I/usr/include/ -lnetcdf')
#os.system('gcc -o misr_exec /home/bsorenson/HighLatitudeStudy/MISR/Updated_Data/netcdf_helper_newest.c -I/usr/include/ -lnetcdf')
#print('gcc -o '+exec_name+' /home/bsorenson/HighLatitudeStudy/MISR/Updated_Data/misr_auto_processor.c -I/usr/include/ -lnetcdf -lm')
#os.system('gcc -o '+exec_name+' /home/bsorenson/HighLatitudeStudy/MISR/Updated_Data/misr_auto_processor.c -I/usr/include/ -lnetcdf -lm')

if(int(sys.argv[1])!=3):
    
    all_years = np.array_split(np.array(['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015']),3)
    months = ['01','02','03','04','05','06','07','08','09','10','11','12']
    #months = ['01','02','03','04','05','06','07','08','09','10','11','12']
    #year='2001'
    for year in all_years[int(sys.argv[1])]:
        for month in months:
            print(base_path+exec_name+' '+year+' '+month)
            os.system(base_path+exec_name+' '+year+' '+month)
else:
    dates = [('2003','09'),('2003','10'),('2003','11'),('2003','12'),         \
             ('2004','01'),('2004','02'),('2004','03'),('2004','04'),         \
             ('2004','05'),('2004','06'),('2004','07'),('2004','08'),         \
             ('2004','09'),('2004','10'),('2004','11'),('2004','12'),         \
             ('2005','01'),('2005','02'),('2005','03')]
    for date in dates:
        print(base_path+exec_name+' '+date[0]+' '+date[1])
        os.system(base_path+exec_name+' '+date[0]+' '+date[1])
   # print('gzip '+base_path+'30to90/aod_no_flags/misr_avg_aod_newData_'+year+month+'_noraw.txt')
   # os.system('gzip '+base_path+'30to90/aod_no_flags/misr_avg_aod_newData_'+year+month+'_noraw.txt')

        

#files = subprocess.check_output('ls '+data_path+temp_date+'/',shell=True).decode('utf-8').strip().split('\n')
#print(data_path+directories[0]+'/'+files[0])

####max_index = int(np.where(directories=='2002.01.01')[0][0])
####for dir in directories[32:max_index]:
####    # Make the average file name based off of the directory name
####    year=dir[:4]    
####    month=dir[5:7]
####    avg_whole_name = avg_base_name+'_'+year+month+'.txt'
####    files = subprocess.check_output('ls '+data_path+dir+'/',shell=True).decode('utf-8').strip().split('\n')
####    for file in files:
####        data_file_name = data_path+dir+'/'+file
####        #print(base_path+exec_name+' '+data_file_name+' '+avg_whole_name)
####        os.system(base_path+exec_name+' '+data_file_name+' '+avg_whole_name)
    

#data_file_name = data_path+directories[0]+'/'+files[0]
#print(exec_name+' '+data_file_name)
#os.system(exec_name+' '+data_file_name)
