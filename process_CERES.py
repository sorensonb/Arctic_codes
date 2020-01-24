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
exec_name = 'autoceres_daily_exec'
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

#for year in all_years[int(sys.argv[1])]:
year='2008'
month='04'
tdays = np.arange(2,31)
days = [str(int(tday)).zfill(2) for tday in tdays]
#days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14',]
for day in days:
    print(base_path+exec_name+' '+year+' '+month+' '+day)
    os.system(base_path+exec_name+' '+year+' '+month+' '+day)

##!#if(int(sys.argv[1])<3):
##!#    
##!#    all_years = np.array_split(np.array(['2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015']),3)
##!#    months = ['01','02','03','04','05','06','07','08','09','10','11','12']
##!#    #months = ['01','02','03','04','05','06','07','08','09','10','11','12']
##!#    #year='2007'
##!#    for year in all_years[int(sys.argv[1])]:
##!#        for month in months:
##!#            print(base_path+exec_name+' '+year+' '+month)
##!#            os.system(base_path+exec_name+' '+year+' '+month)
##!#else:
##!#    if(int(sys.argv[1])==3):
##!#        dates = [('2007','06'),('2007','07'),('2007','08'),('2007','09'),         \
##!#                 ('2007','10'),('2007','11'),('2007','12'),                       \
##!#                 ('2008','01'),('2008','02'),('2008','03'),('2008','04'),         \
##!#                 ('2008','05'),('2008','06'),('2008','07'),('2008','08'),         \
##!#                 ('2008','09'),('2008','10'),('2008','11'),('2008','12'),         \
##!#                 ('2009','01'),('2009','02'),('2009','03'),('2009','04'),         \
##!#                 ('2009','05'),('2009','06'),('2009','07'),('2009','08'),         \
##!#                 ('2009','09'),('2009','10'),('2009','11'),('2009','12'),         \
##!#                 ('2010','01'),('2010','02'),('2010','03'),('2010','04'),         \
##!#                 ('2010','05'),('2010','06'),('2010','07'),('2010','08'),         \
##!#                 ('2010','09'),('2010','10'),('2010','11'),('2010','12')]
##!#    elif(int(sys.argv[1])==4):
##!#        dates = [('2012','07'),('2012','08'),('2012','09'),('2012','10'),         \
##!#                 ('2012','11'),('2012','12'),                                     \
##!#                 ('2013','01'),('2013','02'),('2013','03'),('2013','04'),         \
##!#                 ('2013','05'),('2013','06'),('2013','07'),('2013','08'),         \
##!#                 ('2013','09'),('2013','10'),('2013','11'),('2013','12'),         \
##!#                 ('2014','01'),('2014','02'),('2014','03'),('2014','04'),         \
##!#                 ('2014','05'),('2014','06'),('2014','07'),('2014','08'),         \
##!#                 ('2014','09'),('2014','10'),('2014','11'),('2014','12'),         \
##!#                 ('2015','01'),('2015','02'),('2015','03'),('2015','04'),         \
##!#                 ('2015','05'),('2015','06'),('2015','07'),('2015','08'),         \
##!#                 ('2015','09'),('2015','10'),('2015','11'),('2015','12')]
##!#    for date in dates:
##!#        print(base_path+exec_name+' '+date[0]+' '+date[1])
##!#        os.system(base_path+exec_name+' '+date[0]+' '+date[1])
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
