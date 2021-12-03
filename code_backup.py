#!/usr/bin/env python
"""


"""
import os
import sys
from datetime import datetime

base_dir = "/home/bsorenson/Research/"
dest_dir = "/home/bsorenson/Arctic_codes/"

# ---------------------------------------------------------------------------- 
# Backup the backup script
# ---------------------------------------------------------------------------- 
print("code_backup.py")
cmnd = "cp /home/bsorenson/code_backup.py "+dest_dir
print(cmnd)
os.system(cmnd)
## ---------------------------------------------------------------------------- 
## Ice analysis
## ---------------------------------------------------------------------------- 
#print("Ice_analysis")
#cmnd = "cp "+base_dir+"Ice_analysis/*.py "+dest_dir
#print(cmnd)
#os.system(cmnd)
## ---------------------------------------------------------------------------- 
## CryoSat-2 analysis
## ---------------------------------------------------------------------------- 
#print("CryoSat2")
#cmnd = "cp "+base_dir+"CryoSat2/*.py "+dest_dir
#print(cmnd)
#os.system(cmnd)
## ---------------------------------------------------------------------------- 
## ICESat-2 analysis
## ---------------------------------------------------------------------------- 
#print("ICESat2")
#cmnd = "cp "+base_dir+"ICESat2/*.py "+dest_dir
#print(cmnd)
#os.system(cmnd)
## ---------------------------------------------------------------------------- 
## CERES raw codes
## ---------------------------------------------------------------------------- 
print("CERES")
#cmnd = "cp "+base_dir+"CERES/*.py "+dest_dir
cmnd = "find "+base_dir+"CERES/ -type f -name \"*.py\" | xargs cp -t "+dest_dir
print(cmnd)
os.system(cmnd)
#cmnd = "cp "+base_dir+"CERES/*.c "+dest_dir
#print(cmnd)
#os.system(cmnd)
#cmnd = "cp "+base_dir+"CERES/Make* "+dest_dir
#print(cmnd)
#os.system(cmnd)
#cmnd = "cp "+base_dir+"CERES/README_ceres "+dest_dir
#print(cmnd)
#os.system(cmnd)
## ---------------------------------------------------------------------------- 
## CERES_Ice_comparison codes
## ---------------------------------------------------------------------------- 
#print("CERES_Ice_comparison")
#cmnd = "cp "+base_dir+"CERES_Ice_comparison/*.py "+dest_dir
#cmnd = "find "+base_dir+"CERES_Ice_comparison/ -type f -name \"*.py\" | xargs cp -t "+dest_dir
#print(cmnd)
#os.system(cmnd)
## ---------------------------------------------------------------------------- 
## PIOMAS analysis
## ---------------------------------------------------------------------------- 
#print("PIOMAS")
#cmnd = "cp "+base_dir+"PIOMAS/*.py "+dest_dir
#print(cmnd)
#os.system(cmnd)
# ---------------------------------------------------------------------------- 
# OMI Codes
# ---------------------------------------------------------------------------- 
print("OMI")
final_dir = dest_dir + 'OMI/'
cmnd = "find "+base_dir+"OMI/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
#print(cmnd)
#os.system(cmnd)
# Get codes from Raindrop
# -----------------------
rain_dir = "bsorenson@raindrop.atmos.und.edu:/home/bsorenson/OMI/"

final_dir = dest_dir + 'OMI/processing/shawn_analysis/count_analysis/'
cmnd = "scp "+rain_dir+"shawn_analysis/count_analysis/*.f90 "+final_dir
print(cmnd)
os.system(cmnd)
cmnd = "scp "+rain_dir+"shawn_analysis/count_analysis/Make* "+final_dir
print(cmnd)
os.system(cmnd)

final_dir = dest_dir + 'OMI/processing/shawn_analysis/climo_analysis/'
cmnd = "scp "+rain_dir+"shawn_analysis/climo_analysis/*.f90 "+final_dir
print(cmnd)
os.system(cmnd)
cmnd = "scp "+rain_dir+"shawn_analysis/climo_analysis/Make* "+final_dir
print(cmnd)
os.system(cmnd)

final_dir = dest_dir + 'OMI/processing/JZ_analysis/climo_analysis/'
cmnd = "scp "+rain_dir+"JZ_analysis/climo_analysis/*.f90 "+final_dir
print(cmnd)
os.system(cmnd)
cmnd = "scp "+rain_dir+"JZ_analysis/climo_analysis/Make* "+final_dir
print(cmnd)
os.system(cmnd)

final_dir = dest_dir + 'OMI/processing/JZ_analysis/count_analysis/'
cmnd = "scp "+rain_dir+"JZ_analysis/count_analysis/*.f90 "+final_dir
print(cmnd)
os.system(cmnd)
cmnd = "scp "+rain_dir+"JZ_analysis/count_analysis/Make* "+final_dir
print(cmnd)
os.system(cmnd)

final_dir = dest_dir + 'OMI/processing/JZ_analysis/drift_analysis/'
cmnd = "scp "+rain_dir+"JZ_analysis/drift_analysis/*.f90 "+final_dir
print(cmnd)
os.system(cmnd)
cmnd = "scp "+rain_dir+"JZ_analysis/drift_analysis/Make* "+final_dir
print(cmnd)
os.system(cmnd)

final_dir = dest_dir + 'OMI/processing/JZ_analysis/row_analysis/'
cmnd = "scp "+rain_dir+"JZ_analysis/row_analysis/*.f90 "+final_dir
print(cmnd)
os.system(cmnd)
cmnd = "scp "+rain_dir+"JZ_analysis/row_analysis/Make* "+final_dir
print(cmnd)
os.system(cmnd)

final_dir = dest_dir + 'OMI/processing/JZ_analysis/JZ_lib/'
cmnd = "scp "+rain_dir+"JZ_analysis/JZ_lib/*.f90 "+final_dir
print(cmnd)
os.system(cmnd)

final_dir = dest_dir + 'OMI/processing/fort_lib/'
cmnd = "scp "+rain_dir+"/fort_lib/*.f90 "+final_dir
print(cmnd)
os.system(cmnd)

#cmnd = "cp "+base_dir+"OMI/OMI_simulation/*.py "+dest_dir
#print(cmnd)
#os.system(cmnd)
#cmnd = "cp "+base_dir+"OMI/README_omi "+dest_dir
#print(cmnd)
#os.system(cmnd)


##!#cmnd = "cp "+base_dir+"OMI/*.pro "+dest_dir
##!#print(cmnd)
##!#os.system(cmnd)
##!#cmnd = "cp "+base_dir+"OMI/*.idl "+dest_dir
##!#print(cmnd)
##!#os.system(cmnd)
##!#cmnd = "cp "+base_dir+"OMI/run_Average_AerosolIndexCalculator "+dest_dir
##!#print(cmnd)
##!#os.system(cmnd)
# ---------------------------------------------------------------------------- 
# MISR analysis
# ---------------------------------------------------------------------------- 
print("MISR")
cmnd = "cp "+base_dir+"MISR/*.py "+dest_dir
print(cmnd)
os.system(cmnd)
cmnd = "cp "+base_dir+"MISR/*.c "+dest_dir
print(cmnd)
os.system(cmnd)

## Get codes from JPSS
#jpss_dir = "blake.sorenson@jpss.atmos.und.edu:/home/blake.sorenson/MISR/"
#cmnd = "scp "+jpss_dir+"*.c "+dest_dir
#print(cmnd)
#os.system(cmnd)
#cmnd = "scp "+jpss_dir+"*.py "+dest_dir
#print(cmnd)
#os.system(cmnd)

# ---------------------------------------------------------------------------- 
# MODIS analysis
# ---------------------------------------------------------------------------- 
print("MODIS")
final_dir = dest_dir + 'MODIS/'
cmnd = "find "+base_dir+"MODIS/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)
# ---------------------------------------------------------------------------- 
# AIRS analysis
# ---------------------------------------------------------------------------- 
print("AIRS")
final_dir = dest_dir + 'AIRS/'
cmnd = "find "+base_dir+"AIRS/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

cmnd = "find /home/bsorenson/data/AIRS/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

# ---------------------------------------------------------------------------- 
# SBDART analysis
# ---------------------------------------------------------------------------- 
print("SBDART")
work_dir = 'MODIS/obs_smoke_forcing/SBDART/'

final_dir = dest_dir + 'SBDART/'
cmnd = "find "+base_dir+ work_dir + \
    " -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

## ---------------------------------------------------------------------------- 
## Siphon/Metpy Codes
## ---------------------------------------------------------------------------- 
#print("Siphon/Metpy")
#smcodes = '/home/bsorenson/Programs/Python/auto_weather_codes/'
#cmnd = "cp "+smcodes+"*.py "+dest_dir
#print(cmnd)
#os.system(cmnd)
### ---------------------------------------------------------------------------- 
### Jianglong/Jeff's Mie stuff
### ---------------------------------------------------------------------------- 
##print("Mie stuff")
##smcodes = '/home/bsorenson/Research/Jianglong_stuff/'
##cmnd = "cp "+smcodes+"*.pro "+dest_dir
##cmnd = "cp "+smcodes+"run_mie_code "+dest_dir
##print(cmnd)
##os.system(cmnd)

# Automate the uploading to Github
# --------------------------------

# Change to Arctic codes directory
os.chdir(dest_dir)
#os.chdir('/home/bsorenson/Arctic_codes/')

# Add new stuff
os.system('git add .')

# Determine today's date for the command
today_str = datetime.today().strftime('%Y/%m/%d')
cmnd = 'git commit -m \"Backup '+today_str + '\"'
os.system(cmnd)

# Push
os.system('git push origin master')
