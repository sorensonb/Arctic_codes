#!/usr/bin/env python
"""


"""
import os
import sys
from datetime import datetime


home_dir = os.environ['HOME']

base_dir = home_dir + "/Research/"
dest_dir = home_dir + "/Arctic_codes/"

copy_raindrop = False
copy_calipso  = False
copy_talon    = False

## Add new stuff
#os.chdir(dest_dir)
#os.system('git add .')
#
## Determine today's date for the command
#today_str = datetime.today().strftime('%Y/%m/%d')
#cmnd = 'git commit -m \"Backup '+today_str + '\"'
#os.system(cmnd)
#
## Push
#os.system('git push origin master')
#sys.exit()

# ---------------------------------------------------------------------------- 
# Backup the backup script
# ---------------------------------------------------------------------------- 
print("Home directory python codes")
final_dir = dest_dir 
cmnd = "find "+ home_dir + "/*.py | xargs cp -t "+final_dir
print(cmnd)
os.system(cmnd)

#print("code_backup.py")
#cmnd = "cp home_dir + /code_backup.py "+dest_dir
#print("python_lib.py")
#cmnd = "cp home_dir + /python_lib.py "+dest_dir
#print("bracketology.py")
#cmnd = "cp home_dir + /bracketology.py "+dest_dir
#print(cmnd)
#os.system(cmnd)

# ---------------------------------------------------------------------------- 
# Backup the job search stuff
# ---------------------------------------------------------------------------- 
print("Job search")
final_dir = dest_dir + 'job_search/'
cmnd = "find "+ home_dir + "/miscellaneous/job_search/*.py | xargs cp -t "+final_dir
print(cmnd)
os.system(cmnd)

cmnd = "find "+ home_dir + "/miscellaneous/job_search/*.csv | xargs cp -t "+final_dir
print(cmnd)
os.system(cmnd)


# ---------------------------------------------------------------------------- 
# Backup CSCI final project stuff
# ---------------------------------------------------------------------------- 
print("CSCI")
final_dir = dest_dir + 'CSCI/CSCI_543/jpss_work/'
cmnd = "find "+ home_dir + "/CSCI/final_project/jpss_work/ -type f -name \"*.py\" | xargs cp -t "+final_dir
print(cmnd)
os.system(cmnd)

final_dir = dest_dir + 'CSCI/CSCI_544/'
cmnd = "cp "+ home_dir + "/CSCI/CSCI_544/final_project/*.py "+final_dir
print(cmnd)
os.system(cmnd)

if(copy_talon):
    talon_dir = "blake.sorenson@134.129.128.241:/home/blake.sorenson/CSCI_544/final_project/"
    #rain_dir = "bsorenson@raindrop.atmos.und.edu:/home/bsorenson/OMI/"

    final_dir = dest_dir + 'CSCI/CSCI_544/jpss_work/'

    cmnd = "scp " + talon_dir + "*.py "+final_dir
    print(cmnd)
    os.system(cmnd)

    cmnd = "scp " + talon_dir + "CSCI/CSCI_544/slurm.sh "+final_dir
    print(cmnd)
    os.system(cmnd)

# ---------------------------------------------------------------------------- 
# Backup USIP stuff
# ---------------------------------------------------------------------------- 
print("USIP")
final_dir = dest_dir + 'USIP/'
cmnd = "find "+base_dir+"thermosonde/ -type f -name \"*.py\" | xargs cp -t "+final_dir
print(cmnd)
os.system(cmnd)
# ---------------------------------------------------------------------------- 
# OMI Codes
# ---------------------------------------------------------------------------- 
print("Arctic_compares")
final_dir = dest_dir + 'Arctic_compares/'
cmnd = "find "+base_dir+"Arctic_compares/ -type f -name \"*.py\" | xargs cp -t "+final_dir
print(cmnd)
os.system(cmnd)
cmnd = "find "+base_dir+"Arctic_compares/ -type f -name \"*.txt\" | xargs cp -t "+final_dir
print(cmnd)
os.system(cmnd)

if(copy_calipso):
    calip_dir = "blake.sorenson@134.129.222.8:/home/blake.sorenson/Research/Arctic_compares/"

    final_dir = dest_dir + 'Arctic_compares/calipso_processing/'

    cmnd = "scp " + calip_dir + "*.py "+final_dir
    print(cmnd)
    os.system(cmnd)
    cmnd = "scp " + calip_dir + "*.txt "+final_dir
    print(cmnd)
    os.system(cmnd)

if(copy_talon):
    talon_dir = "blake.sorenson@134.129.128.241:/home/blake.sorenson/OMI/force_efficiency_calc/"
    #rain_dir = "bsorenson@raindrop.atmos.und.edu:/home/bsorenson/OMI/"

    final_dir = dest_dir + 'Arctic_compares/talon_processing/force_efficiency_calc/'

    cmnd = "scp " + talon_dir + "*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)

    cmnd = "scp " + talon_dir + "slurm.sh "+final_dir
    print(cmnd)
    os.system(cmnd)

    cmnd = "scp " + talon_dir + "Make* "+final_dir
    print(cmnd)
    os.system(cmnd)

    talon_dir = "blake.sorenson@134.129.128.241:/home/blake.sorenson/OMI/arctic_comp/"
    final_dir = dest_dir + 'Arctic_compares/talon_processing/arctic_comp/'
    cmnd = "scp " + talon_dir + "*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)

    cmnd = "scp " + talon_dir + "Make* "+final_dir
    print(cmnd)
    os.system(cmnd)

    talon_dir = "blake.sorenson@134.129.128.241:/home/blake.sorenson/OMI/colocate_lib/"

    final_dir = dest_dir + 'Arctic_compares/talon_processing/colocate_lib/'
    cmnd = "scp " + talon_dir + "*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)

    cmnd = "scp " + talon_dir + "Make* "+final_dir
    print(cmnd)
    os.system(cmnd)


    talon_dir = "blake.sorenson@134.129.128.241:/home/blake.sorenson/OMI/ai_force_eff_test/"
    #rain_dir = "bsorenson@raindrop.atmos.und.edu:/home/bsorenson/OMI/"

    final_dir = dest_dir + 'Arctic_compares/talon_processing/ai_force_efficiency/'

    cmnd = "scp " + talon_dir + "*.py "+final_dir
    print(cmnd)
    os.system(cmnd)

    cmnd = "scp " + talon_dir + "tested_models/*.keras "+final_dir
    print(cmnd)
    os.system(cmnd)

    cmnd = "scp " + talon_dir + "tested_models/min_max*.json "+final_dir
    print(cmnd)
    os.system(cmnd)

    cmnd = "scp " + talon_dir + "slurm.sh "+final_dir
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
# CATS analysis
# ---------------------------------------------------------------------------- 
print("CATS")
final_dir = dest_dir + 'CATS/'
cmnd = "find "+base_dir+"CATS/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)
# ---------------------------------------------------------------------------- 
# CRIS analysis
# ---------------------------------------------------------------------------- 
print("CrIS")
final_dir = dest_dir + 'CrIS/'
cmnd = "find "+base_dir+"CrIS/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)
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

cmnd = "cp "+base_dir+"OMI/README* "+dest_dir
print(cmnd)
os.system(cmnd)

# Get codes from Raindrop
# -----------------------

if(copy_raindrop):
    rain_dir = "bsorenson@134.129.222.68:/home/bsorenson/OMI/"
    #rain_dir = "bsorenson@raindrop.atmos.und.edu:/home/bsorenson/OMI/"

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
    
    final_dir = dest_dir + 'OMI/processing/shawn_analysis/daily_gridder/'
    cmnd = "scp "+rain_dir+"shawn_analysis/daily_gridder/*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)
    cmnd = "scp "+rain_dir+"shawn_analysis/daily_gridder/Make* "+final_dir
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
    
    final_dir = dest_dir + 'OMI/processing/JZ_analysis/daily_gridder/'
    cmnd = "scp "+rain_dir+"JZ_analysis/daily_gridder/*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)
    cmnd = "scp "+rain_dir+"JZ_analysis/daily_gridder/Make* "+final_dir
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
    
    final_dir = dest_dir + 'Arctic_compares/processing/'
    cmnd = "scp "+rain_dir+"arctic_comp/*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)
    cmnd = "scp "+rain_dir+"arctic_comp/Make* "+final_dir
    print(cmnd)
    os.system(cmnd)
    cmnd = "scp "+rain_dir+"arctic_comp/auto_process.sh "+final_dir
    print(cmnd)
    os.system(cmnd)
    # Get colocate_vars stuff
    final_dir = dest_dir + 'Arctic_compares/processing/colocate_vars/'
    cmnd = "scp "+rain_dir+"colocate_lib/*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)
    
    # comp grid climo stuff
    final_dir = dest_dir + 'Arctic_compares/processing/comp_grid_climo/'
    cmnd = "scp "+rain_dir+"comp_grid_climo/*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)
    cmnd = "scp "+rain_dir+"comp_grid_climo/Make* "+final_dir
    print(cmnd)
    os.system(cmnd)
    
    # type_analysis stuff
    final_dir = dest_dir + 'Arctic_compares/processing/type_analysis/'
    cmnd = "scp "+rain_dir+"type_analysis/*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)
    cmnd = "scp "+rain_dir+"type_analysis/Make* "+final_dir
    print(cmnd)
    os.system(cmnd)
    #cmnd = "scp "+rain_dir+"arctic_comp/auto_process.sh "+final_dir
    #print(cmnd)
    #os.system(cmnd)
    
    rain_dir = "bsorenson@134.129.222.68:/home/bsorenson/CERES/"
    #rain_dir = "bsorenson@raindrop.atmos.und.edu:/home/bsorenson/CERES/"
    
    final_dir = dest_dir + 'CERES/processing/'
    cmnd = "scp "+rain_dir+"*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)
    cmnd = "scp "+rain_dir+"Make* "+final_dir
    print(cmnd)
    os.system(cmnd)

#cmnd = "cp "+base_dir+"OMI/OMI_simulation/*.py "+dest_dir
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
# FuLiou analysis
# ---------------------------------------------------------------------------- 
print("FuLiou")
final_dir = dest_dir + 'FuLiou/'
cmnd = "find "+base_dir+"FuLiou/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

cmnd = "cp "+base_dir+"FuLiou/Ed4_LaRC_FuLiou/src/simple/*.f90 "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

cmnd = "cp "+base_dir+"FuLiou/Ed4_LaRC_FuLiou/src/simple/REAMDE_FuLiou "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

# ---------------------------------------------------------------------------- 
# GOES analysis
# ---------------------------------------------------------------------------- 
print("GOES")
final_dir = dest_dir + 'GOES/'
cmnd = "find "+base_dir+"GOES/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

# ---------------------------------------------------------------------------- 
# GISSTEMP analysis
# ---------------------------------------------------------------------------- 
print("GISSTEMP")
final_dir = dest_dir + 'GISSTEMP/'
cmnd = "find "+base_dir+"GISSTEMP/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

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
# NAAPS analysis
# ---------------------------------------------------------------------------- 
print("NAAPS")
final_dir = dest_dir + 'NAAPS/'
cmnd = "find "+base_dir+"NAAPS/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)
# ---------------------------------------------------------------------------- 
# NCEP analysis
# ---------------------------------------------------------------------------- 
print("NCEP")
final_dir = dest_dir + 'NCEP/'
cmnd = "find "+base_dir+"NCEP/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)
# ---------------------------------------------------------------------------- 
# NEXRAD analysis
# ---------------------------------------------------------------------------- 
print("NEXRAD")
final_dir = dest_dir + 'NEXRAD/'
cmnd = "find "+base_dir+"NEXRAD/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)
# ---------------------------------------------------------------------------- 
# NSIDC analysis
# ---------------------------------------------------------------------------- 
print("NSIDC")
final_dir = dest_dir + 'NSIDC/'
cmnd = "find "+base_dir+"NSIDC/ -type f -name \"*.py\" | xargs cp -t "+final_dir
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

cmnd = "find " + base_dir + "/data/AIRS/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

# ---------------------------------------------------------------------------- 
# VIIRS analysis
# ---------------------------------------------------------------------------- 
print("VIIRS")
final_dir = dest_dir + 'VIIRS/'
cmnd = "find "+base_dir+"VIIRS/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

# ---------------------------------------------------------------------------- 
# WRF analysis
# ---------------------------------------------------------------------------- 
print("WRF")
final_dir = dest_dir + 'WRF/'
cmnd = "find "+base_dir+"WRF/ -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

# ---------------------------------------------------------------------------- 
# SBDART analysis
# ---------------------------------------------------------------------------- 
print("SBDART")
work_dir = 'SBDART/'

final_dir = dest_dir + 'SBDART/'
cmnd = "find "+base_dir+ work_dir + \
    " -type f -name \"*.py\" | xargs cp -t "+final_dir
#cmnd = "cp "+base_dir+"OMI/*.py "+dest_dir
print(cmnd)
os.system(cmnd)

# ---------------------------------------------------------------------------- 
# TROPOMI analysis
# ---------------------------------------------------------------------------- 
print("TROPOMI")
final_dir = dest_dir + 'TROPOMI/'
cmnd = "find "+base_dir+"TROPOMI/ -type f -name \"*.py\" | xargs cp -t "+final_dir
print(cmnd)
os.system(cmnd)

# Get codes from Raindrop
# -----------------------
if(copy_raindrop):
    rain_dir = "bsorenson@134.129.222.68:/home/bsorenson/OMI/tropomi_colocate/"
    #rain_dir = "bsorenson@raindrop.atmos.und.edu:/home/bsorenson/OMI/tropomi_colocate/"
    
    final_dir = dest_dir + 'TROPOMI/processing/'
    cmnd = "scp "+rain_dir+"*.f90 "+final_dir
    print(cmnd)
    os.system(cmnd)
    cmnd = "scp "+rain_dir+"Make* "+final_dir
    print(cmnd)
    os.system(cmnd)
    
    cmnd = "scp "+rain_dir+"/auto_trop_process.sh "+final_dir
    print(cmnd)
    os.system(cmnd)

# ---------------------------------------------------------------------------- 
# Siphon/Metpy Codes
# ---------------------------------------------------------------------------- 
print("Siphon/Metpy")
smcodes = home_dir + '/programs/'
cmnd = "cp "+smcodes+"*.py "+dest_dir
print(cmnd)
os.system(cmnd)
### ---------------------------------------------------------------------------- 
### Jianglong/Jeff's Mie stuff
### ---------------------------------------------------------------------------- 
##print("Mie stuff")
##smcodes = 'home_dir + /Research/Jianglong_stuff/'
##cmnd = "cp "+smcodes+"*.pro "+dest_dir
##cmnd = "cp "+smcodes+"run_mie_code "+dest_dir
##print(cmnd)
##os.system(cmnd)

# Automate the uploading to Github
# --------------------------------

# Change to Arctic codes directory
#os.chdir(dest_dir)
##os.chdir('home_dir + /Arctic_codes/')
#
## Add new stuff
#os.system('git add .')
#
## Determine today's date for the command
#today_str = datetime.today().strftime('%Y/%m/%d')
#cmnd = 'git commit -m \"Backup '+today_str + '\"'
#os.system(cmnd)
#
## Push
#os.system('git push origin master')
