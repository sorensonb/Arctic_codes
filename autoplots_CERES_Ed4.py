#!/usr/bin/env python

"""
  NAME:
    autoplots_CERES_Ed4.py

  PURPOSE:
    Automate the generation of single-swath CERES SWF and LWF 
    figures for new days of CERES data. Used for comparing to 
    single-swath OMI data for Arctic aerosol radiative impacts study.

  SYNTAX:
    $ python autoplots_CERES_Ed4.py date
    (date is of format YYYYMMDD)

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>  : 2020/03/30
      Written

"""

import os
import sys

if(len(sys.argv) != 2):
    print("SYNTAX: ./autoplots_CERES_Ed4.py date")
    print("         date: YYYYMMDD")
    sys.exit()

in_date = sys.argv[1]

base_path = "/home/bsorenson/Research/CERES/"
os.chdir(base_path)

# First make the single_swath directory
# -------------------------------------
os.system('mkdir single_swath_' + in_date)

# Then, make the subdirectories
# -------------------------------------
os.system('mkdir single_swath_' + in_date + '/LWF')
os.system('mkdir single_swath_' + in_date + '/SWF')

# Then, change into the LWF directory
# -------------------------------------
os.chdir('single_swath_' + in_date + '/LWF')

# Make daily average and single-swath figures
# -------------------------------------
os.system('../../plot_single_CERES_Ed4.py ' + in_date + ' LWF')
for ii in range(24):
    os.system('../../plot_single_CERES_Ed4.py ' + in_date + str(ii).zfill(2) + ' LWF')

# Then, go to SWf and do the same
# -------------------------------------
os.chdir('../SWF')
os.system('../../plot_single_CERES_Ed4.py ' + in_date + ' SWF')
for ii in range(24):
    os.system('../../plot_single_CERES_Ed4.py ' + in_date + str(ii).zfill(2) + ' SWF')
