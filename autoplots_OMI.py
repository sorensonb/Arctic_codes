#!/usr/bin/env python

"""
  NAME:
    autoplots_OMI.py

  PURPOSE:
    Automate the generation of single-swath OMI figures for new days
    of OMI data. Used for comparing to single-swath CERES data
    for Arctic aerosol radiative impacts study.

  SYNTAX:
    $ python autoplots_OMI.py date
    (date is of format YYYYMMDD)

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>  : 2020/03/30
      Written (modification of autoplots_CERES_Ed4.py)

"""

import os
import sys
from subprocess import check_output

if(len(sys.argv) != 2):
    print("SYNTAX: ./autoplots_OMI.py date")
    print("         date: YYYYMMDD")
    sys.exit()

in_date = sys.argv[1]

base_path = "/home/bsorenson/Research/OMI/"
os.chdir(base_path)

# First make the single_swath directory
# -------------------------------------
os.system('mkdir single_swath_' + in_date)

# Then, change into the single_swath directory
# -------------------------------------
os.chdir('single_swath_' + in_date + '/')

# Make daily average figure
# -------------------------------------
os.system('../../plot_single_OMI.py ' + in_date + ' LWF')

# Pull the data file names and make single-swath figures
# ------------------------------------------------------
files = check_output('ls /home/bsorenson/data/OMI/H5_files/OMI*' + \
    in_date[:4] + 'm' + in_date[4:8] + '*.he5',\
    shell=True).decode('utf-8').strip().split('\n')
for infile in files:
    os.system('../plot_single_OMI.py ' + in_date + \
    infile.split('_')[-2][10:14] + ' UVAerosolIndex 60')

