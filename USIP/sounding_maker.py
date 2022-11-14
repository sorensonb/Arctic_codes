#!/usr/bin/env python

"""
NAME:
  sounding_maker.py

PURPOSE:
  Creates basic soundings for any number of readsounding-compatible soundings

SYNTAX:
  ./sounding_maker.py any_number_of_soundings

MODIFICATIONS:
  Blake Sorenson  - 180119: Written 

"""

# Import functions
import os
import sys

if(len(sys.argv) < 2):
    print "SYNTAX: ./sounding_maker.py any_number_of_soundings"
    sys.exit(1)

args = sys.argv[1:]

path='/nas/home/bsorenson/prg/plotsounding/'
# Add check for compare_cn2.py, plot_cn2.py, and readsounding.py
if (os.path.isfile(path+'plotsounding.py') is False):
    print "ERROR: "+path+"plotsounding.py not found"
    print "       Please ensure the necessary program is in the correct location"
    sys.exit(1)

# Go through and make all the image files from the selected .txt files
for arg in args:
    os.system('$ADPAA_DIR/src/python_lib/plotsounding.py '+arg+' save')
