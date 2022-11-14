#!/usr/bin/env python

"""
NAME:
  cn2_profile_maker.py

PURPOSE:
  Creates cn2 profile images from any number of readsounding-compatible 
  soundings.

SYNTAX:
  ./cn2_profile_maker.py any_number_of_soundings

MODIFICATIONS:
  Blake Sorenson  - 180119: Written 

"""

# Import functions
import os
import sys

if(len(sys.argv) < 2):
    print "SYNTAX: ./cn2_profile_maker.py any_number_of_soundings"
    sys.exit(1)

args = sys.argv[1:]

# Add check for compare_cn2.py, plot_cn2.py, and readsounding.py
if (os.path.isfile('/nas/home/bsorenson/prg/thermosonde1617/plot_cn2.py') is False):
    print "File not found: /nas/home/bsorenson/prg/thermosonde1617/plot_cn2.py"
    print "Please ensure plot_cn2.py is in the correct directory before proceeding"
    sys.exit(1)

# Go through and make all the image files from the selected .txt files
for arg in args:
    os.system('/nas/home/bsorenson/prg/thermosonde1617/plot_cn2.py '+arg+' save')
