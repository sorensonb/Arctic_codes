#!/usr/bin/env python

"""
  NAME:
    add_folder.py
  
  PURPOSE:
    Automatically create the directory tree for a new year or location in
    the sounding archive.
  
  SYNTAX:
    ./add_folder.py year location [time] [type]
  
  EXAMPLE:
    To add the year 2012 with all locations to the tree:
      $$ ./add_folder.py 2012 all

    To add only Bismarck to the 2013 directory tree:
      $$ ./add_folder.py 2013 BIS
  
  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>  - 2018/03/22:
      Written
  
"""

# Import functions
import glob
import os
import sys
import commands
import numpy
from time import strptime

def syntax():
    print "SYNTAX: ./add_folder.py year location [time] [type] [clean]"

if(len(sys.argv) < 3):
    syntax()
    sys.exit(1)

CLEAN=False
for arg in sys.argv:
    if((arg=='help') | (arg=='-h')):
        syntax()
        sys.exit(1)
    if(arg=='clean'):
        CLEAN=True

base_path = '/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/sounding_archive/'

all_locations = ['BIS','JAX','TUS']
all_times = ['00Z','12Z']
plottypes = ['cn2','temperature']

years = [sys.argv[1]]
locations=[sys.argv[2]]
if(locations[0]=='all'):
    locations = all_locations
times = all_times
### plottypes = [sys.argv[4]]
### if(plottypes[0]=='all'):
###     plottypes=all_types

# Open all the text files

# Process years
for year in years:
    try:
        os.makedirs(base_path+year)
        print "Created "+base_path+year
    except OSError as exc:  
        print base_path+year+'  already exists'
    # If the user wants all locations for the current year, automatically find
    # the location directory names within the current year directory
    if(locations[0]=='all'):
        locations = all_locations
    # Process locations
    for loc in locations:
        try:
            os.makedirs(base_path+year+'/'+loc)
            print "Created "+base_path+year+'/'+loc
        except OSError as exc:
            print base_path+year+'/'+loc + '  already exists'
        try:
            os.makedirs(base_path+year+'/'+loc+'/images')
            print "Created "+base_path+year+'/'+loc+'/images'
        except OSError as exc:
            print base_path+year+'/'+loc+'/images' + '  already exists'

        try:
            os.makedirs(base_path+year+'/'+loc+'/soundings')
            print "Created "+base_path+year+'/'+loc+'/soundings'
        except OSError as exc:
            print base_path+year+'/'+loc+'/soundings' + '  already exists'

        try:
            os.makedirs(base_path+year+'/'+loc+'/images/boxplots')
            print "Created "+base_path+year+'/'+loc+'/images/boxplots'
        except OSError as exc:
            print base_path+year+'/'+loc+'/images/boxplots' + '  already exists'
            
        if(times[0]=='all'):
            times=all_times                   
 
        for time in times:
            try:
                os.makedirs(base_path+year+'/'+loc+'/images/boxplots/'+time)
                print "Created "+base_path+year+'/'+loc+'/images/boxplots/'+time
            except OSError as exc:
                print base_path+year+'/'+loc+'/images/boxplots/'+time + '  already exists'

            for ptype in plottypes:
                print ptype
                try:
                    os.makedirs(base_path+year+'/'+loc+'/images/boxplots/'+time+'/'+ptype)
                    print "Created "+base_path+year+'/'+loc+'/images/boxplots/'+time+'/'+ptype
                except OSError as exc:
                    print base_path+year+'/'+loc+'/images/boxplots/'+time+'/'+ptype + '  already exists'
