#!/usr/bin/env python

"""
  NAME:
    autoplots.py
  
  PURPOSE:
  
  SYNTAX:
      ./autoplots.py year location time type [clean]"
          year      - 2015, 2016, 2017, or all"
          location  - BIS, JAX, TUS, or all"
          time      - 00, 12, or all"
          type      - cn2, temperature, or all "
          clean     - removes extracted sounding files from sounding"
                      archive when done"
  
  EXAMPLE:
  
  MODIFICATIONS:
    Blake Sorenson  - 180129: Written 
  
"""

# Import functions
import glob
import os
import sys
import commands
import numpy
from time import strptime

def syntax():
    print "SYNTAX: ./autoplots.py year location time type [clean]"
    print "        year      - 2015, 2016, 2017, or all"
    print "        location  - BIS, JAX, TUS, or all"
    print "        time      - 00, 12, or all"
    print "        type      - cn2, temperature, or all "
    print "        clean     - removes extracted sounding files from sounding"
    print "                    archive when done"


if(len(sys.argv) < 5):
    syntax()
    sys.exit(1)

CLEAN=False
for arg in sys.argv:
    if((arg=='help') | (arg=='-h')):
        syntax()
        sys.exit(1)
    if(arg=='clean'):
        CLEAN=True

base_path = '/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/'

status,start_dir = commands.getstatusoutput('pwd') 
years = [sys.argv[1]]
locations=[sys.argv[2]]
times = [sys.argv[3]]
plottypes = [sys.argv[4]]
if(years[0]=='all'):
    status, output = commands.getstatusoutput('ls '+base_path+'sounding_archive/')
    years = output.split('\n')
if(plottypes[0]=='all'):
    plottypes=['cn2','temperature']

# Open all the text files

# Process years
for year in years:
    # If the user wants all locations for the current year, automatically find
    # the location directory names within the current year directory
    if(locations[0]=='all'):
        status, output = commands.getstatusoutput('ls '+base_path+'sounding_archive/'+year+'/')
        locations = output.split('\n')
    # Process locations
    for loc in locations:
        # Get a list of all the archive files for the current location
        status,archives = commands.getstatusoutput('ls '+base_path+'sounding_archive/'+year+'/'+loc+'/*.txt') 
        archives = archives.split('\n')
        # Sort the archive files in order from January to December
        try:
            archives = sorted(archives, key=lambda f : strptime(f.split('/')[-1].split('_')[0][:3],"%b").tm_mon)
        except ValueError as exc:
            print "ERROR: Files not located in "+base_path+'sounding_archive/'+year+'/'+loc
            print "Please ensure the correct files are in the directory before proceeding"
            continue
        # If the user wants all times for the current location, automatically 
        # find the time directory names within the current location directory
        if(times[0]=='all'):
            status, output = commands.getstatusoutput('ls '+base_path+'sounding_archive/'+year+'/'+loc+'/images/boxplots/')
            times = output.split('\n')
            times = numpy.array(times)
            # Remove any unwanted directory names
            times = numpy.delete(times,numpy.where(times=='both'))
                    
        # Extract individual soundings from archive files? 
        for time in times:
            # Ignore the extra 'Z' on the end of the time
            time = str(time[:2])
            # If the user wants all image types for the current time, 
            # automatically find the time directory names within the 
            # current location directory
            if(plottypes[0]=='all'):
                status, output = commands.getstatusoutput('ls '+base_path+'sounding_archive/'+year+'/'+loc+'/images/boxplots/'+time+'Z/')
                plottypes = output.split('\n')
            for ptype in plottypes:
                image_path = base_path+'sounding_archive/'+year+'/'+loc+'/images/boxplots/'+time+'Z/'+ptype+'/'
                # Loop over each of the archive files for the current location
                for afile in archives:
                    # Get the number of the current month to determine which 
                    # files to remove later
                    month = afile.split('/')[-1].split('_')[0]
                    month_no = strptime(month[:3],'%b').tm_mon 
                    str_monthno = str(month_no)
                    if(month_no<10):
                        str_monthno = '0'+str_monthno
                    # Use one of the monthly analyzer programs to generate the
                    # appropriate boxplot 
                    os.system(base_path+'/monthly_'+ptype+'_analyzer.py '+time+' '+afile)
                print "Moving images to "+ image_path
                os.system('mv *_'+year+'_'+loc+'_'+time+'Z_*.png '+image_path)
        if CLEAN is True:
            # Remove the extra archive text files for the current location 
            print "Removing extracted sounding files from" + base_path+'sounding_archive/'+year+'/'+loc+'/soundings/'
            os.system('rm '+base_path+'sounding_archive/'+year+'/'+loc+'/soundings/'+year[2:]+'*'+loc+'_UWYO.txt')
