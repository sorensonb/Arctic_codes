#!/usr/bin/env python

"""
  NAME:
    convert_ARL.py

  PURPOSE:
    Read in an archived sounding file from the Air Research Laboratory 
    (www.arl.noaa.gov/...) and output the file in GSD format. The output
    is printed to the command line and piped to an output file.

    Input file must be of format
      <model>_YYYYMMDDHH_ANALYSIS_<location>_GSD
      YYMMDD_HHMMSS_<model>_<location>_ARL.txt
      160501_000000_BIS_ARL.txt:
    - Compatible Models:
      - HRRR
      - NAM


    WEBSITE: https://www.ready.noaa.gov/READYamet.php
      This contains the tool for accessing archived HRRR soundings. Only
      goes back about 4 years.

  SYNTAX:
    ./convert_ARl.py archived_sounding_file > output_file_name

  EXAMPLE:
    To convert an ARL file named May2016.txt to a GSD-formatted file named
    NAM_2016050100_ANALYSIS_BIS_ARL

    $ ./convert_ARl.py 160501_000000_NAM_BIS_ARL.txt > NAM_2016050100_ANALYSIS_BIS_ARL

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu> -   2018/04/30: Written    
"""

import numpy as np
import sys
import re
import math
from datetime import date

# Check to make sure there are enough arguments
if(len(sys.argv)!=2):
    print "SYNTAX: ./convert_ARl.py sounding_file"
    print "Input file must be of format"
    print "      YYMMDD_HH0000_<model>_<location>_ARL.txt"
    print "Compatible models:"
    print " - HRRR"
    print " - NAM"
    print "EXAMPLE: 171026_120000_NAM_GFK_ARL.txt"
    sys.exit()


pressIndex = 0 
altIndex = 1 
tempIndex = 2 
dpIndex = 3
spdIndex = 5 
dirIndex = 4

press = np.array([])
alt = np.array([])
temp = np.array([])
dp = np.array([])
spd = np.array([])
dct = np.array([])


# UNIVERSAL READER
count=0
cline=''
model=''
with open(sys.argv[1]) as f:
    inlines = np.array(f.readlines())
    xindex=0
    for i in range(0,len(inlines)):
        if(re.sub(' +',' ',inlines[i].strip()).split(' ')[0]=='YR:'):
            xindex=i
        if(re.sub(' +',' ',inlines[i].strip()).split(' ')[0]=='PRESS'):
            startIndex=i+4
            break
    cline = inlines[xindex]
    for line in inlines[startIndex:]:
        # Remove any tabs to account for GRAW files
        templine = line.strip()
        templine = re.sub(" +"," ", templine.strip()).split(" ")
        if(len(templine)==6):
            press = np.append(press[:],int(float(templine[pressIndex][:-1]))*10)
            alt = np.append(alt[:],int(templine[altIndex][:-2]))
            temp = np.append(temp[:],int(float(templine[tempIndex])*10))
            dp = np.append(dp[:],int(float(templine[dpIndex])*10))
            spd = np.append(spd[:],int( float(templine[spdIndex][:-1])*1.94384 ))
            dct = np.append(dct[:],int(float(templine[dirIndex])))

tline =re.sub(' +',' ',cline.strip()).split(' ') 
lat = tline[-3]
lon = tline[-1][:-1]
year = tline[1]
month = tline[3]
day = tline[5]
hour = tline[7]
lat = abs(float(tline[13]))
lon = abs(float(tline[15]))
model = sys.argv[1].split("_")[-3]
site = sys.argv[1].split("_")[-2]

final_month = date(int(year),int(month),1).strftime('%b')

# Make output file name
outname = model+"_"+year+month+day+hour+"_ANALYSIS_"+site+"_ARL"

# FIX: - Update model name
#      - Update date
#      - Figure out what to do with the station identifier
#      - Add correct lat and lon to line 4
with open(outname,'w') as opf:
    print >> opf, model+" analysis valid for grid point 10.2 nm / 28 deg from "+site+":"
    print >> opf, model[:3]+"{:>10} {:>6}      ".format(int(hour),int(day))+final_month+"    "+year
    print >> opf, "   CAPE      0    CIN     -0  Helic  99999     PW     12"
    print >> opf, "      1  23062      0  "+str(lat)+"  "+str(lon)+"    273  99999"
    print >> opf, "      2  99999  99999  99999 {:6}  99999  99999".format(len(press)-2+6)
    print >> opf, "      3           "+site+"                   12     kt"
    print >> opf, '      9 {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}'.format(int(press[0]),int(alt[0]),int(temp[0]),int(dp[0]),int(dct[0]),int(spd[0]))
    for i in range(1,len(press)):
        print >> opf, '      4 {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}'.format(int(press[i]),int(alt[i]),int(temp[i]),int(dp[i]),int(dct[i]),int(spd[i]))
print "Saved file: "+outname 
