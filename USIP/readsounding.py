#!/usr/bin/env python

"""
  NAME:
    readsounding.py

  PURPOSE:
    Reads in multiple sounding types and returns a dictionary containing
    all of the data from the sounding file. 

    Eight keys in the dictionary,
        - PRESS
        - ALT
        - TEMP
        - DP
        - SPD
        - DP
        - UWIND
        - VWIND
    which correspond to the eight "standard" sounding data, are in all caps.
    These eight entries will always be present in the output dictionary 
    no matter what the sounding format is. The other data in the dictonary will
    use the same variable names as those given in the sounding file.

  CALLS:
    Built-in Python modules np, re, sys, math,  netCDF4

  PYTHON VERSION:
    2.7.5

  IMPORTING:
    To import all functions from readsounding, use this syntax:
      >>> from readsounding import *

    or
      >>> import readsounding as rs

    Then, to use the functions, use this syntax:
      >>> a = readsounding('sounding file')
      >>> spd = kts2ms(a['SPD'])
      >>> kts_spd = ms2kts(spd)

    or
      >>> a = rs.readsounding('sounding file')
      >>> spd = rs.kts2ms(a['SPD'])
      >>> kts_spd = rs.ms2kts(spd)

    To import single functions, use this syntax:
      >>> from readsounding import readsounding
      >>> from readsounding import ms2kts
      >>> from readsounding import kts2ms

  USAGE:
    Depending on how the functions were imported (see IMPORTING), sounding
    data can be read either by
      >>> S = readsounding('input_filename', [clat=optional latitude], [clon=optional longitude], [allGraw=True/False], [keepMissing=True/False])
      
    or by 
      >>> S = rs.readsounding('input_filename')

    In the example above:
        - 'input_filename' is a sounding text file in a compatible format 
        - clat and clon are optional coordinates used for extracting a profile 
          from a WRF output file centered on those coordinates,
        - allGraw is a boolean (set to False by default) that determines whether
          or not to read all
        - keepMissing is a boolean (set to False by default) that allows the
          use to keep any missing values in the data. 


    This returns a dictionary containing all the data from the sounding file.
    The standard keys in the dictionary that are present no matter what the
    sounding type is are:
      S['NAME']  : file title containing the year, month, day, hour, and type of
                   sounding
      S['TITLE'] : Plot title containing the location, year, month, day, and UTC
                   hour of the sounding
      S['LABEL'] : Plot label containing type, date, location, and hour of
                   sounding
      S['UNITS'] : Dictionary containing the units of each variable in the file
      S['PRESS'] : pressure (in hPa)
      S['ALT']   : altitude (in m)
      S['TEMP']  : temperature (in degC)
      S['DP']    : dew point temperature (in degC)
      S['SPD']   : wind speed (in kts)
      S['DIR']   : wind direction (in degrees)
      S['UWIND'] : zonal wind speed (in kts)
      S['VWIND'] : meridional wind speed (in kts)

    To see what the other keys in the dictionary are, either look at the 
    header in the sounding text file or type
    >>> print S.keys()

  EXAMPLE:
    Read in sounding data from MCR_8_Aug_1000a_1.txt:
      >>> S = readsounding('MCR_8_Aug_1000a_1.txt:')  

    Access data from the returned dictionary:
      >>> S['TEMP']:
      array([ 28. ,  27.5,  27.2, ..., -37.8, -37.5, -37.4])
      >>> S['PRESS']:
      array([ 1012.5 ,  1010.5 ,  1007.09, ...,    11.56,    11.5 ,    11.45])   

    Convert zonal wind speeds from knots to meters/second:
      >>> S['UWIND'] = kts2ms(S['UWIND'])

    Convert zonal wind speeds from meter/second to knots:
      >>> S['UWIND'] = ms2kts(S['UWIND'])

    Read all data from a Graw sounding file (both the ascent and descent) 
    titled /nas/Radiosonde/20170503/170503_204757_GFK_GRAW.txt:
      >>> G = readsounding('/nas/Radiosonde/20170503/170503_204757_GFK_GRAW.txt',allGraw=True)

    The returned dictionary G will be formatted the same as any other
    readsounding dictionary, but it will contain any and all data after balloon
    burst in addition to all ascent data. 

  NOTES:
    readsounding works with these sounding types:
    - GRAW RTS and ProfileDataTables   
                       (RTS: /nas/Radiosonde/20170503/170503_000000_GFK_GRAW.txt
                        PDT: /nas/Radiosonde/Flights/140923_000000_GFK_GRAW.txt)
    - NWS soundings (UWyo upper-air database)           
                 (ex: /nas/Radiosonde/UWYO_Soundings/140721_180000_BIS_UWYO.txt)
        - Format:
          YYMMDD_HHMMSS_<location>_UWYO.txt
    - Model soundings from https://rucsoundings.noaa.gov
      - Files must have this name format:
        <model>_YYYYMMDDHH_<type>_<location>_GSD    
      - Compatible Models:
        - RAP.Op40
        - RAP.Bak40
        - NAM    
        - GFS
        - FIM
          ex: RAP.Op40_2017052412_ANALYSIS_BIS_GSD
              RAP analysis sounding for Bismarck for May 24, 2017 at 12Z
          ex: GFS_2017052612_F012_BIS_GSD
              GFS 12-hour forecast sounding for Bismarck at 12Z on May 26, 2017
    - Model soundings from the Air Resources Laboratory Archive 
        (https://ready.arl.noaa.gov/READYamet.php)
      - Files must have this name format (similar to rucsoundings soundings):
        <model>_YYYYMMDDHH_<type>_<location>_GSD    
      - Compatible Models:
        - HRRR
        - NAM
        ex: HRRR_2018050121_ANALYSIS_GFK_ARL
            HRRR analysis sounding for Grand Forks for May 1, 2018 at 21Z
    - Model soundings from weather.cod.edu
      - Files must have this name format (similar to rucsoundings soundings):
          <model>_YYYYMMDDHH_<type>_<location>_COD
        where <location> is the nearest METAR site. To download a sounding 
        from COD, open up a model and click a point on the map. Click 
        "Raw Sounding Text" just below the plot and copy and paste everything
        that comes up into a file.
      - Compatible models:
        - RUC
        - NAM
        - NAMNST
        - GFS
    - Plymouth State archived soundings
        - Files must have similar format to NWS/UWYO soundings:
          YYMMDD_HHMMSS_<location>_PSU.txt
        - Soundings found at http://vortex.plymouth.edu/myo/upa/raobplt-a.html 
    - Model soundings from weather.admin.niu.edu/machine/textfcstsound.html
      - Files must have this name format:
          <model>_YYYYMMDDHH_<type>_<location>_NIU
      - Compatible models:
        - RUC
        - WRF
        - GFS
    - WRF output
        - When using WRF output, add latitude and longitude to the function
          call to specify the location of the profile. If no latitude and
          longitude are passed to the function, it uses the middle of the
          model area to find the profile.
          SYNTAX: a = readsounding(<wrf file>, <lat>, <lon>)
    - COAMPS forecast soundings 
                                (ex: /nas/und/Florida/2015/Analysis/ModelFiles/)
    - Snow White Radiosonde soundings from the CAPE2015 project
                              (ex: /nas/und/Florida/2015/Analysis/BalloonFiles/)
          
    https://rucsoundings.noaa.gov/ (forecast soundings)
    https://ruc.noaa.gov/raobs/    (FSL soundings) 

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>      - 2017/06/27: 
        Written
    Blake Sorenson <blake.sorenson@und.edu>      - 2017/06/28: 
        Added UTC time for Graw soundings
    Blake Sorenson <blake.sorenson@und.edu>      - 2017/07/18:
        Added plot label to dictionary
    Blake Sorenson <blake.sorenson@und.edu>      - 2017/08/02:
        Added WRF output option
    Blake Sorenson <blake.sorenson@und.edu>      - 2017/08/08:
        Added PSU archived soundings             
    Blake Sorenson <blake.sorenson@und.edu>      - 2017/10/12:
        Added readWRF to extract WRF profile following GRAW radiosonde data
    Blake Sorenson <blake.sorenson@und.edu>      - 2018/01/17:
        Added readAllGraw function
    Blake Sorenson <blake.sorenson@und.edu>      - 2018/02/12:
        Updated GRAW data section
    Blake Sorenson <blake.sorenson@und.edu>      - 2018/02/26:
        Rewrote to make the output dictionary contain all data from the
        sounding file 
    Blake Sorenson <blake.sorenson@und.edu>      - 2018/03/26:
        Fixed UWYO data reading problem
    Blake Sorenson <blake.sorenson@und.edu>      - 2018/05/02:
        Updated known missing values to work with new Graw RTS tables. 
        Added archived model soundings from Air Resources Laboratory archive.
"""

import numpy as np
import re
import sys
import math
import datetime
from netCDF4 import Dataset

# kts2ms converts wind speeds (U, V, or whole) from knots to meters/second
def kts2ms(speeds):
    speeds *= 0.514444
    return speeds 

# ms2kts converts wind speeds (U, V, or whole) from meters/second to knots
def ms2kts(speeds):
    speeds *= 1.94384
    return speeds 

###############################################################################
#
# readsounding is the main function that reads the sounding data from the file
# and creates the output dictionary.
#
# If a WRF profile is passed, a center lat and lon must be provided
#
###############################################################################
def readsounding(inputfile, clat=None, clon=None, allGraw=False, keepMissing=False): 
    pi    = 3.14159265
    graw  = False  # graw is used as a boolean to determine whether or not to
                   # add 'TIME' and 'UTC' keys to the output dictionary. Graw 
                   # sounding files are the only sounding files that contain
                   # time
    inlines = np.array(open(inputfile,'r').readlines())
    outdict = dict()
    keys = []
    units = [] 
    COAMPS=False
    MCR=False
    GRAW=False
    GSD=False
    winds_to_ms=False
    winds_to_kts=True
    #--------------------------------------------------------------------------
    #
    # Use the sounding type to determine which line in the sounding to pull
    # variable names and units from. Find the lines to begin and end reading
    # data.
    #
    # startline  - the index of the first line from which data is read
    # endline    - the index of the last line from which data is red
    # keys       - an array that contains the variable names extracted from
    #              the file
    # units      - an array that contains the units that correspond to the 
    #              variables that were extracted from the file. If the file
    #              does not have any units, the units array is hard coded
    #
    #--------------------------------------------------------------------------
    #COAMPS
    if(inputfile.split("/")[-1].split("_")[0]=='coamps.prof'):
        COAMPS=True
        # Split the data file name to get date, time, and location
        tempname = inputfile.split('/')[-1].split('_')
        model    = "COAMPS"
        datetemp = tempname[4]
        _type    = tempname[5]
        location = tempname[6]+tempname[7]
        year     = datetemp[2:4]
        month    = datetemp[4:6]
        day      = datetemp[6:8]
        thour     = datetemp[8:10]
        minute   = "00"
        second   = "00"
        # Calculate the sounding hour by adding the starting hour to the
        # forecast hour
        hour  = str(int(tempname[5][2:])+int(thour) )    
        saver    = model+" "+_type+" "+location
        outputstart = year+"_"+month+"_"+day+"_"+hour+"_"+minute+"_"+second+\
                      "."+location+"."+_type+"."+model

        startline = 1
        endline = len(inlines)
        keys = re.sub(' +',' ',inlines[0].strip()).split(' ')
        units = ['None','None','None','None','None','None','None','hPa','m','degC','degC','kts','deg','m/s','m/s']
        pressIndex = 7 
        altIndex = 8 
        tempIndex = 9 
        dpIndex = 10
        spdIndex = 11 
        dirIndex = 12

    # NOAA model sounding (GFS, NAM, FIM, RAP.Op40, RAP.Bak40) (rucsoundings.noaa.gov) (GSD)
    # Air Reserach Laboratory Archive model data (ready.arl.noaa.gov/) (ARL)
    elif((inputfile.split('/')[-1].split('_')[-1]=='GSD') | \
         (inputfile.split('/')[-1].split('_')[-1]=='ARL')):
        GSD=True
        winds_to_kts=False
        # Split the data file name to get date, time, location, and model
        tempname = inputfile.split('/')[-1].split('_')
        model    = tempname[0]
        datetemp = tempname[1]
        _type    = tempname[2]
        location = tempname[3]
        year     = datetemp[2:4]
        month    = datetemp[4:6]
        day      = datetemp[6:8]
        hour     = datetemp[8:10]
        minute   = "00"
        second   = "00"
        saver    = model+" "+_type+" "+location
        outputstart = year+"_"+month+"_"+day+"_"+hour+"_"+minute+"_"+second+\
                      "_"+location+"_"+_type+"_"+model

        startline = 6
        endline = len(inlines)
        keys = ['level','PRESS','ALT','TEMP','DP','DIR','SPD']
        units = ['None','hPa','m','degC','degC','deg','kts']
        pressIndex = 1 
        altIndex = 2 
        tempIndex = 3 
        dpIndex = 4
        spdIndex = 6 
        dirIndex = 5
    # GRAW
    elif(inputfile.split("/")[-1].split("_")[-1]=="GRAW.txt"):
        GRAW=True
        saver="GRAW RAOB"
        # Split the data file name to get date, time, location, and model
        tempname=inputfile.split('/')[-1].split('_')
        year   = tempname[0][:2]
        month  = tempname[0][2:4]
        day    = tempname[0][4:6]
        hour   = tempname[1][:2]
        minute = tempname[1][2:4]
        second = tempname[1][4:6]
        outputstart = year+"_"+month+"_"+day+"_"+hour+"_"+minute+"_"+second+\
                      "_GRAW"
        if((outputstart=='14_09_23_00_00_00_GRAW') | (outputstart=='17_02_14_00_00_00_GRAW')):
            winds_to_kts=False
        ghead = re.sub(' +','',inlines[0].strip()).split('\t')
        # Startline and endline tell the reader below where to read lines
        startline = 1
        endline = len(inlines)
        for x in ghead:
            tempvar = x.split('[')
            keys.append(tempvar[0])
            if(len(tempvar)==2):
                units.append(str(tempvar[1].split(']')[0]))
            else:
                units.append('None')
        keys = np.array(keys) 
        units = np.array(units) 
        pressIndex = np.where((keys=='Pressure') | (keys=='P'))
        # Add in a case to use the less precise 'Geopotential height' instead of higher
        # resolution 'geometric height'
        # Shut off for tethered test
        ###if(outputstart=='17_11_17_23_02_35_GRAW'):
        ###    altIndex=np.where(keys=='GeoPot')
        ###else:
        altIndex = np.where((keys=='GeometricHeight') | (keys=='Altitude'))
        tempIndex = np.where(keys=='T') 
        dpIndex = np.where(keys=='Dew')
        spdIndex = np.where(keys=='Wsp')
        dirIndex = np.where(keys=='Wdir')
    # MCR
    elif(inputfile.split('/')[-1].split('_')[0]=='MCR'):
        MCR=True
        winds_to_kts = False
        # Split the filename to get information
        tempname = inputfile.split("/")[-1].split("_")
        saver="MCR RAOB"
        year="15"
        tmonth=tempname[2]
        if((tmonth=="Aug") | (tmonth=="AUG")):
            month="08"
        else:
            month="07"
        tday=tempname[1]
        if(int(tday)<10):
            day="0"+tday
        else:
            day=tday
        ttime=tempname[3]
        # The hour in the filename is in local time, so convert to UTC
        hour=str(int(ttime[:2])+4)
        minute=ttime[2:4]
        second="00"
        outputstart = year+"_"+month+"_"+day+"_"+hour+"_"+minute+"_"+second+\
                      "_MCR"

        startline = 9
        endline = np.where(inlines=='\r\n')[0][4]-2 
        keys  = np.array(re.sub(' +',' ',inlines[6].strip()).split(' '))
        units = ['m','deg','kts','/sec','degC','degC','hPa','%','g/m3','g/m3','N','kts','hPa','mm']
        pressIndex = 6 
        altIndex = 0 
        tempIndex = 4 
        dpIndex = 5
        spdIndex = 2 
        dirIndex = 1
        
    # weather.cod.edu or weather.admin.niu.edu sounding (see header)
    # (WRF, NAM, NAMNST, RAP, RUC, GFS) 
    elif((inputfile.split('/')[-1].split('_')[-1]=='COD') | \
         (inputfile.split('/')[-1].split('_')[-1]=='NIU')):
        winds_to_kts = False
        # Split the filename to get information
        tempname = inputfile.split('/')[-1].split('_')
        model    = tempname[0]
        datetemp = tempname[1]
        _type    = tempname[2]
        location = tempname[3]
        year     = datetemp[2:4]
        month    = datetemp[4:6]
        day      = datetemp[6:8]
        hour     = datetemp[8:10]
        minute   = "00"
        second   = "00"
        saver    = model+" "+_type+" "+location
        outputstart = year+"_"+month+"_"+day+"_"+hour+"_"+minute+"_"+second+\
                      "."+location+"."+_type+"."+model

        headline = np.where(inlines=='-------------------------------------------------------------------------------\n')[0][0]+1
        unitline = np.where(inlines=='-------------------------------------------------------------------------------\n')[0][0]+2
        startline = np.where(inlines=='-------------------------------------------------------------------------------\n')[0][1]+1
        endline = np.where(inlines=='\n')[0][1]-3
        keys  = np.array(re.sub(' +',' ',inlines[headline].strip()).split(' '))
        units = np.array(re.sub(' +',' ',inlines[unitline].strip()).split(' '))
        pressIndex = 1 
        altIndex = 2 
        tempIndex = 3 
        dpIndex = 4
        spdIndex = 9 
        dirIndex = 8

    # Plymouth State archived sounding
    elif(inputfile.split('/')[-1].split('_')[-1]=='PSU.txt'):
        winds_to_kts = False
        # Split the filename to get information
        tempname = inputfile.split("/")[len(inputfile.split("/"))-1].split("_")
        datetemp=tempname[0]
        timetemp=tempname[1]
        location=tempname[2]
        saver=location+" RAOB"
        # Extract YY,MM,DD,HH,MM,SS from the filename
        year  = datetemp[:2]
        month = datetemp[:4][2:]
        day   = datetemp[:6][4:]
        hour  = timetemp[:2]
        minute= timetemp[2:4]
        second= timetemp[4:6]
        outputstart = year+"_"+month+"_"+day+"_"+hour+"_"+minute+"_"+second+\
                      "_"+location+"_PSU"

        headline = np.where(inlines=='-------------------------------------------------------------------------------\n')[0][0]+1
        unitline = np.where(inlines=='-------------------------------------------------------------------------------\n')[0][0]+2
        startline = np.where(inlines=='-------------------------------------------------------------------------------\n')[0][1]+1
        endline = len(inlines)-2 
        keys  = np.array(re.sub(' +',' ',inlines[headline].strip()).split(' '))
        units = np.array(re.sub(' +',' ',inlines[unitline].strip()).split(' '))
        pressIndex = 1 
        altIndex = 2 
        tempIndex = 3 
        dpIndex = 4
        spdIndex = 9 
        dirIndex = 8
        
    # UWYO
    elif(inputfile[len(inputfile)-8:]=="UWYO.txt"):
        winds_to_kts = False
        # NWS/UWYO sounding
        # Split the filename to get information
        tempname = inputfile.split("/")[len(inputfile.split("/"))-1].split("_")
        datetemp=tempname[0]
        timetemp=tempname[1]
        location=tempname[2]
        saver=location+" RAOB"
        # Extract YY,MM,DD,HH,MM,SS from the filename
        year  = datetemp[:2]
        month = datetemp[:4][2:]
        day   = datetemp[:6][4:]
        hour  = timetemp[:2]
        minute= timetemp[2:4]
        second= timetemp[4:6]

        outputstart = year+"_"+month+"_"+day+"_"+hour+"_"+minute+"_"+second+\
                      "_"+location+"_UWYO"

        headnum = 3 # Index of the header line
        unitnum = 4 # Index of the units line
        # Create keys in the dictonary for each variable in the file
        bhead = inlines[3]
        keys = re.sub(' +',' ',bhead.strip()).split(' ')
        units = re.sub(' +',' ',inlines[4].strip()).split(' ')
        # Assume that the second new-line character in the UWYO file separates
        # the data from the metadata.
        # Also assume that the second to last line of data for all UWYO files
        # is missing wind data, so ignore that line 
        startline = 6
        if(len(re.sub(' +',' ',inlines[6].strip()).split(' ')) != 11):
            startline += 1 
        endline = np.where(inlines=='\n')[0][1]-1
        # Occasionally, a UWYO sounding has two levels at the end without wind
        # speed, so keep decreasing endline until the last line with full data
        # is found.
        while(len(re.sub(' +',' ',inlines[endline].strip()).split(' '))!=11):
            endline-=1
        endline+=1
        pressIndex = 0 
        altIndex = 1 
        tempIndex = 2 
        dpIndex = 3
        spdIndex = 7 
        dirIndex = 6
    else:
        print "ERROR: " + inputfile + " is not compatible with readsounding.py"
        print "\tSee readsounding.py header text for compatible sounding types"
        print "\tExiting"
        sys.exit(1)
 
    # Set the six standard key names
    keys[pressIndex]='PRESS'
    keys[altIndex]='ALT'
    keys[tempIndex]='TEMP'
    keys[dpIndex]='DP'
    keys[spdIndex]='SPD'
    keys[dirIndex]='DIR'

    # Convert the keys and units to numpy arrays 
    keys = np.array(keys)
    units = np.array(units)
 
    # After finding the variable names, create the dictionary
    # Initialize each data array in the dictionary to be full of 
    # MVCs
    for key in keys:
        outdict[key] = np.full(endline-startline,999999.9999)
     
    # UNIVERSAL READER
    count=0
    for line in inlines[startline:endline]:
        # Remove any tabs to account for GRAW files
        templine = line.strip().replace('\t',' ')
        templine = re.sub(" +"," ", templine.strip()).split(" ")
        for i in range(0,len(templine)):
            if((templine[i]!='////./') & (templine[i]!='///') & \
               (templine[i]!='//.//') & (templine[i]!='-----') & \
               (templine[i]!='///./') & (templine[i]!='/////') & \
               (templine[i]!='/////./') & (templine[i]!='99999') & \
               (templine[i]!='SFC') & (templine[i]!='TRP')):
                outdict[keys[i]][count] = templine[i]
        count+=1

    # Create the units dictionary
    unitdict = dict()
    for x, y in zip(keys, units):
        unitdict[x]=y

    # Reformat the units
    for key in unitdict:
        if((unitdict[key]=='\xb0') | (unitdict[key]=='DEG')):
            unitdict[key]='deg'
        elif(unitdict[key]=='\xb0C'):
            unitdict[key]='degC'
        elif(unitdict[key]=='KTS'):
            unitdict[key] = unitdict[key].lower()

    # Find missing values in the six standard data keys and remove the data in
    # each key that corresponds to the location of those values
    # NOTE: Only keep missing data if working with a Graw RTS table. This is
    # mainly for working on the NASA USIP thermosonde project in which good
    # data was being deleted because an obscure column was missing.
    #if(GRAW is not True):
    if(keepMissing is False):
        bad_indices = []
        for i in range(0,len(outdict['TEMP'])):
            if( (outdict['PRESS'][i]==999999.9999) | (outdict['ALT'][i]==999999.9999) | \
                (outdict['TEMP'][i]==999999.9999) | (outdict['DP'][i]==999999.9999) | \
                (outdict['SPD'][i]==999999.9999) | (outdict['DIR'][i]==999999.9999) ):
                bad_indices.append(i)

        for key in outdict:
            outdict[key] = np.delete(outdict[key],bad_indices) 


    # For each sounding type, perform any necessary final edits
    if(COAMPS is True):
        # Reverse data
        for key in outdict:
            outdict[key] = outdict[key][::-1]
        outdict['SPD']*=1.94384
        outdict['UTRU']*=1.94384
        outdict['VTRU']*=1.94384
    elif(MCR is True):
        # Convert altitude from feet to meters
        outdict['ALT']*=0.3048
    elif(allGraw is False):
        # Remove descent data
        # Find the indices where the data is not ascending (includes constant
        # height data)
        indices = np.where(np.diff(outdict['ALT'])<=0)[0]+1
        for key in outdict:
            outdict[key] = np.delete(outdict[key],indices)
    if(GSD is True):
        # Divide press, temp, and dew point by 10
        outdict['PRESS'] = outdict['PRESS']/10.
        outdict['TEMP'] = outdict['TEMP']/10.
        outdict['DP'] = outdict['DP']/10.
    if((inputfile.split('/')[-1].split('_')[-1]=='COD') | \
         (inputfile.split('/')[-1].split('_')[-1]=='NIU')):
        # Remove the 'level' key and data
        outdict.pop('LEV')
        unitdict.pop('LEV')
    if(winds_to_kts is True):
        outdict['SPD'] = ms2kts(outdict['SPD'])
        unitdict['SPD'] = 'kts'

    # Calculate u and v wind components
    UWind = np.array([-tspd*math.sin((pi/180)*tdir) for tspd, tdir in zip(outdict['SPD'], outdict['DIR'])])
    VWind = np.array([-tspd*math.cos((pi/180)*tdir) for tspd, tdir in zip(outdict['SPD'], outdict['DIR'])])
    unitdict['UWIND'] = 'kts'
    unitdict['VWIND'] = 'kts'

    # Get the string version of the month
    months = ["Jan","Feb","Mar","Apr","May","June","July","Aug","Sept","Oct",\
              "Nov","Dec"]
    imonth= int(month)
    str_month = months[imonth-1]

    if(int(year)<70):
        year="20"+year
    else:
        year="19"+year

    title=saver+" "+str_month+" "+day+" "+year+" "+hour+"Z"
    date=day+'-'+str_month+'-'+year
    label= saver+' '+date+' '+hour+'Z'

    # Add the final entries to the dictionary
    outdict['UNITS']=unitdict
    outdict['UWIND'] = UWind
    outdict['VWIND'] = VWind
    outdict['NAME'] = outputstart
    outdict['TITLE']=title
    outdict['LABEL']=label

    return outdict

###############################################################################
#
# readWRF extracts a WRF profile that follows the path of a sounding through
# the model volume. This assumes that the sounding path goes through the model
# volume. Otherwise, the model profile comes from the corner of the volume that
# is closest to the radiosonde's flight path. 
#
# The function requires a WRF output file and a Graw sounding as arguments.
# A dictionary containing the model profile data, as well as a plot title and
# labels, is returned.
#
###############################################################################
def readWRF(wrffile, sounding):
    DIR   = np.array([])
    SPD   = np.array([])
    TEMP  = np.array([])
    DP    = np.array([]) 
    PRESS = np.array([])
    ALT   = np.array([])
    U     = np.array([])
    V     = np.array([])
    TIME  = np.array([])
    UTC   = np.array([])
    LAT   = np.array([])
    LON   = np.array([])

    # WRF sounding
    saver="WRF" 
    tempname=wrffile.split('/')[-1].split('_')
    year =tempname[2][2:4]
    month=tempname[2][4:6]
    day  =tempname[2][6:8]
    hour =tempname[2][8:10]
    minute="00"
    second="00"

    outputstart = year+"_"+month+"_"+day+"_"+hour+"_"+minute+"_"+second+\
                  "_WRF"

    # Read the netCDF file
    data = Dataset(wrffile, "r")

    # Read in the GRAW sounding data with readsounding
    radio = readsounding(sounding)


    # Convert the 3 dimensional netCDF data to np arrays
    T = data.variables['T'][0,:,:,:]
    P = data.variables['P'][0,:,:,:]  # 1-d array of pressure profile
    PB = data.variables['PB'][0,:,:,:]  # 1-d array of pressure profile
    u   = data.variables['U'][0,:,:,:]
    v   = data.variables['V'][0,:,:,:]
    # Perturbation and base-state geopotential are used to find geopotential height
    PH = data.variables['PH'][0,:,:,:]    # Perturbation geopotential
    PHB = data.variables['PHB'][0,:,:,:]  # Base-state geopotential
    xlat = data.variables['XLAT'][0,:]
    xlon = data.variables['XLONG'][0,:]
    xlat_u = data.variables['XLAT_U'][:]
    xlon_u = data.variables['XLONG_U'][:]
    xlat_v = data.variables['XLAT_V'][:]
    xlon_v = data.variables['XLONG_V'][:]
    
    # Calculate T, P, Z, U, and V for the 3-D arrays
    # Calculate geopotential height from geopotential
    hgt = (PH+PHB)/9.81
    # Calculate pressure by adding pressure perturbation to base-state pressure
    p=(PB+P)
    # Convert perturbation potential temperature to temperature
    TH=T+300
    temp = TH*((p/1e5)**(2.0/7.0))

    # Convert data to arrays
    hgt=np.asarray(hgt)
    press=np.asarray(p/100.) 
    temp=np.asarray(temp-273.15) 
    DP  =np.asarray(temp-273.15) # NOTE: change when profiles are made
    # SPD=np.asarray(spd) 
    # DIR=np.asarray(dir)
    u=np.asarray(u*1.94384) 
    v=np.asarray(v*1.94384) 
    
    # Define lists to hold the model data indices for the matching data. These
    # indices will be used to extract data following the path of the radiosonde
    # through the model volume.
    first_dimension_indices    = np.arange(59)
    second_dimension_indices   = [0]*len(press[:,0,0])
    third_dimension_indices    = [0]*len(press[:,0,0])
    u_second_dimension_indices = [0]*len(press[:,0,0])
    u_third_dimension_indices  = [0]*len(press[:,0,0])
    v_second_dimension_indices = [0]*len(press[:,0,0])
    v_third_dimension_indices  = [0]*len(press[:,0,0])
    for i in range(0,len(p[:,0,0])):
        # Find the ballpark model pressure for the given level to look for in the radiosonde data
        model_press = np.mean(press[i,:,:])
        # Find the coordinates in the radiosonde data at the point where the pressure is closest to
        # the ballpark model pressure
        rlat = radio['LAT'][:].flat[np.abs(radio['PRESS'][:]-model_press).argmin()]
        rlon = radio['LON'][:].flat[np.abs(radio['PRESS'][:]-model_press).argmin()]
        # Get the model coordinates that are closest to the radiosonde coordinates
        mlat = xlat[:].flat[np.abs(xlat[:]-rlat).argmin()]
        mlon = xlon[:].flat[np.abs(xlon[:]-rlon).argmin()]

        # Find the index of the close lat/lon in the data
        mlat_ind = int(np.where(xlat[:]==mlat)[0][0])
        mlon_ind = int(np.where(xlon[:]==mlon)[0][0])
        # An older wrf output file has a different format, so account for that
        if(wrffile.strip().split('/')[-1]=='wrfout_d01_2008033100.nc'):
            u_mlat = xlat_u[:].flat[np.abs(xlat_u[:]-rlat).argmin()]
            u_mlon = xlon_u[:].flat[np.abs(xlon_u[:]-rlon).argmin()]
            v_mlat = xlat_v[:].flat[np.abs(xlat_v[:]-rlat).argmin()]
            v_mlon = xlon_v[:].flat[np.abs(xlon_v[:]-rlon).argmin()]
            u_mlat_ind = int(np.where(xlat_u[:]==u_mlat)[0][0])
            u_mlon_ind = int(np.where(xlon_u[:]==u_mlon)[0][0])
            v_mlat_ind = int(np.where(xlat_v[:]==v_mlat)[0][0])
            v_mlon_ind = int(np.where(xlon_v[:]==v_mlon)[0][0])
        # For other files, use the original lat/lon indices
        else:
            u_mlat_ind = v_lat_ind = mlat_ind
            u_mlon_ind = v_lon_ind = mlon_ind
        second_dimension_indices[i]=mlat_ind
        third_dimension_indices[i]=mlon_ind
        u_second_dimension_indices[i]=u_mlat_ind
        u_third_dimension_indices[i]=u_mlon_ind
        v_second_dimension_indices[i]=v_mlat_ind
        v_third_dimension_indices[i]=v_mlon_ind

    # Use the lat and lon indices to extract the profile from the model data   
    ALT = hgt[first_dimension_indices,second_dimension_indices,third_dimension_indices]
    PRESS = press[first_dimension_indices,second_dimension_indices,third_dimension_indices]
    TEMP = temp[first_dimension_indices,second_dimension_indices,third_dimension_indices]
    DP = temp[first_dimension_indices,second_dimension_indices,third_dimension_indices]
    U = u[first_dimension_indices,u_second_dimension_indices,u_third_dimension_indices]
    V = v[first_dimension_indices,v_second_dimension_indices,v_third_dimension_indices]
    DP  =np.asarray(temp-273.15) # NOTE: change when profiles are made
    # SPD=np.asarray(spd) 
    # DIR=np.asarray(dir)

    # Get the string version of the month
    months = ["Jan","Feb","Mar","Apr","May","June","July","Aug","Sept","Oct",\
              "Nov","Dec"]
    imonth= int(month)
    str_month = months[imonth-1]

    if(int(year)<70):
        year="20"+year
    else:
        year="19"+year

    title=saver+" "+str_month+" "+day+" "+year+" "+hour+"Z"
    date=day+'-'+str_month+'-'+year
    label= saver+' '+date+' '+hour+'Z'

    units=dict({'PRESS':'hPa','ALT':'m','TEMP':'degC','DP':'degC','SPD':'kts',\
                'DIR':'deg','U':'kts','V':'kts'})
    # Create a dictionary to hold all the data
    outdict=dict({'NAME': outputstart,'TITLE': title,'LABEL':label,\
                  'UNITS': units,'PRESS': PRESS,'ALT': ALT,'TEMP': TEMP,\
                  'DP': DP,'SPD': SPD,'DIR': DIR,'U': U,'V': V})

    # Make the output dictionary
    outdict=dict({'NAME': 'WRFPROFILE.txt','TITLE': 'WRF PROFILE','LABEL':'no',\
                  'PRESS': PRESS,'ALT': ALT,'TEMP': TEMP,\
                  'DP': DP,'U': U,'V': V})

    return outdict 


# Reads in all data from a GRAW RTS table, not just where the altitude is 
# increasing
def readAllGraw(inputfile):

    readsounding(inputfile,allGraw=True) 
