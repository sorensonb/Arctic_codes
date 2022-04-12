#!/usr/bin/env python

"""
  NAME:
    compare_cn2.py
  
  PURPOSE:
    Calculate the refractive index structure parameters for radiosonde, 
    model, and thermosonde soundings and compare the profiles. 
  
  CALLS:
    Built-in Python modules np, sys, re, matplotlib.pyplot, os
    Custom modules ADPAA and readsounding

  PYTHON VERSION:
    2.7.5 

  IMPORTING:
    In addition to being a normal python script, compare_cn2.py has the ability
    to be imported as a function and be run inside other python scripts. 
    Functions may be imported individually or all at once. 

    Importing all at once:
      To import all three functions from compare_cn2 at once:
        >>> import compare_cn2
        or
        >>> from compare_cn2 import *

      To import it under a different name:
        >>> import compare_cn2 as cn2

      To run the entire program as a function:
        >>> cn2.compare_cn2(<radiosonde file>,<model file>)

      To use a different function (for example, to find the tropopause
      pressure):
        >>> trop = cn2.trop_calc(<pressure array>,<temperature array>,<altitude array>)
    
    Importing individually:
      To import a single function from compare_cn2.py, use this syntax:
        >>> from compare_cn2 import <function>

  USAGE:
    There are two ways to use compare_cn2.py. The first is to run it as a script. 
    To run as a script, use this syntax:
 
      /nas/home/bsorenson/prg/thermosonde1617/compare_cn2.py <radiosonde file> <model file> 

    The other way is to run it as a function imported from a module. 
    See IMPORTING for details

  EXAMPLE:
    /nas/home/bsorenson/prg/thermosonde1617/compare_cn2.py 170605_120000_BIS_UWYO.txt RAP.Op40_2017060512_ANALYSIS_BIS_GSD

    This will calculate the refractive index structure parameter for each file
    and plot each file's cn2 profiles in log units with altitude. The average
    difference and standard deviation of the differences between the log10 cn2
    is printed.
 
  NOTES:
    As of the latest modification, the thermosonde data has not been generated,
    so this program only works with radiosonde and model data.

    This program makes use of readsounding.py to read in the soundings. As 
    described in the comments of readsounding.py, this program is compatible
    with these sounding types:
        - GRAW Profile Data Tables
        - NWS radiosonde soundings (UWYO database)
        - GSD model soundings, including RAP.Op40, RAP.Bak40, NAM, GFS, and FIM
          (rucsoundings.noaa.gov)
        - COAMPS model soundings
        - MCR Snow White radiosonde soundings (CAPE2015)
        - Model soundings from weather.cod.edu, including RAP, NAM, NAMNST, GFS

  MODIFICATIONS:
    Blake Sorenson  - 170712: Written 
    Blake Sorenson <blake.sorenson@und.edu>         - 2018/02/26:
        Updated dictionaries to work with new version of readsounding
  
"""
# Import modules
import numpy as np 
import numpy.ma as ma
import sys
import matplotlib.pyplot as plt
import os
import warnings
import datetime as dt
from adpaa import ADPAA
from readsounding import *
from cn2_lib import *

# trop_calc finds the tropopause for a sounding and returns the tropopause
# pressure
def trop_calc(press, temp, alt):

    var1=5
    length=len(press)
    maxlength=length-var1
    # Skip the lower third of the atmosphere to avoid low level inversions
    starter=int(maxlength/3)
    # By default, the tropopause is first set as the highest sounding level 
    # pressure so that if the definition of a tropopause is not found, 
    # something can still be returned.
    trop_pres=press[-1]
    for i in range(starter, maxlength):
        t1=temp[i+var1] 
        t0=temp[i-var1]
        alt1=alt[i+var1]
        alt0=alt[i-var1]
        tdtdz=(t1-t0)/(alt1-alt0)
        # Definition of tropopause is lowest level at which the lapse rate is
        # -0.002 degC/m or higher
        if((tdtdz>-0.002)):
        #    trop_pres = i
            return press[i]
    return press[-1]

###############################################################################
#
# smooth puts an 11-point smoother on the t, u, and v data and returns
# the smoothed t, u, and v as well as pressures and altitudes that
# correspond to the averaged data. 
#
# There are two ways to call the function. The first way is to pass a 
# dictionary containing pressure, temperature, altitude, u wind, and v
# wind as the only argument. This dictionary can be easily created using 
# the readsounding function from $ADPAA_DIR/src/python_lib/readsounding.py.
# The other way is to pass five separate arrays (pressure, temperature, 
# altitude, u wind, and v wind) to the function. 
#
# The function returns a dictionary containing the smoothed data. If a 
# dictionary was passed to the function, the data arrays are updated with the 
# smoothed data and the same dictionary is returned, containing any other data
# the original dictionary had. If individual arrays were passed, a dictionary
# containing only the smoothed data will be passed back. 
#
###############################################################################
def smooth(first, temp=None, alt=None, u=None, v=None):
    # Determine if a single dictionary was passed
    if(type(first) is dict):
        press = first['PRESS']
        temp = first['TEMP']
        alt = first['ALT']
        u = first['UWIND']
        v = first['VWIND']
    # If the first argument is not a dictionary, then there must be four other
    # arguments. If not, the function cannot continue.
    elif((temp is not None) & (alt is not None) & (u is not None) & \
         (v is not None)):
        press = first
    else:
        print("Invalid function call. See comments")
        return

    # Declare arrays to hold smoothed data
    avg_t = np.array([])
    avg_u = np.array([])
    avg_v = np.array([])
    mid_alt = np.array([])
    mid_press = np.array([])

    # Set the totals equal to the first elements of the t, u, and v. When
    # the for loop started at 0, an extra 0 was showing up at the beginning
    # of the averaged data arrays. Starting the loop at 1 removes the 0, 
    # but requires the first elements of the data arrays to be added
    # to the total before the start.
    total_t=temp[0]
    total_u=u[0]
    total_v=v[0]
    # Loop through the t, u, and v data
    for i in range(1, len(temp)):
        # If 11 elements have been summed, average the current total and 
        # append the averages to arrays. 
        if(i%11==0):
            avg_t = np.append(avg_t[:],total_t/11)
            avg_u = np.append(avg_u[:],total_u/11)
            avg_v = np.append(avg_v[:],total_v/11)
            mid_alt = np.append(mid_alt[:],alt[i-6])
            mid_press = np.append(mid_press[:],press[i-6])
            total_t=0
            total_u=0
            total_v=0
        # Add the current data to the totals
        total_t+=temp[i]
        total_u+=u[i]
        total_v+=v[i]

    # If a dictionary was passed as the only argument, update the data in 
    # the dictionary with the new smoothed data. Make a copy of the dictionary
    # so that the returned version has the edited data, but the original
    # dictionary remains the same.
    if(type(first) is dict):
        outdict = dict(first)
        outdict['PRESS'] = mid_press
        outdict['TEMP'] = avg_t
        outdict['ALT'] = mid_alt
        outdict['UWIND'] = avg_u
        outdict['VWIND'] = avg_v
    # If no dictionary was passed, create a dictionary to hold the smoothed
    # data
    else:
        outdict = dict({'PRESS':mid_press,'TEMP':avg_t,'ALT':mid_alt,'UWIND':avg_u,\
                      'VWIND':avg_v})

    # Return the dictionary
    return outdict 

################################################################################
#
# cn2_calc_thermo calculates the refractive index structure parameter for
# thermosonde data. This function needs the temperature difference data
# measured by the thermosonde as well as associated pressure and temperature
# data from the radiosonde in the thermosonde package. 
#
# There are two ways to call the function. The first way is to pass an array
# of temperature differences and a dictionary containing pressure and 
# temperature as the second argument. While the dictionary can contain other 
# data as well, it must contain at least pressure and temperature. This 
# dictionary can be easily created using the readsounding function from 
# $ADPAA_DIR/src/python_lib/readsounding.py. The other way is to pass the 
# temperature difference array and two separate arrays (pressure and 
# temperature) to the function. 
#
# The function returns a dictionary containing the calculated cn2 data. If a 
# dictionary was passed to the function, the dictionary is updated to contain
# the cn2 and temp diff data. If individual arrays were passed, a dictionary 
# containing the pressure, temperature, temperature difference, and cn2 is 
# returned.  
#
################################################################################
def cn2_calc_thermo(tempdiff, first, temp=None):
    # Determine if a dictionary was passed
    if(type(first) is dict):
        press = first['PRESS']
        temp = first['TEMP']
    # If first is not a dictionary, temp must be a temperature array
    elif((temp is not None)):
        press = first
    else:
        print("Invalid function call. See comments")
        return
    length=len(press)
    cn2 = np.array([]) 
    r=1
    # Loop over the data and calculate cn2
    for i in range(0, length):
        if tempdiff[i] is ma.masked:
            cn2 = np.append(cn2[:],tempdiff[i])
        else:
            p    = press[i]*100
            t    = temp[i]+273.15
            diff = (tempdiff[i])**2

            ########################################################################
            # The actual equation to use is:
            # cn2(h) = (((76*10**-8)*(p/(t**2)))**2)*ct2(h), where
            # ct2(h) = (tempdiff**2)*r**-0.6666667
            #
            ########################################################################
 
                # ct2 = ((tempdiff)**2)*(1**0.6666667)
#                ct2 = ((t1-t0)**2)*z

            # Calculate refractive index structure parameter
            #ct2 = (diff**2)*r**-(2/3)
            ct2 = diff
            x =(((76*10**-8)*(p/(t**2)))**2)*ct2 
            cn2   = np.append(cn2[:], x)

    # If a dictionary was passed, update the dictionary to contain the 
    # temperature differences and cn2
    if(type(first) is dict):
        first['TEMPDIFF'] = tempdiff
        first['CN2T'] = cn2
    # If not, create a dictionary to hold all the data
    else:
        first = dict({'PRESS':press,'TEMP':temp,'TEMPDIFF':tempdiff,'CN2T':cn2})

    # Return the dictionary
    return first

################################################################################
#
# cn2_calc calculates the refractive index structure parameter (cn2) for
# a given set of sounding data. This function may be called by passing
# a dictionary containing pressure, temperature, altitude, u wind, and v wind 
# as the only argument or by passing five individual arrays containing the data
# listed above. The dictionary can be easily created using the readsounding 
# function in $ADPAA_DIR/src/python_lib/readsounding.py. 
# 
# By default, this function assumes that a radiosonde file is being compared 
# with a model file, so the levels immediately above and below the current level
# are used to calculate lapse rate and wind shear. By adding another argument 
# to the function call, method="thermo", the function will use the levels 150 
# meters above and below the current level to calcualte lapse rates and wind 
# shear. However, this requires high resolution soundings, such as GRAW sounding
# files. 
#
# If a dictionary was passed to the function, that dictionary is updated to 
# hold the cn2 as well as the trimmed sounding data (to account for the cn2 
# calculation, the first and last elements of each sounding data array is 
# deleted). If individual arrays were passed, a new dictionary is created
# to hold the data.
#
# EQUATIONS:
#
#   The AFGL model for estimating Cn2 from radiosonde data (Dewan et al, 1993)
#   is used in this function. The equation is: 
#
#     cn2_r = 2.8*(((79*10**-6)*(p/(t**2)))**2)*((dtdz+(9.8*(10**-3)))**2)*l0
#
#   where p is pressure in hPa, t is temperature in Kelvin, dtdz is  
#
#
################################################################################
def cn2_calc(first, temp=None, alt=None, u=None, v=None, method="model"):
    # Determine if the first argument is a dictionary
    if(type(first) is dict):
        press = first['PRESS']
        temp = first['TEMP']
        alt = first['ALT']
        u = first['UWIND']
        v = first['VWIND']
        
        # If the u/v winds are in kts, convert them to meters per second
        if(str(first['UNITS']['UWIND']) == 'kts'):
            u = kts2ms(u)
            first['UNITS']['UWIND']='m/s'
        if((first['UNITS']['VWIND']) == 'kts'):
            v = kts2ms(v)
            first['UNITS']['VWIND']='m/s'
    # If not, the next four must be arrays
    elif((temp is not None) & (alt is not None) & (u is not None) & \
         (v is not None)):
        press = first
    else:
        print("Invalid function call. See comments")

    # Find the tropopause pressure using trop_calc
    trop_pres = trop_calc(press, temp, alt)

    # Subtract one from the length of the data in order to get lapse rates
    # and wind shear
    length=len(press)-1
    cn2 = np.array([]) 
    for i in range(1, length):
        p    = press[i]
        t    = temp[i]+273.15
        # If the model method is being used, the levels directly above and 
        # below the current level are used to calculate lapse rates and     
        # wind shear. 
        if(method=="model"):
            t1   = temp[i+1]+273.15
            t0   = temp[i-1]+273.15
            alt1 = alt[i+1]
            alt0 = alt[i-1]
            u1   = u[i+1]
            u0   = u[i-1]
            v1   = v[i+1]
            v0   = v[i-1]
        # If the thermo method is being used, the levels that are closest to
        # 150 meters above and below the current level are used to calculate
        # lapse rates and wind shear.
        else:
            mid = alt[i]
            # Add and subtract 150 meters to the current altitude 
            upper = mid+150
            lower = mid-150
            # Find the index in the altitude array where the altitude is 
            # closest to 150 meters above the current altitude 
            upperindex = np.where(alt==(alt[:].flat[np.abs(alt[:]-upper).argmin()]))[0][0]
            # Find the index in the altitude array where the altitude is 
            # closest to 150 meters below the current altitude 
            lowerindex = np.where(alt==(alt[:].flat[np.abs(alt[:]-lower).argmin()]))[0][0]
            # If this new index is the same as the current loop index, add
            # one to the index to avoid problems when trying to calculate lapse
            # rates.
            if(upperindex==i):
                upperindex+=1
            # If this new index is the same as the current loop index, subtract
            # one to the index to avoid later problems 
            if(lowerindex==i):
                lowerindex-=1
            # Find the sounding data at the new indices
            t1   = temp[upperindex]+273.15
            t0   = temp[lowerindex]+273.15
            alt1 = alt[upperindex]
            alt0 = alt[lowerindex]
            u1   = u[upperindex]
            u0   = u[lowerindex]
            v1   = v[upperindex]
            v0   = v[lowerindex]

        # Calculate lapse rate
        # Ignore RuntimeWarnings in this step
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            dtdz=(t1-t0)/(alt1-alt0)
            # Calculate wind shear
            S=((((u1-u0)/(alt1-alt0))**2)+(((v1-v0)/(alt1-alt0))**2))**0.5
            if(S>0.04):
                S=0.04
            # Outer limit (l0) calcualation varies slightly depending on 
            # whether or not the sounding is above the tropopause. 
            # If troposphere:
            #if(i<trop_pres):
            if(p>trop_pres):
                l0=(0.1**(4/3))*10**(1.64+42*S)
                if(dtdz>0):
                    dtdz=0
            else:
                l0=(0.1**(4/3))*10**(0.506+50*S)
        # Calculate refractive index structure parameter 
        ########################################################################
        #
        # NOTE: Dry adiabatic lapse rate positive or negative?
        #       According to Dewan et al, set to positive
        #
        ########################################################################
        cn2_r = 2.8*(((79*10**-6)*(p/(t**2)))**2)*((dtdz+(9.8*(10**-3)))**2)*l0
        cn2 = np.append(cn2[:], cn2_r)

    # Since the above loop started at index 1 and stopped one short of the max
    # length of the data, the first and last elements of the sounding data 
    # arrays must be deleted to maintain the same size between all the data
    ##press = np.delete(press,[len(press)-1,0])
    ##temp  = np.delete(temp,[len(temp)-1,0])
    ##alt   = np.delete(alt,[len(alt)-1,0])
    ##u     = np.delete(u,[len(u)-1,0])
    ##v     = np.delete(v,[len(v)-1,0])

    # If a dictionary was passed to the function, update the data in the 
    # dictionary and add the cn2
    if(type(first) is dict):
        for key in first:
            if type(first[key]) is np.ndarray :
                first[key] = np.delete(first[key],[len(first[key])-1,0])
        #first['PRESS'] = press
        #first['TEMP'] = temp
        #first['ALT'] = alt
        #first['UWIND'] = u
        #first['VWIND'] = v
        first['CN2'] = cn2
    # If not, create a dictionary to hold all the data
    else:
        press = np.delete(press,[len(press)-1,0])
        temp  = np.delete(temp,[len(temp)-1,0])
        alt   = np.delete(alt,[len(alt)-1,0])
        u     = np.delete(u,[len(u)-1,0])
        v     = np.delete(v,[len(v)-1,0])
        first = dict({'PRESS':press,'TEMP':temp,'ALT':alt,'UWIND':u,'VWIND':v,'CN2':cn2})
    
    # Return the data
    return first 

################################################################################
#
# plot_cn2 creates a figure for plotting cn2 profiles. Due to the large
# range of the cn2 data (10**-13 near the surface decreasing to 10**-19 in the
# stratosphere), the x axis is plotted in logarithmic units.
#
# Optional arguments may be passed to the function. Any arguments are assumed
# to be dictionaries that contain cn2 data. If arguments have been passed,
# the data will be plotted on the graph and the graph will be displayed.
#
################################################################################
def plot_cn2(*args):
    import matplotlib.pyplot as plt
    plt.figure(1)
    ax = plt.figure(1).add_subplot(111)
    plt.ylim(0, 8)
    leg_loc='upper right'
    # To view cn2 data, axis must be logarithmic
    #plt.xscale('log')
    plt.xlim(-20, -13)
    plt.title('$C_{n}^{2}$ Profile Comparison')
    plt.ylabel('Altitude [km]')
    plt.xlabel("log$_{10}$ [$C_{n}^{2}$ [m$^{-2/3}$]]")

    save_file = False
    filename = ''
    if('save' in args):
        save_file = True
        filename = args[0]['NAME']+'_CN2.png'
        args = np.delete(args,args.index('save'))
    # If arguments have been passed, plot the data from the arguments
    if(len(args) != 0):
        count=0
        colors=['black','red','cyan','green','olive','blue','purple','yellow']
        for arg in args:
            plotalt = arg['ALT']
            if str(arg['UNITS']['ALT']) == 'm':
                plotalt = plotalt/1000.
            try:
                plt.plot(np.log10(arg['CN2']),plotalt,color=colors[count%8],label=arg['LABEL']) 
            except:
                if 'CN2T' in arg.keys():
                    plt.plot(np.log10(arg['CN2T']),plotalt,color=colors[count%8],label='Thermosonde') 
                elif 'CN2' in arg.keys():
                    plt.plot(np.log10(arg['CN2']),plotalt,color=colors[count%8],label=arg['LABEL']) 
                
            count+=1
        plt.legend(loc='upper right', fontsize=12)
        if save_file is True:
            plt.savefig(filename,dpi=300)
            print("Saved image: "+filename)
        plt.show()

#def compare_cn2(radio_file,model_file,thermo_file):
def compare_cn2(radio_file,model_file,thermo_file,mlat=None,mlon=None):

    ###########################################################################
    #
    # NOTE: The final version will take in three files: radiosonde, thermosonde
    #       and analysis (model) profiles. For now, radiosonde and model 
    #       profiles will be used to test the code
    #
    ###########################################################################
    
    RADIO=True
    THERMO=False
    mlat=None
    mlon=None
    #if(len(sys.argv)<3):
    #    print "SYNTAX: ./compare_calc.py <radiosonde profile> <model profile>"
    #    sys.exit(1)
    #if(len(sys.argv)==5):
    #    mlat = sys.argv[-2]
    #    mlon = sys.argv[-1]
    # Read in radiosonde, model, and thermosonde data
    sounding_1 = readsounding(radio_file) # radio_file
    #sounding_1 = readsounding(sys.argv[1]) # radio_file
    if((mlat is None) and (mlon is None)):
        sounding_2 = readsounding(model_file) # model_file
        #sounding_2 = readsounding(sys.argv[2]) # model_file
    else:
        sounding_2 = readsounding(model_file, float(mlat),float(mlon))# model_file
        #sounding_2 = readsounding(sys.argv[2],float(mlat),float(mlon))# model_file

    # Declare Arrays
    out_cn2r     = np.array([])
    in_cn2r      = np.array([])
    out_temp     = np.array([])
    in_temp      = np.array([])
    out_press    = np.array([])
    in_press     = np.array([])
    in_u         = np.array([])
    in_v         = np.array([])
    in_alt       = np.array([])
    tempdiff     = np.array([])
    tm_tempdiff  = np.array([])
    percent_diff = np.array([])
    diff = np.array([])
    
    # Extract date and time information from the radiosonde's readsounding name
    rname   = sounding_1['NAME'].split("_")
    ryear   = rname[0]
    rmonth  = rname[1]
    rday    = rname[2]
    rhour   = rname[3]
    rminute = rname[4]
    rsecond = rname[5].split(".")[0]
    # Convert hh:mm:ss time to seconds from midnight
    rtime   = int(dt.timedelta(hours=int(rhour), minutes=int(rminute), seconds=int(rsecond)).total_seconds()) 
    
    # Extract date and time information from the model's readsounding name
    mname   = sounding_2['NAME'].split("_")
    myear   = mname[0]
    mmonth  = mname[1]
    mday    = mname[2]
    mhour   = mname[3]
    mminute = mname[4]
    msecond = mname[5].split(".")[0]
    # Convert hh:mm:ss time to seconds from midnight
    mtime   = int(dt.timedelta(hours=int(mhour), minutes=int(mminute), seconds=int(msecond)).total_seconds()) 
    
    ############################################################################
    #  
    # COMPARISON: Model vs. Radiosonde
    #
    ############################################################################
    
    ############################################################################
    #
    # NOTE: In most cases, the radiosonde profiles will have more data than the 
    #       model profiles, so the model soundings will be in the model loop most 
    #       of the time. However, this will not be hard coded to account for 
    #       the slim possibility that this does not occur
    #
    ############################################################################
    # Determine which dataset is larger. The larger dataset will be used
    # in the inside loop for determining similar pressure levels
    length1=len(sounding_1['PRESS'])
    length2=len(sounding_2['PRESS'])
    # Convert wind speeds from knots to meters per second
    sounding_1['UWIND'] = kts2ms(sounding_1['UWIND'])
    sounding_1['VWIND'] = kts2ms(sounding_1['VWIND'])
    sounding_2['UWIND'] = kts2ms(sounding_2['UWIND'])
    sounding_2['VWIND'] = kts2ms(sounding_2['VWIND'])
    if(length1>length2):
        radio = dict(sounding_1)
        model = dict(sounding_2)
    else:
        model = dict(sounding_1)
        radio = dict(sounding_2)
    """

    # If radiosonde UTC times don't match up with thermo times, do
    # thermo[-1]-radio[-1] to find offset of end times. Add offset
    # to radio times to get matching times.


    # THERMO: Use time from radiosonde file to match thermosonde data to 
    # radiosonde data
    # Get UTC thermosonde time (either from data itself or from thermosonde
    # data file start time and add on each file time
    starttime = hours+minutes+seconds
    thermo.data['Time'] += starttime
    # Use UTC times to compare with radiosonde data
    #closetime = thermo.data['Time'][:].flat[np.abs(thermo.data['Time'][:]-\
    #                        time).argmin()]

    # Account for the offset in thermosonde and radiosonde UTC time
    time_offset = thermo.data['Time'][-1]-radio['UTC'][-1]
    radio['UTC']+=time_offset
    tempdiff = np.array([])
    voltdiff = np.array([])
    for time in radio['UTC']:
        closetime=thermo.data['Time'][:].flat[np.abs(thermo.data['Time'][:]-\
                  time).argmin()]
        tempdiff=np.append(tempdiff[:],thermo.data['TempDiff'][np.where(\
                    thermo.data['Time']==closetime)])
    """

    """
    fig12 = plt.figure()
    plt.plot(model['ALT'],model['TEMP'],color='black')
    plt.plot(radio['ALT'],radio['TEMP'],color='red')
    plt.title('TEMP')
    fig13 = plt.figure()
    plt.plot(model['ALT'],model['UWIND'],color='black')
    plt.plot(radio['ALT'],radio['UWIND'],color='red')   
    plt.title('UWIND')
    fig14 = plt.figure()
    plt.plot(model['ALT'],model['VWIND'],color='black')
    plt.plot(radio['ALT'],radio['VWIND'],color='red')   
    plt.title('VWIND')
    """

    # Find the tropopause pressure for later. Assume that the soundings are
    # similar, so only the radiosonde tropopause will be used.
    trop_pres = trop_calc(radio['PRESS'],radio['TEMP'],radio['ALT'])
    
    # If one profile has data that goes above the other profile, delete the
    # data from the higher profile that is above the lower profile
    model['TEMP']  = np.delete(model['TEMP'],np.argwhere(model['PRESS']<radio['PRESS'][-1]))
    model['ALT'] = np.delete(model['ALT'],np.argwhere(model['PRESS']<radio['PRESS'][-1]))
    model['UWIND'] = np.delete(model['UWIND'],np.argwhere(model['PRESS']<radio['PRESS'][-1]))
    model['VWIND'] = np.delete(model['VWIND'],np.argwhere(model['PRESS']<radio['PRESS'][-1]))
    model['PRESS'] = np.delete(model['PRESS'],np.argwhere(model['PRESS']<radio['PRESS'][-1]))
    
    # Find similar pressures between the model and radiosonde profiles
    for i in range (len(model['PRESS'])):
        temp_outrpres = 0
        # Use the first pressure in the radio profile as the closest pressure
        # to start with
        temp_innrpres =  radio['PRESS'][0]
        temp_innralt  =  radio['ALT'][0]
        # Loop through radio pressures and compare each pressure to the current
        # model pressure
        for j in range(0, len(radio['PRESS'])):
            # If the current radio pressure is closer to the current model pressure
            # declare the current radio pressure to be the closest pressure
            if(abs(model['PRESS'][i] - radio['PRESS'][j]) < abs(model['PRESS'][i] - temp_innrpres)):
                temp_outrpres = model['PRESS'][i] 
                temp_innrpres = radio['PRESS'][j]
                temp_innralt  = radio['ALT'][j] 
            else:
                continue
    
        # Append the close pressures and altitudes to arrays
        out_press = np.append(out_press[:], temp_outrpres)
        in_press = np.append(in_press[:], temp_innrpres)
        in_alt = np.append(in_alt[:], temp_innralt)
    
    
    # Average the temperature, u wind, and v wind around each pressure level for 
    # the radiosonde data
    inloop = len(in_press)-1
    for i in range(1, inloop):
        # Find where in the radiosonde data the pressure equals the current
        # pressure, the pressure above the current pressure, and the pressure
        # below
        top=int(np.where(radio['PRESS']==in_press[i+1])[0][0])
        mid=int(np.where(radio['PRESS']==in_press[i])[0][0])
        bot=int(np.where(radio['PRESS']==in_press[i-1])[0][0])
        if(top==bot):
            top+=1
        topmid=int((top+mid)/2)
        midbot=int((mid+bot)/2)
        if(midbot==topmid):
            topmid+=1
   
 
        # Average t, u, and v between the matching pressures above and below
        # the current level

        #######################################################################
        #
        # NOTE: Commented out version uses the matching layers above and below 
        # each current layer to find the averaged values. New version uses the 
        # observed levels above and below
        #
        #######################################################################

        # Old
        #addtemp = np.average(radio['TEMP'][midbot:topmid]) 
        #addu = np.average(radio['UWIND'][midbot:topmid]) 
        #addv = np.average(radio['VWIND'][midbot:topmid]) 
        # New
        addtemp = np.average(radio['TEMP'][mid-1:mid+1]) 
        addu = np.average(radio['UWIND'][mid-1:mid+1]) 
        addv = np.average(radio['VWIND'][mid-1:mid+1]) 
        # Once temperature difference measurements are matched with radiosonde
        # data, average the thermosonde measurements for comparison with model
        #adddiff = np.average(thermo.data['TempDiff'][bot:top])
        in_temp = np.append(in_temp[:], addtemp)
        in_u = np.append(in_u[:], addu)
        in_v = np.append(in_v[:], addv)
    
    radio['PRESS'] = in_press
    radio['TEMP']  = in_temp
    radio['ALT']   = in_alt
    radio['UWIND']     = in_u
    radio['VWIND']     = in_v
    
    # Remove the first and last elements of each of these arrays to account for
    # the averaging done above
    model['PRESS'] = np.delete(model['PRESS'],[len(model['PRESS'])-1,0])
    model['TEMP']  = np.delete(model['TEMP'],[len(model['TEMP'])-1,0])
    model['ALT']   = np.delete(model['ALT'],[len(model['ALT'])-1,0])
    model['UWIND']     = np.delete(model['UWIND'],[len(model['UWIND'])-1,0])
    model['VWIND']     = np.delete(model['VWIND'],[len(model['VWIND'])-1,0])
    radio['PRESS'] = np.delete(radio['PRESS'],[len(radio['PRESS'])-1,0])
    radio['ALT']   = np.delete(radio['ALT'],[len(radio['ALT'])-1,0])
    
    """ 
    fig2 = plt.figure()
    plt.plot(model['ALT'],model['TEMP'],color='black')
    plt.plot(radio['ALT'],radio['TEMP'],color='red')
    plt.title('TEMP AVG')
    fig3 = plt.figure()
    plt.plot(model['ALT'],model['UWIND'],color='black')
    plt.plot(radio['ALT'],radio['UWIND'],color='red')   
    plt.title('U AVG')
    fig4 = plt.figure()
    plt.plot(model['ALT'],model['VWIND'],color='black')
    plt.plot(radio['ALT'],radio['VWIND'],color='red')   
    plt.title('V AVG')
    plt.show()
    """
 
    old=False
    if(old==True):
        # Use cn2_calc function to calculate cn2 for each profile
         model = cn2_calc(model)
         radio = cn2_calc(radio)
    else:
        #print "model"
        if(len(out_press)>600):
            model = cn2_calc(model,method="thermo")
        else:
            model = cn2_calc(model)
        #print "radio"
        if(len(in_press)>600):
            radio = cn2_calc(radio,method="thermo")
        else:
            radio = cn2_calc(radio)
    
    # Since cn2 is only being evaluated above the PBL, remove data below 3000 m
    radio['PRESS'] = np.delete(radio['PRESS'], np.argwhere(model['ALT']<3000))
    radio['CN2']   = np.delete(radio['CN2'], np.argwhere(model['ALT']<3000))
    radio['ALT']   = np.delete(radio['ALT'], np.argwhere(model['ALT']<3000))
    model['PRESS'] = np.delete(model['PRESS'], np.argwhere(model['ALT']<3000))
    model['CN2']   = np.delete(model['CN2'], np.argwhere(model['ALT']<3000))
    model['ALT']   = np.delete(model['ALT'], np.argwhere(model['ALT']<3000))
    
    ############################################################################
    #
    # Statistical calculations taken from Frehlich et al, 2010
    #
    ############################################################################
    
    # Compute log10 cn2 arrays
    # Use the new arrays to calculate an average difference and standard devation
    model_log = np.log10(model['CN2'])
    radio_log = np.log10(radio['CN2'])
    rmlog_diff = model_log-radio_log
    rmlog_diff_avg = abs(np.average(rmlog_diff))
    rmlog_diff_std = np.std(rmlog_diff)
   
    # Calculate the average values of the differences for the troposphere 
    # and stratosphere
    trop_cn2 = rmlog_diff[np.where(radio['PRESS']>=trop_pres)]
    trop_cn2_avg = abs(np.average(trop_cn2))
    trop_cn2_std = np.std(trop_cn2) 
    strat_cn2 = rmlog_diff[np.where(radio['PRESS']<trop_pres)]
    strat_cn2_avg = abs(np.average(strat_cn2))
    strat_cn2_std = np.std(strat_cn2) 

    # Calculate percent difference between the two profiles
    for i in range(0, len(model['CN2'])):
        t_dif=abs(model['CN2'][i]-radio['CN2'][i])
        prcnt=(t_dif/max(abs(model['CN2'][i]), abs(radio['CN2'][i])))*100.
    
        diff = np.append(diff[:], t_dif)
        percent_diff = np.append(percent_diff[:], prcnt)
    
    
    # Plot pressure or altitude on the x axis
    alt=True
    if(alt==True):
        out_plotAxis=model['ALT']
        in_plotAxis=radio['ALT']
    else:
        out_plotAxis=list(reversed(model['PRESS']))
        in_plotAxis=list(reversed(radio['PRESS']))
        out_cn2r=list(reversed(model['CN2']))
        in_cn2r=list(reversed(radio['CN2']))
        percent_diff=list(reversed(percent_diff))
        diff=list(reversed(diff))
        plt.xlim(1, 1050)
        leg_loc='upper left'
        axis="Pressure [hPa]"

    # Generate plot
    plot_cn2(radio,model)
    
    # Use the radiosonde and model titles to generate plot labels
    #olabel=model['TITLE'].split(" ")[0]+" "+model['TITLE'].split(" ")[1]+" "+\
    #       model['TITLE'].split(" ")[-1]
    #olabel = model['TITLE'].split(" ")[0]+" Model Analysis"
    #ilabel=radio['TITLE'].split(" ")[0]+" "+radio['TITLE'].split(" ")[1]+" "+\
    #       radio['TITLE'].split(" ")[-1]
    #ilabel = radio['TITLE'].split(" ")[0]+" Radiosonde Analysis"
    # Plot data
    ###plt.plot(model['CN2'], (out_plotAxis/1000.), color='black', \
    ###         label=model['LABEL'])
    ###plt.plot(radio['CN2'], (in_plotAxis/1000.), color='red', \
    ###         label=radio['LABEL'])
    final_prcntdif = np.delete(percent_diff, \
             np.where(np.isnan(percent_diff)))
    print(radio['TITLE']+" vs. "+model['TITLE']+"\t  Log Diff Avg: "+\
          str(round(rmlog_diff_avg, 5)) + "\t  STD: "+str(round(rmlog_diff_std, 5))+\
          "\tAVR PRCNT DIFF: " + str(round(np.average(final_prcntdif), 5)))
    print("Troposphere:  avg_diff = ",round(trop_cn2_avg,5),\
        "  std_dev = ",round(trop_cn2_std,5))
    print("Stratosphere: avg_diff = ",round(strat_cn2_avg,5),\
        "  std_dev = ",round(strat_cn2_std,5))
 
    
    PLOT=False
    if(PLOT==True):
        ax2 = ax.twinx()
        ax2.plot(out_plotAxis, percent_diff, color='blue', label='% diff')
        ax2.set_ylabel("Percent Difference")
        ax2.set_ylim(0, 200)
    #plt.title("$C_{n}^{2}$ Profiles for Bismarck, ND at 00Z 12 June 2017")
    #plt.title("June 12 2017 00Z Bismarck $C_{n}^{2}$ Radiosonde and Model Profiles")
    ###plt.legend(loc='upper right', fontsize=12)
    ###plt.show()
 
    ###########################################################################
    #   
    # COMPARISON: Thermosonde vs. Radiosonde
    #
    ###########################################################################
    if(thermo_file is not None):
        THERMO=True
        thermo=ADPAA()
        thermo.ReadFile(thermo_file)
        thermo.mask_MVC()
   
        # Reset the radio dictionary 
        radio = readsounding(radio_file,allGraw=True,keepMissing=True)
    
        from cn2_lib import read_temp_diffs
        thermo_scn2 = read_temp_diffs(radio_file, thermo_file)
        ##!################################################################################
        ##!##
        ##!## TRANSPLANTED FROM thermo_test.py (and then from comparison_checker.py)
        ##!##
        ##!## Stretch the thermosonde times to match the radiosonde's times. 
        ##!##
        ##!################################################################################
        ##!#
        ##!## Delete radiosonde data where the altitudes are missing (only at the end of the file)
        ##!#radio['Time'] = np.delete(radio['Time'],np.where(radio['ALT']==999999.9999))
        ##!#radio['UTCTime'] = np.delete(radio['UTCTime'],np.where(radio['ALT']==999999.9999))
        ##!#radio['TEMP'] = np.delete(radio['TEMP'],np.where(radio['ALT']==999999.9999))
        ##!#radio['PRESS'] = np.delete(radio['PRESS'],np.where(radio['ALT']==999999.9999))
        ##!#radio['UWIND'] = np.delete(radio['UWIND'],np.where(radio['ALT']==999999.9999))
        ##!#radio['VWIND'] = np.delete(radio['VWIND'],np.where(radio['ALT']==999999.9999))
        ##!#radio['SPD'] = np.delete(radio['SPD'],np.where(radio['ALT']==999999.9999))
        ##!#radio['DIR'] = np.delete(radio['DIR'],np.where(radio['ALT']==999999.9999))
        ##!#radio['ALT'] = np.delete(radio['ALT'],np.where(radio['ALT']==999999.9999))
        ##!#
       
        ##!##----------------------------------------------------------------------
        ##!##
        ##!## NOTE: These starting and stopping times will vary depending on the
        ##!##       launch, so change for each launch.
        ##!##
        ##!##       The values used here are for the first full launch on 
        ##!##       May 4th and 5th, 2018.
        ##!##---------------------------------------------------------------------- 
        ##!###!## Thermo data (Values for the 2019/05/04 launch)
        ##!###!#start_time = 13152.0
        ##!###!#stop_time = 14220.0 
        ##!###!## Radio data
        ##!###!#radio_ascent_start = 13334.7
        ##!###!#radio_ascent_stop  = 14410.7
        ##!###!#radio_last_contact = 14410.7
        ##!###!##radio_last_index = 7734 # last contact
        ##!###!#radio_last_index = 1076 # burst
        ##!###!##thermo_last_index = 10226 # last contact 
        ##!###!#thermo_start_index= 2622
        ##!###!#thermo_last_index = 3690
        ##!###!#thermo_final_last_index = 3918  # last index in the file
        ##!###!##radio_last_contact = 26383.0

        ##!###!## Values for the 2018/05/05 launch
        ##!###!##### Thermo data 
        ##!###!####start_time = 23331.0
        ##!###!####stop_time = 27380.0 
        ##!###!##### Radio data
        ##!###!####radio_ascent_start = 20154.5
        ##!###!####radio_ascent_stop  = 24628.5
        ##!###!####radio_last_contact = 26398.5
        ##!###!#####radio_last_index = 7734 # last contact
        ##!###!####radio_last_index = 5964 # burst
        ##!###!#####thermo_last_index = 10226 # last contact 
        ##!###!####thermo_last_index = 8675
        ##!###!####thermo_final_last_index = 10232  # last index in the file
        ##!###!#####radio_last_contact = 26383.0

        ##!#if(radio_file.strip().split('/')[-1][:2] == '19'):
        ##!#    pre2019 = False
        ##!#else:
        ##!#    pre2019 = True
        ##!#if(pre2019 == True):
        ##!#    # Thermo data
        ##!#    start_time = 23331.0
        ##!#    stop_time = 27380.0 
        ##!#    # Radio data
        ##!#    radio_ascent_start = 20154.5
        ##!#    radio_ascent_stop  = 24628.5
        ##!#    radio_last_contact = 26398.5
        ##!#     
        ##!#    radio_last_index = 5964 # burst
        ##!#
        ##!#    thermo_start_index = 4626
        ##!#    thermo_last_index = 8675
        ##!#    thermo_final_last_index = 10232  # last index in the file
        ##!#else:
        ##!#    # Thermo data
        ##!#    start_time = 13152.0
        ##!#    stop_time = 14220.0 
        ##!#    # After 2019/05/03
        ##!#    radio_ascent_start = 13334.7
        ##!#    radio_ascent_stop  = 14410.7
        ##!#    radio_last_contact = 14410.7
        ##!#
        ##!#    # After 2019/05/03
        ##!#    radio_last_index = 1076 # burst
        ##!#
        ##!#    thermo_start_index= 2622
        ##!#    thermo_last_index = 3690
        ##!#    thermo_final_last_index = 3918  # last index in the file


        ##!#radio_start_time  = radio_ascent_start 
        ##!#radio_last_time  = radio_last_contact 
        ##!##radio_diff = (radio_last_time-radio_start_time)/(7054.0-1490.0)
        ##!#radio_diff = (radio_last_time-radio_start_time)/(np.where(radio['UTCTime']==radio_last_time)[0][0]-np.where(radio['UTCTime']==radio_start_time)[0][0])
        ##!#
        ##!#thermo_diff = 1.0/radio_diff
        ##!#
        ##!#diff = start_time-radio_start_time # first launch 2018 05 05 
        ##!##plotTime = radio['UTCTime']+diff
        ##!#
        ##!## Match the launch times between the thermosonde and radiosonde
        ##!#plotTime = thermo.data['Time']-diff
        ##!#thermo.data['Time'] = thermo.data['Time']-diff
        ##!#
        ##!#ascent_rtime = radio['UTCTime'][0:radio_last_index]
        ##!#ascent_ralt = radio['ALT'][0:radio_last_index]
        ##!#ascent_ttime = plotTime[thermo_start_index:thermo_last_index]
        ##!#ascent_talt = thermo.data['Alt'][thermo_start_index:thermo_last_index]
        ##!## First full launch
        ##!##ascent_rtime = radio['UTCTime'][1490:radio_last_index]
        ##!##ascent_ralt = radio['ALT'][1490:radio_last_index]
        ##!##ascent_ttime = plotTime[4626:thermo_last_index]
        ##!##ascent_talt = thermo.data['Alt'][4626:thermo_last_index]
        ##!#
        ##!####descent_rtime = radio['UTCTime'][radio_last_index:len(radio['UTCTime'])]
        ##!####descent_ralt = radio['ALT'][radio_last_index:len(radio['UTCTime'])]
        ##!####descent_ttime = plotTime[thermo_last_index:thermo_final_last_index]
        ##!#
        ##!## Find the ratio of the number of radiosonde ascent times to the number of
        ##!## thermosonde ascent times
        ##!#ascent_diff_ratio = float(len(ascent_rtime))/float(len(ascent_ttime))
        ##!## Use that ratio to create new thermosonde ascent times that match with the
        ##!## radiosonde's ascent times. Now, the radiosonde and thermosonde ascent times
        ##!## have launches at the same time and bursts at the same time, although the 
        ##!## burst altitudes differ by a few dozen meters.
        ##!#ascent_ttime = np.arange(ascent_ttime[0],ascent_rtime[-1]+1,ascent_diff_ratio)
        ##!#
        ##!###!#descent_rtime = radio['UTCTime'][radio_last_index:len(radio['UTCTime'])]
        ##!###!#descent_ralt = radio['ALT'][radio_last_index:len(radio['UTCTime'])]
        ##!#descent_ttime = plotTime[thermo_last_index:thermo_final_last_index]
        ##!## Repeat the process for the descent
        ##!#descent_rtime = radio['UTCTime'][radio_last_index:len(radio['UTCTime'])]
        ##!#descent_ralt = radio['ALT'][radio_last_index:len(radio['UTCTime'])]
        ##!##descent_ttime = plotTime[thermo_last_index:len(thermo.data['Time'])]
        ##!#descent_talt = thermo.data['Alt'][thermo_last_index:thermo_final_last_index]
        ##!#
        ##!#descent_diff_ratio = float(len(descent_rtime))/float(len(descent_ttime)-1.0)
        ##!#descent_ttime = np.arange(ascent_ttime[-1]+1,descent_rtime[-1]+1,descent_diff_ratio)
        ##!#
        ##!#rplot = np.arange(0,len(ascent_rtime))
        ##!#tplot = np.arange(0,len(ascent_ttime))

        ##!#if(pre2019 == True):
        ##!#    # Repeat the process for the descent
        ##!#    descent_rtime = radio['UTCTime'][radio_last_index:len(radio['UTCTime'])]
        ##!#    descent_ralt = radio['ALT'][radio_last_index:len(radio['UTCTime'])]
        ##!#    #descent_ttime = plotTime[thermo_last_index:len(thermo.data['Time'])]
        ##!#    descent_talt = thermo.data['Alt'][thermo_last_index:thermo_final_last_index]
        ##!#    
        ##!#    descent_diff_ratio = float(len(descent_rtime))/float(len(descent_ttime)-1.0)
        ##!#    descent_ttime = np.arange(ascent_ttime[-1]+1,descent_rtime[-1]+1,descent_diff_ratio)
        ##!#    # Combine the new thermosonde times into a single array
        ##!#    ad_ttime = np.concatenate([ascent_ttime,descent_ttime])
        ##!#    ad_talt  = np.concatenate([ascent_talt,descent_talt])
        ##!#
        ##!#    # Before 2019/05/03
        ##!#    thermo.data['Time'][thermo_start_index:thermo_final_last_index] = ad_ttime
        ##!#    thermo.data['Alt'][thermo_start_index:thermo_final_last_index] = ad_talt
        ##!#else:
        ##!#    # After 2019/05/03
        ##!#    thermo.data['Time'][thermo_start_index:thermo_last_index] = ascent_ttime
        ##!#    thermo.data['Alt'][thermo_start_index:thermo_last_index] = ascent_talt
        ##!#
        ##!##### Combine the new thermosonde times into a single array
        ##!####ad_ttime = np.concatenate([ascent_ttime,descent_ttime])
        ##!####ad_talt  = np.concatenate([ascent_talt,descent_talt])
        ##!#
        ##!###!#thermo.data['Time'][thermo_start_index:thermo_last_index] = ascent_ttime
        ##!###!#thermo.data['Alt'][thermo_start_index:thermo_last_index] = ascent_talt
        ##!##thermo.data['Time'][4626:10232] = ad_ttime
        ##!##thermo.data['Alt'][4626:10232] = ad_talt
        ##!#thermo.mask_MVC()
        ##!#
        ##!################################################################################
        ##!##
        ##!## END OF comparison_checker.py TRANSPLANT
        ##!##
        ##!################################################################################
        ##!#
        ##!## Account for the fact that the 
        ##!####time_offset = thermo.data['Time'][-1]-radio['UTCTime'][-1]
        ##!##thermo_start_time = 86332.0
        ##!##thermo_stop_time = 86606.0   # these three are for the tethered test
        ##!##radio_start_time  = 83093.0
        ##!#thermo_start_time = 23331.0
        ##!#thermo_stop_time = 27380.0   # these three are for the tethered test
        ##!#radio_start_time  = 20154.0
        ##!#
        ##!#combined_start_time = 13334.7
        ##!#combined_stop_time = 14410.7
        ##!## First launch values
        ##!###combined_start_time = 20154.5
        ##!###combined_stop_time = 24629.5
        ##!#
        ##!## This part is made obselete by the comparison_checker.py transplant
        ##!#####time_offset = thermo_start_time-radio_start_time-23.0 # tethered test
        ##!####time_offset = thermo_start_time-radio_start_time
        ##!#####radio['UTCTime']+=time_offset
        ##!####thermo.data['Time']-=time_offset
        ##!#
        ##!## Grab the "matching" times to plot a time series of the data later
        ##!#closetime_thermo = np.array([])
        ##!#matchalt_thermo = np.array([])
        ##!#closetime_radio = np.array([])
        ##!#
        ##!## Get rid of data that is outside the ascent time
        ##!##closeindex_radio = np.where((radio['UTCTime']>=86322.0) & (radio['UTCTime']<=86606.0))[0] # tethered test
        ##!#closeindex_radio = np.where((radio['UTCTime']>=combined_start_time) & (radio['UTCTime']<=combined_stop_time))[0]
        ##!#for key in radio.keys():
        ##!#    if((key != 'UNITS') & (type(radio[key]) is not str)):
        ##!#        radio[key] = radio[key][closeindex_radio]
        ##!#
        ##!#tempdiff = np.array([])
        ##!#for time in radio['UTCTime']:
        ##!#    closetime=thermo.data['Time'][:].flat[np.abs(thermo.data['Time'][:]-time).argmin()]
        ##!#    closetime_thermo = np.append(closetime_thermo[:],closetime)
        ##!#    close_index = np.where(thermo.data['Time']==closetime)[0]
        ##!#    matchalt_thermo = np.append(matchalt_thermo[:],thermo.data['Alt'][close_index])
        ##!##    print radio['ALT'][np.where(radio['UTCTime']==time)],thermo.data['Alt'][close_index]
        ##!##    print thermo.data['TempDiff'][np.where(thermo.data['Time']==closetime)]
        ##!#    tempdiff = np.append(tempdiff[:],thermo.data['TempDiff'][np.where(thermo.data['Time']==closetime)])
        ##!#
        ##!#
        ##!#
        ##!#thermo_cn2 = dict(radio)
        ##!#thermo_cn2 = cn2_calc_thermo(tempdiff,thermo_cn2)
        ##!## Adding the method='thermo' flag causes the radiosonde Cn2 to be several
        ##!## of magnitude higher than without
        sradio = dict(radio)
        sradio = smooth(sradio)
        sradio = cn2_calc(sradio,method='thermo')
        radio = cn2_calc(radio)
        ##!#
        ##!## Ignore the missing thermosonde_data
        ##!#masked_tempdiff = ma.masked_values(tempdiff,1e6)
        ##!#masked_indices = np.where(masked_tempdiff!=ma.masked)
        ##!#for key in thermo_cn2:
        ##!#    if type(thermo_cn2[key]) is np.ndarray:
        ##!#        thermo_cn2[key] = thermo_cn2[key][masked_indices]
        ##!#
        ##!## SMOOTHER
        ##!## Put an 11-point smoother on the thermosonde data
        ##!## Declare arrays to hold smoothed data
        ##!#avg_t = np.array([])
        ##!#avg_u = np.array([])
        ##!#avg_v = np.array([])
        ##!#mid_alt = np.array([])
        ##!#mid_press = np.array([])
        ##!#avg_cn2t = np.array([])
        ##!#temp_cn2t = np.zeros(11)
        ##!#j=0
        ##!## Set the totals equal to the first elements of the t, u, and v. When
        ##!## the for loop started at 0, an extra 0 was showing up at the beginning
        ##!## of the averaged data arrays. Starting the loop at 1 removes the 0, 
        ##!## but requires the first elements of the data arrays to be added
        ##!## to the total before the start.
        ##!#
        ##!## Convert thermo_cn2['CN2T'] to logarithmic for smoothing
        ##!#thermo_scn2 = dict(thermo_cn2)
        ##!#thermo_scn2['LABEL'] = 'Smoothed Thermosonde'
        ##!#thermo_scn2['CN2T'] = np.log10(thermo_scn2['CN2T'])
        ##!#total_t=thermo_scn2['TEMP'][0]
        ##!#total_cn2t=thermo_scn2['CN2T'][0]
        ##!#total_u=thermo_scn2['UWIND'][0]
        ##!#total_v=thermo_scn2['VWIND'][0]
        ##!## Loop through the t, u, and v data
        ##!#for i in range(1, len(thermo_cn2['CN2T'])):
        ##!#    # If 11 elements have been summed, average the current total and 
        ##!#    # append the averages to arrays. 
        ##!#    if(i%11==0):
        ##!#        avg_t = np.append(avg_t[:],total_t/11)
        ##!#        avg_u = np.append(avg_u[:],total_u/11)
        ##!#        avg_v = np.append(avg_v[:],total_v/11)
        ##!#        mid_alt = np.append(mid_alt[:],thermo_cn2['ALT'][i-6])
        ##!#        mid_press = np.append(mid_press[:],thermo_cn2['PRESS'][i-6])
        ##!#        #avg_cn2t = np.append(avg_cn2t[:],total_cn2t/11)
        ##!#        avg_cn2t = np.append(avg_cn2t[:],np.average(temp_cn2t))
        ##!#        j=0
        ##!#        total_t=0
        ##!#        total_u=0
        ##!#        total_v=0
        ##!#    # Add the current data to the totals
        ##!#    total_t+=thermo_scn2['TEMP'][i]
        ##!#    #total_cn2t+=thermo_cn2['CN2T'][i]
        ##!#    temp_cn2t[j] = thermo_scn2['CN2T'][i]
        ##!#    j+=1
        ##!#    total_u+=thermo_scn2['UWIND'][i]
        ##!#    total_v+=thermo_scn2['VWIND'][i]
        ##!#
        ##!## REMOVE to prevent resetting the data with the smoothed data
        ##!#thermo_scn2['CN2T']=avg_cn2t
        ##!#thermo_scn2['ALT']=mid_alt
        ##!#thermo_scn2['PRESS'] = mid_press
        ##!#thermo_scn2['TEMP'] = avg_t
        ##!#thermo_scn2['UWIND`'] = avg_u
        ##!#thermo_scn2['VWIND`'] = avg_v
        ##!#
        ##!## Convert thermo_cn2['CN2T'] back to linear
        ##!#thermo_scn2['CN2T'] = 10.**(thermo_scn2['CN2T'])
       
        """
        #######################################################################
        #
        # Normal comparison section for radiosonde vs thermosonde
        #
        #######################################################################
        # For each smoothed thermo value, find the corresponding altitude in
        # the radiosonde data and average the radiosonde data around the 
        # thermosonde altitude.
        close_indices = [] 
        for talt in thermo_scn2['ALT']:
            closealt=radio['ALT'][:].flat[np.abs(radio['ALT'][:]-talt).argmin()]
            close_indices.append(int(np.where(radio['ALT']==closealt)[0][0]))


        # Average the temperature, u wind, and v wind around each pressure level for 
        # the radiosonde data
        inloop = len(close_indices)-1
        in_temp = np.array([])
        in_press = np.array([])
        in_alt = np.array([])
        in_u = np.array([])
        in_v = np.array([])
        for i in range(1, inloop):
            # Find where in the radiosonde data the pressure equals the current
            # pressure, the pressure above the current pressure, and the pressure
            # below
            ###top=int(np.where(radio['PRESS']==in_press[i+1])[0][0])
            ###mid=int(np.where(radio['PRESS']==in_press[i])[0][0])
            ###bot=int(np.where(radio['PRESS']==in_press[i-1])[0][0])
            ###if(top==bot):
            ###    top+=1
            ###topmid=int((top+mid)/2)
            ###midbot=int((mid+bot)/2)
            ###if(midbot==topmid):
            ###    topmid+=1
   
 
            # Average t, u, and v between the matching pressures above and below
            # the current level

            #######################################################################
            #
            # NOTE: Commented out version uses the matching layers above and below 
            # each current layer to find the averaged values. New version uses the 
            # observed levels above and below
            #
            #######################################################################

            # Old
            #addtemp = np.average(radio['TEMP'][midbot:topmid]) 
            #addu = np.average(radio['UWIND'][midbot:topmid]) 
            #addv = np.average(radio['VWIND'][midbot:topmid]) 
            # New
            addtemp = np.average(radio['TEMP'][close_indices[i-1]:close_indices[i+1]]) 
            addu = np.average(radio['UWIND'][close_indices[i-1]:close_indices[i+1]]) 
            addv = np.average(radio['VWIND'][close_indices[i-1]:close_indices[i+1]]) 
            # Once temperature difference measurements are matched with radiosonde
            # data, average the thermosonde measurements for comparison with model
            #adddiff = np.average(thermo.data['TempDiff'][bot:top])
            in_temp = np.append(in_temp[:], addtemp)
            in_u = np.append(in_u[:], addu)
            in_v = np.append(in_v[:], addv)
            in_press = np.append(in_press[:], radio['PRESS'][close_indices[i]])
            in_alt = np.append(in_alt[:], radio['ALT'][close_indices[i]])
        
        radio['PRESS'] = in_press
        radio['TEMP']  = in_temp
        radio['ALT']   = in_alt
        radio['UWIND']     = in_u
        radio['VWIND']     = in_v


        # Remove the first and last elements of each of these arrays to account for
        # the averaging done above
        thermo_scn2['PRESS'] = np.delete(thermo_scn2['PRESS'],[len(thermo_scn2['PRESS'])-1,len(thermo_scn2['PRESS'])-2,0,1])
        thermo_scn2['TEMP']  = np.delete(thermo_scn2['TEMP'], [len(thermo_scn2['TEMP'])-1,len(thermo_scn2['TEMP'])-2,0,1])
        thermo_scn2['CN2T']  = np.delete(thermo_scn2['CN2T'], [len(thermo_scn2['CN2T'])-1,len(thermo_scn2['CN2T'])-2,0,1])
        thermo_scn2['ALT']   = np.delete(thermo_scn2['ALT'],  [len(thermo_scn2['ALT'])-1,len(thermo_scn2['ALT'])-2,0,1])
        thermo_scn2['UWIND']     = np.delete(thermo_scn2['UWIND'],  [len(thermo_scn2['UWIND'])-1,len(thermo_scn2['UWIND'])-2,0,1])
        thermo_scn2['VWIND']     = np.delete(thermo_scn2['VWIND'],  [len(thermo_scn2['VWIND'])-1,len(thermo_scn2['VWIND'])-2,0,1])

        radio = cn2_calc(radio,method='thermo')

        # Since cn2 is only being evaluated above the PBL, remove data below 3000 m
        radio['PRESS'] = np.delete(radio['PRESS'], np.argwhere(radio['ALT']<3000))
        radio['CN2']   = np.delete(radio['CN2'], np.argwhere(radio['ALT']<3000))
        radio['ALT']   = np.delete(radio['ALT'], np.argwhere(radio['ALT']<3000))
        thermo_scn2['PRESS'] = np.delete(thermo_scn2['PRESS'], np.argwhere(thermo_scn2['ALT']<3000))
        thermo_scn2['CN2T']   = np.delete(thermo_scn2['CN2T'], np.argwhere(thermo_scn2['ALT']<3000))
        thermo_scn2['ALT']   = np.delete(thermo_scn2['ALT'], np.argwhere(thermo_scn2['ALT']<3000))
      

        """ 
        #----------------------------------------------------------------------
        #
        # The comparison methods taken for the radiosonde-model section
        # do not yield good results, so experimenting with using the normally-
        # smoothed radiosonde data and finding the closest altitudes
        #
        #----------------------------------------------------------------------
        # Since cn2 is only being evaluated above the PBL, remove data below 3000 m
        #sradio['PRESS'] = np.delete(sradio['PRESS'], np.argwhere(sradio['ALT']<3000))
        #sradio['CN2']   = np.delete(sradio['CN2'], np.argwhere(sradio['ALT']<3000))
        #sradio['ALT']   = np.delete(sradio['ALT'], np.argwhere(sradio['ALT']<3000))
        #thermo_scn2['PRESS'] = np.delete(thermo_scn2['PRESS'], np.argwhere(thermo_scn2['ALT']<3000))
        #thermo_scn2['CN2T']   = np.delete(thermo_scn2['CN2T'], np.argwhere(thermo_scn2['ALT']<3000))
        #thermo_scn2['ALT']   = np.delete(thermo_scn2['ALT'], np.argwhere(thermo_scn2['ALT']<3000))
        
        close_indices = [] 
        for talt in thermo_scn2['ALT']:
            closealt=sradio['ALT'][:].flat[np.abs(sradio['ALT'][:]-talt).argmin()]
            close_indices.append(int(np.where(sradio['ALT']==closealt)[0][0]))

        sradioCn2 = sradio['CN2'][close_indices] 
        sradioAlt = sradio['ALT'][close_indices] 
        sradio['CN2'] = sradioCn2
        sradio['ALT'] = sradioAlt
   
 
        ##!#########################################################################
        ##!##
        ##!## Statistical calculations taken from Frehlich et al, 2010
        ##!##
        ##!#########################################################################
        ##!#
        ##!## Compute log10 cn2 arrays
        ##!## Use the new arrays to calculate an average difference and standard
        ##!## devation
        ##!#thermo_log = np.log10(thermo_scn2['CN2T'])
        ##!#radio_log = np.log10(sradio['CN2'])
        ##!##radio_log = np.log10(radio['CN2'])
        ##!#rtlog_diff = thermo_log-radio_log
        ##!#rtlog_diff_avg = abs(np.average(rtlog_diff))
        ##!#rtlog_diff_std = np.std(rtlog_diff)

        ##!## Calculate the average values of the differences for the troposphere 
        ##!## and stratosphere
        ##!#trop_cn2 = rtlog_diff[np.where(thermo_scn2['PRESS']>=trop_pres)]
        ##!#num_trop_cn2 = len(trop_cn2)
        ##!#trop_cn2_avg = abs(np.average(trop_cn2))
        ##!#trop_cn2_std = np.std(trop_cn2) 
        ##!#strat_cn2 = rtlog_diff[np.where(thermo_scn2['PRESS']<trop_pres)]
        ##!#num_strat_cn2 = len(strat_cn2)
        ##!#strat_cn2_avg = abs(np.average(strat_cn2))
        ##!#strat_cn2_std = np.std(strat_cn2) 
       
        ##!#diff = np.array([]) 
        ##!#percent_diff = np.array([]) 
        ##!## Calculate percent difference between the two profiles
        ##!#for i in range(0, len(thermo_scn2['CN2T'])):
        ##!#    t_dif=abs(thermo_scn2['CN2T'][i]-sradio['CN2'][i])
        ##!#    #t_dif=abs(thermo_scn2['CN2T'][i]-radio['CN2'][i])
        ##!#    prcnt=(t_dif/max(abs(thermo_scn2['CN2T'][i]), abs(sradio['CN2'][i])))*100.
        ##!#    #prcnt=(t_dif/max(abs(thermo_scn2['CN2T'][i]), abs(radio['CN2'][i])))*100.
        ##!#
        ##!#    diff = np.append(diff[:], t_dif)
        ##!#    percent_diff = np.append(percent_diff[:], prcnt)

        ##!#final_prcntdif = np.delete(percent_diff, \
        ##!#                 np.where(np.isnan(percent_diff)))
        ##!#print(radio['TITLE']+" vs. Thermosonde"+"\t  Log Diff Avg: "+\
        ##!#      str(round(rtlog_diff_avg, 5)) + "\t  STD: "+str(round(rtlog_diff_std, 5))+\
        ##!#      "\tAVR PRCNT DIFF: " + str(round(np.average(final_prcntdif), 5)))
        ##!#print("Troposphere:  avg_diff = ",round(trop_cn2_avg,5),\
        ##!#    "  std_dev = ",round(trop_cn2_std,5))
        ##!#print("Stratosphere: avg_diff = ",round(strat_cn2_avg,5),\
        ##!#    "  std_dev = ",round(strat_cn2_std,5))
        ##!#print("Num_obs troposphere=",num_trop_cn2)
        ##!#print("Num_obs stratosphere=",num_strat_cn)

        # Plot the thermosonde Cn2 and radiosonde Cn2 on a graph
        plot_cn2(thermo_scn2,sradio)
        #plot_cn2(thermo_scn2,thermo_scn2)

        #plot_cn2()
        #plt.plot(np.log10(radio['CN2']),radio['ALT']/1000.)
        #plt.plot(np.log10(sradio['CN2']),sradio['ALT']/1000.)
        #plt.plot(np.log10(thermo_cn2['CN2T']),thermo_cn2['ALT']/1000.)


    """
    
    # Smooth thermosonde measurements?? 
    smooth_diff = np.array([]) 
    total_diff=tempdiff[0]
    for i in range(1, len(tempdiff)):
        # If 11 elements have been summed, average the current total and 
        # append the totals to arrays. 
        if(i%11==0):
            smooth_diff = np.append(avg_t[:],total_diff/11)
            total_diff=0
        # Add the current data to the totals
        total_t+=tempdiff[i]
    
    # 11 point smoother for T, u, v
    smooth_r = smooth(radio)
    #smooth_r['PRESS'], smooth_r['TEMP'], smooth_r['ALT'], su, sv = smooth(\
    #                               radio['PRESS'],radio['TEMP'],\
    #                               radio['ALT'],radio['UWIND'],radio['VWIND'])
    
    # Calculate cn2 for thermosonde and radiosonde data using P, T, u, v, z    
    smooth_r = cn2_calc(smooth_r,method="thermo")
    smooth_r = cn2_calc_thermo(smooth_diff,smooth_r)
    #smooth_r['ALT'] = np.delete(smooth_r['ALT'],-1)
    #smooth_r['ALT'] = np.delete(smooth_r['ALT'],0)
    smooth_r['CN2'] = np.delete(smooth_r['CN2'], \
                      np.where(smooth_r['ALT']<3000))
    smooth_r['CN2T'] = np.delete(smooth_r,['CN2T'] np.where(smooth_r['ALT']<3000))
    smooth_r['ALT']         = np.delete(smooth_r['ALT'], np.where(smooth_r['ALT']<3000))
    
    # Calculate average log difference and standard deviation of log difference
    smooth_r_log = np.log10(smooth_cn2_r)
    smooth_t_log = np.log10(smooth_cn2_t)
    rtlog_diff = smooth_r_log-smooth_radio_log
    rtlog_diff_avg = abs(np.average(rtlog_diff))
    rtlog_diff_std = np.std(rtlog_diff)
    ###########################################################################
    
    # This section is only to test the smoothing function 
    
    # Smooth the high-resolution sounding data
    tempdict = dict(radio)
    rawdict = dict(radio)
    smooth_rm = smooth(tempdict)
    
    # Calculate raw cn2 using the thermosonde approach 
    raw = cn2_calc(rawdict,method="thermo")
    
    # If the sounding file is very large (Graw high-density sounding), use the
    # thermosonde approach, in which the smoothed levels closest to 150 
    # meters above and below the current level are used to calculate
    # lapse rates and wind shear 
    if(len(radio['PRESS'])>1500):
        smooth_rm = cn2_calc(smooth_rm,method="thermo") 
    else:
        smooth_rm = cn2_calc(smooth_rm) 
    smooth_rm['CN2'] = np.delete(smooth_rm['CN2'], np.where(smooth_rm['ALT']<3000))
    smooth_rm['ALT'] = np.delete(smooth_rm['ALT'], np.where(smooth_rm['ALT']<3000))
    raw['CN2'] = np.delete(raw['CN2'], np.where(raw['ALT']<3000))
    raw['ALT'] = np.delete(raw['ALT'], np.where(raw['ALT']<3000))
    
    plot_cn2()
    plt.plot(raw['CN2'],(raw['ALT']/1000.),color='black',label="new raw")
    plt.plot(smooth_rm['CN2'],(smooth_rm['ALT']/1000.),color='red',label="new smooth")
    plt.legend()
    #plt.show()
    
    
    ###########################################################################
    #   
    # COMPARISON: Thermosonde vs. Model
    #
    ###########################################################################
    
    # Average thermosonde data around pressure points??
    
    tempdiff = thermo.data['TEMPDIFF'][match_indices]
    
    # 11 point smoother for T, u, v
    smooth_t = smooth(radio)
    
    # Calculate cn2 for thermosonde and radiosonde data using P, T, u, v, z    
    smooth_t = cn2_calc(smooth_t,"thermo")
    smooth_t = cn2_calc_thermo(tempdiff,smooth_t)
    #smooth_t['ALT'] = np.delete(smooth_t['ALT'],-1)
    #smooth_t['ALT'] = np.delete(smooth_t['ALT'],0)
    
    # Calculate average log difference and standard deviation of log difference
    smooth_r_log = np.log10(smooth_t['CN2'])
    smooth_t_log = np.log10(smooth_t['CN2T'])
    rtlog_diff = smooth_r_log-smooth_radio_log
    rtlog_diff_avg = abs(np.average(rtlog_diff))
    rtlog_diff_std = np.std(rtlog_diff)
    ###########################################################################
"""

# This section is necessary to allow the option to either run this program as
# a script or to import functions from this program (trop_calc, cn2_calc)
# for use in other programs
if __name__ == "__main__":
    # Check arguments
    mlat=None
    mlon=None
    radio=None
    model=None
    thermo=None
    if(len(sys.argv)<3):
        print("SYNTAX: ./compare_calc.py radiosonde_profile model_profile [thermosonde_tempdiff_file]")
        sys.exit(0)
    if(len(sys.argv)==4):
        radio=sys.argv[1]
        model=sys.argv[2]
        thermo=sys.argv[3]
    elif(len(sys.argv)==5):
        mlat=sys.argv[-2]
        mlon=sys.argv[-1]
    compare_cn2(radio,model,thermo,mlat,mlon)
