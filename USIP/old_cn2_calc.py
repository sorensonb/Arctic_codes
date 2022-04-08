#!/usr/bin/env python

"""
  NAME:
    cn2_calc.py
  
  PURPOSE:
    Calculate the refractive index structure parameters for radiosonde, 
    model, and thermosonde soundings and compare the profiles. 
  
  CALLS:
    Built-in Python modules numpy, sys, re, matplotlib.pyplot, os
    Custom modules ADPAA and readsounding
 
  USAGE:
    There are two ways to use cn2_calc.py. The first is to run it as a script. 
    To run cn2_calc.py as a script, use this syntax:
 
      /nas/home/bsorenson/prg/thermosonde_1617/cn2_calc.py <radiosonde file> <model file> 

    
  
  EXAMPLE:
    /nas/home/bsorenson/prg/thermosonde_1617/cn2_calc.py 170605_120000_BIS_UWYO.txt RAP.Op40_2017060512_ANALYSIS_BIS_GSD
  
    This will create a UND NASA-formatted file with time, pressure, altitude,
    temperature, and the refractive index structure parameter. It will also
    create a plot of the optical turbulence with altitude on the x axis and
    the title "GRAW Optical Turbulence: September 23, 2014".  

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

  11-point smoother possibility
  https://stackoverflow.com/questions/28536191/how-to-filter-smooth-with-scipy-numpy


  
  MODIFICATIONS:
    Blake Sorenson  - 170605: Written 
  
"""
# Import modules
import numpy
import sys
import matplotlib.pyplot as plt
import os
import warnings
import datetime as dt
from adpaa import ADPAA
from readsounding import readsounding

# trop_calc finds the tropopause for a sounding and returns the index in the 
# sounding data of the tropopause. 
def trop_calc(press, temp, alt):
    length=len(press)

    # Find tropopause pressures
    var1=5
    maxlength=length-var1
    # Skip the lower third of the atmosphere to avoid low level inversions
    starter=int(maxlength/3)
    # By default, the tropopause is first set as the highest sounding level
    trop_pres=length
    for i in range(starter, maxlength):
        t1=temp[i+var1] 
        t0=temp[i-var1]
        alt1=alt[i+var1]
        alt0=alt[i-var1]
        tdtdz=(t1-t0)/(alt1-alt0)
        # Definition of tropopause is lowest level at which the lapse rate is
        # -0.002 degC/km or higher
        if((tdtdz>-0.002)):
            trop_pres = i
            return trop_pres

#def cn2_calc(press, temp, alt, u, v, trop_pres):
# cn2_calc calculates the refractive index structure parameter (cn2) for
# a given set of sounding data. 
def cn2_calc(press, temp, alt, u, v, trop_pres, chooser):
    RADIO=True
    # Subtract one from the length of the data in order to get lapse rates
    # and wind shear
    length=len(press)-1
    cn2 = numpy.array([]) 
    for i in range(1, length):
        p    = press[i]*100
        t1   = temp[i+1]+273.15
        t    = temp[i]+273.15
        t0   = temp[i-1]+273.15
        alt1 = alt[i+1]
        alt0 = alt[i-1]
        u1   = u[i+1]
        u0   = u[i-1]
        v1   = v[i+1]
        v0   = v[i-1]

 
        ########################################################################
        # IMPORTANT:
        # 
        # This setup is only for temporary, testing purposes.
        # The real code must use the temperature difference measured by the
        # thermosonde in place of t1-t0 and the distance between the probes
        # in place of alt1-alt0, which will end up being 1 meter.

        # The actual equation to use is:
        # cn2(h) = (((76*10**-8)*(p/(t**2)))**2)*ct2(h), where
        # ct2(h) = (tempdiff**2)*r**-0.6666667
        #
        ########################################################################
 
#        if(THERMO==True):
        if(chooser=="thermo"):
#            z = (alt1-alt0)**(-0.666667)
            # ct2 = ((tempdiff)**2)*(1**0.6666667)
#            ct2 = ((t1-t0)**2)*z

            # Calculate refractive index structure parameter
            ct2 = (tempdiff**2)*r**-(2/3)
            x =(((76*10**-8)*(p/(t1**2)))**2)*ct2     # using alternative method
            cn2   = numpy.append(cn2[:], x)

        ########################################################################
        #									      
        # IMPORTANT: 							      
        #   The above equations are for use with the thermosonde measurements,    
        #   where there are temperature difference measurements.                  
        #
        #   For other soundings, use these equations
        #
        ########################################################################
#        elif(RADIO==True):
        else:
            # Calculate lapse rate
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
                if(i<trop_pres):
                    l0=(0.1**(4/3))*10**(1.64+42*S)
                    if(dtdz>0):
                        dtdz=0
                else:
                    l0=(0.1**(4/3))*10**(0.506+50*S)
            # Calculate refractive index structure parameter 
            ###################################################################
            #
            # NOTE: Dry adiabatic lapse rate positive or negative?
            #       According to Dewan et al, set to positive
            #
            ###################################################################
            cn2_r = 2.8*(((79*10**-8)*(p/(t**2)))**2)*((dtdz+(9.8*(10**-3)))**2)*l0
            cn2 = numpy.append(cn2[:], cn2_r)
       # if((alt1-alt0)<400):
       #     print "BAD"
       # print "P: " + str(p) + "  T1: " + str(t1)+"  T: "+str(t)+"  T0: "+str(t0)+\
       #       "  ALT1: "+str(alt1)+"  ALT0: "+str(alt0)+"  CN2: "+str(cn2_r)
       # print "P: " + str(p) + "  U1: " + str(u1)+"  U0: "+str(u0)+"  V1: "+str(v1)+\
       #       "  V0: "+str(v0)+"  ALT1: "+str(alt1)+"  CN2: "+str(cn2_r)
    return cn2

    
###############################################################################
#
# NOTE: The final version will take in three files: radiosonde, thermosonde,
#       and analysis (model) profiles. For now, radiosonde and model profiles
#       will be used to test the code
#
###############################################################################

# Check arguments
if(len(sys.argv)!=3):
    print "SYNTAX: ./compare_calc.py <radiosonde profile> <model profile>"
#    print "SYNTAX: ./cn2_calc.py <sounding> <x axis variable (alt, press,  or time)> <plot title> <verbose mode optional (-v)>"
    sys.exit(0)

RADIO=True
THERMO=False
# Read in radiosonde, model, and thermosonde data
radio = readsounding(sys.argv[1])
model = readsounding(sys.argv[2])
#thermo=ADPAA(sys.argv[3])

# Declare Arrays
out_cn2r  = numpy.array([])
in_cn2r   = numpy.array([])
out_temp  = numpy.array([])
in_temp   = numpy.array([])
out_press = numpy.array([])
in_press  = numpy.array([])
in_u      = numpy.array([])
in_v      = numpy.array([])
in_alt    = numpy.array([])
percent_diff = numpy.array([])
diff = numpy.array([])

# THERMO: Use time from radiosonde file to match thermosonde data to 
# thermosonde data
# Extract date and time information from the radiosonde's readsounding name
rname   = radio['NAME'].split("_")
ryear   = rname[0]
rmonth  = rname[1]
rday    = rname[2]
rhour   = rname[3]
rminute = rname[4]
rsecond = rname[5].split(".")[0]
# Convert hh:mm:ss time to seconds from midnight
rtime   = int(dt.timedelta(hours=int(rhour), minutes=int(rminute), seconds=int(rsecond)).total_seconds()) 

# Extract date and time information from the model's readsounding name
mname   = model['NAME'].split("_")
myear   = mname[0]
mmonth  = mname[1]
mday    = mname[2]
mhour   = mname[3]
mminute = mname[4]
msecond = mname[5].split(".")[0]
# Convert hh:mm:ss time to seconds from midnight
mtime   = int(dt.timedelta(hours=int(mhour), minutes=int(mminute), seconds=int(msecond)).total_seconds()) 

# Get thermosonde start time from file name
#for rtime in radio['TIME']:
 

###############################################################################
#
# NOTE: In most cases, the radiosonde profiles will have more data than the 
#       model profiles, so the model soundings will be in the outer loop most 
#       of the time. However, this will not be hard coded to account for 
#       the slim possibility that this does not occur
#
###############################################################################

# Determine which dataset is larger. The larger dataset will be used
# in the inside loop for determining similar pressure levels
length1=len(radio['PRESS'])
length2=len(model['PRESS'])
if(length1>length2):
    outermax=length2
    outrpres=model['PRESS']  
    out_temp=model['TEMP']  
    out_alt=model['ALT']
    # U and V wind components must be converted from knots to m/s
    out_u=model['U']*0.514444
    out_v=model['V']*0.514444
    innermax=length1
    innrpres=radio['PRESS']
    innrtemp=radio['TEMP']
    innralt=radio['ALT']
    # U and V wind components must be converted from knots to m/s
    innru=radio['U']*0.514444
    innrv=radio['V']*0.514444
else:
    outermax=length1
    outrpres=radio['PRESS']
    out_temp=radio['TEMP']
    out_alt =radio['ALT']
    # U and V wind components must be converted to knots to m/s
    out_u=radio['U']*0.514444
    out_v=radio['V']*0.514444
    innermax=length2
    innrpres=model['PRESS']
    innrtemp=model['TEMP']
    innralt=model['ALT']
    # U and V wind components must be converted to knots to m/s
    innru=model['U']*0.514444
    innrv=model['V']*0.514444

# Find similar pressures between the two profiles
for i in range (0, outermax):
    temp_outrpres = 0
    # Use the first pressure in the inner profile as the closest pressure
    # to start with
    temp_innrpres = innrpres[0] 
    temp_innralt  = innralt[0]
    # Loop through inner pressures and compare each pressure to the current
    # outer pressure
    for j in range(0, innermax):
        # If the current inner pressure is closer to the current outer pressure
        # declare the current inner pressure to be the closest pressure
        if(abs(outrpres[i] - innrpres[j]) < abs(outrpres[i] - temp_innrpres)):
            temp_outrpres = outrpres[i]
            temp_innrpres = innrpres[j]
            temp_innralt  = innralt[j]
        else:
            continue

    # Append the close pressures and altitudes to arrays
    out_press = numpy.append(out_press[:], temp_outrpres)
    in_press = numpy.append(in_press[:], temp_innrpres)
    in_alt = numpy.append(in_alt[:], temp_innralt)

outloop = len(out_press)
inloop = len(in_press)-1
#for i in range(0, outloop):
#    out_temp = numpy.append(out_temp[:], outrtemp[int(numpy.where(outrpres==out_press[i])[0][0])])
#    out_temp = numpy.append(out_temp[:], outrtemp[i])

# Average the temperature, u wind, and v wind around each pressure level for 
# the radiosonde data
for i in range(1, inloop):
    top=int(numpy.where(innrpres==in_press[i+1])[0][0])
    mid=int(numpy.where(innrpres==in_press[i])[0][0])
    bot=int(numpy.where(innrpres==in_press[i-1])[0][0])
    topmid=int((top+mid)/2)
    midbot=int((mid+bot)/2)
    if(midbot==topmid):
        topmid+=1

#    addtemp = numpy.average(innrtemp[midbot:topmid]) 
#    addu = numpy.average(innru[midbot:topmid]) 
#    addv = numpy.average(innrv[midbot:topmid]) 
    addtemp = numpy.average(innrtemp[bot:top]) 
    addu = numpy.average(innru[bot:top]) 
    addv = numpy.average(innrv[bot:top]) 
    in_temp = numpy.append(in_temp[:], addtemp)
    in_u = numpy.append(in_u[:], addu)
    in_v = numpy.append(in_v[:], addv)

#print "OLD"
#print "NEW"

# Remove the first and last elements of each of these arrays to account for
# the averaging done above
out_press = numpy.delete(out_press, 0)
out_press = numpy.delete(out_press, -1)
out_temp  = numpy.delete(out_temp, 0)
out_temp  = numpy.delete(out_temp, -1)
out_u     = numpy.delete(out_u, 0)
out_u     = numpy.delete(out_u, -1)
out_v     = numpy.delete(out_v, 0)
out_v     = numpy.delete(out_v, -1)
out_alt   = numpy.delete(out_alt, 0)
out_alt   = numpy.delete(out_alt, -1)
in_press  = numpy.delete(in_press, 0)
in_press  = numpy.delete(in_press, -1)
in_alt    = numpy.delete(in_alt, 0)
in_alt    = numpy.delete(in_alt, -1)

# Calculate cn2 and append to array
#   Convert pressure to Pascals
#   Convert temperature to Kelvin
out_trop_pres = trop_calc(out_press, out_temp, out_alt)
in_trop_pres = trop_calc(in_press, in_temp, in_alt)

# Use cn2_calc function to calculate cn2 for each profile
out_cn2r = cn2_calc(out_press, out_temp, out_alt, out_u, out_v, out_trop_pres,\
                    "radio")
in_cn2r = cn2_calc(in_press, in_temp, in_alt, in_u, in_v, in_trop_pres, "radio")
#thrm_cn2 = cn2_calc(out_press, out_temp, out_alt, out_u, out_v, out_trop_pres,\
#                    "thermo")

# Remove the first and last elements of each of these arrays to account for 
# the cn2 calculation
out_alt     = numpy.delete(out_alt,-1)
out_alt     = numpy.delete(out_alt,0)
out_press   = numpy.delete(out_press,-1)
out_press   = numpy.delete(out_press,0)
in_alt      = numpy.delete(in_alt,-1)
in_alt      = numpy.delete(in_alt,0)
in_press    = numpy.delete(in_press,-1)
in_press    = numpy.delete(in_press,0)

# Since cn2 is only being evaluated above the PBL, remove data below 3000 m
in_press  = numpy.delete(in_press, numpy.argwhere(out_alt<3000))
in_cn2r   = numpy.delete(in_cn2r, numpy.argwhere(out_alt<3000))
in_alt    = numpy.delete(in_alt, numpy.argwhere(out_alt<3000))
out_press = numpy.delete(out_press, numpy.argwhere(out_alt<3000))
out_cn2r  = numpy.delete(out_cn2r, numpy.argwhere(out_alt<3000))
out_alt   = numpy.delete(out_alt, numpy.argwhere(out_alt<3000))

###############################################################################
#
# Statistical calculations taken from Frehlich et al, 2010
#
###############################################################################

# Compute log10 cn2 arrays
# Use the new arrays to calculate an average difference and standard devation
out_log = numpy.log10(out_cn2r)
in_log = numpy.log10(in_cn2r)
log_diff = out_log-in_log
log_diff_avg = abs(numpy.average(log_diff))
log_diff_std = numpy.std(log_diff)

# Calculate percent difference between the two profiles
for i in range(0, len(out_cn2r)):
    t_dif=abs(out_cn2r[i]-in_cn2r[i])
    prcnt=(t_dif/max(abs(out_cn2r[i]), abs(in_cn2r[i])))*100.
    #prcnt=(t_dif/(abs((out_cn2r[i]+in_cn2r[i])/2)))*100.
    #prcnt=(t_dif/min(abs(out_cn2r[i]), abs(in_cn2r[i])))

    diff = numpy.append(diff[:], t_dif)
    percent_diff = numpy.append(percent_diff[:], prcnt)

# Generate plot
fig1 = plt.figure()
ax = fig1.add_subplot(111)
# Plot pressure or altitude on the x axis
alt=True
if(alt==True):
    out_plotAxis=out_alt
    in_plotAxis=in_alt
    axis="Altitude [m]"
    plt.xlim(1, 35000)
    leg_loc='upper right'
else:
    out_plotAxis=list(reversed(out_press))
    in_plotAxis=list(reversed(in_press))
    out_cn2r=list(reversed(out_cn2r))
    in_cn2r=list(reversed(in_cn2r))
    percent_diff=list(reversed(percent_diff))
    diff=list(reversed(diff))
    plt.xlim(1, 1050)
    leg_loc='upper left'
    axis="Pressure [hPa]"

# Use the radiosonde and model titles to generate plot labels
olabel=model['TITLE'].split(" ")[0]+" "+model['TITLE'].split(" ")[1]+" "+\
       model['TITLE'].split(" ")[-1]
#olabel = model['TITLE'].split(" ")[0]+" Model Analysis"
ilabel=radio['TITLE'].split(" ")[0]+" "+radio['TITLE'].split(" ")[1]+" "+\
       radio['TITLE'].split(" ")[-1]
#ilabel = radio['TITLE'].split(" ")[0]+" Radiosonde Analysis"
# Plot data
ax.plot(out_cn2r, out_plotAxis, color='black', label=olabel)
ax.plot(in_cn2r, in_plotAxis, color='red', label=ilabel)
# To view cn2 data, axis must be logarithmic
plt.xscale('log')
# Old limits: 10**-24, 10**-10
plt.xlim(10**-21, 10**-13)
plt.ylim(0,30000)
plt.ylabel(axis)
#imagename=outputstart+"_cn2_"+str(plotter)+".png" 
plt.xlabel("log$_{10}$ [$C_{n}^{2}$ [m$^{-2/3}$]]")
plt.legend(loc=leg_loc, fontsize=12)
final_prcntdif = numpy.delete(percent_diff, numpy.where(numpy.isnan(percent_diff)))
print radio['TITLE']+" vs. "+model['TITLE']+"\t  Log Diff Avg: "+str(round(log_diff_avg, 5)) + "\t  STD: "+str(round(log_diff_std, 5))+"\tAVR PRCNT DIFF: " + str(round(numpy.average(final_prcntdif), 5))
# plt.show()

PLOT=False
if(PLOT==True):
#    fig2 = plt.figure()
    ax2 = ax.twinx()
    ax2.plot(out_plotAxis, percent_diff, color='blue', label='% diff')
#    ax2.set_yscale('log')
#    plt.ylim(10**-24, 10**-10)
    ax2.set_ylabel("Percent Difference")
    ax2.set_ylim(0, 200)
#plt.title("$C_{n}^{2}$ Profiles for Bismarck, ND at 00Z 12 June 2017")
#plt.title("June 12 2017 00Z Bismarck $C_{n}^{2}$ Radiosonde and Model Profiles")
plt.title(radio['TITLE']+" vs. "+model['TITLE'])
plt.show()
 
#    plt.savefig(imagename, dpi=300)
#    print "------- Saved image: " + imagename


"""
sys.exit(1)

# ORDER: Time, press, alt, temp, cn2
mvc_arr = ['999999.9999','999999.9999','999999.9999','9.99999e-09']
vscal   = [     '1.0000',     '1.0000',     '1.0000','1.0000']
vmiss   = mvc_arr
vname   = ['Air pressure [hPa]','Altitude [m]',
           'Air temp [degC]',
           'Refractive index structure parameter [m^-2/3]']
vdesc   = ['Time','Pressure','Altitude','Air_Temp','Refrct_prm']
vunits  = ['second','hPa','m','degC','m^-2/3']

# Add to dictionary
dicty = {}
dicty['Time']=time
dicty['Pressure']=press
dicty['Altitude']=alt
dicty['Air_Temp']=temp
if(RADIO==True):
    dicty['Refrct_prm']=cn2r
else:
    dicty['Refrct_prm']=cn2

# Make an ADPAA object
adpaa = ADPAA()

adpaa.NLHEAD  = '22'
adpaa.FFI     = '1001'
adpaa.ONAME   = 'Delene, David'
adpaa.ORG     = 'University of North Dakota'
adpaa.SNAME   = 'NASA USIP Thermosonde'
adpaa.MNAME   = 'NASAUSIP2017'
adpaa.IVOL    = 1
adpaa.NVOL    = 1
adpaa.DATE    = '20'+year+' '+month+' '+day
adpaa.RDATE   = re.sub("-"," ",str(datetime.date.today()))+""
adpaa.DX      = 1.0
adpaa.XNAME   = 'Time [second]; Time in UT seconds from midnight on day \
                aircraft flight started'
adpaa.NV      = len(mvc_arr)
adpaa.VSCAL   = vscal
adpaa.VMISS   = vmiss
adpaa.VNAME   = vname
adpaa.NSCOML  = 0
adpaa.NNCOML  = 4
adpaa.DTYPE   = 'Final Data'
adpaa.VFREQ   = '1 Hz Data'
adpaa.VDESC   = vdesc
adpaa.VUNITS  = vunits
adpaa    = dicty
# Create format array
formats = adpaa.create_format_arr(mvc_arr)
adpaa.formats=formats

# Save file
if(verbose==True):
    print "----- Saving to UND-NASA file"
outname=outputstart+".cn2.1Hz"
adpaa.WriteFile(outname)
print "------- Saved file " + outname

if(verbose==True):
    print "Move output files to /nas/home/bsorenson/prg/thermosonde1617/output_files/"
    print "Copy images to ~/Research/thermosonde/"
    print "Move images to /nas/home/bsorenson/prg/thermosonde1617/images"
os.system('mv *.cn2.1Hz /nas/home/bsorenson/prg/thermosonde1617/output_files/')
#os.system('gthumb '+imagename)
os.system('cp *.png ~/Research/thermosonde/')
os.system('mv '+imagename+' /nas/home/bsorenson/prg/thermosonde1617/images/')
    """

