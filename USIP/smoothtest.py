#!/usr/bin/env python

"""
  NAME:
    smoothtest.py
  
  PURPOSE:
    Test the 11-point smoother for the radiosonde/thermosonde data comparisons
  
  NOTES:

  11-point smoother possibility
  https://stackoverflow.com/questions/28536191/how-to-filter-smooth-with-scipy-numpy
  
  MODIFICATIONS:
    Blake Sorenson  - 170623: Written 
  
"""

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
            #trop_pres = i
            return press[i]
            #return trop_pres

#def cn2_calc(press, temp, alt, u, v, trop_pres):
# cn2_calc calculates the refractive index structure parameter (cn2) for
# a given set of sounding data. 
def cn2_calc(press, temp, alt, u, v, chooser):
    RADIO=True
    trop_pres=trop_calc(press,temp,alt)
    # Subtract one from the length of the data in order to get lapse rates
    # and wind shear
    length=len(press)-1
    cn2 = numpy.array([]) 
    for i in range(1, length):
        t    = temp[i]+273.15
        p    = press[i]*100
        if(chooser=="new"):
            mid = alt[i] 
            upper = mid+150
            lower = mid-150
            upperindex = numpy.where(alt==(alt[:].flat[numpy.abs(alt[:]-upper).argmin()]))[0][0]
            lowerindex = numpy.where(alt==(alt[:].flat[numpy.abs(alt[:]-lower).argmin()]))[0][0]
            if(upperindex==i):
                upperindex+=1
            if(lowerindex==i):
                lowerindex-=1
            t1   = temp[upperindex]+273.15
            t0   = temp[lowerindex]+273.15
            alt1 = alt[upperindex]
            alt0 = alt[lowerindex]
            u1   = u[upperindex]
            u0   = u[lowerindex]
            v1   = v[upperindex]
            v0   = v[lowerindex]
        else:
            t1   = temp[i+1]+273.15
            t0   = temp[i-1]+273.15
            alt1 = alt[i+1]
            alt0 = alt[i-1]
            u1   = u[i+1]
            u0   = u[i-1]
            v1   = v[i+1]
            v0   = v[i-1]
        #alt1 = alt[:].flat[numpy.abs((alt[:]-alt[i])<150).argmax()]
        #alt0 = alt[:].flat[numpy.abs((alt[i]-alt[:])<150).argmax()]

 
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
#        if(chooser=="thermo"):
#            z = (alt1-alt0)**(-0.666667)
            # ct2 = ((tempdiff)**2)*(1**0.6666667)
#            ct2 = ((t1-t0)**2)*z

            # Calculate refractive index structure parameter
#            ct2 = (tempdiff**2)*r**-(2/3)
#            x =(((76*10**-8)*(p/(t1**2)))**2)*ct2     # using alternative method
#            cn2   = numpy.append(cn2[:], x)

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
#        else:
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
            #if(i<trop_pres):
            if(p>trop_pres*100):
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

# Import modules
import numpy
import sys
import matplotlib.pyplot as plt
import os
import warnings
import datetime as dt
from adpaa import ADPAA
from metpy.plots import SkewT 
from metpy.units import units
from scipy.interpolate import interp1d
import readsounding as rs

###############################################################################
#
# NOTE: The final version will take in three files: radiosonde, thermosonde,
#       and analysis (model) profiles. For now, radiosonde and model profiles
#       will be used to test the code
#
###############################################################################

# Check arguments
if(len(sys.argv)!=3):
    print "SYNTAX: ./cn2_calc.py <radiosonde profile> <model profile>"
#    print "SYNTAX: ./cn2_calc.py <sounding> <x axis variable (alt, press,  or time)> <plot title> <verbose mode optional (-v)>"
    sys.exit(0)

RADIO=True
THERMO=False
# Read in radiosonde, model, and thermosonde data
radio = rs.readsounding(sys.argv[1])
#thermo=ADPAA(sys.argv[3])

# Declare Arrays
out_cn2r  = numpy.array([])
out_temp  = numpy.array([])
out_press = numpy.array([])
in_press  = numpy.array([])
percent_diff = numpy.array([])
diff = numpy.array([])

# THERMO: Use time from radiosonde file to match radiosonde data to 
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


# 11-point T, u, v smoother
avg_t = numpy.array([])
avg_u = numpy.array([])
avg_v = numpy.array([])
indices = numpy.array([])
mid_alt = numpy.array([])
mid_press = numpy.array([])
radio['U']*=0.514444
radio['V']*=0.514444
total_t=radio['TEMP'][0]
total_u=radio['U'][0]
total_v=radio['V'][0]
for i in range(1, len(radio['TEMP'])):
    if(i%11==0):
#        print str(i)+" Average T: "+str(total_t/11) + "\tAverage U: "+\
#                str(total_u/11)+"\tAverage V: "+str(total_v/11)+\
#        "\tAround: "+str(i-6)
        avg_t = numpy.append(avg_t[:],total_t/11)
        avg_u = numpy.append(avg_u[:],total_u/11)
        avg_v = numpy.append(avg_v[:],total_v/11)
        indices = numpy.append(indices[:], i-6)
        mid_alt = numpy.append(mid_alt[:],radio['ALT'][i-6])
        mid_press = numpy.append(mid_press[:],radio['PRESS'][i-6])
        total_t=0
        total_u=0
        total_v=0
#    print str(i)+"\tADDING T: " + str(radio['TEMP'][i]) + " to total: "+str(total_t)
#    print str(i)+"\tADDING U: " + str(radio['U'][i]) + " to total: "+str(total_u)
#    print str(i)+"\tADDING V: " + str(radio['V'][i]) + " to total: "+str(total_v)
    total_t+=radio['TEMP'][i]
    total_u+=radio['U'][i]
    total_v+=radio['V'][i]
# Average whatever is left
#avg_t.append(total_t/(len(radio['TEMP'])%11))
#avg_u.append(total_u/(len(radio['TEMP'])%11))
#avg_v.append(total_v/(len(radio['TEMP'])%11))

# Calculate raw cn2
old_raw_cn2 = cn2_calc(radio['PRESS'],radio['TEMP'],radio['ALT'],radio['U'],radio['V'],"old")
new_raw_cn2 = cn2_calc(radio['PRESS'],radio['TEMP'],radio['ALT'],radio['U'],radio['V'],"new")
old_smooth_cn2 = cn2_calc(mid_press,avg_t,mid_alt,avg_u,avg_v,"old")
new_smooth_cn2 = cn2_calc(mid_press,avg_t,mid_alt,avg_u,avg_v,"new")
ralt = numpy.delete(radio['ALT'],-1)
ralt = numpy.delete(ralt,0)
salt = numpy.delete(mid_alt,-1)
salt = numpy.delete(salt,0)
old_raw_cn2 = numpy.delete(old_raw_cn2, numpy.where(ralt<3000))
new_raw_cn2 = numpy.delete(new_raw_cn2, numpy.where(ralt<3000))
old_smooth_cn2 = numpy.delete(old_smooth_cn2, numpy.where(salt<3000))
new_smooth_cn2 = numpy.delete(new_smooth_cn2, numpy.where(salt<3000))
ralt = numpy.delete(ralt, numpy.where(ralt<3000))
#smooth_cn2 = numpy.delete(smooth_cn2, numpy.where(salt<3000))
salt = numpy.delete(salt, numpy.where(salt<3000))

###############################################################################
#
# Old=level immediately above and below current diagnosis level
# New=levels +/- 150 meters above and below current diagnosis level
#
# With the new lapse rate calculation, the raw and smoothed data are nearly
# identical. Old method is close but not as exact as new method. New method
# smooths data almost exactly like smoother
#
###############################################################################
fig1 = plt.figure()
ax = fig1.add_subplot(111)
plt.xscale('log')
plt.xlim(10**-21,10**-13)
#plt.ylim(1,30000)
plt.plot(old_smooth_cn2,salt,color='black',label='old smooth')
plt.plot(new_smooth_cn2,salt,color='red',label='new smooth')
plt.legend(loc='upper right')
#plt.savefig("140923_smooth1.png",dpi=300)
plt.title("Smooth")
fig2 = plt.figure()
ax = fig1.add_subplot(111)
plt.xscale('log')
plt.xlim(10**-21,10**-13)
#plt.ylim(1,35000)
#plt.ylim(1,30000)
plt.plot(old_raw_cn2,ralt,color='black',label='old raw')
plt.plot(new_raw_cn2,ralt,color='red',label='new raw')
plt.title("Raw")
plt.legend(loc='upper right')
#plt.savefig("140923_smooth1.png",dpi=300)
fig3 = plt.figure()
ax = fig1.add_subplot(111)
plt.xscale('log')
plt.xlim(10**-21,10**-13)
#plt.ylim(1,35000)
#plt.ylim(1,30000)
plt.plot(old_raw_cn2,ralt,color='black',label='old raw')
plt.plot(old_smooth_cn2,salt,color='red',label='old smooth')
plt.title("Old")
plt.legend(loc='upper right')
#plt.savefig("140923_smooth1.png",dpi=300)
###############################################################################
#
# Create interpolated cubic function for smoothed cn2. Could be used to make
# smoother smoothed plots on raw data
#
###############################################################################

testralt = numpy.delete(ralt[:],numpy.where(ralt[:]>salt[-1]))
testralt = numpy.delete(testralt,numpy.where(ralt[:]<salt[0]))
short_raw_cn2 = numpy.delete(new_raw_cn2[:],numpy.where(ralt[:]>salt[-1]))
short_raw_cn2 = numpy.delete(short_raw_cn2[:],numpy.where(ralt[:]<salt[0]))
f1 = interp1d(salt,new_smooth_cn2,kind='cubic')
interp_cn2 = f1(testralt)

cn2_diff = numpy.array([])
#cn2_diff = ((interp_cn2-short_raw_cn2)/interp_cn2)*100.
for interp, raw in zip(interp_cn2, short_raw_cn2):  
    prcnt_diff = ((interp-raw)/max(interp,raw))*100.
    cn2_diff = numpy.append(cn2_diff[:],prcnt_diff)

average = numpy.average(cn2_diff)
std = numpy.std(cn2_diff)
zeros = numpy.zeros(len(testralt))
mu = numpy.empty(len(testralt))
nstd = numpy.empty(len(testralt))
pstd = numpy.empty(len(testralt))
mu.fill(average)
nstd.fill(average-std)
pstd.fill(average+std)

fig4 = plt.figure()
ax = fig4.add_subplot(111)
plt.xscale('log')
plt.xlim(10**-21,10**-13)
#plt.ylim(1,35000)
#plt.ylim(1,30000)
plt.plot(new_raw_cn2,ralt,color='black',label='new raw')
plt.plot(new_smooth_cn2,salt,color='red',label='new smooth')
plt.plot(interp_cn2,testralt,color='green',label='interpolated smooth')
plt.title("New")
plt.legend(loc='upper right')
#plt.savefig("140923_smooth1.png",dpi=300)
fig5 = plt.figure()
ax = fig1.add_subplot(111)
plt.xscale('log')
plt.xlim(10**-21,10**-13)
#plt.ylim(1,35000)
#plt.ylim(1,30000)
plt.plot(new_raw_cn2,ralt,color='black',label='new raw')
plt.plot(old_smooth_cn2,salt,color='red',label='old smooth')
plt.title("New")
plt.legend(loc='upper right')
#plt.savefig("140923_smooth1.png",dpi=300)
plt.show()
"""
# Percent difference plot
fig6 = plt.figure()
ax = fig6.add_subplot(111)
plt.xlim(-100,100)
plt.scatter(cn2_diff,testralt,color='black')
plt.plot(zeros,testralt,color='red')
plt.plot(mu,testralt,color='green')
plt.plot(nstd,testralt,color='blue')
plt.plot(pstd,testralt,color='blue')
print "AVG: "+str(average)+"\tSTD: "+str(std)
plt.title(sys.argv[1])
#plt.savefig("140923_smooth1.png",dpi=300)
plt.show()

"""
"""
plt.rcParams['figure.figsize'] = (9, 9)
fig1 = plt.figure()

# Create skew-T plot
skew = SkewT(fig1, rotation=45)
skew.plot(radio['PRESS']*units.mbar, radio['TEMP']*units.degC, 'red')
skew.plot(radio['PRESS']*units.mbar, radio['DP']*units.degC,'green')
skew.plot(mid_press, avg_t, 'black')
# Generate wind barbs at every 25 millibars
plotP = []
my_interval=numpy.arange(100, 1000, 25) * units('mbar')
for center in my_interval:
    index=(numpy.abs((radio['PRESS']*units.mbar)-center)).argmin()
    if index not in plotP:
        plotP.append(index) 
# Plot the wind barbs using the same color as the temp/dp profiles
skew.plot_barbs(radio['PRESS'][plotP], radio['U'][plotP], radio['V'][plotP], color='red')

# Generate wind barbs at every 25 millibars
plotP2 = []
my_interval=numpy.arange(100, 1000, 25) * units('mbar')
for center in my_interval:
    index=(numpy.abs((mid_press*units.mbar)-center)).argmin()
    if index not in plotP2:
        plotP2.append(index) 
# Plot the wind barbs using the same color as the temp/dp profiles
#skew.plot_barbs(mid_press[plotP2], avg_u[plotP2], avg_v[plotP2], color='black')
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
skew.ax.set_ylim(1050, 50)
skew.ax.set_xlim(-50,50)

#plt.title(sounding1['TITLE']+" vs "+sounding2['TITLE'])
plt.legend(loc='upper left', fontsize=14)
plt.show()
"""

sys.exit(1)


# Use cn2_calc function to calculate cn2 for each profile
out_cn2r = cn2_calc(radio['PRESS'],radio['TEMP'],radio['ALT'],rs.kts2ms(radio['U']),\
                    rs.kts2ms(radio['V']),"radio")

# Remove the first and last elements of each of these arrays to account for 
# the cn2 calculation
out_alt     = numpy.delete(radio['ALT'],-1)
out_alt     = numpy.delete(out_alt,0)
out_press   = numpy.delete(radio['PRESS'],-1)
out_press   = numpy.delete(out_press,0)

# Since cn2 is only being evaluated above the PBL, remove data below 3000 m
out_press = numpy.delete(out_press, numpy.argwhere(out_alt<3000))
out_cn2r  = numpy.delete(out_cn2r, numpy.argwhere(out_alt<3000))
out_alt   = numpy.delete(out_alt, numpy.argwhere(out_alt<3000))

###############################################################################
#
# Statistical calculations taken from Frehlich et al, 2010
#
###############################################################################


# Generate plot
fig1 = plt.figure()
ax = fig1.add_subplot(111)


# Plot data
ax.plot(out_cn2r, out_alt, color='black')
ax.plot(cn2_smooth, out_alt, color='red')
# To view cn2 data, axis must be logarithmic
plt.xscale('log')
# Old limits: 10**-24, 10**-10
plt.xlim(10**-21, 10**-13)
plt.ylim(0,30000)
plt.ylabel(axis)


#imagename=outputstart+"_cn2_"+str(plotter)+".png" 
plt.xlabel("log$_{10}$ [$C_{n}^{2}$ [m$^{-2/3}$]]")
plt.legend(loc=leg_loc, fontsize=12)
# plt.show()

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
