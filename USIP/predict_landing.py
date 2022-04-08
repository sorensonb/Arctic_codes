#!/usr/bin/env python

"""
  NAME:
    predict_landing.py

  PURPOSE:
    Calculate the predicted landing point of a radiosonde package. 

  CALLS:
    - Python modules numpy, sys, re, matplotlib.pyplot, matplotlib.gridspec,
      mpl_toolkits.basemap, scipy.interpolate
      (Requires scipy version 0.17.0 or greater)
    - Custom module readsounding

  SYNTAX:
    /nas/home/bsorenson/prg/thermosonde1617/predict_landing.py <GRAW RTS file> <whole or end>
    
    This program needs a Graw RTS table (/nas/Radiosonde/20170503/170503_204757_GFK_GRAW.txt) 
    
  MODIFICATIONS:
    Blake Sorenson  - 170616: Written

"""

import numpy
import sys
import re
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import commands
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
from readsounding import readsounding

if(len(sys.argv)!=3):
    print "SYNTAX: ./predict_landing.py GRAW_RTS_file whole/final > <outputfilename>.txt"
    print "         Whole: Predicts the entire balloon path starting with the 100th point"
    print "                and using the u and v winds to calculate the path"
    print "         Final: Predicts the final few minutes of the balloon path starting "
    print "                at the last data point received from the radiosonde."
    sys.exit(1)

status, output = commands.getstatusoutput("wc -l "+sys.argv[1])
length = int(output.split(' ')[0])

TIME = numpy.zeros([length-1])
SPD = numpy.zeros([length-1])
DIR = numpy.zeros([length-1])
LON = numpy.zeros([length-1])
LAT = numpy.zeros([length-1])
ALT = numpy.zeros([length-1])

fx = numpy.array([])
fy = numpy.array([])
flat = numpy.array([])
flon = numpy.array([])
#newx = numpy.array([])
#newy = numpy.array([])

# Check the file name to determine how to read in the data from the soundings
filechooser = sys.argv[1].split("/")[-1]
if((filechooser=="170214_000000_SIM_GRAW.txt") | (filechooser=="140923_000000_GFK_GRAW.txt")):
    spdindex = 4
    spdmult  = 0.514444
    dirindex = 5
    lonindex = 6
    latindex = 7
    altindex = 8
else:
    spdindex = 16
    spdmult  = 1.0
    dirindex = 17
    lonindex = 21
    latindex = 22
    altindex = 23

count=0
with open(sys.argv[1]) as f:
    next(f)
    for i, line in enumerate(f):    
        templine=line.strip().replace("\t"," ")
        finline =re.sub(" +"," ",templine).split(" ")
        badData=False
        for entry in finline:
            if(entry=='-----'):
                badData=True
        # If the current line has any missing data, insert missing values
        # into each array in the current spot
        if(badData==True):
            TIME[count] = 999999.9999
            SPD[count]  = 999999.9999
            DIR[count]  = 999999.9999
            LON[count]  = 999999.9999
            LAT[count]  = 999999.9999
            ALT[count]  = 999999.9999
        else:
            TIME[count] = float(finline[0])
            SPD[count]  = float(finline[spdindex])*spdmult
            DIR[count]  = float(finline[dirindex])
            LON[count]  = float(finline[lonindex])
            LAT[count]  = float(finline[latindex])
            ALT[count]  = float(finline[altindex])
        count+=1

TIME = numpy.delete(TIME[:], numpy.where(TIME==999999.9999))
SPD = numpy.delete(SPD[:], numpy.where(SPD==999999.9999))
DIR = numpy.delete(DIR[:], numpy.where(DIR==999999.9999))
LON = numpy.delete(LON[:], numpy.where(LON==999999.9999))
LAT = numpy.delete(LAT[:], numpy.where(LAT==999999.9999))
ALT = numpy.delete(ALT[:], numpy.where(ALT==999999.9999))

# Extrapolation

# To predict landing location, extrapolate altitude to surface. Assume a
# constant rate of descent 
tx = TIME
ty = ALT
tlat = LAT
tlon = LON
altf = interpolate.interp1d(tx,ty,fill_value='extrapolate')
i = int(tx[1])

print altf(10)

# Append the extrapolated altitudes and times to arrays
while(altf(i)>ALT[0]):
    fx = numpy.append(fx[:],float(i))
    fy = numpy.append(fy[:], altf(i))
    i+=1 
#latf = interpolate.interp1d(tx,tlat,fill_value='extrapolate')
#lonf = interpolate.interp1d(tx,tlon,fill_value='extrapolate')
#i = int(tx[0])

# Generate basemap
fig5 = plt.figure(figsize=(10,8))
m = Basemap(llcrnrlon=-98., llcrnrlat=45., urcrnrlon=-90., urcrnrlat=50.,\
           resolution='f',projection='cass',lon_0=-95.,lat_0=45.)
#m.etopo()
#m.shadedrelief()
m.drawcoastlines()
m.drawcountries()
m.drawstates()
m.drawcounties()
# Convert balloon latitude and longitude to map projection coordinates
# (in meters)
plotx, ploty = m(LON,LAT)

# Predict balloon path from point 100
# Use readsounding as a shortcut to get the U and V wind components 
a = readsounding(sys.argv[1])
a['UWIND']*=0.514444
a['VWIND']*=0.514444

count=0
if(sys.argv[2].lower()=="whole"):
    # Predict the entire balloon path starting at the 100th point
    #closealt = a['ALT'][:].flat[numpy.abs(a['ALT'][:]-fy[i]).argmin()]
    # Proof of concept: predict the entire flight path
    startx = plotx[100]
    starty = ploty[100]
    newx = numpy.zeros(len(fy)-100)
    newy = numpy.zeros(len(fx)-100)
    for i in range(100, len(a['UWIND'])):
        startx+=a['UWIND'][i]
        starty+=a['VWIND'][i]
        newx[count] = startx 
        newy[count] = starty 
    for i in range(len(a['UWIND']),len(fy)):
        closealt = a['ALT'][:].flat[numpy.abs(a['ALT'][:]-fy[i]).argmin()]
        caindex = numpy.where(a['ALT']==closealt)     
        startx+=a['UWIND'][caindex]
        starty+=a['VWIND'][caindex]
        newx[count] = startx 
        newy[count] = starty 
    count+=1
else:
    # Predict the unknown portion of the balloon path starting at the last
    # known point. To do this, for each extrapolated altitude, a comparible
    # altitude will be found from the ascent. The u and v wind components
    # from that altitude on the ascent will be used to find the next 
    # balloon position. Assumes a similar wind pattern on the descent 
    startx = plotx[-1]
    starty = ploty[-1]
   
    newx = numpy.zeros(len(fy)-len(LAT))
    newy = numpy.zeros(len(fy)-len(LAT))
    
    for i in range(len(LAT),len(fy)):
        # Find the altitude from the ascent profile that is closest to 
        # the current extrapolated descent altitude
        closealt = a['ALT'][:].flat[numpy.abs(a['ALT'][:]-fy[i]).argmin()]
        # Find the index of that ascent altitude 
        caindex = numpy.where(a['ALT']==closealt)     
        # Use the index of the ascent altitude to find the u and v
        # components at that altitude.
    
        ###########################################################################
        # 
        # Since it is assumed that the winds on the descent are similar
        # to the winds on the ascent and the movement of the balloon is only due
        # to wind, the next x and y coordinates can be found by adding the current 
        # u and v wind components to the current x and y positions, respectively.  
        #
        # To find the next balloon position, these equations are used:
        #   new_x_position = old_x_position + (current_u_speed*1s)
        #   new_y_position = old_y_position + (current_v_speed*1s)
        #
        # Because the wind speeds are in meters per second and the projection
        # coordinates are in meters, the new x and y coordinates can be found
        # by first multiplying the u and v winds by 1 second to find the change
        # in the x and y directions caused by the u and v winds (the 1 Hz data from
        # the radiosonde allows this). Adding the change in x and y coordinates
        # to the current x and y coordinates gives the new x and y coordinates.
        #
        ###########################################################################
        startx+=a['UWIND'][caindex]
        starty+=a['VWIND'][caindex]
        newx[count] = startx 
        newy[count] = starty 
        count+=1
    
# Convert the predicted x and y projection coordinates to lat/lon
printlon, printlat = m(newx, newy, inverse=True)
print "Latitude,Longitude"
for x, y, in zip(printlat, printlon):
    print str(x)+","+str(y)

#m.drawrivers()
#plotfx, plotfy = m(flon,flat)
#m.plot(plotfx, plotfy, color='red')
m.plot(plotx, ploty)
m.plot(newx,newy)
#fig5.show()
#plt.show()
"""
x = np.arange(0,10)
y = np.exp(-x/3.0)
f = interpolate.interp1d(x, y, fill_value='extrapolate')

print f(9)
print f(11)
"""
