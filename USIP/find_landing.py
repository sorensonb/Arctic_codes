#!/usr/bin/env python

"""
  NAME:
    compare_cn2.py


"""

import numpy
import sys
import re
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
from readsounding import readsounding

if(len(sys.argv)!=2):
    print "SYNTAX: ./find_landing.py <sounding file>"
    sys.exit(1)

TIME = numpy.array([])
SPD = numpy.array([])
DIR = numpy.array([])
LON = numpy.array([])
LAT = numpy.array([])
ALT = numpy.array([])

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


with open(sys.argv[1]) as f:
    next(f)
    for i, line in enumerate(f):    
        templine=line.strip().replace("\t"," ")
        finline =re.sub(" +"," ",templine).split(" ")
        skip=False
        for entry in finline:
            if(entry=='-----'):
                skip=True
        if(skip==False):
            TIME = numpy.append(TIME[:], float(finline[0]))
            SPD  = numpy.append(SPD[:], float(finline[spdindex])*spdmult)
            DIR  = numpy.append(DIR[:], float(finline[dirindex]))
            LON  = numpy.append(LON[:], float(finline[lonindex]))
            LAT  = numpy.append(LAT[:], float(finline[latindex]))
            ALT  = numpy.append(ALT[:], float(finline[altindex]))

#print "Latitude,Longitude"
#for lat, lon in zip(LAT, LON):
#    print str(lat)+","+str(lon)
fig1 = plt.figure()
plt.plot(TIME, ALT)

#plt.show()


# Produces map of northern Minnesota

#m = Basemap(llcrnrlon=-98., llcrnrlat=45., urcrnrlon=-90., urcrnrlat=50.,\
#            epsg=5520)
#m.arcgisimage(service='ESRI_Imagery_World_2D',xpixels=1500,verbose=True)

#plt.show()

# Extrapolation

fx = numpy.array([])
fy = numpy.array([])
flat = numpy.array([])
flon = numpy.array([])
flatx = numpy.array([])
tx = TIME[:-10]
ty = ALT[:-10]
tlat = LAT[:-10]
tlon = LON[:-10]
altf = interpolate.interp1d(tx,ty,fill_value='extrapolate')
i = int(tx[1])
while(altf(i)>0):
    fx = numpy.append(fx[:],float(i))
    fy = numpy.append(fy[:], altf(i))
    i+=1 
latf = interpolate.interp1d(tx,tlat,fill_value='extrapolate')
lonf = interpolate.interp1d(tx,tlon,fill_value='extrapolate')
i = int(tx[0])
print "Latitude,Longitude"
for i in range(int(tx[0]),int(fx[-1])+1):
    flatx = numpy.append(flatx[:], i)
    flat = numpy.append(flat[:], latf(i))
    flon = numpy.append(flon[:], lonf(i))
    print str(latf(i))+","+str(lonf(i))


fig2 = plt.figure()
plt.plot(fx, fy, color='red')
plt.plot(tx, ty, color='black')

fig3 = plt.figure()
plt.plot(flatx, flat, color='red')
plt.plot(tx, tlat, color='black')

fig4 = plt.figure()
plt.plot(flatx, flon, color='red')
plt.plot(tx, tlon, color='black')

fig5 = plt.figure(figsize=(10,8))
m = Basemap(llcrnrlon=-98., llcrnrlat=45., urcrnrlon=-90., urcrnrlat=50.,\
           resolution='f',projection='cass',lon_0=-95.,lat_0=45.)
#m.etopo()
#m.shadedrelief()
m.drawcoastlines()
m.drawcountries()
m.drawstates()
m.drawcounties()
#m.drawrivers()
plotfx, plotfy = m(flon,flat)
m.plot(plotfx, plotfy, color='red')
plotx, ploty = m(LON,LAT)
m.plot(plotx, ploty)
fig5.show()
plt.show()
"""
x = np.arange(0,10)
y = np.exp(-x/3.0)
f = interpolate.interp1d(x, y, fill_value='extrapolate')

print f(9)
print f(11)
"""
