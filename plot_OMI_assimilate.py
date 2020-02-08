#!/usr/bin/env python
"""
  NAME:

  PURPOSE:

  SYNTAX:

  EXAMPLE:

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2018/09/07:

"""

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as color
##import matplotlib.colors as colors
from matplotlib.colors import rgb2hex,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase

if(len(sys.argv)<2):
    print("SYNTAX: python plot_OMI_assimilate.py assimilation_output_file")
    sys.exit()

infile = sys.argv[1]
plot_date = infile.strip().split('/')[-1].split('_')[0]
#plot_date = infile.strip().split('/')[-1].split('_')[-1].split('.')[0]
if(plot_date=='flat'):
    plot_date = infile.strip().split('/')[-1].split('_')[0]


##LAT[i,j],LON[i,j],AI[i,j],0.5,SZA[i,j],VZA[i,j],RAZ[i,j], \
##ALBEDO[i,j,0],ALBEDO[i,j,1],REFLECTANCE[i,j,0],\
##REFLECTANCE[i,j,1]))

# Set up grid arrays
n_p = 1440
nl = 720

UVAI     = np.zeros(shape=(n_p,nl))
#realUVAI = np.zeros(shape=(n_p,nl))
count      = np.zeros(shape=(n_p,nl))
#real_count = np.zeros(shape=(n_p,nl))

# Open the input file
with open(infile,'r') as f:
    # Read in all lines from the file
    flines = f.readlines()
    # Use the number at the beginning of the file to determine how large to
    # make the data arrays
    file_len = len(flines)

    # For now, only plot AI, so only extract LAT, LON, and AI from the
    # file.
    LON     = np.zeros(file_len)
    LAT     = np.zeros(file_len)
    baseAI  = np.zeros(file_len)
    #realAI  = np.zeros(file_len)

    # Loop over the data file and insert the AI values into the 0.25x0.25 grid
    for i, line in enumerate(flines[1:]):
        ###LAT[i] = float(lines.strip()[0])
        ###LON[i] = float(lines.strip()[1])
        ###AI[i]  = float(lines.strip()[2])
        tlon      = float(line.strip().split()[0])
        tlat      = float(line.strip().split()[1])
        tbase_ai  = float(line.strip().split()[6])
        #tbase_ai  = float(line.strip().split()[2])
        #treal_ai  = float(line.strip().split()[2])

        index1 = int(np.floor(tlat*4 + 360.))
        index2 = int(np.floor(tlon*4 + 720.))
        #index1 = int(np.floor(plotLAT[i,j]*4 + 360.))
        #index2 = int(np.floor(plotLON[i,j]*4 + 720.))
        
        if(index1 < 0): index1 = 0
        if(index1 > 719): index1 = 719
        if(index2 < 0): index2 = 0                                                                                            
        if(index2 > 1439): index2 = 1439
        
        UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + tbase_ai)/(count[index2,index1]+1)
        count[index2, index1] = count[index2,index1] + 1


lonmin = -180
lonmax = 180

latmax =  90
latmin =  -90
# Set up the polar stereographic projection map
fig1 = plt.figure(figsize=(8,8))
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')
fig = plt.gcf()
m.drawcoastlines()
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))

#define the lat which will be gridded into 0.25 x 0.25 bins
LATalt = np.arange(-90,90,0.25)
LONalt = np.arange(-180,180,0.25)


# Set up the color bar
cmap = plt.get_cmap('jet')
v_min = 0.000  # AOD
v_max = 1.000
#v_min = -1.000  # AI
#v_max = 3.000

norm = Normalize(vmin=v_min,vmax=v_max)
mapper = ScalarMappable(norm=norm,cmap=cmap)
nodatacolor="black"

# Loop over the grid and plot the data
for ii in range(0,n_p-1):
    for jj in range(0,nl-1):
        if(count[ii,jj]>0):
            colors = mapper.to_rgba(UVAI[ii,jj])
            
            lon0 = LONalt[ii]
            lat0 = LATalt[jj]
            lat1 = LATalt[jj+1]
            lon1 = LONalt[ii+1]
       
            if(lat1>latmin):
     
                y = [lat0,lat1,lat1,lat0]
                x = [lon0,lon0,lon1,lon1]
                mx, my = m(x,y)
                mxy = zip(mx,my)
                pair1 = (mx[0],my[0])
                pair2 = (mx[1],my[1])
                pair3 = (mx[2],my[2])
                pair4 = (mx[3],my[3])
                # Plot the box on the map using color
                color2 = rgb2hex(colors)
                poly = Polygon([pair1,pair2,pair3,pair4],facecolor=color2,edgecolor=color2)
                plt.gca().add_patch(poly)
        

plt.title('OMI Assimilated Aerosol Optical Depth '+plot_date)
#plt.title('OMI Assimilated Aerosol Index '+plot_date)
cax = fig.add_axes([0.16,0.075,0.7,0.025])
cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
#cb.ax.set_xlabel('Aerosol Optical Depth')
cb.ax.set_xlabel('Aerosol Index')
#out_name = 'omi_assimilate_aod_'+plot_date+'.png'       
out_name = 'omi_assimilate_aod_'+plot_date+'_10282019.png'       
plt.savefig(out_name)
print('Saved image '+out_name)
