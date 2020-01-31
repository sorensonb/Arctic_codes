#!/usr/bin/env python
"""


"""

import numpy as np
import numpy.ma as ma
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy
import cartopy.crs as ccrs
import glob
import subprocess
from IceLib import read_ice,ice_trendCalc


data_loc = '/home/bsorenson/HighLatitudeStudy/Ice_Analysis/data/'

if(len(sys.argv)!=2):
    print("SYNTAX: ./ice_trend.py season")
    print("        season = 'spring','summer','autumn','winter','all'")
    sys.exit()

thielSen=False
season=sys.argv[1]

season_adder = ' '+season
file_adder = '_'+season
if(season=='all'):
    season_adder = ''
    file_adder = ''

ice_data = read_ice(season)

"""
spring=False
summer=False
autumn=False
winter=False
season_adder = ' '+season
file_adder = '_'+season
if(season=='spring'):
    spring=True
    ls_check = ['03','04','05']
elif(season=='summer'):
    summer=True
    ls_check = ['06','07','08']
elif(season=='autumn'):
    autumn=True
    ls_check = ['09','10','11']
elif(season=='winter'):
    winter=True
    ls_check = ['12','01','02']
else:
    season_adder = ''
    file_adder = ''
    ls_check=['01','02','03','04','05','06','07','08','09','10','11','12']


# Grab all the ice files
#file_names = glob.glob(data_loc+'*.bin')
# Using this syntax, ignores 2000
cmnd = "ls /home/bsorenson/Research/Ice_analysis/data/nt_*.bin"
#status,output = commands.getstatusoutput(cmnd)
#file_initial = output.strip().split('\n')
#
#file_names = []
#for fname in file_initial:
#    if(fname[50:52] in ls_check):
#        file_names.append(fname)


file_names = subprocess.check_output(cmnd,shell=True).decode('utf-8').strip().split('\n')

# Read in the latitude and longitude data
latfileo = open('/home/bsorenson/Research/Ice_analysis/psn25lats_v3.dat','r')
lats = np.reshape(np.fromfile(latfileo,dtype=np.uint32)/100000.,(448,304))
lonfileo = open('/home/bsorenson/Research/Ice_analysis/psn25lons_v3.dat','r')
tlons = np.fromfile(lonfileo,dtype=np.uint32)/100000.
tlons[tlons>180.] = tlons[tlons>180.]-42949.67296
lons = np.reshape(tlons,(448,304))
#lons = np.reshape(np.fromfile(lonfileo,dtype=np.uint32)/100000.,(448,304))

num_files = len(file_names)
ice_data = {}
ice_data['data'] = np.full((num_files,448,304),-9.)
ice_data['lat']  = lats
ice_data['lon']  = lons
ice_data['trends'] = np.full((448,304),-9.)
ice_data['land_trends'] = np.full((448,304),-9.)
ice_data['titles'] = []

count = 0
for fname in file_names:

    total_name = fname
    #total_name = data_loc+fname
    print(total_name)
    fileo = open(total_name,'rb')
    test = bytearray(fileo.read())
    
    ## To print the header
    #print(test[0:299])
    
    # psg = polar stereographic grid
    
    missing_val             = int(test[0:5].decode('utf-8'))
    ncol_psg                = int(test[6:11].decode('utf-8'))
    nrow_psg                = int(test[12:17].decode('utf-8'))
    lat_enclosed_psg        = float(test[24:29].decode('utf-8'))
    greenwich_orient_psg    = float(test[30:35].decode('utf-8'))
    j_coord_pole            = float(test[42:47].decode('utf-8'))
    i_coord_pole            = float(test[48:53].decode('utf-8'))
    inst_desc               = test[54:59].decode('utf-8')
    data_desc               = test[60:65].decode('utf-8')
    julian_start_day        = int(test[66:71].decode('utf-8'))
    start_hour              = int(test[72:77].decode('utf-8'))
    start_min               = int(test[78:83].decode('utf-8'))
    julian_end_day          = int(test[84:89].decode('utf-8'))
    end_hour                = int(test[90:95].decode('utf-8'))
    end_min                 = int(test[96:101].decode('utf-8'))
    year                    = int(test[102:107].decode('utf-8'))
    julian_day              = int(test[108:113].decode('utf-8'))
    channel_desc            = int(test[114:119].decode('utf-8'))
    scaling_factor          = float(test[120:125].decode('utf-8'))
    file_name               = test[126:149].decode('utf-8')
    image_title             = test[150:229].decode('utf-8')
    data_info               = test[230:299].decode('utf-8')
    
    # Move the data into a np array
    data = np.reshape(np.array(test[300:]),(448,304))
    
    #for i in range(len(data)):
    #    templine = data[i,:].astype(np.float)
    #    templine[templine<251] = (templine[templine<251]/scaling_factor)*100.
    #    templine = templine.astype(int)
    #    print(templine)
    
    #data[np.where(data>=251)]=data[np.where(data>=251)]*0.
    data[np.where(data<251)]=  \
        (data[np.where(data<251)]/scaling_factor)*100.

    ice_data['data'][count,:,:] = data[:,:]
    ice_data['titles'].append(file_name)
#    total_data[:,:,count] = data[:,:]
    count+=1
"""

ice_data = ice_trendCalc(ice_data,thielSen=thielSen)

"""
# Loop over the data and calculate trends
for i in range(448):
    print(i)
    max_trend = -99.
    min_trend = 99.
    for j in range(304):
        ice_data['trends'][i,j],temp_pcnt = trend_calc(ice_data,i,j,thielsen=thielSen)
        ice_data['land_trends'][i,j],temp_pcnt = trend_calc(ice_data,i,j,thielsen=thielSen)
        temp_trend = ice_data['trends'][i,j]
        if(temp_trend>max_trend):
            max_trend = temp_trend
        if(temp_trend<min_trend):
            min_trend = temp_trend
    print("  max trend = ",max_trend)
    print("  min trend = ",min_trend)
    # Deal with land masks
    good_indices = np.where(ice_data['data'][0,i,:]<251.)
    land_indices = np.where(ice_data['data'][0,i,:]>=251.)
    ice_data['trends'][i,land_indices] = np.nan
    ice_data['land_trends'][i,good_indices] = np.nan
"""

#data = ice_data['trends']
plot_good_data = ma.masked_invalid(ice_data['trends'].T)
plot_land_data = ma.masked_invalid(ice_data['land_trends'].T)
fig1 = plt.figure(figsize=(9,5))
plt.pcolormesh(plot_good_data,cmap=plt.cm.bwr,vmin=-50,vmax=50)
plt.colorbar(label='Percent Ice Concentration Trend (%)')
plt.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.title('NSIDC Sea Ice Concentration'+season_adder.title()+' Trends\nJan 2001 to Dec 2018')
plt.savefig('ice_trend_200101_201812'+file_adder+'.png',dpi=300)
print("Saved image ice_trend_200101_201812"+file_adder+".png")
plt.show()

#fileo.close()
