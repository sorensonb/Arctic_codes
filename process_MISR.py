#!/home/blake.sorenson/anaconda3/bin/python

"""


"""

from subprocess import check_output
from netCDF4 import Dataset
import numpy as np
import sys

if(len(sys.argv) != 3):
    print("SYNTAX: ./process_MISR.py year month")
    sys.exit()

in_year = sys.argv[1]
in_month = sys.argv[2]

latmin = 45.

lat_ranges = np.arange(latmin,90,1.0)
lon_ranges = np.arange(-180.,180.,1.0)

# Declare arrays
misr_aod = np.zeros((len(lat_ranges),len(lon_ranges)))
misr_counts = np.zeros((len(lat_ranges),len(lon_ranges)))

# List directories in the main MISR directory
data_path = "/data/MISR/l5ftl01.larc.nasa.gov/misrl2l3/MISR/MIL2ASAE.003/"
dirs_list = check_output('ls '+data_path,shell=True).decode('utf-8').strip().split('\n')

for tdir in dirs_list:  
    dir_split = tdir.split('.')
    if(len(dir_split) > 1):
        if((dir_split[0] == in_year) & (dir_split[1] == in_month)):
            # Loop over the netCDF files in this directory
            print(tdir)
            new_list = check_output('ls '+data_path + tdir + '/*.nc',shell=True).decode('utf-8').strip().split('\n')

            # Loop over each file here, read data, calculate averages within each grid box
            for tfile in new_list:
                print(tfile)
                data = Dataset(tfile,'r')
                
                LAT = data['4.4_KM_PRODUCTS']['Latitude'][:]
                LON = data['4.4_KM_PRODUCTS']['Longitude'][:]
                AOD = data['4.4_KM_PRODUCTS']['Aerosol_Optical_Depth'][:]
                #AOD = data['4.4_KM_PRODUCTS/AUXILIARY']['Aerosol_Optical_Depth_Raw'][:]
  
                for xi in range(len(lat_ranges)-1):
                    #print(lat_ranges[xi])
                    check_lon = LON[np.where((LAT >= lat_ranges[xi]) & (LAT < lat_ranges[xi+1]))]

                    if(len(check_lon) > 0):
    
                        # Find the lon index that matches with the minimum longitude here
                        min_lon_idx = np.where(lon_ranges[:] > np.min(check_lon))[0][0]
                        max_lon_idx = np.where(lon_ranges[:] < np.max(check_lon))[0][-1]
                        for yj in range(min_lon_idx,max_lon_idx):
                        #for yj in range(len(lon_ranges)-1):
                            good_AOD = AOD[np.where((AOD > -900) & \
                                        ((LAT > lat_ranges[xi]) & (LAT < lat_ranges[xi+1])) & \
                                        ((LON > lon_ranges[yj]) & (LON < lon_ranges[yj+1])))]
                            if(len(good_AOD) > 0):
                                if(misr_counts[xi,yj] == 0):
                                    misr_aod[xi,yj] = np.nanmean(good_AOD) 
                                    misr_counts[xi,yj] = len(good_AOD)
                                else:
                                    misr_aod[xi,yj] = ((misr_aod[xi,yj]*misr_counts[xi,yj])+np.nanmean(good_AOD))/\
                                                        (misr_counts[xi,yj]+len(good_AOD))
                                    misr_counts[xi,yj] += len(good_AOD)

# Save data to an output file
in_year = sys.argv[1]
in_month = sys.argv[2]

num_lat = len(lat_ranges)
num_lon = len(lon_ranges)

outname = 'misr_gridded_aod_'+in_year+in_month+'.nc'

nc = Dataset(outname,'w',format='NETCDF4')
n_lon = nc.createDimension('dlon',num_lon)
n_lat = nc.createDimension('dlat',num_lat)

LAT = nc.createVariable('Latitude','i2',('dlat','dlon'))
LAT.description = 'Latitude'
LAT.units = 'Degrees'
LON = nc.createVariable('Longitude','i2',('dlat','dlon'))
LON.description = 'Longitude'
LON.units = 'Degrees'
AOD = nc.createVariable('AOD','f4',('dlat','dlon'))
AOD.description('Monthly Averaged Aerosol Optical Depth for '+in_year + in_month)
OB_COUNT = nc.createVariable('OB_COUNT','i2',('dlat','dlon'))
OB_COUNT.description = '# of MISR AOD measurements used in each monthly average'

for xi in range(num_lat):
    for yj in range(num_lon):
        LAT[xi,yj] = lat_ranges[xi]
        LON[xi,yj] = lon_ranges[yj]

for xi in range(num_lat):
    for yj in range(num_lon):
        if(misr_counts[xi,yj] == 0.0):
            AOD[xi,yj] = -999.9
            OB_COUNT[xi,yj] = -99
        else:
            AOD[xi,yj] = misr_aod[xi,yj]
            OB_COUNT[xi,yj] = misr_counts[xi,yj]

nc.close()

print("Saved file",outname)
