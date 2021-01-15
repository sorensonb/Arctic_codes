#!/usr/bin/env python
"""
  NAME:

  PURPOSE:

  SYNTAX:

  EXAMPLE:

  MODIFICATIONS:
    2020/09/10

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
import cartopy.crs as ccrs

if(len(sys.argv)<3):
    print("SYNTAX: python plot_ICESat2.py file variable")
    print("        variable options:")
    print("        - sea_ice_height")
    sys.exit()

var_dict = {
    'sea_ice_height': 'gt2r/sea_ice_segments/heights/height_segment_height',
    'sea_ice_conc':   'gt2r/sea_ice_segments/stats/ice_conc' 
}

max_dict = {
    'sea_ice_height': 5.0,
    'sea_ice_conc': 102.
}

infile = sys.argv[1]
variable = sys.argv[2]

#new_time = '20130611t1322_new'
data_path = '/home/bsorenson/data/ICESat2/'
new_base_path = '/home/shared/OMAERUV/OMAERUV_Parameters/new_files/'

# Read in the data
data = h5py.File(infile,'r')

# Pull out lat and lon
latitude  = data['gt2r/sea_ice_segments/latitude'][:]
longitude = data['gt2r/sea_ice_segments/longitude'][:]

plot_data = data[var_dict[variable]][:]


# mask bad data
divider = 10
mask_data = np.ma.masked_where(plot_data > 1000, plot_data)
split_data = np.array_split(mask_data,len(mask_data)/divider)
avg_split = np.full(int(len(plot_data)/divider),-99.)

for i,avg_arr in enumerate(split_data):
    if(avg_arr.count() > 0):
        avg_split[i] = np.nanmean(avg_arr)

#[np.nanmean(avg_arr) for avg_arr in np.array_split(test_mask,len(test_mask)/10)]


# Make the figure
plt.close()
plt.figure()
#plt.scatter(np.arange(len(plot_data)),plot_data,s=4)
plt.scatter(np.arange(len(avg_split)),avg_split,s=4)
plt.ylim(np.min(plot_data)*0.98,max_dict[variable])
plt.title(variable)
plt.grid()
plt.show()

data.close()
