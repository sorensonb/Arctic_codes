#!/usr/bin/env python
"""


"""

import numpy as np
import sys
import matplotlib.pyplot as plt

data_loc = '/home/bsorenson/HighLatitudeStudy/Ice_Analysis/data/'

if(len(sys.argv)!=2):
    print("SYNTAX: ./read_ice_bin.py ice_bin_file")
    sys.exit()

infile = sys.argv[1]
fileo = open(infile,'rb')
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

data[np.where(data>=251)]=data[np.where(data>=251)]*0.
data[np.where(data<251)]=(data[np.where(data<251)]/scaling_factor)*100.

data[0,0] = 100.

plt.pcolormesh(data.T)
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.title(file_name)
plt.show()

fileo.close()
