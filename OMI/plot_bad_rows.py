#!/usr/bin/env python
"""
  NAME:
    plot_bad_rows.py    

  PURPOSE:

  SYNTAX:

  EXAMPLE:

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2021/07/22:
      Written
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import pandas as pd

if(len(sys.argv)<2):
    print("SYNTAX: python plot_OMI_fort_out.py bad_row_file [xtrack_row_file]")
    sys.exit()

file_name = sys.argv[1]
do_xtrack = False
bad_color = 'jet'
xtrack_add = ''
if(len(sys.argv) == 3):
    bad_color = 'bwr_r'
    xtrack_color = 'jet'
    xtrack_add = '_xtrack'
    do_xtrack = True
    xtrack_name = sys.argv[2]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Extract date information from the file name
name_split = file_name.strip().split('/')[-1].split('_')
start_year = name_split[-2]
end_year   = name_split[-1].split('.')[0]

bad_dict = {}
max_num = 0

str_dates = []

# Read bad row file
# -----------------
with open(file_name,'r') as f_in:
    flines = f_in.readlines()
    total_array = np.zeros((len(flines),61))
    #int_dates = np.zeros(len(flines))
    for ii, line in enumerate(flines):
    #for line in f_in:
        templine = line.strip().split()
        #int_dates[ii] = int(templine[0])
        str_dates.append(templine[0][:4])
        #int_dates[ii] = int(templine[0])
        #bad_dict[templine[0]] = np.zeros(int(templine[1]))
        if(int(templine[1]) > 0):
            for jj in templine[2:]:
                total_array[ii,int(jj)] = 1
            #bad_dict[templine[0]][ii] = int(templine[2+ii])

if(do_xtrack):
    str_dates_X = []
    
    # Read xtrack row file
    # -----------------
    with open(xtrack_name,'r') as f_in:
        flines = f_in.readlines()
        total_array_X = np.zeros((len(flines),61))
        #int_dates = np.zeros(len(flines))
        for ii, line in enumerate(flines):
        #for line in f_in:
            templine = line.strip().split()
            #int_dates[ii] = int(templine[0])
            str_dates_X.append(templine[0][:4])
            #int_dates[ii] = int(templine[0])
            #bad_dict[templine[0]] = np.zeros(int(templine[1]))
            if(int(templine[1]) > 0):
                for jj in templine[2:]:
                    total_array_X[ii,int(jj)] = 1
                #bad_dict[templine[0]][ii] = int(templine[2+ii])
    
    x_vals_X = np.arange(total_array_X.shape[0])
    y_vals_X = np.arange(total_array_X.shape[1])
    mask_array_X = np.ma.masked_where(total_array_X == 0,total_array_X)


#index = 0
#for key in bad_dict:
#    total_array[index,:(len(bad_dict[key]))] = bad_dict[key]
#    index += 1

x_vals = np.arange(total_array.shape[0])
y_vals = np.arange(total_array.shape[1])
mask_array = np.ma.masked_where(total_array == 0,total_array)

fig1, ax = plt.subplots(figsize=(10,4))
ax.pcolormesh(x_vals,y_vals,mask_array.T,shading='auto',cmap=bad_color)
if(do_xtrack):
    ax.pcolormesh(x_vals_X,y_vals_X,mask_array_X.T,shading='auto',cmap=xtrack_color)
ax.set_yticks(y_vals,minor=True)
ax.set_xticks(x_vals[::183])
ax.set_xticklabels(str_dates[::183])
ax.grid(which='minor',axis='y')
ax.grid(which='major',axis='y',linewidth=1.0,color='black')
ax.set_xlabel('Year')
ax.set_ylabel('Row Number')

outname = 'bad_rows_'+start_year + '_'+end_year+xtrack_add + '_2.png'
plt.savefig(outname,dpi=300)
print("Saved image",outname)
plt.show()
