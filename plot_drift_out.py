#!/usr/bin/env python
"""
  NAME:
    plot_drift_out.py

  PURPOSE:

  SYNTAX:

  EXAMPLE:

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2021/04/01:
      Written (modification of AerosolIndexPlotting_Procedure_August112014t0046.pro by rcontreras) 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import pandas as pd

if(len(sys.argv)<2):
    print("SYNTAX: python plot_drift_out.py out_file")
    sys.exit()
     
file_name = sys.argv[1]
if(len(sys.argv) == 3):
    file_name2 = sys.argv[2]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Extract date information from the file name
name_split = file_name.strip().split('/')[-1].split('_')
version = name_split[-1].split('.')[0]
if(version == 'jz2'):
    title_string = 'Only good rows'
elif(version == 'jz27'):
    title_string = 'Only rows 49, 50, 53, and 56-60'
elif(version == 'jz28'):
    title_string = 'Only rows 49, 50, and 56-60'
elif(version == 'jz210'):
    title_string = 'Only good rows, no 43 for Jul/Aug 2019'
elif(version == 'bs2'):
    title_string = 'Only rows 1 - 21'
else:
    title_string = ''
if(len(sys.argv) == 3):
    name_split2 = file_name2.strip().split('/')[-1].split('_')
    version2 = name_split2[-1].split('.')[0]
    if(version2 == 'jz2'):
        title_string2 = 'Only good rows'
    elif(version2 == 'jz27'):
        title_string2 = 'Only rows 49, 50, 53, and 56-60'
    elif(version2 == 'jz28'):
        title_string2 = 'Only rows 49, 50, and 56-60'
    elif(version2 == 'jz210'):
        title_string2 = 'Only good rows, no 43 for Jul/Aug 2019'
    elif(version2 == 'bs2'):
        title_string2 = 'Only rows 1 - 21'
    else:
        title_string2 = ''

# Read data from the first file
# -----------------------------
in_data = pd.read_csv(file_name, delim_whitespace=True)
dates  = in_data['Date'].values
dt_dates = [datetime.strptime(str(tmpdate),"%Y%m") \
    for tmpdate in dates]
avg_ai  = in_data['Avg'].values
counts  = in_data['#_obs'].values
#count70  = in_data['Cnt70'].values
#count75  = in_data['Cnt75'].values
#count80  = in_data['Cnt80'].values
#count85  = in_data['Cnt85'].values
x_range = np.arange(len(avg_ai))

if(len(sys.argv) == 3):
    # Read data from the first file
    # -----------------------------
    in_data2 = pd.read_csv(file_name2, delim_whitespace=True)
    dates2  = in_data2['Date'].values
    dt_dates2 = [datetime.strptime(str(tmpdate),"%Y%m") \
        for tmpdate in dates2]
    avg_ai2  = in_data2['Avg'].values
    counts2  = in_data2['#_obs'].values
    #count70  = in_data['Cnt70'].values
    #count75  = in_data['Cnt75'].values
    #count80  = in_data['Cnt80'].values
    #count85  = in_data['Cnt85'].values
    x_range2 = np.arange(len(avg_ai2))

if(len(sys.argv) == 3):
    version = version + 'd'+ version2
    fig1, (ax1, ax2) = plt.subplots(2, figsize=(8,8))
else:
    fig1, ax1 = plt.subplots(figsize=(8,4))
    #fig1.set_figsize_inches(6,3)
ln11 = ax1.plot(dt_dates,avg_ai,label='AI')
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
ax12 = ax1.twinx()
ln12 = ax12.plot(dt_dates,counts,color='tab:orange',label='counts')
lns1 = ln11 + ln12
labs1 = [l.get_label() for l in lns1]
ax1.set_ylabel('Average AI',color='tab:blue')
ax12.set_ylabel('Ob counts',color='tab:orange')
ax1.set_title(title_string)
ax1.grid()
if(len(sys.argv) == 3):
    ax1.tick_params(axis='x',bottom=False,labelbottom=False)
    ln21 = ax2.plot(dt_dates2,avg_ai2,label='AI')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax22 = ax2.twinx()
    ln22 = ax22.plot(dt_dates2,counts2,color='tab:orange',label='counts')
    lns2 = ln21 + ln22
    labs2 = [l.get_label() for l in lns2]
    ax2.legend(lns2,labs2,loc='upper center', bbox_to_anchor=(0.5,-0.15),\
        ncol=2)
    ax2.text(1,ax2.get_ylim()[1]*0.9,'B')
    ax2.set_ylabel('Average AI',color='tab:blue')
    ax22.set_ylabel('Ob counts',color='tab:orange')
    ax2.set_title(title_string2)
    ax2.grid()
    ax2.set_xlabel('Year')
else:
    ax1.legend(lns1,labs1)
    ax1.set_xlabel('Year')

plt.subplots_adjust(left=0.09,right=0.88)
plt.savefig('omi_'+version+'_drift.png',dpi=300)
print("Saved image omi_"+version+"_drift.png")
plt.show()
#fig2 = plt.figure()
#plt.plot(xrange_00,count_00,label='00Z')
#plt.plot(xrange_06,count_06,label='06Z')
#plt.plot(xrange_12,count_12,label='12Z')
#plt.plot(xrange_18,count_18,label='18Z')
#plt.legend()
#plt.show()
