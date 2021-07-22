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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Extract date information from the file name
name_split = file_name.strip().split('/')[-1].split('_')
version = name_split[-1].split('.')[0]
#start_year = name_split[2]
#end_year   = name_split[3]
#ai_thresh  = float(int(name_split[4])/100)
#dtype      = name_split[5].split('.')[0]

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

fig1, ax = plt.subplots(figsize=(8,4))
#fig1.set_figsize_inches(6,3)
ln1 = ax.plot(dt_dates,avg_ai,label='AI')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
ax2 = ax.twinx()
ln2 = ax2.plot(dt_dates,counts,color='tab:orange',label='counts')
lns = ln1 + ln2
labs = [l.get_label() for l in lns]
ax.legend(lns,labs)
ax.set_ylabel('Average AI',color='tab:blue')
ax2.set_ylabel('Ob counts',color='tab:orange')
ax.set_title('AI Counts in Remote Ocean')
ax.grid()
plt.subplots_adjust(left=0.09,right=0.88)
plt.savefig('omi_'+version+'_drift.png',dpi=300)
print("Saved image omi_"+version+"_drift.png")

#fig2 = plt.figure()
#plt.plot(xrange_00,count_00,label='00Z')
#plt.plot(xrange_06,count_06,label='06Z')
#plt.plot(xrange_12,count_12,label='12Z')
#plt.plot(xrange_18,count_18,label='18Z')
#plt.legend()
#plt.show()
