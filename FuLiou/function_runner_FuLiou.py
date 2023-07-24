#!/usr/bin/env python

"""


"""

import importlib, FuLiouLib
from FuLiouLib import *
import sys

#data = read_fuliou_output()
#data01 = read_fuliou_cloud_output(aertype = 1, spec_type = 'loftplume')
#data11 = read_fuliou_cloud_output(aertype = 11, spec_type = 'loftplume')
data01 = read_fuliou_cloud_output(aertype = 1)
data11 = read_fuliou_cloud_output(aertype = 11)

fig = plt.figure(figsize = (8,4))
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

sfc_idx = 1
zen_idx = 0
c_idx = 5

fu_var = 'SWF_TOTAL'

# dims: ['ALB','SolarZenith','CloudOptDepth','CloudFrac','AOD']
ax1.plot(data01['AOD'], data01[fu_var][sfc_idx,zen_idx,0,c_idx,:], label =  str(data01['CloudOptDepth'][0]))
ax1.plot(data01['AOD'], data01[fu_var][sfc_idx,zen_idx,2,c_idx,:], label =  str(data01['CloudOptDepth'][2]))
ax1.plot(data01['AOD'], data01[fu_var][sfc_idx,zen_idx,4,c_idx,:], label =  str(data01['CloudOptDepth'][4]))
ax1.plot(data01['AOD'], data01[fu_var][sfc_idx,zen_idx,6,c_idx,:], label =  str(data01['CloudOptDepth'][6]))
ax1.plot(data01['AOD'], data01[fu_var][sfc_idx,zen_idx,8,c_idx,:], label =  str(data01['CloudOptDepth'][8]))
ax1.legend()

ax1.set_ylabel('TOA SWF [W/m2]')
ax1.set_xlabel('AOD (0.55 μm)')
ax1.set_title('Continental Aerosol\nIncreasing COD')

ax2.plot(data11['AOD'], data11[fu_var][sfc_idx,zen_idx,0,c_idx,:], label =  str(data11['CloudOptDepth'][0]))
ax2.plot(data11['AOD'], data11[fu_var][sfc_idx,zen_idx,2,c_idx,:], label =  str(data11['CloudOptDepth'][2]))
ax2.plot(data11['AOD'], data11[fu_var][sfc_idx,zen_idx,4,c_idx,:], label =  str(data11['CloudOptDepth'][4]))
ax2.plot(data11['AOD'], data11[fu_var][sfc_idx,zen_idx,6,c_idx,:], label =  str(data11['CloudOptDepth'][6]))
ax2.plot(data11['AOD'], data11[fu_var][sfc_idx,zen_idx,8,c_idx,:], label =  str(data11['CloudOptDepth'][8]))
ax2.legend()

ax2.set_ylabel('TOA SWF [W/m2]')
ax2.set_xlabel('AOD (0.55 μm)')
ax2.set_title('Soot (BC) Aerosol\nIncreasing COD')

plt.suptitle('Cloud Faction = ' + str(data01['CloudFrac'][c_idx]))

fig.tight_layout()

plt.show()

sys.exit()
plot_FuLiou_boxes(save = True)
