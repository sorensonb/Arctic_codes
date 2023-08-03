#!/usr/bin/env python

"""


"""

import importlib, FuLiouLib
from FuLiouLib import *
import sys

file_list = 'fuliouout_file_list.txt'
fu_data = read_fuliou_cloudaer_output(file_list)

sys.exit()

infile = 'Ed4_LaRC_FuLiou/src/simple/fuliou_profile_aertype11_alb1_aertop22_aerbot24_cldhgt28.txt'

data = pd.read_csv(infile, delim_whitespace = True, header = 0)

z1s   = np.asarray(data['Z1'].values)
z2s   = np.asarray(data['Z2'].values)
aprof = np.asarray(data['aprofs'].values)
cld   = np.asarray(data['cldpres'].values)

cloud_idx = np.where(cld == 1)

#ax5.axhspan(CrIS_data_day_smoke['press'][tmp_idxs[ii][0]], CrIS_data_day_smoke['press'][tmp_idxs[ii][-1]], color = 'purple', alpha = 0.5)
fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
ax1.plot(aprof, z1s)
ax2.plot(np.cumsum(aprof), z1s)
ax1.set_ylim([0, 30])
ax2.set_ylim([0, 30])
ax1.set_ylabel('Height [km]')
ax2.set_ylabel('Height [km]')
ax1.set_xlabel('Aerosol profile')
ax2.set_xlabel('Cumulative aer. prof.')
ax1.axhspan(z1s[cloud_idx], z2s[cloud_idx], color = 'purple', alpha = 0.5)
fig.tight_layout()
plt.show()

sys.exit()

nv = 36
fill_val = 0.5
initial_aprof = np.full(nv, fill_val)
final_total = 100.
aeridx_top = 22
aeridx_bot = 28
num_idx    = (aeridx_bot - aeridx_top) + 1
init_total_without_aeridx = final_total - (fill_val * (nv - num_idx))
prof_val_even = init_total_without_aeridx / num_idx
final_aprof = initial_aprof
final_aprof[aeridx_top:aeridx_bot+1] = prof_val_even
xvals = np.arange(len(final_aprof))

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
ax1.plot(final_aprof, xvals)
ax2.plot(np.cumsum(final_aprof), xvals)
ax1.set_title('Prof value = '+str(np.round(prof_val_even,2)))
fig.tight_layout()
plt.show()

sys.exit()

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
