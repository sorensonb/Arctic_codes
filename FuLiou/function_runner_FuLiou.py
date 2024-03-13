#!/usr/bin/env python

"""


"""

import importlib, FuLiouLib
from FuLiouLib import *
import sys

dates = ['20220909','20220915','20220922']
for date in dates:
    old = h5py.File('fuliou_output_naaps_' + date + '.h5')
    new = h5py.File('fuliou_output_naaps_' + date + '_2stream.h5')

    print(date)
    print("NAAPS:")
    print(old['net_heating_rate'][10,:])
    print(new['net_heating_rate'][10,:])

    old.close()
    new.close()

    old = h5py.File('fuliou_output_lidar_' + date + '_highres3.h5')
    new = h5py.File('fuliou_output_lidar_' + date + '_2stream.h5')

    print("NAAPS:")
    print(old['net_heating_rate'][10,:])
    print(new['net_heating_rate'][10,:])

    old.close()
    new.close()


sys.exit()

data = h5py.File('fuliou_output_lidar_20220922.h5')

time      = data['time'][:]
mask_heat = data['net_heating_rate'][:,:]
mask_ext  = data['dust_ext'][:,:]
mask_ext = np.ma.masked_where(mask_ext <= 0, mask_ext)
mask_ext = np.log10(mask_ext)
#mask_heat = np.ma.masked_invalid(mask_heat)

xvals = np.arange(mask_heat.shape[0])
yvals = data['press_lvl'][0,:-1]
yvals2 = data['press_lvl'][0,:]

data.close()

fig = plt.figure(figsize = (9, 6))
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
mesh = ax1.pcolormesh(xvals, yvals, mask_heat.T, shading = 'auto', vmin = -10, vmax = 10, cmap = 'bwr')
cbar = fig.colorbar(mesh, ax = ax1, label = 'Net Heating Rate')
ax1.invert_yaxis()
ax1.set_title('Heating Rate')

mesh = ax2.pcolormesh(xvals, yvals2, mask_ext.T, shading = 'auto', vmin = -6, vmax = -3.0)
cbar = fig.colorbar(mesh, ax = ax2, label = 'Log10 of Dust Extinction')
ax2.invert_yaxis()
ax2.set_title('Dust extinction')
fig.tight_layout()
plt.show()


sys.exit()


date_str = '20220712081447'
lat = 12.55566
lon = -17.51049
calc_solar_zenith_angle(date_str, lat, lon, use_pysolar = True)

sys.exit()

file_list = 'fuliouout_file_list.txt'
fu_data = read_fuliou_cloudaer_output(file_list)

##plot_FuLiou_cloudaer_multivar(fu_data, 'SolarZenith', 'CloudOptDepth', fu_var = 'SWF_TOTAL', \
#plot_FuLiou_cloudaer_multivar(fu_data, 'CloudOptDepth','CloudFrac',  fu_var = 'SWF_TOTAL', \
#    aerhgt_idx = 0, cldhgt_idx = 0, zen_idx = 0, cod_idx = 5, \
#    cldfrac_idx = 5, divide_by_aod = True, save = True)
#plot_FuLiou_cloudaer_multivar_mesh(fu_data, 'SolarZenith', 'CloudOptDepth', fu_var = 'SWF_TOTAL', \
plot_FuLiou_cloudaer_multivar_mesh(fu_data, 'CloudOptDepth', 'CloudFrac', fu_var = 'SWF_TOTAL', \
    aerhgt_idx = 0, cldhgt_idx = 0, zen_idx = 10, cod_idx = 5, \
    cldfrac_idx = 5, divide_by_aod = True, save = False)


sys.exit()

sfc_types = ['Ocean','Land','Ice']

aer_titles = ['Continental Aerosol','Soot (BC) Aerosol']

aertype_idx = 0
sfc_idx = 2
cod_idx = 5
aerhgt_idx = 0
cldhgt_idx = 0
zen_idx = 0
cldfrac_idx = 5

#fu_var = 'SWF_TOTAL_NO_AER'
fu_var = 'SWF_TOTAL'

#delt_swf_delt_aod = fu_data[fu_var][1,sfc_idx,aerhgt_idx,cldhgt_idx,:,cod_idx,cldfrac_idx,0] - \
#                    fu_data[fu_var][1,sfc_idx,aerhgt_idx,cldhgt_idx,:,cod_idx,cldfrac_idx,-1]

# new_dims: ['aertype','ALB','aertop_bot','cldhgt','SolarZenith','CloudOptDepth','CloudFrac','AOD']
xvar = 'CloudOptDepth'
zvar = 'aertop_bot'

fig = plt.figure(figsize = (11, 7))
axs = fig.subplots(2,3)
for aertype_idx in range(2):
    local_title = aer_titles[aertype_idx]
    #for sfc_idx, ax in enumerate(axs):
    for sfc_idx in range(3):
        local_sfctype = sfc_types[sfc_idx]
        #for cldfrac_idx in range(fu_data['CloudFrac'].shape[0]):
        for cldfrac_idx in range(fu_data['CloudFrac'].shape[0]):
            delt_swf_delt_aod = (fu_data[fu_var][aertype_idx,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,:,cldfrac_idx,-1] - \
                                fu_data[fu_var][aertype_idx,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,:,cldfrac_idx,0])
            #delt_swf_delt_aod = (fu_data[fu_var][aertype_idx,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,:,cldfrac_idx,-1] - \
            #                    fu_data[fu_var][aertype_idx,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,:,cldfrac_idx,0]) / \
            #                    (fu_data['AOD'][-1] - fu_data['AOD'][0])
            #ax.plot(fu_data['SolarZenith'], delt_swf_delt_aod)
            axs[aertype_idx,sfc_idx].plot(fu_data['CloudOptDepth'], delt_swf_delt_aod, \
                label = 'CldFrac = ' + str(fu_data['CloudFrac'][cldfrac_idx]))

        #ax.set_xlabel('SZA')
        axs[aertype_idx,sfc_idx].set_xlabel('COD')
        axs[aertype_idx,sfc_idx].set_ylabel('ΔSWF/ΔAOD')
        #axs[aertype_idx,sfc_idx].set_title(sfc_types[sfc_idx])
        axs[aertype_idx,sfc_idx].set_title(local_title + '\n' + local_sfctype)

axs[0,2].legend()

fig.tight_layout()

plt.show()
sys.exit()

fig = plt.figure(figsize = (11,7))
axs = fig.subplots(2,3)
#ax1 = fig.add_subplot(2,3,1)
#ax2 = fig.add_subplot(2,3,2)
#ax3 = fig.add_subplot(2,3,3)
#ax4 = fig.add_subplot(2,3,4)
#ax5 = fig.add_subplot(2,3,5)
#ax6 = fig.add_subplot(2,3,6)

for aertype_idx in range(2):
    local_title = aer_titles[aertype_idx]
    for sfc_idx in range(3):
        local_sfctype = sfc_types[sfc_idx]

        #for cod_idx in range(fu_data['CloudOptDepth'].shape[0]):
        for cod_idx in range(0,10,2):
            axs[aertype_idx,sfc_idx].plot(fu_data['AOD'], \
                fu_data[fu_var][aertype_idx,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,cod_idx,cldfrac_idx,:], \
                label =  str(fu_data['CloudOptDepth'][cod_idx]))

        axs[aertype_idx,sfc_idx].set_ylabel('TOA SWF [W/m2]')
        axs[aertype_idx,sfc_idx].set_xlabel('AOD (0.55 μm)')
        axs[aertype_idx,sfc_idx].set_title(local_title + '\n' + local_sfctype)

fig.tight_layout()
plt.show()

sys.exit()

# new_dims: ['aertype','ALB','aertop_bot','cldhgt','SolarZenith','CloudOptDepth','CloudFrac','AOD']
ax1.plot(fu_data['AOD'], fu_data[fu_var][0,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,0,cldfrac_idx,:], label =  str(fu_data['CloudOptDepth'][0]))
ax1.plot(fu_data['AOD'], fu_data[fu_var][0,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,2,cldfrac_idx,:], label =  str(fu_data['CloudOptDepth'][2]))
ax1.plot(fu_data['AOD'], fu_data[fu_var][0,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,4,cldfrac_idx,:], label =  str(fu_data['CloudOptDepth'][4]))
ax1.plot(fu_data['AOD'], fu_data[fu_var][0,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,6,cldfrac_idx,:], label =  str(fu_data['CloudOptDepth'][6]))
ax1.plot(fu_data['AOD'], fu_data[fu_var][0,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,8,cldfrac_idx,:], label =  str(fu_data['CloudOptDepth'][8]))
ax1.legend()

ax1.set_ylabel('TOA SWF [W/m2]')
ax1.set_xlabel('AOD (0.55 μm)')
ax1.set_title('Continental Aerosol\nIncreasing COD')

ax2.plot(fu_data['AOD'], fu_data[fu_var][1,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,0,cldfrac_idx,:], label =  str(fu_data['CloudOptDepth'][0]))
ax2.plot(fu_data['AOD'], fu_data[fu_var][1,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,2,cldfrac_idx,:], label =  str(fu_data['CloudOptDepth'][2]))
ax2.plot(fu_data['AOD'], fu_data[fu_var][1,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,4,cldfrac_idx,:], label =  str(fu_data['CloudOptDepth'][4]))
ax2.plot(fu_data['AOD'], fu_data[fu_var][1,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,6,cldfrac_idx,:], label =  str(fu_data['CloudOptDepth'][6]))
ax2.plot(fu_data['AOD'], fu_data[fu_var][1,sfc_idx,aerhgt_idx,cldhgt_idx,zen_idx,8,cldfrac_idx,:], label =  str(fu_data['CloudOptDepth'][8]))
ax2.legend()

ax2.set_ylabel('TOA SWF [W/m2]')
ax2.set_xlabel('AOD (0.55 μm)')
ax2.set_title('Soot (BC) Aerosol\nIncreasing COD')

plt.suptitle('Cloud Faction = ' + str(fu_data['CloudFrac'][cldfrac_idx]))

fig.tight_layout()

plt.show()

sys.exit()
plot_FuLiou_boxes(save = True)





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
