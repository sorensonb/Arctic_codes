#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *
import random
#from sklearn.model_selection import train_test_split

if(sys.argv[1] == 'download'):
    l_download_shawn = True
else:
    l_download_shawn = False

#swaths_needed = [
#    #'200607250251',\
#    #'200607270239',\
#    #'200607270557',\
#    #'200607270736',\
#    #'200607272047',\
#    #'201408121758',\
#    #'201507071754',\
#    ##'201507071933',\ DIDN'T WORK
#    #'201507090748',\
#    #'201507091245',\
#    #'201708160014',\
#    ##'201708161504',\ DIDN'T WORK
#    '201708162139',\
#    '201708171050',\
#    '201708171547',\
#    '201708191713',\
#    '201808260158',\
#    '201808260337',\
#    '201808260655',
#]

#days_needed = [
#    '20060725',\
#    '20060727',\
#    '20140812',\
#    '20150707',\
#    '20150709',\
#    '20170816',\
#    '20170817',\
#    '20180826'
#]

# NEED TO FIX 20060723, 20170821, 20180811
# NO ICE DATA FOR 20180807, 20180808, 20180809
#run_list = ['20060728']

#run_list = ['20180709',\
#            '20180710',\
#run_list = ['20180711',\
#run_list = ['20180712',\
#            '20180713',\
#            '20180714',\
#            '20180715',\
#            '20180716',\
#            '20180717',\
#            '20180725',\
#            '20180726',\
#            '20180727',\
#            '20180728',\
#            '20180729',\
#            '20180730',\
#            '20180731']

#run_list = ['20150701',\
#run_list = ['20150702',\
#            '20150714',\
#            '20150715',\
#            '20150716',\
#            '20150717',\
#            '20150718',\
#            '20150719',\
#]

# REDO 20080418, 20080419, 20180630, 20190808. They didn't process many files.

run_list = ['20180630',
        ]

#run_list = ['20060731',\
#            '20140816',\
#            '20170810',\
#            '20170811',\
#            '20170821',\
#            ]
#            #'20180822',\
#            #'20180823',\
#            #'20180824',\
#            #'20180825',\
#            #'20180827',\
#            #'20180828',\
#            #'20180829',\
#            #'20180830',\

#run_list = ['20060723',\
#            '20060728',\
#            '20060731',\
#            '20080418',\
#            '20080419',\
#            '20140809',\
#            '20140810',\
#            '20140815',\
#            '20140816',\
#            '20150625',\
#            '20150626',\
#            '20170810',\
#            '20170811',\
#            '20170815']
#run_list = ['20180630',\
#            '20180702',\
#            '20180709',\
#            '20180719',\
#            '20180720',\
#            '20180807',\
#            '20180808',\
#            '20180809',\
#            '20180811',\
#            #'20180822',\
#            #'20180823',\
#            #'20180824',\
#            #'20180825',\
#            #'20180827',\
#            #'20180828',\
#            #'20180829',\
#            #'20180830',\
#            '20190807',\
#            '20190808',\
#            '20190814']


#run_list = ['20060720',\
#            '20060801',\
#            '20080418',\
#            '20080426',\
#            '20140807',\
#            '20140815',\
#            '20150623',\
#            '20170812',\
#            '20170823',\
#            '20180806',\
#            '20180818',\
#            '20190806',\
#            '20190815',\
#]




# PROCESS THESE
#run_list = ['20190811',\
#            '20190812','20190813','20190814']
#run_list = ['20190813','20190814']
#run_list = ['20180724']
#run_list = ['20180701','20180702','20180703','20180704','20180705','20180706','20180707','20180708']
#run_list = ['20060721','20060722','20060724','20060725']

#run_list = ['20060721','20060722','20060723','20060724','20060725',\
#run_list = ['20060726','20060727','20060728','20060729','20060730']
#run_list = ['20080419','20080420','20080421','20080422','20080423','20080424','20080425']
#run_list = ['20140808','20140809','20140810','20140811','20140812','20140813','20140814']
#run_list = ['20150624','20150625','20150626','20150627','20150628','20150629','20150630',\
#            '20150703','20150704','20150705','20150706','20150707','20150708','20150709',\
#            '20150710','20150711','20150712','20150713',\
#            '20170813','20170814','20170815','20170816','20170817','20170818','20170819',\
#            '20170820',\
#run_list = ['20170822',\
#            '20180701','20180702','20180703','20180704','20180705','20180706','20180707',\
#            '20180708',\
#            '20180718','20180719','20180720','20180721','20180722','20180723','20180724',\
#run_list = ['20180810','20180812','20180813',\
#            '20180814','20180815','20180816','20180817',\
#run_list = ['20190807','20190808','20190809','20190810','20190811','20190812','20190813',\
#            '20190814']

shawn_path = '/data1/OMI/shawn_files/ltc3/'

#run_list = days_needed
#
#for ttime in swaths_needed:
#    automate_all_preprocess(ttime, download = False, \
#        images = False, process = True, \
#        omi_dtype = 'ltc3', include_tropomi = False,\
#        copy_to_raindrop = False, \
#        minlat = 70., \
#        remove_empty_scans = True, \
#        reprocess_only_omi = False, \
#        remove_ch2_file = False)
#
#sys.exit()

# TESTING ON 2024/07/31: Setting minlat to 65
# Extended the minimum lat values in the OMI and MODIS files down to ~58 oN
# When set to 70, the min lats in OMI and MODIS files are ~61.5. 
# Setting back to 70 to reduce the sizes of coloc files

if(l_download_shawn):
    for ddate in run_list:
        dt_date_str = datetime.strptime(ddate, '%Y%m%d')
        shawn_files = glob(dt_date_str.strftime(shawn_path + '%Y%m%d*'))
        if(len(shawn_files) < 9):
        #if(len(shawn_files) == 0):
            cmnd = dt_date_str.strftime(\
                'scp -r bsorenson@134.129.222.68:/Research/OMI/'+\
                'out_files-ltc3/new/%Y%m%d* ') + shawn_path
            print(cmnd)
            os.system(cmnd)
    
            shawn_files = glob(dt_date_str.strftime(shawn_path + '%Y%m%d*'))

    final_list = entire_wrapper(min_AI = -2.0, minlat = 70., download = True, \
        images = False, process = False, run_list = run_list, copy_to_raindrop = False, \
        #include_tropomi = False, remove_ch2_file = False, reprocess_only_omi = False)
        include_tropomi = False, remove_ch2_file = False, reprocess_only_omi = False)
        #skiptimes = \
        #['200607270736','200607280958','200607281951','200607291359','200607301800',\
        #'200804221524','200804231607','201408140931','201506250915','201507071933',\
        #'201507120819','201507120501','201507121812','201507132352','201708211840',\
        #'201708151103','201708160153','201708162139','201708160828','201708160014',\
        #'201708161504','201708172231','201708180140','201708191713','201708190720','201708192349',\
        #'201708190402','201708211840','201807190914','201808110601','201808110919'])

else:
    final_list = entire_wrapper(min_AI = -2.0, minlat = 70., download = False, \
        images = False, process = True, run_list = run_list, copy_to_raindrop = False, \
        include_tropomi = False, remove_ch2_file = False, reprocess_only_omi = False, \
        skiptimes = \
        ['201807120054','201507020424','200607201000']
        )
        #['200607220359','200607230800','200607231436','200607232250','200607250927',\
        # '200607251920','200607270736','200607280004','200607280958','200607281951',\
        # '200607280819','200607281633','200607282309','200607291359','200607301800',\
        #'200804221524','200804231607','201408140931','201506250915','201507071933',\
        #'201507120819','201507120501','201507121812','201507132352','201708211840',\
        #'201708151103','201708160153','201708162139','201708160828','201708160014',\
        #'201708161504','201708180140','201708191713','201708190720','201708192349',\
        #'201708190402','201708211840','201807190914','201808110601','201808110919',\
        #'201808151708','201908081949','201908130835','201908130517','201908131510'])
        ##include_tropomi = False, remove_ch2_file = False, reprocess_only_omi = False)

sys.exit()


base_file = '/home/bsorenson/Research/Arctic_compares/comp_data/with_cldpres/colocated_subset_201807052213.hdf5'
calc_data = 'test_calc_out.hdf5'

in_base = h5py.File(base_file)
in_calc = h5py.File(calc_data)

mask_orig = np.ma.masked_where((in_base['ceres_swf'][:,:] == -999.) | \
                               (in_base['ceres_swf'][:,:] > 3000), in_base['ceres_swf'][:,:])
mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])

diff_calc = mask_orig - mask_calc

fig = plt.figure(figsize = (12, 4))
ax1 = fig.add_subplot(1,3,1, projection = ccrs.NorthPolarStereo())
ax2 = fig.add_subplot(1,3,2, projection = ccrs.NorthPolarStereo())
ax3 = fig.add_subplot(1,3,3, projection = ccrs.NorthPolarStereo())

mesh = ax1.pcolormesh(in_base['omi_lon'][:,:], in_base['omi_lat'][:,:], mask_orig, \
    transform = ccrs.PlateCarree(), shading = 'auto')
cbar = fig.colorbar(mesh, ax = ax1, label = 'Original SWF')
ax1.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
ax1.coastlines()


mesh = ax2.pcolormesh(in_base['omi_lon'][:,:], in_base['omi_lat'][:,:], mask_calc, \
    transform = ccrs.PlateCarree(), shading = 'auto')
cbar = fig.colorbar(mesh, ax = ax2, label = 'Calculated SWF')
ax2.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
ax2.coastlines()

mesh = ax3.pcolormesh(in_base['omi_lon'][:,:], in_base['omi_lat'][:,:], diff_calc, \
    transform = ccrs.PlateCarree(), shading = 'auto', vmin = -40, vmax = 40, cmap = 'bwr')
cbar = fig.colorbar(mesh, ax = ax3, label = 'ΔSWF')
ax3.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
ax3.coastlines()

in_base.close()
in_calc.close()

fig.tight_layout()
plt.show()

sys.exit()

"""

orig_path = '/home/bsorenson/Research/Arctic_compares/comp_data/'
new_path  = '/home/bsorenson/Research/Arctic_compares/comp_data/with_cldpres/'

date_str = '201807041633'

var = 'modis_cld'

data1 = h5py.File(orig_path + 'colocated_subset_' + date_str + '.hdf5')
data2 = h5py.File(new_path + 'colocated_subset_' + date_str + '.hdf5')

mask_ceres1 = np.ma.masked_where( (data1[var][:,:] == -999.) | (data1[var][:,:] > 3000), data1[var][:,:])
mask_ceres2 = np.ma.masked_where( (data2[var][:,:] == -999.) | (data2[var][:,:] > 3000), data2[var][:,:])
mask_ceres3 = np.ma.masked_where( (data2['modis_cld_top_pres'][:,:] == -999.) | \
    (data2['modis_cld_top_pres'][:,:] > 3000), data2['modis_cld_top_pres'][:,:])

print(np.max(mask_ceres3), np.min(mask_ceres3))

fig = plt.figure()
mapcrs = ccrs.NorthPolarStereo()
datacrs = ccrs.PlateCarree()
ax1 = fig.add_subplot(1,3,1, projection = mapcrs)
ax2 = fig.add_subplot(1,3,2, projection = mapcrs)
ax3 = fig.add_subplot(1,3,3, projection = mapcrs)

ax1.pcolormesh(data1['omi_lon'][:,:], data1['omi_lat'][:,:], mask_ceres1, shading = 'auto', transform = datacrs)
ax2.pcolormesh(data2['omi_lon'][:,:], data2['omi_lat'][:,:], mask_ceres2, shading = 'auto', transform = datacrs)
ax3.pcolormesh(data2['omi_lon'][:,:], data2['omi_lat'][:,:], mask_ceres3, shading = 'auto', transform = datacrs)
ax1.coastlines()
ax2.coastlines()
ax3.coastlines()
ax1.set_extent([-180,180,65,90], datacrs)
ax2.set_extent([-180,180,65,90], datacrs)
ax3.set_extent([-180,180,65,90], datacrs)

data1.close()
data2.close()
fig.tight_layout()
plt.show()

sys.exit()
"""

date_strs = [\
        '201807040005',\
        '201807040144',\
        '201807040322',\
        '201807041633',\
        '201807041812',\
        '201807041951',\
        '201807042130',\
        '201807042309',\
        '201807050048',\
        '201807050227',\
        '201807051538',\
        '201807051717',\
        '201807051856',\
        '201807052034',\
        '201807052213',\
        '201807052352',\
        '201807210047',\
        '201807211358',\
        '201807211537',\
        '201807211716',\
        '201807211855',\
        '201807212034'
    ]

automate_all_preprocess(date_strs, download = False, images = False, process = True,\
    omi_dtype = 'ltc3', minlat = 70., copy_to_raindrop = False, include_tropomi = False)

sys.exit()



"""
data = h5py.File('force_effic_values.h5')

fig = plt.figure()
ax1 = fig.add_subplot(1,4,1)
ax2 = fig.add_subplot(1,4,2)
ax3 = fig.add_subplot(1,4,3)
ax4 = fig.add_subplot(1,4,4)

ax1.pcolormesh(data['regress_slopes'][0,:,1,:].T, cmap = 'bwr', shading = 'auto', vmin = -30, vmax = 30)
ax1.set_title('Ocean')
ax2.pcolormesh(data['regress_slopes'][3,:,1,:].T, cmap = 'bwr', shading = 'auto', vmin = -30, vmax = 30)
ax2.set_title('Mix')
ax3.pcolormesh(data['regress_slopes'][1,:,1,:].T, cmap = 'bwr', shading = 'auto', vmin = -30, vmax = 30)
ax3.set_title('Ice')
ax4.pcolormesh(data['regress_slopes'][2,:,1,:].T, cmap = 'bwr', shading = 'auto', vmin = -30, vmax = 30)
ax4.set_title('Land')

data.close()

fig.tight_layout()
plt.show()

sys.exit()

#data = pd.read_csv('test_out_file_nocod.txt') 
data = pd.read_csv('test_out_file.txt') 

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.scatter(data['OMI'], data['CERES'], s = 6, color = 'k')
plot_trend_line(ax, data['OMI'], data['CERES'], color='tab:red', linestyle = '-', \
    slope = 'linregress')
ax.set_xlabel('OMI AI')
ax.set_ylabel('CERES SWF')
plt.show()

sys.exit()

data = h5py.File('force_effic_values.h5')

fig = plt.figure(figsize = (10, 4))
ax1 = fig.add_subplot(1,4,1)
ax2 = fig.add_subplot(1,4,2)
ax3 = fig.add_subplot(1,4,3)
ax4 = fig.add_subplot(1,4,4)

ax1.plot(-data['regress_slopes'][0,:,1], label = 'cloud')
ax1.plot(-data['regress_slopes'][0,:,0], label = 'clear')
ax1.axhline(0, color = 'grey', linestyle = '--')

ax2.plot(-data['regress_slopes'][3,:,1], label = 'cloud')
ax2.plot(-data['regress_slopes'][3,:,0], label = 'clear')
ax2.axhline(0, color = 'grey', linestyle = '--')

ax3.plot(-data['regress_slopes'][1,:,1], label = 'cloud')
ax3.plot(-data['regress_slopes'][1,:,0], label = 'clear')
ax3.axhline(0, color = 'grey', linestyle = '--')

ax4.plot(-data['regress_slopes'][2,:,1], label = 'cloud')
ax4.plot(-data['regress_slopes'][2,:,0], label = 'clear')
ax4.axhline(0, color = 'grey', linestyle = '--')

fig.tight_layout()
plt.show()

sys.exit()
"""


##date_str = '201908110033'
#date_str = '201807041951'
##date_str = '201807042130'
##date_str = '201807042309'
#date_str = '201807052213'
#plot_compare_OMI_CERES_MODIS_NSIDC(date_str, 7, \
#    omi_dtype = 'ltc3', minlat = 65., zoom = True, save = True)
#sys.exit()

"""
date_str = '201807052213'
plot_compare_OMI_MODIS_v2(date_str, 7, \
    omi_dtype = 'ltc3', minlat = 65., zoom = True, save = False)
sys.exit()
"""

"""
date_str = '201408112211'
#date_str = '201807052213'
plot_compare_OMI_MODIS_NSIDC_v2(date_str, 'true_color', \
    omi_dtype = 'ltc3', minlat = 65., zoom = True, save = False)
sys.exit()
"""

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Steps for daily Arctic comp analysis
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

begin_date = '200504'
end_date   = '202009'
##!#season     = 'sunlight'
##!#minlat = 70.
##!#maxlat = 87.
##!##NSIDC_data = readNSIDC_monthly_grid_all(begin_date, end_date, \
##!##    season, calc_month = True, minlat = minlat, maxlat = maxlat)
##!#
##!## HERE: Calculate the daily-averaged monthly averages, can use as a substitute
##!##       for the old monthly data
##!#shawn_file = home_dir + '/Research/OMI/omi_shawn_daily_2005_2020.hdf5'
##!#jz_file    = home_dir + '/Research/OMI/omi_VJZ211_daily_2005_2020.hdf5'
##!##daily_VSJ4 = h5py.File(shawn_file, 'r')
##!#daily_VSJ4 = readOMI_daily_HDF5(shawn_file, minlat = minlat, maxlat = maxlat)
##!##OMI_daily_VSJ4  = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -0.10, minlat = minlat, maxlat = maxlat)
##!##OMI_daily_VJZ211 = calcOMI_MonthAvg_FromDaily(jz_file, min_AI = -0.10, minlat = minlat, maxlat = maxlat)
##!#
##!## BEGIN SINGLE-DAY FUNCTION HERE
##!#
##!#date_str = '20180704'
##!##tidx = 100
##!#
##!## Load in the single-day MODIS cloud fraction
##!#MYD08_data = read_MODIS_MYD08_single(date_str, minlat = minlat, \
##!#    maxlat = maxlat)
##!#
##!## Load in the single-day NSIDC ice concentration
##!#NSIDC_data =  readNSIDC_daily(date_str, maxlat = maxlat)
##!#NSIDC_data = grid_data_conc(NSIDC_data, minlat = minlat, maxlat = maxlat)
##!#
##!## Figure out which OMI day index matches the date_str
##!#file_strs = np.array([str(tval) for tval in daily_VSJ4['day_values']])
##!#match_idx = np.where(date_str == file_strs)[0][0]
##!#
##!#local_shawn = np.ma.masked_where(daily_VSJ4['count_AI'][match_idx,:,:] == 0, \
##!#    daily_VSJ4['grid_AI'][match_idx,:,:])
##!#
##!#min_idx = np.where(MYD08_data['lat'][:] >= minlat)[0][0]
##!#max_idx = np.where(MYD08_data['lat'][:] <= maxlat)[0][-1] + 1
##!#
##!#mask_ice = np.ma.masked_where(NSIDC_data['grid_ice_conc'][:,:] < 0, \
##!#    NSIDC_data['grid_ice_conc'][:,:])
##!#
##!#
##!#
##!#
##!#
##!#fig = plt.figure()
##!#ax1 = fig.add_subplot(2,2,1, projection = mapcrs)
##!#ax2 = fig.add_subplot(2,2,2, projection = mapcrs)
##!#ax3 = fig.add_subplot(2,2,3, projection = mapcrs)
##!#
##!#ax1.pcolormesh(daily_VSJ4['lon_values'][:], daily_VSJ4['lat_values'][:], \
##!#    local_shawn, transform = datacrs, \
##!#    shading = 'auto', cmap = 'jet')
##!#ax1.coastlines()
##!#ax1.set_extent([-180,180,65,90], datacrs)
##!#
##!#ax2.pcolormesh(MYD08_data['lon'][:], MYD08_data['lat'][min_idx:max_idx], \
##!#    MYD08_data['cld_frac_mean'][:,:], transform = datacrs, \
##!#    shading = 'auto')
##!#ax2.coastlines()
##!#ax2.set_extent([-180,180,65,90], datacrs)
##!#
##!#ax3.pcolormesh(NSIDC_data['grid_lon'][:,:], NSIDC_data['grid_lat'][:,:], \
##!#    mask_ice, transform = datacrs, \
##!#    shading = 'auto')
##!#ax3.coastlines()
##!#ax3.set_extent([-180,180,65,90], datacrs)
##!#
##!#plt.suptitle(file_strs[match_idx])
##!#
##!#fig.tight_layout()
##!#plt.show()
##!#
##!#sys.exit()



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# End daily Arctic comp analysis
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Read in the OMI data
# ---------------------

#testfile = 'grid_coloc_test_res050.hdf5'
testfile = 'grid_coloc_test_res100.hdf5'
data = h5py.File(testfile)

#idx_dict, lats, lons = match_aeronet_to_grid_AI(data, aeronet_file = 'aeronet_site_info.txt', \
#    min_ai = 1.5)
#
#sys.exit()


#filename = 'comp_grid_climo_v1.hdf5'
#filename = 'comp_grid_climo_v2.hdf5'
#filename = 'comp_grid_climo_v3.hdf5'
#filename = 'comp_grid_climo_v4.hdf5'
#filename = 'comp_grid_climo_v5.hdf5'
filename1 = 'comp_grid_climo_v6.hdf5'
#filename = 'comp_grid_climo_v7.hdf5'
#filename1 = 'comp_grid_climo_v8.hdf5'
#filename1 = 'comp_grid_climo_v9.hdf5'
#filename1 = 'comp_grid_climo_v10.hdf5'
#filename2 = 'comp_grid_climo_v11.hdf5'
filename2 = 'comp_grid_climo_v12.hdf5'
filename3 = 'comp_grid_climo_v14.hdf5'
comp_grid_data_v6 = read_comp_grid_climo(filename1)
comp_grid_data_v11 = read_comp_grid_climo(filename2)
comp_grid_data_v14 = read_comp_grid_climo(filename3)


#comp_grid_data_v7 = read_comp_grid_climo(filename)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# NOTE: Beginning with comp_grid_climo version 10, the new comp_grid_climo
#       versions use 'modis_cod' as a binning variable. Since the
#       comp_grid_climo dictionary is used for binning the raw data, 
#       use of the new comp_grid_climo versions causes the raw data binning
#       to be switched to 'modis_cod'. 
#
#       Thus, to go back to the old way of binning by 'modis_ch7' rather
#       than 'modis_cod', run this code using version 6.
#
#       As a note, reworking this code to specifically use separate
#       bins (CH7 vs COD) apart from relying on the comp_grid_climo
#       dictionary would simplify this process. 
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#files = glob(home_dir + \
#    '/Research/Arctic_compares/comp_data/colocated_subset_20*.hdf5')

data_path = home_dir + \
    '/Research/Arctic_compares/comp_data/colocated_subset_'

# NOTE: These dates correspond to the dates that are in comp_grid_climo_v4
dates = [
        '200607240029',\
        '200607240208',\
        '200607240347',\
        '200607240526',\
        '200607242016',\
        '200607242155',\
        '200607242334',\
        '200607250112',\
        '200607250251',\
        '200607250430',\
        '200607252238',\
        '200607260017',\
        '200607260156',\
        '200607260335',\
        '200607260513',\
        '200607260652',\
        '200607260831',\
        '200607262003',\
        '200607262142',\
        '200607262321',\
        '200607270100',\
        '200607270239',\
        '200607270418',\
        '200607270557',\
        '200607270736',\
        '200607270914',\
        '200607272047',\
        '200607272226',\
        '200804221841',\
        '200804222159',\
        '201408110046',\
        '201408110404',\
        '201408111853',\
        '201408112032',\
        '201408112211',\
        '201408120308',\
        '201408121758',\
        '201408121937',\
        '201408122115',\
        '201408122254',\
        '201506271220',\
        '201506271359',\
        '201506271538',\
        '201506271717',\
        '201506271856',\
        '201507061353',\
        '201507061532',\
        '201507061711',\
        '201507061850',\
        '201507062028',\
        '201507062207',\
        '201507071615',\
        '201507071754',\
        '201507071933',\
        '201507081837',\
        '201507082016',\
        '201507082155',\
        '201507090113',\
        '201507090748',\
        '201507090927',\
        '201507091245',\
        '201507091424',\
        '201507091603',\
        '201507100514',\
        '201507100653',\
        '201708160014',\
        '201708161146',\
        '201708161325',\
        '201708161504',\
        '201708161643',\
        '201708161821',\
        '201708162000',\
        '201708162139',\
        '201708171050',\
        '201708171229',\
        '201708171547',\
        '201708171726',\
        '201708172222',\
        '201708181133',\
        '201708181312',\
        '201708181451',\
        '201708181630',\
        '201708181809',\
        '201708181948',\
        '201708191038',\
        '201708191217',\
        '201708191355',\
        '201708191534',\
        '201708191713',\
        '201708191852',\
        '201807040005',\
        '201807040144',\
        '201807040322',\
        '201807041633',\
        '201807041812',\
        '201807041951',\
        '201807042130',\
        '201807042309',\
        '201807050048',\
        '201807050227',\
        '201807051538',\
        '201807051717',\
        '201807051856',\
        '201807052034',\
        '201807052213',\
        '201807052352',\
        '201807210047',\
        '201807211358',\
        '201807211537',\
        '201807211716',\
        '201807211855',\
        '201807212034',\
        '201808100200',\
        '201808140135',\
        '201808141804',\
        '201808141942',\
        '201808142121',\
        '201808142300',\
        '201808260158',\
        '201808260337',\
        '201808260655',\
        '201908100129',\
        '201908100308',\
        '201908101936',\
        '201908102115',\
        '201908102254',\
        '201908110033',\
        '201908110212',\
        '201908110351',\
        '201908110708',\
        '201908111523',\
        '201908111702',\
        '201908111841',\
        '201908112019',\
        '201908112158',\
        '201908112337',\
        #'202107040232',\
        #'202107041722',\
        #'202107051627',\
        #'202107112225',\
        #'202107121454',\
        #'202107121633',\
        #'202107121812',\
        #'202107121950',\
        #'202107122129',\
        #'202107122308',\
        #'202107131219',\
        #'202107131537',\
        #'202107290046',\
        #'202107290225',\
        #'202107292033',\
        #'202107292212',\
        #'202107292351',\
        #'202107300129',\
        #'202107300308',\
        #'202107300447',\
        #'202107300626',\
        #'202107301937',\
        #'202107302116',\
        #'202107302255',\
        # NOTE: All July 2021 times are added for comp_grid_climo_v8
        #'202108010117',\
        #'202108010256',\
        #'202108010435',\
        #'202108010614',\
        #'202108011607',\
        #'202108011925',\
        #'202108012103',\
        #'202108012242',\
    ]

## NOTE: These dates correspond to the dates that are in comp_grid_climo_v3
#dates = [
#         '200607240029',
#         '200607240208',
#         '200607240347',
#         '200607240526',
#         '200607240705',
#         '200607240844',
#         '200607242016',
#         '200607242155',
#         '200607242334',
#         '200607250112',
#         '200607250251',
#         '200607250430',
#         '200607250609',
#         '200607250748',
#         '200607252238',
#         '200607270100',
#         '200607270239',
#         '200607270418',
#         '200607270557',
#         '200607270736',
#         '200607270914',
#         '200607272047',
#         '200607272226',
#         '200804221841',
#         '200804222020',
#         '200804222159',
#         '201408110046',
#         '201408110404',
#         '201408111853',
#         '201408112032',
#         '201408112211',
#         '201408112350',
#         '201408120129',
#         '201408120308',
#         '201408121758',
#         '201408121937',
#         '201408122115',
#         '201408122254',
#         '201506271042',
#         '201506271220',
#         '201506271359',
#         '201506271538',
#         '201506271717',
#         '201506271856',
#         '201506272214',
#         '201507061035',
#         '201507061353',
#         '201507061532',
#         '201507061711',
#         '201507061850',
#         '201507062028',
#         '201507062207',
#         '201507062346',
#         '201507070125',
#         '201507071118',
#         '201507071257',
#         '201507071436',
#         '201507071615',
#         '201507071754',
#         '201507071933',
#         '201507072112',
#         '201507080347',
#         '201507080526',
#         '201507080705',
#         '201507081023',
#         '201507081202',
#         '201507081340',
#         '201507081519',
#         '201507081658',
#         '201507081837',
#         '201507082016',
#         '201507082155',
#         '201507082334',
#         '201507090113',
#         '201507090252',
#         '201507090430',
#         '201507090609',
#         '201507090748',
#         '201507090927',
#         '201507091106',
#         '201507091245',
#         '201507091424',
#         '201507091603',
#         '201507091741',
#         '201507091920',
#         '201507092059',
#         '201507092238',
#         '201507100017',
#         '201507100156',
#         '201507100335',
#         '201507100514',
#         '201507100653',
#         '201507101010',
#         '201507101328',
#         '201708160014',
#         '201708161146',
#         '201708161325',
#         '201708161504',
#         '201708161643',
#         '201708161821',
#         '201708162000',
#         '201708162139',
#         '201807040005',
#         '201807040144',
#         '201807040322',
#         '201807041633',
#         '201807041812',
#         '201807041951',
#         '201807042130',
#         '201807042309',
#         '201807050048',
#         '201807050227',
#         '201807051538',
#         '201807051717',
#         '201807051856',
#         '201807052034',
#         '201807052213',
#         '201807052352',
#         '201807210047',
#         '201807211358',
#         '201807211537',
#         '201807211716',
#         '201807211855',
#         '201807212034',
#         '201808100200',
#         '201808140135',
#         '201808141804',
#         '201808141942',
#         '201808142121',
#         '201808142300',
#         '201808260158',
#         '201808260337',
#         '201808260655',
#         '201908100129',
#         '201908100308',
#         '201908101936',
#         '201908102115',
#         '201908102254',
#         '201908110033',
#         '201908110212',
#         '201908110351',
#         '201908110708',
#         '201908111523',
#         '201908111702',
#         '201908111841',
#         '201908112019',
#         '201908112158',
#         '201908112337',
#         '202108010117',
#         '202108010256',
#         '202108010435',
#         '202108010614',
#         '202108011607',
#         '202108011925',
#         '202108012103',
#         '202108012242',
#        ]

#calc_pcnt_aerosol_over_type(dates, 1.5, minlat = 70., dtype = 'PERT', ax = None)
#sys.exit()
#plot_aerosol_over_types(dates[102], min_AI = 2.0, ai_val = 'TROP_AI', save = False)
#plot_aerosol_over_types(dates[125], min_AI = 2.0, ai_val = 'TROP_AI', save = False)


#plot_aerosol_over_type_combined(data, dates, min_ai = 1.5, save = False, plot_map = True)
#sys.exit()

##!#fig = plt.figure(figsize = (9, 6))
##!#ax1 = fig.add_subplot(2,1,1)
##!#ax2 = fig.add_subplot(2,1,2)
##!##ax3 = fig.add_subplot(3,1,3)
##!#calc_pcnt_aerosol_over_type(dates, 1.5, ax = ax1)
##!#calc_pcnt_aerosol_over_type_dayavgs(data, 1.5, ax = ax2, area_calc = True, hatch_cloud = True, plot_map = False)
##!##calc_pcnt_aerosol_over_type_dayavgs(data, 1.5, ax = ax3, area_calc = True, hatch_cloud = True)
##!#
##!#ax1.set_ylabel('Pcnt of OMI Pixels')
##!##ax2.set_ylabel('Pcnt of 0.5 deg. Grid Boxes')
##!#ax2.set_ylabel('Pcnt oF Aerosol Area')
##!#
##!#fig.tight_layout()
##!#
##!#plt.show()
#sys.exit()

files = [data_path + date + '.hdf5' for date in dates]

# Figure out the total size to insert the data
# ---------------------------------------------
#minlat = 70.
minlat = 65.    # 2024/01/10: changed to 65.
total_size = 0
for ff in files:
    data = h5py.File(ff,'r')
    local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
    local_data = np.ma.masked_where(abs(local_data) > 12, local_data)
    local_data = np.ma.masked_where(data['omi_lat'][:,:] < minlat, \
        local_data) 
    local_data = np.ma.masked_where((data['ceres_swf'][:,:] < -200.) | \
        (data['ceres_swf'][:,:] > 3000), \
        local_data) 
    local_size = local_data.compressed().shape[0]
    print(ff, local_size)
    total_size += local_size

    data.close()


# Set up the data structure to hold all the data
# ----------------------------------------------
combined_data = {}
combined_data['omi_uvai_pert'] = np.full(total_size, np.nan)
combined_data['omi_uvai_raw']  = np.full(total_size, np.nan)
combined_data['modis_cld']     = np.full(total_size, np.nan)
combined_data['modis_cod']     = np.full(total_size, np.nan)
combined_data['ceres_swf']     = np.full(total_size, np.nan)
combined_data['modis_ch7']     = np.full(total_size, np.nan)
combined_data['omi_sza']       = np.full(total_size, np.nan)
combined_data['omi_lat']       = np.full(total_size, np.nan)
combined_data['nsidc_ice']     = np.full(total_size, np.nan)

print("Loading data")

# Loop back over the files and insert the data into the structure
# ---------------------------------------------------------------
total_size = 0
beg_idx = 0
end_idx = 0
for ff in files:

    data = h5py.File(ff,'r')
    # NOTE: Changed the omi variable here from "pert" to "raw" on 20230623.
    #       This move should allow for coloc data to be read after 2020
    local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
    local_data = np.ma.masked_where(abs(local_data) > 12, local_data)
    local_data = np.ma.masked_where(data['omi_lat'][:,:] < minlat, \
        local_data) 
    local_data = np.ma.masked_where((data['ceres_swf'][:,:] < -200.) | \
        (data['ceres_swf'][:,:] > 3000), \
        local_data) 
    local_size = local_data.compressed().shape[0]

    beg_idx = end_idx
    end_idx = beg_idx + local_size

    for tkey in combined_data.keys():
        combined_data[tkey][beg_idx:end_idx] = \
            data[tkey][~local_data.mask]

    print(local_size)
    total_size += local_size

    data.close()


print(np.min(combined_data['omi_uvai_pert']), np.max(combined_data['omi_uvai_pert']))
print(np.min(combined_data['ceres_swf']), np.max(combined_data['ceres_swf']))
print(np.min(combined_data['modis_cod']), np.max(combined_data['modis_cod']))
print(np.min(combined_data['omi_sza']), np.max(combined_data['omi_sza']))
print(np.min(combined_data['nsidc_ice']), np.max(combined_data['nsidc_ice']))

combined_data['nsidc_ice'][:] = \
    np.where(combined_data['nsidc_ice'][:] == 254., 101., combined_data['nsidc_ice'][:])

min_max_dict = {}

key_variables = ['omi_uvai_pert', 'omi_sza', 'modis_cod', 'nsidc_ice']

for key in key_variables:
    min_max_dict[key] = {}
    min_max_dict[key]['min'] = np.min(combined_data[key])
    min_max_dict[key]['max'] = np.max(combined_data[key])

    drange = min_max_dict[key]['max'] - min_max_dict[key]['min']
    combined_data[key] = (combined_data[key] - min_max_dict[key]['min']) / drange

print(np.min(combined_data['omi_uvai_pert']), np.max(combined_data['omi_uvai_pert']))
print(np.min(combined_data['ceres_swf']), np.max(combined_data['ceres_swf']))
print(np.min(combined_data['modis_cod']), np.max(combined_data['modis_cod']))
print(np.min(combined_data['omi_sza']), np.max(combined_data['omi_sza']))
print(np.min(combined_data['nsidc_ice']), np.max(combined_data['nsidc_ice']))

#num_test = 100
pcnt_test = 0.25
num_test = int(combined_data['omi_uvai_pert'].shape[0] * pcnt_test)
num_train = combined_data['omi_uvai_pert'].shape[0] - num_test
ranges = np.arange(0, combined_data['omi_uvai_pert'].shape[0])

train_idxs, test_idxs = train_test_split(ranges, test_size = num_test)

#test_idxs = random.sample(range(0,combined_data['omi_uvai_pert'].shape[0]) , num_test)
#tf_no_idxs = [num not in test_idxs for num in ranges]
#train_idxs = ranges[tf_no_idxs]

# CONTINUE HERE

sys.exit()

#mask_data = np.ma.masked_invalid(OMI_base['UVAI_raw'])
#mask_data = np.ma.masked_invalid(data['uvai_raw'])
#
#mask_dims = np.array([ (False in mask_data[ii,:].mask) \
#    for ii in range(mask_data.shape[0])])
#
#keep_idxs = np.where(mask_dims == True)[0]
#
#dset.create_dataset('VARIABLE', data = data[key][keep_idxs,:])

sys.exit()

ai_min  = 2
ai_max  = None
sza_min = 50
sza_max = 55
ice_min = None
ice_max = None
ch7_min = None
ch7_max = None
cld_min = None
cld_max = None

#trend_type = 'theil-sen'
trend_type = 'linregress'
##!#lin_raw_dict_v6 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v6, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        )
##!#lin_smth_dict_v6 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v6, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smooth', sizer = 1)
lin_smth2_dict_v6 = calc_raw_grid_slopes(\
        combined_data, comp_grid_data_v6, \
        ai_min = ai_min, ai_max = ai_max, \
        trend_type = trend_type, \
        # 2024-01-10: decided to use the perturbed OMI data, since that is 
        #             what's used in the daily averages.
        xval = 'omi_uvai_pert', \
        smoother = 'smoother', sizer = 1)
lin_smth2_dict_v6['minlat'] = minlat
##!#lin_raw_dict_v11 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v11, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type)
##!#lin_smth_dict_v11 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v11, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smooth', sizer = 1)
##!#lin_smth2_dict_v11 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v11, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smoother', sizer = 1)
##!#lin_raw_dict_v14 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v14, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type)
##!#lin_smth_dict_v14 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v14, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smooth', sizer = 1)
##!#lin_smth2_dict_v14 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v14, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smoother', sizer = 1)
##!#
##!#trend_type = 'theil-sen'
##!#thl_raw_dict_v6 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v6, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type)
##!#thl_smth_dict_v6 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v6, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smooth', sizer = 1)
##!#thl_smth2_dict_v6 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v6, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smoother', sizer = 1)
##!#thl_raw_dict_v11 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v11, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type)
##!#thl_smth_dict_v11 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v11, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smooth', sizer = 1)
##!#thl_smth2_dict_v11 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v11, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smoother', sizer = 1)
##!#thl_raw_dict_v14 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v14, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type)
##!#thl_smth_dict_v14 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v14, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smooth', sizer = 1)
##!#thl_smth2_dict_v14 = calc_raw_grid_slopes(\
##!#        combined_data, comp_grid_data_v14, \
##!#        ai_min = ai_min, ai_max = ai_max, \
##!#        trend_type = trend_type, \
##!#        smoother = 'smoother', sizer = 1)


#return_dict = \
#    plot_compare_slopes_scatter(thl_raw_dict_v14, combined_data, comp_grid_data_v14, \
#    5, 3, dtype = 'raw', ice_idx = 0, ai_min = 2, \
#    ai_max = None, show_trend = False, save = False)


min_cloud = 0.95
maxerr = 2
data_type = 'raw'
# data_type: 'raw' or 'grid'
ocean_slopes = calc_slope_clear_clean_sfctype(lin_smth2_dict_v6, 0, 0, \
    maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)
ice_slopes   = calc_slope_clear_clean_sfctype(lin_smth2_dict_v6, 1, 0, \
    maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)
land_slopes  = calc_slope_clear_clean_sfctype(lin_smth2_dict_v6, 2, 0, \
    maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)

combined_slope_dict = {
    'ocean': ocean_slopes, \
    'ice': ice_slopes, \
    'land': land_slopes
}

# Calculate solar declination angles
# ----------------------------------

begin_date = '200504'
end_date   = '202009'
season     = 'sunlight'
#minlat = 65.
maxlat = 87.
NSIDC_data = readNSIDC_monthly_grid_all(begin_date, end_date, \
    season, calc_month = True, minlat = minlat, maxlat = maxlat)

# HERE: Calculate the daily-averaged monthly averages, can use as a substitute
#       for the old monthly data
#shawn_file = home_dir + '/Research/OMI/omi_shawn_daily_2005_2020.hdf5'
shawn_file = home_dir + '/Research/OMI/omi_shawn_daily_2005_2020_v2.hdf5'
jz_file    = home_dir + '/Research/OMI/omi_VJZ211_daily_2005_2020.hdf5'
OMI_daily_VSJ4  = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -0.10, minlat = minlat, maxlat = maxlat)
OMI_daily_VJZ211 = calcOMI_MonthAvg_FromDaily(jz_file, min_AI = -0.10, minlat = minlat, maxlat = maxlat)

daily_VSJ4 = readOMI_daily_HDF5(shawn_file, minlat = minlat, maxlat = maxlat)

OMI_VSJ4   = readOMI_NCDF(infile = \
    '/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', \
    minlat = minlat, maxlat = maxlat - 1)
OMI_VJZ211 = readOMI_NCDF(infile = \
    '/home/bsorenson/Research/OMI/omi_ai_VJZ211_2005_2020.nc', \
    minlat = minlat, maxlat = maxlat - 1)

MYD08_data = read_MODIS_MYD08_monthrange(begin_date,end_date,\
    minlat = minlat, maxlat = maxlat, calc_month = False)

#plot_test_forcing_v2(OMI_VSJ4, NSIDC_data, MYD08_data, lin_smth2_dict_v6, \
#    tidx, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
#    cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
#    save = False)

##!#def plot_test_data(OMI_data, tidx, min_ai):
##!#    local_data = np.ma.masked_where(OMI_data['AI'][tidx,:,:] < min_ai, OMI_data['AI'][tidx,:,:])
##!#
##!#    fig = plt.figure(figsize = (7, 3))
##!#    ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
##!#    ax2 = fig.add_subplot(1,2,2)
##!#    mesh = ax1.pcolormesh(OMI_data['LON'], OMI_data['LAT'], local_data, \
##!#        transform = datacrs, shading = 'auto', cmap = 'jet', vmin = -0.75, vmax = 1.)
##!#    cbar = fig.colorbar(mesh, ax = ax1)
##!#    ax1.coastlines()
##!#    ax1.set_extent([-180,180,65,90], datacrs)
##!#    ax2.hist(local_data.flatten(), bins = 'auto')
##!#    plt.suptitle(OMI_data['DATES'][tidx])
##!#    fig.tight_layout()
##!#    plt.show()

     

##!## This uses method 1: trend
##!## -------------------------
##!#plot_type_forcing_all_months(OMI_data, NSIDC_data, 'clear', minlat = minlat, \
##!#    maxlat = maxlat, use_szas = False, save = False, coloc_dict = lin_smth2_dict_v6)

# This uses method 2: individual month-based
# ------------------------------------------

ai_thresh = 0.05
#ai_thresh = -0.15

#plot_test_forcing_v2(OMI_daily_VSJ4, NSIDC_data, MYD08_data, lin_smth2_dict_v6, \
#    81, minlat = 70., maxlat = 87., ai_thresh = ai_thresh, \
#    cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
#    save = False)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#                      BEGIN FORCING VERSION 3 STUFF
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#                             ERROR ANALYSIS
#
# First method:
#     Daily-gridded monthly averages of forcing are calculated using 
#     the following two scenarios:
#    
#         - '_adderror': Using the standard deviation of the different
#               AI-SWF slopes that are grouped into the "clear" and 
#               "cloudy" types in determining the forcing efficiency,
#               the standard error of the mean SZA-mean forcing efficiency
#               estimate is derived. In the "add error" scenario,
#               the standard error of each mean SZA-mean forcing efficiency
#               is added to the calculated mean SZA-mean forcing efficiency,
#               following:
#       
#                   δSWF/δAI = δSWF/δAI_mean + δSWF/δAI_std_error
#
#               and these new values are then used to calculate the daily
#               forcing estimates, which are then averaged into monthly
#               values.
#         - '_suberror': Same as in '_adderror', but the standard error
#               of the mean SZA-mean forcing efficiency is subtracted
#               from the calculated mean SZA-mean forcing efficiency,
#               following:
#
#                   δSWF/δAI = δSWF/δAI_mean + δSWF/δAI_std_error
#
#               Similarly, these new forcing efficiency values are 
#               used to find monthly forcing estimates. 
#
#     Note that in this simplified error analysis method, the 
#     standard errors of the AI/SWF slopes from each individual
#     bin are NOT accounted for. Also, the error in the MODIS daily
#     cloud fraction is not accounted for either. Both of these
#     must be accounted for for a more robust analysis.
#         
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#                        PLOTTING/ANALYSIS FUNCTIONS
#    
# VERSION 3: Individual forcing values are calculated for each Apr-Sep day
#            from 2005 - 2020, using the daily-gridded shawn OMI data, the
#            subsetted daily MODIS MYD08 cloud product, and the daily
#            NSIDC ice concentration data. 
#
# plot_slopes_cloud_types_szamean: Plots the four-panel forcing
#            efficiency estimates as a function of SZA and averaged
#            along the CH7 bins.
#
# plot_test_forcing_v3: calculates the single-day forcing estimate and
#            plots the estimate for that day.
#
# calculate_type_forcing_v3_monthly: calculates the single-day forcing 
#            estimates and averages them into monthly averages.
#
# plot_test_forcing_v3_monthly: Takes the "all_month_vals" daily-gridded
#           monthy averaged forcing estimates as input and plots a single
#           month of the data.
#
# plot_test_forcing_v3_daily_monthly: Combination of 
#           "plot_test_forcing_v3" and "plot_test_forcing_v3_monthly", 
#           plots the daily-averaged AI, daily estimated forcing values,
#           and daily-averaged single-month forcing value for a YYYYMMDD
#           date string.
#
# plot_test_forcing_v3_all_months: Using the "all_month_vals" daily-gridded
#           monthy averaged forcing estimates, calculates the trends over
#           the monthly averaged forcing estimates for each month range
#           and plots the results. 
#
# plot_type_forcing_v3_all_months_arctic_avg_manyrefice: using each of 
#           the daily-gridded monthly averaged forcing estimates calculated
#           with reference ices from 2005 - 2020, plot each of the 
#           ref-ice "simulations" on one graph. 
#           
# calc_forcing_slopes_v3_all_months_arctic_avg_manyrefice:
#           This function returns the ΔFlux (slope times # years) value
#           if "return_slope = True" is included as an argument to the 
#           function call. There is one slope for each month and for 
#           each year (and by year, this means the "reference ice/cloud"
#           year. Uses the same slope calculations as the above 
#           function.
#
# calc_print_forcing_slope_error_v3(all_month_vals:
#           This function gathers the refice and refcld simulation results,
#           calculates the trends from each refice and refcld simulation,
#           and calculates the mean and standard error of those
#           ref ice/cld simulations. The results are printed in a table

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Plot the raw or grid slopes and the associated slope errors
# -----------------------------------------------------------
#plot_compare_grid_climo_stderr(lin_smth2_dict_v6, 'raw', 2, save = False)

# Plot the 2-d forcing efficiencies
# ----------------------------------
#plot_slopes_cloud_types(lin_raw_dict, save = False, vmin = -15, vmax = 15)

# Plot the four-panel sza-meaned forcing efficiency estimates
# -----------------------------------------------------------
plot_slopes_cloud_types_szamean(lin_smth2_dict_v6,cld_idx = 0, maxerr = 2, \
    data_type = 'raw', remove_high_error = True, hatch_cloud = False, \
    min_cloud = 0.95, save = False)
sys.exit()
# Plot the calculated forcing values for a single day
# ---------------------------------------------------
#plot_test_forcing_v3(daily_VSJ4, OMI_daily_VSJ4, '20080422', \
#    coloc_dict, minlat = 65., maxlat = 87., ai_thresh = 0.7, \
#    cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
#    save = False, filter_bad_vals = False)

# Calculate daily forcing values and average the daily values into
# monthly forcing estimates
# ----------------------------------------------------------------
maxerr = 1.5
ai_thresh = 0.7
all_month_vals = calculate_type_forcing_v3_monthly(daily_VSJ4, \
    OMI_daily_VSJ4, lin_smth2_dict_v6, 'all', ai_thresh = ai_thresh, \
    maxerr = maxerr, minlat = minlat)
#write_daily_month_force_to_HDF5(all_month_vals, OMI_daily_VSJ4, \
#    name_add = '_dayaithresh07_aipert_dataminlat70')
write_daily_month_force_to_HDF5(all_month_vals, OMI_daily_VSJ4, \
    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
    dtype = 'pert', 
    name_add = '_dayaithresh07_v3')

# Calculate daily forcing values and average the daily values into
# monthly forcing estimates, but by adding or subtracting the 
# slope error from the SZA-mean AI/SWF slopes. 
# ----------------------------------------------------------------
all_month_vals_adderror = calculate_type_forcing_v3_monthly(daily_VSJ4, \
    OMI_daily_VSJ4, lin_smth2_dict_v6, 'all', ai_thresh = ai_thresh, minlat = minlat, \
    maxerr = maxerr, mod_slopes = 'add')
write_daily_month_force_to_HDF5(all_month_vals_adderror, OMI_daily_VSJ4, \
    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
    dtype = 'pert', 
    name_add = '_dayaithresh07_v3_adderror')
all_month_vals_suberror = calculate_type_forcing_v3_monthly(daily_VSJ4, \
    OMI_daily_VSJ4, lin_smth2_dict_v6, 'all', ai_thresh = ai_thresh, minlat = minlat, \
    maxerr = maxerr, mod_slopes = 'subtract')
write_daily_month_force_to_HDF5(all_month_vals_suberror, OMI_daily_VSJ4, \
    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
    dtype = 'pert', 
    name_add = '_dayaithresh07_v3_suberror')

# Calculate daily forcing values and average the daily values into
# monthly forcing estimates, but here using the daily ice concentration
# values from 2005 as a reference. Am using this to try to see how
# the change in sea ice affects the aerosol forcing. May need to use
# the average of the first three years of sea ice rather than one
# year...
# ----------------------------------------------------------------
#all_month_vals = calculate_type_forcing_v3_monthly(daily_VSJ4, \
#    OMI_daily_VSJ4, lin_smth2_dict_v6, 'all', ai_thresh = 0.7, minlat = 65., \
#    reference_ice = '2005')

ref_ice_vals = np.arange(2005,2021)
for ref_ice in ref_ice_vals:
    all_month_vals_ice = calculate_type_forcing_v3_monthly(daily_VSJ4, \
        OMI_daily_VSJ4, lin_smth2_dict_v6, 'all', ai_thresh = ai_thresh, minlat = minlat, \
        maxerr = maxerr, reference_ice = str(ref_ice))
    write_daily_month_force_to_HDF5(all_month_vals_ice, OMI_daily_VSJ4, \
        name_add = '_dayaithresh07_v3_refice' + str(ref_ice))

ref_cld_vals = np.arange(2005,2021)
for ref_cld in ref_cld_vals:
    all_month_vals_cld = calculate_type_forcing_v3_monthly(daily_VSJ4, \
        OMI_daily_VSJ4, lin_smth2_dict_v6, 'all', ai_thresh = ai_thresh, minlat = minlat, \
        maxerr = maxerr, reference_cld = str(ref_cld))
    write_daily_month_force_to_HDF5(all_month_vals_cld, OMI_daily_VSJ4, \
        name_add = '_dayaithresh07_v3_refcld' + str(ref_cld))
    

# Write the all_month_vals to an HDF5 file
# ----------------------------------------
#write_daily_month_force_to_HDF5(all_month_vals, OMI_daily_VSJ4, \
#    name_add = '_dayaithresh07')

sys.exit()

# Read the all_month_vals to an HDF5 file
# NOTE: If running any of the following after reading in the all_month_vals
#       from the HDF5 file, must change the 'all_month_vals' to 
#       all_month_dict['FORCE_EST']
#
# NOTE: Unless specified, minlat is 70 for all files
# -------------------------------------------------------------------------
infile = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07.hdf5'
#infile_pert65 = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_aipert_dataminlat65.hdf5'
infile_pert65 = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v3.hdf5'
#infile_min70_pert = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_aipert_dataminlat70.hdf5'
#infile_2005 = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_refice2005.hdf5'
#infile_2020 = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_refice2020.hdf5'
infile_add = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v3_adderror.hdf5'
infile_sub = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v3_suberror.hdf5'
#infile_add = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_adderror.hdf5'
#infile_sub = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_suberror.hdf5'
all_month_dict = read_daily_month_force_HDF5(infile)
#all_month_dict_ref05  = read_daily_month_force_HDF5(infile_2005)
#all_month_dict_pert65  = read_daily_month_force_HDF5(infile_pert65)
#all_month_dict_pert70  = read_daily_month_force_HDF5(infile_min70_pert)
#all_month_dict_ref20  = read_daily_month_force_HDF5(infile_2020)
#all_month_dict_adderr = read_daily_month_force_HDF5(infile_add)
#all_month_dict_suberr = read_daily_month_force_HDF5(infile_sub)

# Plot the daily-gridded monthly forcing values for a single month
# ----------------------------------------------------------------
#plot_test_forcing_v3_monthly(all_month_vals, OMI_daily_VSJ4, '201807', \
#    minlat = 65, save = False)

# Plot the daily gridded OMI AI, daily estimated forcing, and daily-averaged
# monthly forcing value for a single month
# --------------------------------------------------------------------------
#date_str = '20180705'
#plot_test_forcing_v3_daily_monthly(date_str, all_month_dict['FORCE_EST'], \
#    daily_VSJ4, OMI_daily_VSJ4, lin_smth2_dict_v6, \
#    minlat = 65., maxlat = 87., ai_thresh = 0.7, \
#    cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
#    save = False, filter_bad_vals = False)

# Plot the trend of daily-gridded monthly forcing values
# ------------------------------------------------------
#plot_type_forcing_v3_all_months(all_month_vals, OMI_daily_VSJ4, \
#           minlat = 65, omi_data_type = 'pert')

# Calculate the Arctic-wide average of the daily-gridded monthly forcing
# values for each month and plot them
# ----------------------------------------------------------------------
#plot_type_forcing_v3_all_months_arctic_avg(all_month_vals, OMI_daily_VSJ4, \
#    minlat = 65, trend_type = 'standard', omi_data_type = 'pert')
 
# Calculate the Arctic-wide average of the daily-gridded monthly forcing
# values for each month and plot them, but also plotting the '_adderror'
# and '_suberror' results for the first look at an error analysis.
# ----------------------------------------------------------------------
#plot_type_forcing_v3_all_months_arctic_avg(all_month_dict['FORCE_EST'], 
#    OMI_daily_VSJ4, minlat = 65, trend_type = 'standard', 
#    month_values2 = all_month_dict_adderr['FORCE_EST'], 
#    month_values3 = all_month_dict_suberr['FORCE_EST'], \
#    omi_data_type = 'pert', labels = ['adderror', 'suberror'])

# Same as plot_type_forcing_v3_all_months_arctic_avg, but it plots the 
# Arctic-averaged results for each of the 2005 - 2020 reference ice
# simulations. Also works for the 'refcld' simulations.a
# ptype: 'forcing', 'error', 'pcnt_error'
# return_slope: True to return the calculated Δflux
# --------------------------------------------------------------------
plot_type_forcing_v3_all_months_arctic_avg_manyrefice(all_month_dict['FORCE_EST'], \
       OMI_daily_VSJ4, minlat = 65, trend_type = 'standard', stype = 'cld', \
       ptype = 'forcing', vtype = 'v3')

# Calculates the slopes of the refice or refcld simulations plotted
# in the "plot_type_forcing_v3_all_months_arctic_avg_manyrefice" function.
# ------------------------------------------------------------------------
#ice_slopes = calc_forcing_slopes_v3_all_months_arctic_avg_manyrefice(all_month_dict['FORCE_EST'], \
#       OMI_daily_VSJ4, minlat = 65, trend_type = 'standard', stype = 'ice', \
#       ptype = 'forcing')

# Print the results of the refice/cld simulation trend comparisons as 
# a table
# -------------------------------------------------------------------
#calc_print_forcing_slope_error_v3(all_month_dict['FORCE_EST'], \
#       OMI_daily_VSJ4, minlat = 65, trend_type = 'standard', \
#       vtype = 'v3', ptype = 'forcing')

sys.exit()
 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#                      END FORCING VERSION 3 STUFF
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

plot_dual_combined_multi_type(comp_grid_data_v6, combined_data, \
    xval = 'ai', \
    #cld_min = None, cld_max = None,\
    ch7_min = 0.05, ch7_max = 0.1,\
    ice_min1 = 0, ice_max1 = 20,\
    ice_min2 = 105, ice_max2 = None,\
    sza_min = 50, sza_max = 55,\
    ai_min = 2,  ai_max = None,\
    save = False, show_trend = False, shade_density = False, \
    trend_type = 'theil-sen')

plot_dual_combined_multi_type(comp_grid_data_v6, combined_data, \
    xval = 'ai', \
    #cld_min = None, cld_max = None,\
    ch7_min1 = 0.15, ch7_max1 = 0.2,\
    ch7_min2 = 0.15, ch7_max2 = 0.2,\
    ice_min1 = 0, ice_max1 = 20,\
    ice_min2 = 105, ice_max2 = None,\
    sza_min = 45, sza_max = 50,\
    ai_min = 2,  ai_max = None,\
    save = False, show_trend = False, shade_density = False, \
    trend_type = 'theil-sen')

sys.exit()
#plot_type_forcing_all_months(OMI_data, NSIDC_data, 'average', \
#        minlat = minlat, maxlat = maxlat, \
#        save = False, coloc_dict = lin_smth2_dict_v6)
min_cloud = 0.95
cld_idx = 0
sfc_type_idx = 2
maxerr = 2
plot_slopes_cloud_types_szamean(lin_smth2_dict_v6,cld_idx = 0, maxerr = 2, \
    data_type = 'raw', remove_high_error = True, hatch_cloud = False, \
    min_cloud = 0.95, save = False)

sys.exit()



sys.exit()

debug_data = pd.read_csv('debug_file_iceidx0_szaidx5_ch7idx3.txt', \
    delim_whitespace = True, names = \
    ['lat','lon','ai','sza','swf','sza_bin','sza_low_edge',\
    'sza_hgh_edge','ch7_bin','ch7_low_edge','ch7_hgh_edge'])

# Figure out which local pixels are in the debug file from raindrop
found_pixels = np.array([ppixel in debug_data['swf'].values \
    for ppixel in return_dict['local_ydata']])

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

ax.scatter(return_dict['local_xdata'][found_pixels], \
    return_dict['local_ydata'][found_pixels], s = 6, color = 'g', \
    label = 'used in comp_grid')
ax.scatter(return_dict['local_xdata'][~found_pixels], \
    return_dict['local_ydata'][~found_pixels], s = 6, color = 'r', \
    label = 'not in comp_grid')
ax.legend()
plt.show()

sys.exit()


sys.exit()

trend_type = 'linregress'
lin_raw_dict = calc_raw_grid_slopes(\
        combined_data, comp_grid_data, \
        ai_min = ai_min, ai_max = ai_max, \
        trend_type = trend_type, \
        )
lin_smth_dict = calc_raw_grid_slopes(\
        combined_data, comp_grid_data, \
        ai_min = ai_min, ai_max = ai_max, \
        trend_type = trend_type, \
        smoother = 'smooth', sizer = 1)
lin_smth2_dict = calc_raw_grid_slopes(\
        combined_data, comp_grid_data, \
        ai_min = ai_min, ai_max = ai_max, \
        trend_type = trend_type, \
        smoother = 'smoother', sizer = 1)


#lin_raw_dict_v7 = calc_raw_grid_slopes(\
#        combined_data, comp_grid_data_v7, \
#        ai_min = ai_min, ai_max = ai_max, \
#        trend_type = trend_type, \
#        )
#lin_smth_dict_v7 = calc_raw_grid_slopes(\
#        combined_data, comp_grid_data_v7, \
#        ai_min = ai_min, ai_max = ai_max, \
#        trend_type = trend_type, \
#        smoother = 'smooth', sizer = 1)
#lin_smth2_dict_v7 = calc_raw_grid_slopes(\
#        combined_data, comp_grid_data_v7, \
#        ai_min = ai_min, ai_max = ai_max, \
#        trend_type = trend_type, \
#        smoother = 'smoother', sizer = 1)

trend_type = 'theil-sen'
thl_raw_dict = calc_raw_grid_slopes(\
        combined_data, comp_grid_data, \
        ai_min = ai_min, ai_max = ai_max, \
        trend_type = trend_type, \
        )
thl_smth_dict = calc_raw_grid_slopes(\
        combined_data, comp_grid_data, \
        ai_min = ai_min, ai_max = ai_max, \
        trend_type = trend_type, \
        smoother = 'smooth', sizer = 1)
thl_smth2_dict = calc_raw_grid_slopes(\
        combined_data, comp_grid_data, \
        ai_min = ai_min, ai_max = ai_max, \
        trend_type = trend_type, \
        smoother = 'smoother', sizer = 1)

#thl_raw_dict_v7 = calc_raw_grid_slopes(\
#        combined_data, comp_grid_data_v7, \
#        ai_min = ai_min, ai_max = ai_max, \
#        trend_type = trend_type, \
#        )
#thl_smth_dict_v7 = calc_raw_grid_slopes(\
#        combined_data, comp_grid_data_v7, \
#        ai_min = ai_min, ai_max = ai_max, \
#        trend_type = trend_type, \
#        smoother = 'smooth', sizer = 1)
#thl_smth2_dict_v7 = calc_raw_grid_slopes(\
#        combined_data, comp_grid_data_v7, \
#        ai_min = ai_min, ai_max = ai_max, \
#        trend_type = trend_type, \
#        smoother = 'smoother', sizer = 1)


min_cloud = 0.95
cld_idx = 0
sfc_type_idx = 2
maxerr = 2
plot_slopes_cloud_types_szamean(lin_smth2_dict_v6,cld_idx = 0, maxerr = 2, \
    data_type = 'raw', remove_high_error = True, hatch_cloud = False, \
    min_cloud = 0.95, save = False)
sys.exit()

plot_compare_grid_climo_counts_3version(\
    thl_smth_dict_v6, thl_smth_dict_v11, thl_smth_dict_v12, 'raw', 1)

sys.exit()


plot_raw_grid_slopes(thl_raw_dict, save = False, vmin = -15, vmax = 15)
sys.exit()

mask_cloud = np.ma.masked_where(\
    lin_smth2_dict['raw_cldvals'][sfc_type_idx,:,:,cld_idx] < -9, \
    lin_smth2_dict['raw_cldvals'][sfc_type_idx,:,:,cld_idx])
hasher = np.ma.masked_where(mask_cloud < min_cloud, \
    mask_cloud)

cloud_slopes = np.ma.masked_where(mask_cloud < min_cloud, \
    lin_smth2_dict['raw_slopes'][sfc_type_idx,:,:])
clear_slopes = np.ma.masked_where(mask_cloud >= min_cloud, \
    lin_smth2_dict['raw_slopes'][sfc_type_idx,:,:])

cloud_slopes = np.ma.masked_where(\
    lin_smth2_dict['raw_stderr'][sfc_type_idx,:,:] > maxerr, \
    cloud_slopes)
clear_slopes = np.ma.masked_where(\
    lin_smth2_dict['raw_stderr'][sfc_type_idx,:,:] > maxerr, \
    clear_slopes)

cloud_means = np.nanmean(cloud_slopes, axis = 0)
clear_means = np.nanmean(clear_slopes, axis = 0)
cloud_std   = np.std(cloud_slopes, axis = 0)
clear_std   = np.std(clear_slopes, axis = 0)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(lin_smth2_dict_v6['sza_mins'], cloud_means, color = 'tab:blue')
ax.plot(lin_smth2_dict_v6['sza_mins'], cloud_means - cloud_std, ':', color = 'tab:blue')
ax.plot(lin_smth2_dict_v6['sza_mins'], cloud_means + cloud_std, ':', color = 'tab:blue')
ax.plot(lin_smth2_dict_v6['sza_mins'], clear_means, color = 'tab:orange')
ax.plot(lin_smth2_dict_v6['sza_mins'], clear_means - clear_std, ':', color = 'tab:orange')
ax.plot(lin_smth2_dict_v6['sza_mins'], clear_means + clear_std, ':', color = 'tab:orange')

plt.show()

plot_slopes_cloud_types(lin_raw_dict, save = False, vmin = -15, vmax = 15)

sys.exit()

plot_compare_grid_climo_stderr(lin_smth_dict, 'raw_land', 2, save = False)

plot_raw_grid_slopes(thl_raw_dict, save = False, vmin = -15, vmax = 15)
plot_raw_grid_slopes(thl_smth_dict, save = False, vmin = -15, vmax = 15)
plot_raw_grid_slopes(thl_smth2_dict, save = False, vmin = -15, vmax = 15)


sys.exit()


plot_combined_scatter(combined_data, ax = ax1, \
    omi_min = ai_min,  omi_max  = ai_max, \
    sza_min = sza_min, sza_max = sza_max, \
    ice_min = ice_min, ice_max = ice_max, \
    ch7_min = ch7_min, ch7_max = ch7_max, \
    trend_type = 'theil-sen', show_trend = False)

plot_comp_grid_scatter(comp_grid_data, ax = ax2, xval = 'ai', \
    ai_min = ai_min,  ai_max  = ai_max, \
    sza_min = sza_min, sza_max = sza_max, \
    ice_min = ice_min, ice_max = ice_max, \
    ch7_min = ch7_min, ch7_max = ch7_max)












data = h5py.File('comp_data/colocated_subset_200607260156.hdf5','r')
#data = h5py.File('comp_data/colocated_subset_201507082016.hdf5','r')

mask_cod = np.ma.masked_where(data['modis_cod'][:,:] <= 0, data['modis_cod'][:,:])
mask_cod = np.ma.masked_invalid(mask_cod)
lat = data['omi_lat'][:,:]
lon = data['omi_lon'][:,:]
mask_ch1 = np.ma.masked_where(data['modis_ch1'][:,:] < 0, data['modis_ch1'][:,:])
mask_omi = np.ma.masked_where(data['omi_uvai_raw'][:,:] < -100, data['omi_uvai_raw'][:,:])

data.close()

fig = plt.figure(figsize = (7,6))
ax1 = fig.add_subplot(2,2,1, projection = mapcrs)
ax2 = fig.add_subplot(2,2,2, projection = mapcrs)
ax3 = fig.add_subplot(2,2,3, projection = mapcrs)
ax4 = fig.add_subplot(2,2,4)

ax1.pcolormesh(lon, lat, mask_omi, transform = datacrs, shading = 'auto', cmap = 'jet')
ax2.pcolormesh(lon, lat, mask_ch1, transform = datacrs, shading = 'auto', cmap = 'Greys_r')
ax3.pcolormesh(lon, lat, mask_cod, transform = datacrs, shading = 'auto', cmap = 'viridis', vmax = 60)

ax1.coastlines()
ax2.coastlines()
ax3.coastlines()

ax1.set_extent([-180,180,65,90], datacrs)
ax2.set_extent([-180,180,65,90], datacrs)
ax3.set_extent([-180,180,65,90], datacrs)

ax4.hist(mask_cod.compressed(), bins = 'auto')
ax4.set_yscale('log')
ax4.set_xlabel('MODIS COD')
ax4.set_ylabel('Counts')

ax1.set_title('OMI UVAI Raw')
ax2.set_title('MODIS CH1 Reflectance')
ax3.set_title('MODIS Cloud Optical Depth')

fig.tight_layout()

print(np.nanmax(mask_cod))

plt.show()

sys.exit()

run_list = ['20060724','20060725','20060726','20060727', \
            '20080422',\
            '20140811','20140812','20150627','20150706','20150707','20150708',\
            '20150709','20150710','20170816']
run_list = ['20170817','20170818']

run_list = ['20170819']
run_list = [
    '20180704','20180705','20180721','20180810','20180814', \
    '20180826','20190810','20190811','20210801']

#run_list = [\

#final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = True, \
#    images = False, process = False, run_list = run_list, copy_to_raindrop = False, \
#    include_tropomi = True, remove_ch2_file = True)
final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = False, \
    images = False, process = True, run_list = run_list, copy_to_raindrop = True, \
    include_tropomi = True, remove_ch2_file = True)

sys.exit()



# NOTE: As of 20230623, dates "20210704" though "20210730" have been
#       downloaded, processed, AND colocated on raindrop.
run_list = [
#    '20210704',
#    '20210705',
#    '20210711',
#    '20210712',
#    '20210713',
#    '20210729',
#    '20210730',
    '20210731',
    '20210801',
    '20210802',
#    '20210803',
#    '20210804',
#    '20210805',
#    '20210806',
#    '20210807',
#    '20210808',
#    '20210809',
#    '20210810',
#    '20210811',
#    '20210812',
#    '20210813',
#    '20210911',
#    '20210912',
#    '20210913',
]

#final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = True, \
#    images = False, process = False, run_list = run_list, copy_to_raindrop = False, \
#    include_tropomi = True, remove_ch2_file = True)
final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = False, \
    images = True, process = True, run_list = run_list, copy_to_raindrop = True, \
    include_tropomi = True, remove_ch2_file = True)

sys.exit()






#coloc_data = '201807050048'
#plot_compare_colocate_cloud(coloc_data, save = False)
#sys.exit()

run_list = ['20060724','20060725','20060726','20060727','20080422',\
    '20140811','20140812','20150627','20150706','20150707','20150708',\
    '20150709','20150710','20170816','20170817','20170818','20170819',\
    '20180704','20180705','20180721','20180810','20180814']

run_list = [\
    '20180826','20190810','20190811','20210801']

# NOTE: RERUN FOR THE 2017 DAYS, BUT TEMPORARILY MOVING THE GIANT
#       CERES SSFL2 FILES THAT WERE USED FOR THE NAAPS ALBEDO
#       STUDY. SLOWING DOWN THE RUNTIME SUBSTANTIALLY

#run_list = ['20060725']
#run_list = ['20200722', '20200723']
#run_list = ['20180721','20180810','20180826','20180827']
#run_list = ['20180721','20180810','20180814','20180826','20180827']
#run_list = ['20180721','20170814','20100731','20100801']
#final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = True, \
#    images = False, process = False, run_list = run_list, copy_to_raindrop = False)
final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = False, \
    images = False, process = True, run_list = run_list, copy_to_raindrop = True, \
    include_tropomi = True, remove_ch2_file = True)

sys.exit()


#cld_mins = np.arange(0.0,91,10.)
#cld_maxs = np.arange(10,101,10.)

#plot_combined_scatter(combined_data, ax = ax1, \
#    omi_min = ai_min,  omi_max  = ai_max, \
#    sza_min = sza_min, sza_max = sza_max, \
#    ice_min = ice_min, ice_max = ice_max, \
#    ch7_min = ch7_min, ch7_max = ch7_max, \
#    trend_type = 'theil-sen', show_trend = False)
#
#plot_comp_grid_scatter(comp_grid_data, ax = ax2, xval = 'ai', \
#    ai_min = ai_min,  ai_max  = ai_max, \
#    sza_min = sza_min, sza_max = sza_max, \
#    ice_min = ice_min, ice_max = ice_max, \
#    ch7_min = ch7_min, ch7_max = ch7_max)

#ch7_min = 0.05
#ch7_max = 0.1
#ice_min = 105
#ice_max = None
#
#plot_dual_combined_grid_climo(comp_grid_data, combined_data, xval = 'ai', 
#    ch7_min = ch7_min, ch7_max = ch7_max, ice_min = ice_min, 
#    ice_max = ice_max, sza_min = sza_min, sza_max = sza_max, 
#    ai_min = ai_min,  ai_max = ai_max, save = False, show_trend = True)
#    
#sys.exit()

for ii in range(len(ice_mins)):
    #for jj in range(len(ch7_mins)):
    for jj in range(len(ch7_mins)):
        print(ice_mins[ii], ice_maxs[ii], ch7_mins[jj], ch7_maxs[jj])
        #print(ice_mins[ii], ice_maxs[ii], ch7_mins[jj], ch7_maxs[jj])
    
        plot_dual_combined_grid_climo(comp_grid_data, combined_data, \
            xval = 'ai', \
            #cld_min = cld_min, cld_max = cld_max,\
            #ch7_min = ch7_min, ch7_max = ch7_max,\
            ch7_min = ch7_mins[jj], ch7_max = ch7_maxs[jj],\
            ice_min = ice_mins[ii], ice_max = ice_maxs[ii],\
            sza_min = sza_min, sza_max = sza_max,\
            ai_min = ai_min,  ai_max = ai_max,\
            save = False, show_trend = True, shade_density = True, \
            trend_type = 'theil-sen')


sys.exit()







#date_str = '201807052213'
#plot_compare_combined_category(date_str, var1 = 'TROP_AI', \
#    var2 = 'CERES_SWF', var3 = None, cat = "ALL", minlat = 65., \
#    xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
#    colorbar = True, trend = 'lin_regress', zoom = True, color = None, \
#    save = False)
#
#sys.exit()
#
#date_strs = ['201807052213']
#automate_all_preprocess(date_strs, download = False, images = False, process = True,\
#    omi_dtype = 'ltc3', copy_to_raindrop = True)
#sys.exit()


run_list = ['20180704', '20180705']
#run_list = ['20200722', '20200723']
#run_list = ['20180721','20180810','20180826','20180827']
#run_list = ['20180721','20180810','20180814','20180826','20180827']
#run_list = ['20180721','20170814','20100731','20100801']
#final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = True, \
#    images = False, process = False, run_list = run_list)
final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = False, \
    images = True, process = True, run_list = run_list, copy_to_raindrop = True)

sys.exit()

date_strs = [\
           '201807210047',
           ###'201807211358', BAD. terminator line
           '201807211537',
           '201807211716',
           '201807211855',
           '201807212034',
           '201808100200',
           '201808140135',
           '201808141804',
           '201808141942',
           '201808142121',
           ###'201808142300', BAD. Error in MODIS and CERES data extent
           '201808260158',
           '201808260337',
           '201808260655',
           '201908100129',
           '201908100308', # GOOD OCEAN COMPARISON
           '201908101936',
           '201908102115',
           '201908102254',
           '201908110033',
           '201908110212',
           '201908110351',
           '201908110708',
           '201908111523',
           '201908111702',
           '201908111841',
           '201908112019',
           ###'201908112158', # BAD
           ###'201908112337', # BAD
           #'202108010117',
           #'202108010256',
           #'202108010435',
           #'202108010614',
           #'202108011607',
           #'202108011925',
           #'202108012103',
           #'202108012242',
    ]


for date_str in date_strs:
    plot_compare_combined_category(date_str, var1 = 'TROP_AI', \
        var2 = 'CERES_SWF', var3 = None, cat = "ALL", minlat = 65., \
        xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
        colorbar = True, trend = 'lin_regress', zoom = True, color = None, \
        save = False)


sys.exit()




run_list = ['20200722', '20200723']
#run_list = ['20180721','20180810','20180826','20180827']
#run_list = ['20180721','20180810','20180814','20180826','20180827']
#run_list = ['20180721','20170814','20100731','20100801']
final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = True, \
    images = False, process = False, run_list = run_list)
#final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = False, \
#    images = True, process = True, run_list = run_list, copy_to_raindrop = True)

sys.exit()



date_str = '201807052213'
plot_compare_AI_combined_category(date_str, var2 = 'CERES_SWF', \
    cat = "ALL", minlat = 65., \
    xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
    colorbar = True, trend = 'theil-sen', zoom = True, color = None, \
    compare_tropomi = True, save = False)
sys.exit()


run_list = ['20190810','20190811']
#final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = True, \
#    images = False, process = False, run_list = run_list, copy_to_raindrop = True)
final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = False, \
    images = True, process = True, run_list = run_list, copy_to_raindrop = True)

sys.exit()


#date_strs = [\
#                     #'200804221841',  # GOOD
#                     #'200804222020',  # GOOD
#                     #'200804222159',  # GOOD
#                     #'201507061711',        
#                     #'201507061850',
#                     #'201507062028',
#                     '201708191713',
#           #'201908110033',  # GOOD
#           #'201507062207',
#           #'201908102254',  # GOOD
#           #'201807052034',
#            ]



#date_str = '201807051856'
#date_str = '201807052213'
#date_str = '201908110033'

for date_str in case_dates:

#date_str = date_strs[0]
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    
    comp_data = h5py.File('comp_data/colocated_subset_' + date_str + '.hdf5','r')
    
    base_path = home_dir + '/data/OMI/H5_files/'
    total_list = subprocess.check_output('ls '+base_path+\
        dt_date_str.strftime('OMI-Aura_L2-OMAERUV_%Ym%m%dt%H%M*.he5'),\
        shell=True).decode('utf-8').strip().split('\n')
    
    #print(total_list[0])
    omi_data = h5py.File(total_list[0],'r')
    #omi_data  = h5py.File(dt_date_str.strftime(\
    #    '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_%Ym%m%dt%H%M*'),'r')
    
    SSA = omi_data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/FinalAerosolSingleScattAlb'][:,:,:]
    
    mask_ssa = np.ma.masked_where(SSA < -2e5, SSA)
    
    mask_uvai = comp_data['omi_uvai_pert'][~mask_ssa[:,:,0].mask]
    mask_uvai = np.ma.masked_invalid(mask_uvai)
    mask_uvai = np.ma.masked_where(mask_uvai < 1.0, mask_uvai)
    
    mask_ssa0  = mask_ssa[:,:,0].compressed()
    mask_ssa0  = mask_ssa0[~mask_uvai.mask]
    
    mask_uvai = mask_uvai.compressed()
   
    if(len(mask_ssa0) > 100): 
        print(date_str, np.nanmean(mask_uvai), np.nanmean(mask_ssa))

        plt.close('all')
        fig = plt.figure(figsize = (12, 9))
        ax1 = fig.add_subplot(2,3,1, projection = mapcrs)
        ax2 = fig.add_subplot(2,3,2, projection = mapcrs)
        ax3 = fig.add_subplot(2,3,3, projection = mapcrs)
        ax4 = fig.add_subplot(2,3,4, projection = mapcrs)
        ax5 = fig.add_subplot(2,3,5)
        
        # Plot AI
        ax1.pcolormesh(comp_data['omi_lon'], comp_data['omi_lat'], comp_data['omi_uvai_pert'], \
            transform = datacrs, shading = 'auto', cmap = 'jet')
        ax1.coastlines()
        ax1.set_extent([-180, 180, 65, 90], datacrs)
        
        # Plot SSA
        ax2.pcolormesh(comp_data['omi_lon'], comp_data['omi_lat'], mask_ssa[:,:,0], \
            transform = datacrs, shading = 'auto')
        ax2.coastlines()
        ax2.set_extent([-180, 180, 65, 90], datacrs)
        
        ax3.pcolormesh(comp_data['omi_lon'], comp_data['omi_lat'], mask_ssa[:,:,1], \
            transform = datacrs, shading = 'auto')
        ax3.coastlines()
        ax3.set_extent([-180, 180, 65, 90], datacrs)
        
        ax4.pcolormesh(comp_data['omi_lon'], comp_data['omi_lat'], mask_ssa[:,:,2], \
            transform = datacrs, shading = 'auto')
        ax4.coastlines()
        ax4.set_extent([-180, 180, 65, 90], datacrs)
        
        # Plot scatter
        ax5.scatter(mask_ssa0, mask_uvai)
        ax5.set_xlabel('SSA')
        ax5.set_ylabel('UVAI')
  
        plt.title(date_str)

        outname = 'comp_ssa_' + date_str + '.png' 
        fig.savefig(outname)
        print('Saved image', outname)
 
    comp_data.close()
    omi_data.close()

#plt.show()



sys.exit()

slope_dict, extract_dict = plot_compare_all_slopes(date_strs = None, save = False, return_slope_dict = True)

sys.exit()

date_strs = [\
                     '201507100017',
                     '201507100156',
                     '201507100335',
                     '201507100514',
                     '201507100653',
                     '201507101010',
            ]

##!#
##!#date_strs = ['201506271856']  # GOOD
for date_str in date_strs:  
 
    #plot_compare_OMI_CERES_MODIS_NSIDC(date_str, 7, \
    #    omi_dtype = 'ltc3', minlat = 65., zoom = True, save = False)

    plot_compare_combined_category(date_str, var1 = 'OMI', \
        var2 = 'CERES_SWF', var3 = None, cat = "ALL", minlat = 65., \
        xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
        colorbar = True, trend = 'lin_regress', zoom = True, color = None, \
        save = True)
sys.exit()



# Uses the slope statistics returned from the other function, 
# along with the trend values of the month idx passed, to determine
# the forcing for each type
calculate_type_forcing(4, trend_type = 'linear', minlat = 65.)

sys.exit()

date_str = '20150707'
plot_compare_combined_category_event(date_str, var1 = 'OMI', \
        var2 = 'CERES_SWF', cat = "ALL", minlat = 65., \
        xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
        colorbar = True, trend = True, zoom = False, color = None, \
        save = False)
sys.exit()

# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------


##!##print(final_list)
##!##        for ttime in good_list:
##!##            if(new_only and (ttime not in out_time_dict.keys())):
##!#
##!#with open('new_found_swaths.txt','r') as fin:
##!#    total_time_dict = json.load(fin)
##!#
##!#just_days = [ttime[:8] for ttime in list(total_time_dict.keys())]
##!#need_data = {}
##!#for ii in range(len(just_days)):
##!#    need_data[just_days[ii]] = True
##!##need_data = [True for ii in range(len(just_days))]
##!#
##!#for ttime in total_time_dict.keys():
##!#    if(total_time_dict[ttime] == True):
##!#        print(ttime, 'is True. Don\'t need data for today')
##!#        need_data[ttime[:8]] = False

# bad days
# 20060728 19-20
# 20120724
# 20120725
# 20120727

# good days:
# 20100629
# 20100731
# 20100801
# 20130803
# 20140811 GOOD
# 20140812 GOOD
# 20140814 GOOD
# 20140816
# 20150627 GOOD
# 20150706

# other good days
#


sys.exit()

date_str = '20140811'
download_OMI_files(date_str, omi_dtype = 'ltc3')
sys.exit()

#date_str = '200804221935'
#date_str = '200804222110'
#date_str = '200804222250'
#date_str = '201807051950'
#date_str = '201807052125'
date_str = '201807052305'
#date_str = '201908110125'
#date_str = '201908110440'

date_strs = ['200607240029', # GOOD
             #'200607240208', # GOOD / CERES mismatch
             '200607240347', # GOOD
             '200607240526', # GOOD
             '200607240844', # GOOD
             '200607242155', # GOOD
             '200607242334', # GOOD
             '200607250112', # GOOD
             '200607250251', # GOOD
             '200607250748', # GOOD?
             '200607252238', # GOOD
             '200607260017', # GOOD
             '200607260156', # GOOD
             '200607260335', # GOOD
             '200607260513', # GOOD?
             '200607260831', # GOOD
             '200607262142', # GOOD
             '200607270100', # GOOD
             '200607270239', # GOOD?
             '200607270418', # GOOD?
             '200607270557', # GOOD?
             '200607270736', # GOOD?
             '200607272226', # GOOD
             '200804221841',  # GOOD
             '200804222020',  # GOOD
             '200804222159',  # GOOD
             '201408110046',
             '201408110404',
             '201408112032',
             '201408112211',
             '201408112350',
             '201408120129',
             '201408120308',
             '201408122115',
             '201408122254',
             '201506271538',
             '201506271717',
             '201506271856',
             '201708161504',  # GOOD
             '201708161643',  # GOOD
             '201708161821',  # GOOD
             '201708171408',  # GOOD
             '201708171547',  # GOOD
             '201708171726',  # GOOD
             '201708171905',  # GOOD
             '201708172043',  # GOOD
             '201708181312',  # GOOD
             '201708181451',  # GOOD
             '201708181630',  # GOOD
             '201708181809',  # GOOD
             '201708181948',  # GOOD
             '201708191355',  # GOOD
             '201708191534',  # GOOD
             '201708191713',  # GOOD
             '201807051856',  # GOOD
             '201807052034',  # GOOD
             '201807052213',  # GOOD
             '201908102115',  # GOOD
             '201908102254',  # GOOD
             '201908110033',  # GOOD
             '201908110351',  # GOOD
            ]
##             ##!#'201605151925',  # MEDIOCRE
##             ##!#'201605152104',  # MEDIOCRE
##             ##!#'201605152243',  # MEDIOCRE
##             ##!#'201605162148',  # MEDIOCRE
##             ##!#'200607260017',  # GOOD
##             ##!#'200607252238',  # GOOD
##             ##!#'200607260156',  # GOOD
##             ##!#'200607260335',  # GOOD
##             ##!#'200607260513',  # GOOD
##             '201808241343',
##            ]

#auto_all_download(date_strs, download = False, rewrite_json = True)
#automate_all_preprocess(date_strs, download = False, images = False, process = True,\
#    omi_dtype = 'ltc3')
#sys.exit()

#date_str = '201908110033'
#date_str = '201708171547'
#for dstr in date_strs:
#    plot_compare_OMI_CERES_MODIS_NSIDC(dstr, 7, \
#        omi_dtype = 'shawn', minlat = 65., zoom = True, save = True)
#sys.exit()

#coloc_data = date_str
#plot_compare_combined_category(coloc_data, var1 = 'OMI', \
#    var2 = 'CERES_SWF', var3 = None, cat = "ALL", minlat = 65., \
#    xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
#    colorbar = True, trend = False, zoom = False, color = None, \
#    save = False)
#sys.exit()
#date_str = '201708161504'
##!#date_str = '201807052034'

#

slope_dict = {}

#coloc_data = '201708161504'
num_points = 2
for date_str in date_strs:
    print(date_str)
    slope_dict[date_str] = event_category_slopes_all(date_str, 'OMI', 'CERES_SWF', var3 = None, \
        cat = "ALL", minlat = 65., xmin = 1.0, xmax = None, ymin = None, \
        ymax = None, trend = False, num_points = num_points, \
        restrict_sza = False, color = None, save = False)


# Extract just the lin regress slopes
# -----------------------------------
extract_dict = {}
dkeys = slope_dict['201708181451'].keys()

for dkey in dkeys:
    extract_dict[dkey] = {}

    lin_slopes    = np.ma.masked_invalid(np.array([slope_dict[tkey][dkey]['Linear'] for tkey in slope_dict.keys()]))
    lin_pvals     = np.ma.masked_invalid(np.array([slope_dict[tkey][dkey]['lin_pval'] for tkey in slope_dict.keys()]))
    thiel_slopes  = np.ma.masked_invalid(np.array([slope_dict[tkey][dkey]['Thiel'] for tkey in slope_dict.keys()]))

    extract_dict[dkey]['Linear']   = np.ma.masked_where((lin_slopes > 500) | (lin_pvals > 0.1), lin_slopes)
    #extract_dict[dkey]['Linear']   = np.ma.masked_where((lin_slopes > 500), lin_slopes)
    extract_dict[dkey]['Thiel']    = np.ma.masked_where(thiel_slopes > 500, thiel_slopes)
    extract_dict[dkey]['lin_pval'] = lin_pvals


fig = plt.figure(figsize = (9, 11))
ax1 = fig.add_subplot(3,2,1)
ax2 = fig.add_subplot(3,2,2)
ax3 = fig.add_subplot(3,2,3)
ax4 = fig.add_subplot(3,2,4)
ax5 = fig.add_subplot(3,2,5)
ax6 = fig.add_subplot(3,2,6)

num_bins = 20
ax1.hist(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Linear']), bins = num_bins, alpha = 0.5, label = 'Linear')
ax1.hist(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Thiel']), bins = num_bins, alpha = 0.5, label = 'Thiel')
print('ICE_CLOUD')
print('\tLinear:',np.nanmean(extract_dict['ICE_CLOUD']['Linear']), np.nanstd(extract_dict['ICE_CLOUD']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['ICE_CLOUD']['Thiel']), np.nanstd(extract_dict['ICE_CLOUD']['Thiel']))
ax2.hist(extract_dict['ICE_CLEAR']['Linear'], bins = num_bins, alpha = 0.5, label = 'Linear')
ax2.hist(extract_dict['ICE_CLEAR']['Thiel'], bins = num_bins, alpha = 0.5, label = 'Thiel')
print('ICE_CLEAR')
print('\tLinear:',np.nanmean(extract_dict['ICE_CLEAR']['Linear']), np.nanstd(extract_dict['ICE_CLEAR']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['ICE_CLEAR']['Thiel']), np.nanstd(extract_dict['ICE_CLEAR']['Thiel']))
ax3.hist(extract_dict['OCEAN_CLOUD']['Linear'], bins = num_bins, alpha = 0.5, label = 'Linear')
ax3.hist(extract_dict['OCEAN_CLOUD']['Thiel'], bins = num_bins, alpha = 0.5, label = 'Thiel')
print('OCEAN_CLOUD')
print('\tLinear:',np.nanmean(extract_dict['OCEAN_CLOUD']['Linear']), np.nanstd(extract_dict['ICE_CLOUD']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['OCEAN_CLOUD']['Thiel']), np.nanstd(extract_dict['ICE_CLOUD']['Thiel']))
ax4.hist(extract_dict['OCEAN_CLEAR']['Linear'], bins = num_bins, alpha = 0.5, label = 'Linear')
ax4.hist(extract_dict['OCEAN_CLEAR']['Thiel'], bins = num_bins, alpha = 0.5, label = 'Thiel')
print('OCEAN_CLEAR')
print('\tLinear:',np.nanmean(extract_dict['OCEAN_CLEAR']['Linear']), np.nanstd(extract_dict['OCEAN_CLEAR']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['OCEAN_CLEAR']['Thiel']), np.nanstd(extract_dict['OCEAN_CLEAR']['Thiel']))
ax5.hist(extract_dict['LAND_CLOUD']['Linear'], bins = num_bins, alpha = 0.5, label = 'Linear')
ax5.hist(extract_dict['LAND_CLOUD']['Thiel'], bins = num_bins, alpha = 0.5, label = 'Thiel')
print('LAND_CLOUD')
print('\tLinear:',np.nanmean(extract_dict['LAND_CLOUD']['Linear']), np.nanstd(extract_dict['ICE_CLOUD']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['LAND_CLOUD']['Thiel']), np.nanstd(extract_dict['ICE_CLOUD']['Thiel']))
ax6.hist(extract_dict['LAND_CLEAR']['Linear'], bins = num_bins, alpha = 0.5, label = 'Linear')
ax6.hist(extract_dict['LAND_CLEAR']['Thiel'], bins = num_bins, alpha = 0.5, label = 'Thiel')
print('LAND_CLEAR')
print('\tLinear:',np.nanmean(extract_dict['LAND_CLEAR']['Linear']), np.nanstd(extract_dict['ICE_CLEAR']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['LAND_CLEAR']['Thiel']), np.nanstd(extract_dict['ICE_CLEAR']['Thiel']))
#ax2.scatter(np.ma.masked_invalid(extract_dict[var]['Linear']), np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Thiel']))

ax1.set_title('ICE_CLOUD')
ax2.set_title('ICE_CLEAR')
ax3.set_title('OCEAN_CLOUD')
ax4.set_title('OCEAN_CLEAR')
ax5.set_title('LAND_CLOUD')
ax6.set_title('LAND_CLEAR')

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
ax5.legend()
ax6.legend()
plt.show()

sys.exit()


#date_str = '20060726'

####var1 = 'OMI'
####var2 = 'CERES_SWF'
######!##plot_compare_OMI_CERES_MODIS_NSIDC(date_str, 7, \
######!##    omi_dtype = 'shawn', minlat = 65., zoom = True, save = False)
######!##sys.exit()
######!#
######!##plot_compare_scatter(date_str, var1, var2, var3 = 'NSIDC_LAND', minlat = 65., \
######!##    xmin = 1, zoom = False, save = False, trend = True)
######!##plot_compare_colocate_spatial(date_str, minlat = 65., zoom = False, \
######!##    save = False)
######!#cat = 'ICE_CLOUD'
######!##cat = 'OCEAN_CLOUD'
######!##cat = 'LAND_CLEAR'
######!##plot_compare_colocate_spatial_category(date_str, cat = cat, minlat = 65., \
######!##    zoom = True, save = False)
####trend = True 
######!#
#####date_str = '20170816'
#####date_str = '20170818'
####data = read_colocated_combined(date_str, zoom = True)
#####data = '201708171547'
####
####fig = plt.figure(figsize = (12,4))
####ax1 = fig.add_subplot(1,3,1)
####ax2 = fig.add_subplot(1,3,2)
####ax3 = fig.add_subplot(1,3,3)
######ax1 = fig.add_subplot(2,3,1)
######ax2 = fig.add_subplot(2,3,4)
######ax3 = fig.add_subplot(2,3,2)
######ax4 = fig.add_subplot(2,3,5)
######ax5 = fig.add_subplot(2,3,3)
######ax6 = fig.add_subplot(2,3,6)
####
####plot_compare_scatter_category(data, var1, var2, var3 = None, \
####    cat = 'ICE_CLOUD', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
####    ax = ax1, colorbar = True, trend = trend, zoom = False, save = False,\
####    color = 'tab:blue')
####plot_compare_scatter_category(data, var1, var2, var3 = None, \
####    cat = 'ICE_CLEAR', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
####    ax = ax1, colorbar = True, trend = trend, zoom = False, save = False,\
####    color = 'tab:orange')
####plot_compare_scatter_category(data, var1, var2, var3 = None, \
####    cat = 'OCEAN_CLOUD', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
####    ax = ax2, colorbar = True, trend = trend, zoom = False, save = False, \
####    color = 'tab:blue')
####plot_compare_scatter_category(data, var1, var2, var3 = None, \
####    cat = 'OCEAN_CLEAR', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
####    ax = ax2, colorbar = True, trend = trend, zoom = False, save = False, \
####    color = 'tab:orange')
####plot_compare_scatter_category(data, var1, var2, var3 = None, \
####    cat = 'LAND_CLOUD', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
####    ax = ax3, colorbar = True, trend = trend, zoom = False, save = False, \
####    color = 'tab:blue')
####plot_compare_scatter_category(data, var1, var2, var3 = None, \
####    cat = 'LAND_CLEAR', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
####    ax = ax3, colorbar = True, trend = trend, zoom = False, save = False, \
####    color = 'tab:orange')
#####plt.suptitle(data['date_str'])
#####plt.suptitle(data)
####fig.tight_layout()
####
####outname = 'arctic_daily_scatter_' + date_str + '.png'
#####fig.savefig(outname, dpi = 300)
####print("Saved image", outname)
####
####plt.show()

sys.exit()


#date_str = '201908110351'
##date_str = '200804222020'
#date_str = '201908110033'
#plot_compare_OMI_CERES_MODIS_NSIDC(date_str, 7, \
#    omi_dtype = 'shawn', minlat = 65., zoom = True, save = False)
#sys.exit()





##!#
##!#out_time_dict, out_file_dict = auto_all_download(date_strs, download = True, rewrite_json = True)
##!#sys.exit()
##!#
##!#

