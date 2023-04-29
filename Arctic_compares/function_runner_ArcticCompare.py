#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *

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


date_strs = ['201807052213']
automate_all_preprocess(date_strs, download = False, images = False, process = True,\
    omi_dtype = 'ltc3', copy_to_raindrop = True)
sys.exit()

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

