#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *


# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------

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
for date_str in date_strs:
    print(date_str)
    slope_dict[date_str] = event_category_slopes_all(date_str, 'OMI', 'CERES_SWF', var3 = None, \
        cat = "ALL", minlat = 65., xmin = 1.0, xmax = None, ymin = None, \
        ymax = None, trend = False, zoom = False, \
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

    extract_dict[dkey]['Linear']   = np.ma.masked_where(lin_slopes > 500, lin_slopes)
    extract_dict[dkey]['Thiel']    = np.ma.masked_where(thiel_slopes > 500, thiel_slopes)
    extract_dict[dkey]['lin_pval'] = lin_pvals


fig = plt.figure(figsize = (9, 11))
ax1 = fig.add_subplot(3,2,1)
ax2 = fig.add_subplot(3,2,2)
ax3 = fig.add_subplot(3,2,3)
ax4 = fig.add_subplot(3,2,4)
ax5 = fig.add_subplot(3,2,5)
ax6 = fig.add_subplot(3,2,6)

ax1.hist(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Linear']), bins = 35, alpha = 0.5, label = 'Linear')
ax1.hist(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Thiel']), bins = 35, alpha = 0.5, label = 'Thiel')
print('ICE_CLOUD')
print('\tLinear:',np.nanmean(extract_dict['ICE_CLOUD']['Linear']), np.nanstd(extract_dict['ICE_CLOUD']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['ICE_CLOUD']['Thiel']), np.nanstd(extract_dict['ICE_CLOUD']['Thiel']))
ax2.hist(extract_dict['ICE_CLEAR']['Linear'], bins = 35, alpha = 0.5, label = 'Linear')
ax2.hist(extract_dict['ICE_CLEAR']['Thiel'], bins = 35, alpha = 0.5, label = 'Thiel')
print('ICE_CLEAR')
print('\tLinear:',np.nanmean(extract_dict['ICE_CLEAR']['Linear']), np.nanstd(extract_dict['ICE_CLEAR']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['ICE_CLEAR']['Thiel']), np.nanstd(extract_dict['ICE_CLEAR']['Thiel']))
ax3.hist(extract_dict['OCEAN_CLOUD']['Linear'], bins = 35, alpha = 0.5, label = 'Linear')
ax3.hist(extract_dict['OCEAN_CLOUD']['Thiel'], bins = 35, alpha = 0.5, label = 'Thiel')
print('OCEAN_CLOUD')
print('\tLinear:',np.nanmean(extract_dict['OCEAN_CLOUD']['Linear']), np.nanstd(extract_dict['ICE_CLOUD']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['OCEAN_CLOUD']['Thiel']), np.nanstd(extract_dict['ICE_CLOUD']['Thiel']))
ax4.hist(extract_dict['OCEAN_CLEAR']['Linear'], bins = 35, alpha = 0.5, label = 'Linear')
ax4.hist(extract_dict['OCEAN_CLEAR']['Thiel'], bins = 35, alpha = 0.5, label = 'Thiel')
print('OCEAN_CLEAR')
print('\tLinear:',np.nanmean(extract_dict['OCEAN_CLEAR']['Linear']), np.nanstd(extract_dict['OCEAN_CLEAR']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['OCEAN_CLEAR']['Thiel']), np.nanstd(extract_dict['OCEAN_CLEAR']['Thiel']))
ax5.hist(extract_dict['LAND_CLOUD']['Linear'], bins = 35, alpha = 0.5, label = 'Linear')
ax5.hist(extract_dict['LAND_CLOUD']['Thiel'], bins = 35, alpha = 0.5, label = 'Thiel')
print('LAND_CLOUD')
print('\tLinear:',np.nanmean(extract_dict['LAND_CLOUD']['Linear']), np.nanstd(extract_dict['ICE_CLOUD']['Linear']))
print('\tThiel: ',np.nanmean(extract_dict['LAND_CLOUD']['Thiel']), np.nanstd(extract_dict['ICE_CLOUD']['Thiel']))
ax6.hist(extract_dict['LAND_CLEAR']['Linear'], bins = 35, alpha = 0.5, label = 'Linear')
ax6.hist(extract_dict['LAND_CLEAR']['Thiel'], bins = 35, alpha = 0.5, label = 'Thiel')
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

date_strs = ['201708161504']  # GOOD
for date_str in date_strs:  
 
##!#data = read_colocated_combined('20180705', zoom = True)
    plot_compare_combined_category(date_str, var1 = 'OMI', \
        var2 = 'CERES_SWF', var3 = None, cat = "ALL", minlat = 65., \
        xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
        colorbar = True, trend = 'lin_regress', zoom = True, color = None, \
        save = False)
sys.exit()

#date_str = '20060726'
date_str = '20170817'
plot_compare_combined_category_event(date_str, var1 = 'OMI', \
        var2 = 'CERES_SWF', cat = "ALL", minlat = 65., \
        xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
        colorbar = True, trend = True, zoom = False, color = None, \
        save = True)

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

