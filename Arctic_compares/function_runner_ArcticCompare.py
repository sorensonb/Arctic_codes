#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *
import random
#from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score


#minlat = 65
#maxlat = 87
#shawn_file = home_dir + '/Research/OMI/omi_shawn_daily_2005_2020_v2.hdf5'
#OMI_daily_VSJ4_high  = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -0.10, max_AI = 20.0, minlat = minlat, maxlat = maxlat)
#OMI_daily_VSJ4_low   = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -20.0, max_AI = 0.7, minlat = minlat, maxlat = maxlat)
#OMI_daily_VSJ4_all   = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -20.0, max_AI = 20.0, minlat = minlat, maxlat = maxlat)
#
#high_vals = np.ma.masked_where(OMI_daily_VSJ4_high['AI'] > 0.7, OMI_daily_VSJ4_high['AI'])
#low_vals  = np.ma.masked_where(OMI_daily_VSJ4_low['AI'] > 0.7, OMI_daily_VSJ4_low['AI'])
#all_vals  = np.ma.masked_where(OMI_daily_VSJ4_all['AI'] > 0.7, OMI_daily_VSJ4_all['AI'])
#
#fig = plt.figure(figsize = (10, 5))
#axs = fig.subplots(2, 3, sharex = True, sharey = True)
#flat_axs = axs.flatten()
#
#for ii in range(6):
#    flat_axs[ii].hist(np.nanmean(high_vals[ii::6,:,:], axis = 0).flatten(), bins = 50, alpha = 0.5)
#    flat_axs[ii].hist(np.nanmean(low_vals[ii::6,:,:], axis = 0).flatten(), bins = 50, alpha = 0.5)
#    flat_axs[ii].hist(np.nanmean(all_vals[ii::6,:,:], axis = 0).flatten(), bins = 50, alpha = 0.5)
#    flat_axs[ii].set_title(str(ii))
#    flat_axs[ii].grid(alpha = 0.5)
#
#plt.suptitle('Comparison of monthly UVAI averages with different min/maxs\n' + \
#    'Blue:   min_AI = -0.10, max_AI = 20.0\n' + \
#    'Orange: min_AI = -20.0, max_AI = 0.7\n' + \
#    'Green:  min_AI = -20.0, max_AI = 20.0')
#
#fig.tight_layout()
#
#outname = 'omi_month_climo_hist_compare.png'
#fig.savefig(outname, dpi = 200)
#print("Saved image", outname)
#
#plt.show()
#sys.exit()


#make_gif('comp_images_20180705/', 'calc_swf_comp_20180705.gif')

#date_str = '201807052213'
#plot_compare_OMI_MODIS_v2(date_str, 7, \
#    omi_dtype = 'ltc3', minlat = 65., zoom = True, save = True)
#sys.exit()


# CODE FOR PLOTTING RAW VARIABLES FROM OMI FILES
"""
data = h5py.File(sys.argv[1])

lat = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
lon = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
hgt = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/FinalAerosolLayerHeight'][:,:]
AI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]

hgt = np.ma.masked_where(hgt < 0, hgt)
AI  = np.ma.masked_where((AI < -2) | (AI > 10), AI)

data.close()

fig = plt.figure(figsize = (11, 4))
ax1 = fig.add_subplot(1,3,1, projection = ccrs.NorthPolarStereo())
ax2 = fig.add_subplot(1,3,2, projection = ccrs.NorthPolarStereo())
ax3 = fig.add_subplot(1,3,3)
ax1.pcolormesh(lon, lat, hgt, transform = datacrs, shading = 'auto')
ax1.coastlines()
ax1.set_extent([-180,180,65,90], datacrs)
ax2.pcolormesh(lon, lat, AI, transform = datacrs, shading = 'auto', cmap = 'jet')
ax2.coastlines()
ax2.set_extent([-180,180,65,90], datacrs)
ax3.hist(hgt.compressed())
fig.tight_layout()
plt.show()

sys.exit()
"""

# Plots the distribution of MYD08 COD standard deviations
# -------------------------------------------------------
#plotMODIS_MYD08_COD_distribution(save=True)
#sys.exit()

"""
#files = glob('comp_data/colocated_subset_20180705*.hdf5')
#files = glob('comp_data/original_prep_data/colocated_subset_20180705*')
#files = glob('comp_data/testing/colocated_subset_2018070*')
files = glob('neuralnet_output/test_calc_out_noland50_2014081*')

min_ai = -2.0
max_ai = 1.0
minlat = 65.    # 2024/01/10: changed to 65.
min_swf = 0.
max_swf = 3000.
max_cod = 70.
min_ice = 0.
max_ice = 500.

for ff in files:

    data = h5py.File(ff)
    dtime = ff.strip().split('/')[-1].split('_')[-1].split('.')[0]

    print(dtime, data.keys())
    local_data = np.ma.masked_invalid(data['omi_uvai_pert'])

    #local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
    #local_data =  np.ma.masked_where( (local_data < min_ai) | \
    #                           (local_data > max_ai) | \
    #                           (data['omi_lat'][:,:] < minlat) | \
    #                           ( (data['ceres_swf'][:,:] < min_swf) | \
    #                           (data['ceres_swf'][:,:] > max_swf) ) | \
    #                           (np.isnan(data['ceres_swf'][:,:]) == True) | \
    #                           (data['modis_cod'][:,:] == -999) | \
    #                           (data['modis_cod'][:,:] > max_cod) | \
    #                           (np.isnan(data['modis_cod'][:,:]) == True) | \
    #                           #(data['modis_cld_top_pres'][:,:] < 0) | \
    #                           #(data['modis_cld_top_pres'][:,:] > 1025) | \
    #                           #(np.isnan(data['modis_cld_top_pres'][:,:]) == True) | \
    #                           (data['nsidc_ice'][:,:] == -999) | \
    #                           (data['nsidc_ice'][:,:] == 251) | \
    #                           #(data['nsidc_ice'][:,:] == 254) | \
    #                           (data['nsidc_ice'][:,:] == 253), \
    #                           local_data)

    #local_omi   = np.ma.masked_where(local_data.mask == False, data['omi_uvai_raw'])
    #local_ceres = np.ma.masked_where(local_data.mask == False, data['ceres_swf'])
    #local_ch7   = np.ma.masked_where(local_data.mask == False, data['modis_ch7'])
    #local_cod   = np.ma.masked_where(local_data.mask == False, data['modis_cod'])
    #local_ctp   = np.ma.masked_where(local_data.mask == False, data['modis_cld_top_pres'])
    #local_ice   = np.ma.masked_where(local_data.mask == False, data['nsidc_ice'])

    #local_omi   = np.ma.masked_where(data['omi_uvai_raw'][:,:] == -999., data['omi_uvai_raw'][:,:])
    local_omi   = np.ma.masked_where(data['omi_uvai_pert'][:,:] == -999., data['omi_uvai_pert'][:,:])
    local_ceres = np.ma.masked_where((data['ceres_swf'][:,:] == -999.) |
                                     (data['ceres_swf'][:,:] > 3000), data['ceres_swf'][:,:])
    #local_ch7   = np.ma.masked_where(data['modis_ch7'][:,:] == -999., data['modis_ch7'][:,:])
    local_cod   = np.ma.masked_where(data['modis_cod'][:,:] == -999., data['modis_cod'][:,:])
    local_cod2  = np.where(np.isnan(data['modis_cod'][:,:]) == True, 0., data['modis_cod'][:,:])
    local_cod2  = np.ma.masked_where(data['modis_cod'][:,:] == -999., local_cod2)
    local_ctp   = np.ma.masked_where(data['modis_cld_top_pres'][:,:] == -999., data['modis_cld_top_pres'][:,:])
    #local_ctp   = np.ma.masked_where(data['modis_cld_top_pres'][:,:] == -999., data['modis_ctp'][:,:])
    local_ctp2   = np.ma.masked_where(local_ctp != 0., local_ctp)
    local_ice   = np.ma.masked_where(data['nsidc_ice'][:,:] == -999., data['nsidc_ice'][:,:])

    local_omi   = np.ma.masked_invalid(local_omi)
    local_ceres = np.ma.masked_invalid(local_ceres)
    #local_ch7   = np.ma.masked_invalid(local_ch7)
    local_cod   = np.ma.masked_invalid(local_cod)
    local_cod2  = np.ma.masked_invalid(local_cod2)
    local_ctp   = np.ma.masked_invalid(local_ctp)
    local_ctp2  = np.ma.masked_invalid(local_ctp2)
    #local_ctp   = np.ma.masked_invalid(local_ctp)
    local_ice   = np.ma.masked_invalid(local_ice)
 
    #combined_data['modis_cld_top_pres'] = \
    # np.where(combined_data['modis_cld_top_pres'] == 0., 1025., \
    # combined_data['modis_cld_top_pres'])
    ## TESTING: Make all NAN values to be "clear-sky"
    #combined_data['modis_cld_top_pres'] = \
    #    np.where(np.isnan(combined_data['modis_cld_top_pres'][:]) == True, 1025., \
    #    combined_data['modis_cld_top_pres'][:])

    fig = plt.figure(figsize = (10, 6))
    ax1 = fig.add_subplot(2,3,1, projection = ccrs.NorthPolarStereo()) # OMI UVAI
    ax2 = fig.add_subplot(2,3,2, projection = ccrs.NorthPolarStereo()) # CERES SWF
    ax3 = fig.add_subplot(2,3,3, projection = ccrs.NorthPolarStereo()) # MODIS CH7
    ax4 = fig.add_subplot(2,3,5, projection = ccrs.NorthPolarStereo()) # MODIS COD
    ax5 = fig.add_subplot(2,3,4, projection = ccrs.NorthPolarStereo()) # MODIS CLD TOP PRES 
    ax6 = fig.add_subplot(2,3,6, projection = ccrs.NorthPolarStereo()) # NSIDC ICE

    ax1.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], local_omi, \
        transform = ccrs.PlateCarree(), shading = 'auto', cmap = 'jet')
    ax1.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax1.coastlines()
    ax1.set_title('OMI AI')

    ax2.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], local_ceres, \
        transform = ccrs.PlateCarree(), shading = 'auto')
    ax2.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax2.coastlines()
    ax2.set_title('CERES SWF')

    #ax3.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], local_ch7, \
    #    transform = ccrs.PlateCarree(), shading = 'auto')
    #ax3.set_extent([-180,180,65,90], ccrs.PlateCarree())
    #ax3.coastlines()
    #ax3.set_title('MODIS CH7')

    ax4.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], local_cod2, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmax = 50)
    ax4.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax4.coastlines()
    ax4.set_title('MODIS COD NAN = 0')

    ax5.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], local_ctp, \
        transform = ccrs.PlateCarree(), shading = 'auto')
    ax5.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax5.coastlines()
    ax5.set_title('MODIS CTP')

    ax6.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], local_cod, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmax = 50)
    ax6.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax6.coastlines()
    ax6.set_title('MODIS COD NO NAN')

    data.close()

    plt.suptitle(dtime)

    fig.tight_layout()
    plt.show()

sys.exit()
    
"""




# Plot the NN architecture
# ------------------------
#plot_NN_architecture(save = True)
#plot_NN_architecture(plot_lines = True, save = True)
#sys.exit()

# Makes a plot of OMI AI, Observed SWF, NN SWF, and Regression
# Designed for showing the validation of the NN under
# aerosol-free conditions
# -----------------------------------------------------------
#plot_compare_NN_output_noaer(sys.argv[1], astrofit = True, save = False)
#sys.exit()

# Plot combined NN error distribution in clear-sky swaths
# -------------------------------------------------------
#plot_NN_error_dist_bulk('noland74', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland100', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland101', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland103', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland105', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland106', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland107', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland108', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland109', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland110', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland112', num_bins = 500, astrofit = True, \
#plot_NN_error_dist_bulk('noland113', num_bins = 500, astrofit = True, \
plot_NN_error_dist_bulk('noland115', num_bins = 500, astrofit = True, \
    use_correct_error_calc = True, xmin = -100, xmax = 100, save = False)
sys.exit()
#plot_NN_error_dist_bulk('noland74', num_bins = 500, astrofit = True, \
#sys.exit()

# Plots more stuff than the NN_output_noaer function
# --------------------------------------------------
#plot_compare_NN_output(sys.argv[1], save = False)
#sys.exit()

# Plot the zoomed in comparison
# -----------------------------
#plot_compare_NN_output_v2(sys.argv[1], auto_zoom = False, save = False)
#sys.exit()

# Compare the OMI, CERES, NN, and NN - CERES for two swaths
# ---------------------------------------------------------
plot_compare_NN_output_double('test_calc_out_noland103_201807082244.hdf5', \
    'neuralnet_output/test_calc_out_noland103_201807052213.hdf5', \
    save = False, include_scatter = False)
#plot_compare_NN_output_double(sys.argv[1], sys.argv[2], save = False, include_scatter = False)
sys.exit()

# CODE TO SEE HOW WELL THE DIFFERENT NN OUTPUTS COMPARE FOR A GIVEN
# SWATH TIME
# -----------------------------------------------------------------
#date_str = '201807052213'
#skip_version = ['noland51','noland52','noland60','noland61']
#compare_nn_version_output(date_str, skip_version = skip_version, save = False)
#
#sys.exit()

# Compare the plume locations between OMI and MODIS data
# ------------------------------------------------------
#plot_compare_NN_output_overlay(sys.argv[1], auto_zoom = True, save = False)
#sys.exit()

#plot_compare_NN_output_overlay_v2(sys.argv[1], auto_zoom = True, save = True)
#sys.exit()




#calc_data = sys.argv[1]
#in_calc = h5py.File(calc_data)
#
#mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
#                               (in_calc['ceres_swf'][:,:] > 3000), in_calc['ceres_swf'][:,:])
#mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
#mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)
#
#both_orig = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_orig).compressed()
#both_calc = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_calc).compressed()
#
#rmse = mean_squared_error(both_orig, both_calc, squared = False)
#mae = np.mean(abs(both_orig - both_calc))
#
#sys.exit()    

## Calculate time offset between OMI and CERES obs
## -----------------------------------------------
#dir_list = [\
#    '200607240844', \
#    '200607242334', \
#    '200607270100', \
#    '200804221841', \
#    '201408112211', \
#    '201708171050', \
#    '201807052213', \
#    '201908110033', \
#]
#
#check_omi_ceres_time_offset(dir_list, num_points_per_file = 10)
#sys.exit()



"""
date_str = '20180705'
minlat = 65.
maxlat = 87.
MYD08_data = read_MODIS_MYD08_single(date_str, minlat = minlat, maxlat = maxlat)

plot_cod = MYD08_data['cod_mean']
print('min COD', np.min(plot_cod))

plot_cod = np.where(MYD08_data['cod_mean'].mask == True, 0, MYD08_data['cod_mean'])
print('min COD', np.min(plot_cod))

fig = plt.figure(figsize = (9, 4))
ax1 = fig.add_subplot(1,2,1, projection = ccrs.NorthPolarStereo())
ax2 = fig.add_subplot(1,2,2, projection = ccrs.NorthPolarStereo())

ax1.pcolormesh(MYD08_data['lon'], MYD08_data['lat'], plot_cod, shading = 'auto', vmax = 50, transform = ccrs.PlateCarree())
ax1.coastlines()
ax1.set_extent([-180,180,65,90], datacrs)

ax2.pcolormesh(MYD08_data['lon'], MYD08_data['lat'], MYD08_data['day_cld_frac_mean'], shading = 'auto', transform = ccrs.PlateCarree())
ax2.coastlines()
ax2.set_extent([-180,180,65,90], datacrs)

fig.tight_layout()
plt.show()

sys.exit()
"""

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#                      BEGIN FORCING VERSION NN STUFF
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#                        PLOTTING/ANALYSIS FUNCTIONS
#
#  select_idxs(test_dict, ai_min, ice_min, ice_max:
#           This function identifies the indices that give the NN data
#           that meet the AI, ICE, COD, and SZA requirements
#
#  combine_NN_data(sim_name):
#           This function selects all the neural network output files 
#           for the given simulation , while also applying some QC criteria,
#           and returns a dictionary containing the observed and NN-estimated
#           SWF, in addition to the binning variables.
#
#  calc_NN_force_slope_intcpt(test_dict, ice_bin_edges, ...):
#           This function calculates the slope and intercepts for the 
#           forcing vs AI values. NOTE: here, forcing is calculated as NN - OBS.
#           
#           NOTE: This is also where the error stuff can be done by 
#                 switching the slope found by Theil-Sen (res[0])
#                 to either the lower bound slope (res[2]) or 
#                 upper bound slope (res[3])
#
#  test_calculate_type_forcing_v4:
#           Calculate the daily observation-estimated forcing
#           given the NN output-derived forcing slopes 
#           (and intercepts, if desired)
#
#  calculate_type_forcing_v4_monthly: calculates the single-day forcing 
#            estimates and averages them into monthly averages.
#
#
# calc_forcing_slopes_v4_all_months_arctic_avg_manyrefice(all_month_vals):
#           Calculates the trends (original and modified) for either
#           the refice or refcld simulations. The results are returned
#           in a dictionary.
#
#
# calc_forcing_slopes_v4_all_months_arctic_avg_uncert(all_month_files):
#           Same as above, but for the specified uncertainty runs.
#
#  plot_NN_architecture(save = False):
#           Plots a graphic showing the NN architecture
#
#  plot_NN_scatter_multiCOD(test_dict, cod_bin_edges,...):
#           This function plots the NN output for given COD bin edges
#           and for the specified AI, ICE, and SZA min/max values. For
#           each COD bin, the selected, scattered forcing 
#           (calculated as NN - OBS) data are plotted as a function of AI.
#           There is one plot for each COD bin.
#
#  plot_NN_bin_slopes(slope_dict):
#           Plot the binned NN/AI slopes for the 4 surface types as well
#           as the bin counts
#
#
#  plot_NN_bin_slopes_codmean(slope_dict, bin_dict, min_ob = 50, save = False):
#           Plots the NN/AI slopes for the 4 surface types averaged ALONG
#           the COD bins. Thus, the meaned forcing slopes are functions
#           of the COD.
#    
#  plot_NN_forcing_daily(...):
#           This function calculates the version 4 (NN) estimated forcing
#           for a single day based on the forcing slopes (and bins)
#           as well as the daily OMI AI, MODIS COD, NSIDC ICE, and
#           the calculated solar zenith angle for the day. 
#
#  plot_type_forcing_v4_all_months_arctic_avg(all_month_files, \)
#           Calculate Arctic-wide average for all of the provided monthly forcings.
#           Meant to be a plot of the Monte Carlo-esque uncertainty analysis.
#           The first file in the list is the control
#
# plot_type_forcing_v3_all_months_arctic_avg_combined(all_month_values, \)
#           Plots three Arctic-averaged monthly results in one image:
#           The first row is the 65 - 87, the second row is 65 - 75, and
#           the third row is 75 - 87.
#
# plot_compare_NN_output(sys.argv[1], save = False):
#           Makes a plot of Observed SWF, NN SWF, and Regression
#           Designed for showing the validation of the NN under
#           aerosol-free conditions
# 
# plot_compare_NN_output_v2(sys.argv[1], auto_zoom = False, save = False):
#           Plots OMI UVAI, MODIS true color, CERES obs, NN output,
#           and forcing (NN - obs) for either the whole Arctic
#           (auto_zoom = False) or just the plume area (auto_zoom = True)
# 
# plot_compare_NN_output_noaer(sys.argv[1], save = False):
#           Makes a plot of Observed SWF, NN SWF, and Regression
#           Designed for showing the validation of the NN under
#           aerosol-free conditions
# 
# plot_compare_NN_output_overlay(sys.argv[1], auto_zoom = True, save = False)
#           Plots zoomed in images of OMI AI, MODIS true color, and 
#           forcing, with hatched regions showing AI greater than xxx.
#           Function written to determine if there was significant drift
#           between the Aura and Aqua overpasses when viewing a smoke plume.
#           
# compare_nn_version_output(date_str, skip_version = skip_version, save = False)
#           This function reads in all NN output files for a given date_str,
#           calculates the mean and standard deviations of both the nn SWF
#           and the forcing (NN - obs), and plots the CERES obs, 
#           NN mean, and standard deviation of the NN values at each
#           swath point.
# 
#                        L2/L3 ERROR FUNCTIONS
#
# plot_L2_validate_regress_all(sim_name, slope_dict_lin, bin_dict,...
#           This function grabs all the NN output files for the 
#           given simualtion name (i.e. noland50), calculates L2-style
#           (NN - CERES) and L3-style (regressions & AI) forcings for
#           each aerosol-containing point in each swath and combines
#           them into arrays. These are also plotted using...
#
# plot_scatter_hist_L2_L3_errors(direct_forcings, calc_forcings,...
#           This function takes the direct (L2-style) and calculated 
#           (L3-style) forcings from plot_L2_validate_regress_all and
#           generates a scatter plot of the two, with error bars, and
#           a histogram of all the errors. 
#
#plot_NN_forcing_daily_L2L3_errors(date_str, daily_VSJ4, \
#    OMI_daily_VSJ4, slope_dict_lin, bin_dict, L2L3_err_mean, L2L3_err_std, \
#           This funtion calculates daily forcing values originally, and then
#           calculates (num_calcs) additional daily forcing values in which
#           each grid point with a forcing value calculated has L2L3 error
#           added that fits in a Gaussian distribution with the mean and 
#           standard deviation given above.
#
#                   STEPS FOR NEW DAILY FORCING ESTIMATION
#
#
# 1. Grab all the NN data (combine_NN_data)
# 2. Calculate forcing slopes and intercepts on the COD, SZA, and ICE bins
# 3. For each day, loop over each grid box
# 4. Use the COD, SZA, and ICE values to select the matching binned forcing 
#    slope and intercept. NOTE: This is where the "refice" and "refcld"
#    simulations can be easily redone.
# 5. Use the regression slope/intercept values to compute the forcing
#    given the daily AI at that grid box.
#
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# check_omi_ceres_time_offset(dir_list, num_points_per_file = 10):
#           This function grabs the OMI and CERES HDF5 files for given DTGs,
#           grabs 10 (or other specified number) of OMI points, finds the closest
#           CERES grid point to that OMI point, and figures out how many minutes
#           are between the OMI and CERES observations at that point.
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Load in all the NN-estimated SWF values (and other values) for
# simulation 'noland48'
# ----------------------------------------------------------------
#test_dict = combine_NN_data('noland48')
#test_dict = combine_NN_data('noland50')
#print("AS OF 2024/08/27, USING NOLAND72")
#print("AS OF 2024/09/08, USING NOLAND50 AGAIN")
#sim_name = 'noland48'
#sim_name = 'noland50'
#sim_name = 'noland72'
#sim_name = 'noland73'
sim_name = 'noland74'
#sim_name = 'noland103'
#sim_name = 'noland75'
print("AS OF 2024/09/09, USING ", sim_name)
test_dict = combine_NN_data(sim_name)


#calc_pcnt_aerosol_over_type_v2(sim_name, 1.5, minlat = 70., dtype = 'PERT', ax = None)
#sys.exit()

# Identify the dictionary array indices that give the NN data
# associated with:
# - AI >= 1
# - 0 <= NSIDC ICE < 20
# - 2 <= MODIS COD < 4
# - 50 <= OMI SZA < 55
# -----------------------------------------------------------
#good_idxs = select_idxs(test_dict, 1, \
#    0, 20, 2, 4, 50, 55)

# Define COD, SZA, and NSIDC ICE bin ranges
# -----------------------------------------
cod_bin_edges = np.array([0,0.5,2,4,8,12,20,30,150])
cod_bin_means = (cod_bin_edges[1:] + cod_bin_edges[:-1]) / 2
sza_bin_edges = np.arange(40, 85, 5)
sza_bin_means = (sza_bin_edges[1:] + sza_bin_edges[:-1]) / 2
#ice_bin_edges = np.array([0, 20, 80, 100.2, 255])
ice_bin_edges = np.array([0, 20, 40, 60, 80, 100.2, 255])
ice_bin_means = (ice_bin_edges[1:] + ice_bin_edges[:-1]) / 2

bin_dict = {
    'cod_bin_edges': cod_bin_edges, \
    'cod_bin_means': cod_bin_means, \
    'sza_bin_edges': sza_bin_edges, \
    'sza_bin_means': sza_bin_means, \
    'ice_bin_edges': ice_bin_edges, \
    'ice_bin_means': ice_bin_means, \
}



# Calculate the regression slopes and intercepts from the NN data
# ---------------------------------------------------------------
min_ob = 50
ai_min_forslopes = 0.0
#slope_dict = calc_NN_force_slope_intcpt(test_dict, ice_bin_edges, \
#        sza_bin_edges, cod_bin_edges, ai_min = 0, min_ob = min_ob, \
#        trend_type = 'theil-sen')
slope_dict_lin = calc_NN_force_slope_intcpt(test_dict, ice_bin_edges, \
        sza_bin_edges, cod_bin_edges, ai_min = ai_min_forslopes, min_ob = min_ob, \
        trend_type = 'linregress')

# Plot the scattered NN output as function of AI for the bins given above
# -----------------------------------------------------------------------
#plot_NN_scatter_multiCOD(test_dict, cod_bin_edges, 0, 101, 255, 50, 55, save = False)
#sys.exit()

sza_min = 40
sza_max = 90
#plot_NN_scatter_combined_alltypes(test_dict, bin_dict, \
#    ai_min_forslopes, sza_min, sza_max, trend_type = 'linregress', \
#    show_specific_cod = 0, min_ai_for_stats = 2.0, \
#    plot_bounds = False, save = False)
#sys.exit()

# See what happens with differnet SZA bins. Is it necessary to bin
# by the SZA?
# -----------------------------------------------------------------
sfc_idx = 0
maxerr = 1.5
#compare_sza_bin_impact_on_slopes(test_dict, bin_dict, sfc_idx, ai_min = ai_min, \
#    min_ob = min_ob, maxerr = maxerr, trend_type = 'linregress', combined_plot = False)
#sys.exit()
#slope_dict2 = calc_NN_force_slope_intcpt(test_dict, ice_bin_edges, \
#        sza_bin_edges, cod_bin_edges, ai_min = 1, min_ob = min_ob, \
#        trend_type = 'theil-sen')

# Plot the binned NN/AI slopes for the 4 surface types
# ----------------------------------------------------
#plot_NN_bin_slopes(slope_dict_lin, bin_dict, min_ob = min_ob, plot_error = True, save = False)
##plot_NN_bin_slopes(slope_dict_lin, bin_dict, min_ob = min_ob, plot_error = True, save = True)
#sys.exit()

# Plot the binned NN/AI slopes for the 4 surface types
# ----------------------------------------------------
#plot_NN_bin_slopes_6types(slope_dict_lin, bin_dict, 'slopes', min_ob = 50, \
#            plot_error = False, save = False)
#sys.exit()

# Run L2 validation
# -----------------
#plot_compare_NN_output_L2_validate(sys.argv[1], slope_dict_lin, bin_dict, \
#    auto_zoom = True, save = False)
#sys.exit()

# Makes a plot of OMI AI, Observed SWF, NN SWF, and Regression
# Designed for showing the validation of the NN under
# aerosol-free conditions
# -----------------------------------------------------------
#plot_compare_NN_output_noaer(sys.argv[1], save = True)

# Plots more stuff than the NN_output_noaer function
# --------------------------------------------------
#plot_compare_NN_output(sys.argv[1], auto_zoom = True, save = False)
#sys.exit()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Calculate L2- and L3-style forcing estimates for all high AI colocation 
# pixels
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#####sim_name = 'noland72'
#infile = 'validate_values_' + sim_name + '.hdf5'
infile = 'validate_values_' + sim_name + '_useintcptTrue.hdf5'
in_data = h5py.File(infile)
direct_forcings = in_data['direct_forcings'][:]
calc_forcings   = in_data['calc_forcings'][:]

#direct_forcings, calc_forcings = \
#    plot_L2_validate_regress_all(sim_name, slope_dict_lin, bin_dict, \
#    ai_thresh = 0.7, mod_slopes = None, mod_intercepts = None, \
#    mod_cod = None, mod_ice = None, use_intercept = True, \
#    min_cod = None, max_cod = None,\
#    min_sza = None, max_sza = None, \
#    min_ice = None, max_ice = None, \
#    save = False, return_values = True, save_values = False)
###write_L2_L3_validation_values(direct_forcings, calc_forcings, sim_name, ai_min, bin_dict)
#plot_scatter_hist_L2_L3_errors(direct_forcings, calc_forcings, \
#    sim_name, num_bins = 200, delta_calc = 20, astrofit = True, \
#    screen_outliers = False, save = False)
#plot_hist_L2_L3_errors(direct_forcings, calc_forcings, sim_name, \
#    xmin = -100, xmax = 100, ax = None, num_bins = 200, astrofit = True, save = False)
#sys.exit()


## From np.mean and np.std on the errors
## -------------------------------------
#errors = calc_forcings - direct_forcings
#mean_error = np.mean(errors) #  3.72
#std_error  = np.std(errors)  # 29.04

## From astropy fit to error histogram.
## -----------------------------------
## noland72?
#mean_error = -3.321
#std_error  = 22.969
## noland50, 6 surface types
#mean_error = -2.571
#std_error  = 22.167


#files = glob('neuralnet_output/test_calc_out_' + sim_name + '*.hdf5')
#
#for tfile in files:
#    #plot_compare_NN_output(tfile, auto_zoom = True, save = True)
#    #plot_compare_NN_output_overlay(tfile, auto_zoom = True, save = True)
#
#    filedate = tfile.strip().split('/')[-1].split('_')[-1][:12]
#    year = int(filedate[:8])
#
#    print(year)
#    if( (year >= 20170708) & (year < 20190810) ):
#
#        plot_compare_NN_output_L2_validate(tfile, slope_dict_lin, bin_dict, \
#            mod_slopes = None, mod_intercepts = None, mod_cod = None, \
#            mod_ice = None, use_intercept = True, auto_zoom = True, \
#            ai_thresh = 0.7, save = True)
#
#        cmnd = 'mv *validate*.png validation_L2_images/'
#        print(cmnd)
#        os.system(cmnd)
#
#sys.exit()




#sys.exit()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# TEST DAILY FORCING ESTIMATION 
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Read in the daily NSIDC, MODIS, and OMI data

begin_date = '200504'
end_date   = '202009'
season     = 'sunlight'
minlat = 65
maxlat = 87
shawn_file = home_dir + '/Research/OMI/omi_shawn_daily_2005_2020_v2.hdf5'
#jz_file    = home_dir + '/Research/OMI/omi_VJZ211_daily_2005_2020.hdf5'

# Is this the "clean-sky" background version?
OMI_daily_VSJ4  = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -20.00, max_AI = 0.7, minlat = minlat, maxlat = maxlat)
# Is this the OMI AI monthly trend version?
#OMI_daily_VSJ4  = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -0.10, max_AI = 20.0, minlat = minlat, maxlat = maxlat)

#OMI_daily_VJZ211 = calcOMI_MonthAvg_FromDaily(jz_file, min_AI = -0.10, minlat = minlat, maxlat = maxlat)
daily_VSJ4 = readOMI_daily_HDF5(shawn_file, minlat = minlat, maxlat = maxlat)

##!#force_vals = test_calculate_type_forcing_v4(daily_VSJ4, OMI_daily_VSJ4, slope_dict_lin, \
##!#        bin_dict, '20170822', minlat = minlat, maxlat = maxlat, ai_thresh = 0.7, \
##!#        maxerr = 2, 
##!#        reference_ice = None, reference_cld = None, mod_slopes = None, \
##!#        mod_intercepts = None, \
##!#        mod_ice = None, ice_err_mean = None, ice_err_std = None, \
##!#        mod_cod = None, cod_err_mean = None, cod_err_std = None, \
##!#        mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
##!#        filter_bad_vals = False, return_modis_nsidc = False, use_intercept = True, \
##!#        debug = False)
##!#
##!#fig = plt.figure()
##!#ax = fig.add_subplot(1,1,1, projection = ccrs.NorthPolarStereo())
##!#ax.pcolormesh(daily_VSJ4['lon_values'], daily_VSJ4['lat_values'], force_vals, vmin = -60, vmax = 60, \
##!#    cmap = 'bwr', transform = ccrs.PlateCarree(), shading = 'auto')
##!#ax.coastlines()
##!#ax.set_extent([-180,180,65,90], ccrs.PlateCarree())
##!#plt.show()
##!#sys.exit()

# Plot spatial averages of OMI data
# ---------------------------------
#plot_grid_OMI_trends_spatial(shawn_file, min_AI = 0.000, \
#    max_AI = 20.0, minlat = 65.5, maxlat = 87.5,  save = False)
#sys.exit()
#plot_grid_OMI_climo_spatial(shawn_file, min_AI = -0.10, \
#    max_AI = None, minlat = 65.5, maxlat = 90.5,  save = False)

# Plot trends in Arctic averaged OMI AI data
# ------------------------------------------
#plot_OMI_trends_arctic_avg(shawn_file, min_AI = -0.10, minlat = 65.5, maxlat = 87.5, save = False);
#plot_grid_OMI_trends_arctic_avg(OMI_daily_VSJ4, min_AI = -0.1, \
#    max_AI = 20.0, minlat = 65.5, maxlat = 90.5,  flat_axs = None, \
#    save = False)
#plot_grid_OMI_trends_arctic_avg_combined(shawn_file, min_AI = -0.1, \
#    max_AI = 20.0, save = False)
#sys.exit()


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Figure out which days have high AI
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#threshold = 1.5
#for ii in range(daily_VSJ4['day_values'].shape[0]):
#    max_daily = np.max(daily_VSJ4['grid_AI'][ii,:,:])
#    if(max_daily > threshold):
#        print(daily_VSJ4['day_values'][ii])
#sys.exit()


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plots daily AI, the clear-sky background, the departure from the background,
# and the histogram of the departure values
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#plot_daily_OMI(daily_VSJ4, OMI_daily_VSJ4, 20050611)

#ai_thresh = 0.05
#ai_thresh = -0.15
#maxerr = 1.5

"""
date_str = '20050403'
ai_thresh = 0.7
values_false = test_calculate_type_forcing_v4(daily_VSJ4, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, date_str, minlat = minlat, maxlat = maxlat, \
    ai_thresh = ai_thresh, maxerr = maxerr,\
    filter_bad_vals = True, \
    reference_ice = None, \
    reference_cld = None, \
    use_intercept = True, 
    mod_slopes = None, \
    return_modis_nsidc = False, debug = False)
values_true = test_calculate_type_forcing_v4(daily_VSJ4, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, date_str, minlat = minlat, maxlat = maxlat, \
    ai_thresh = ai_thresh, maxerr = maxerr,\
    filter_bad_vals = True, \
    reference_ice = None, \
    reference_cld = None, \
    use_intercept = True, 
    mod_slopes = None, \
    return_modis_nsidc = False, debug = False)

pmax = np.max([np.max(values_false), np.max(values_true)])
pmin = np.min([np.min(values_false), np.min(values_true)])

vmin = np.min([pmin, -abs(pmax)])
vmax = np.max([abs(pmin), pmax])

plt.close('all')
fig = plt.figure()
ax1 = fig.add_subplot(1,2,1, projection = ccrs.NorthPolarStereo())
ax2 = fig.add_subplot(1,2,2, projection = ccrs.NorthPolarStereo())
mesh = ax1.pcolormesh(OMI_daily_VSJ4['LON'][:,:], OMI_daily_VSJ4['LAT'][:,:], \
    values_false, transform = ccrs.PlateCarree(), shading = 'auto', \
    cmap = 'bwr', vmin = vmin, vmax = vmax)
cbar = fig.colorbar(mesh, ax = ax1)
ax1.coastlines()
ax1.set_extent([-180,180,65,90], ccrs.PlateCarree())
ax1.set_title('With bad vals')
mesh = ax2.pcolormesh(OMI_daily_VSJ4['LON'][:,:], OMI_daily_VSJ4['LAT'][:,:], \
    values_true, transform = ccrs.PlateCarree(), shading = 'auto', \
    cmap = 'bwr', vmin = vmin, vmax = vmax)
cbar = fig.colorbar(mesh, ax = ax2)
ax2.coastlines()
ax2.set_extent([-180,180,65,90], ccrs.PlateCarree())
ax2.set_title('Without bad vals')
fig.tight_layout()
plt.show()
sys.exit()
"""


# Pre-calculate the daily SZA values for 1 April to 30 September
# at each latitude in the grid. Resulting shape is (num_days, num_latitudes)
# --------------------------------------------------------------
#apr01_idx = int(datetime(2005,4,1).strftime('%j'))
#sep30_idx = int(datetime(2005,9,30).strftime('%j')) + 1
#days = np.arange(1, 366)
#days = np.arange(apr01_idx, sep30_idx)
#del_angles = -23.45 * np.cos( np.radians((360 / 365) * (days + 10)))
#latitudes = np.arange(minlat + 0.5, maxlat + 0.5, 1.0)
#lat_szas = np.array([tlat - del_angles for tlat in latitudes]).T

#  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Calculate the daily forcing estimate values for a date_string
#date_str = '20180705'
#plot_NN_forcing_daily(date_str, daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, minlat = 65., maxlat = 87., \
#    ai_thresh = 0.7, maxerr = maxerr, filter_bad_vals = True, \
#    save = False, use_intercept = True)
#sys.exit()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Test the L2L3_err stuff for a single day
# 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#date_str = '20180705'
#L2L3_err_mean = -3.321
#L2L3_err_std  = 22.969
#num_calcs = 10
#
#
## This funtion calculates daily forcing values originally, and then
## calculates (num_calcs) additional daily forcing values in which
## each grid point with a forcing value calculated has L2L3 error
## added that fits in a Gaussian distribution with the mean and 
## standard deviation given above.
## ------------------------------------------------------------------
#error_type = 'L2L3'
#plot_NN_forcing_daily_L2L3_errors(date_str, daily_VSJ4, \
#    OMI_daily_VSJ4, slope_dict_lin, bin_dict, L2L3_err_mean, L2L3_err_std, error_type, \
#    num_calcs, minlat = 65., maxlat = 87., \
#    use_intercept = True, filter_bad_vals = False, \
#    mod_L2_L3_error = None, \
#    ai_thresh = 0.7, maxerr = maxerr, save = True)



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Test the L2L3_err stuff for a whole study period
#
# NOTE: RUNNING WITH filter_bad_vals SET TO False
# 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
ai_thresh = 0.7
maxerr = 1.5

# Control daily values
# --------------------
#daily_filename = 'arctic_daily_est_forcing_v1.hdf5'        
daily_filename = 'arctic_daily_est_forcing_numsfcbins6.hdf5' # noland72
daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v1.hdf5' # noland50
daily_filename = 'arctic_daily_est_forcing_numsfcbins4.hdf5' # noland50
daily_filename = 'arctic_daily_est_forcing_numsfcbins4_v1.hdf5' # noland72
daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v2.hdf5' # noland74
# Daily values with ice modifiations
# ----------------------------------
#ice_filename = 'arctic_daily_est_forcing_iceerr_v1.hdf5'
ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr.hdf5' # noland72
ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_v1.hdf5' # noland50
ice_filename = 'arctic_daily_est_forcing_numsfcbins4_iceerr.hdf5' # noland50
ice_filename = 'arctic_daily_est_forcing_numsfcbins4_iceerr_v1.hdf5' # noland72
ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_v2.hdf5' # noland74
# Daily values with COD modifiations
# ----------------------------------
#cod_filename = 'arctic_daily_est_forcing_coderr.hdf5'
cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr.hdf5' # std = 5, noland72
cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_v1.hdf5' # std = 5, noland50
cod_filename = 'arctic_daily_est_forcing_numsfcbins4_coderr.hdf5' # std = 5, noland50
cod_filename = 'arctic_daily_est_forcing_numsfcbins4_coderr_v1.hdf5' # std = 5, noland72
cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_v2.hdf5' # std = 5, noland74

#daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)
#ice_dict   = read_daily_month_force_L2L3_error_from_HDF5(ice_filename)
#cod_dict   = read_daily_month_force_L2L3_error_from_HDF5(cod_filename)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Functions for 4-panel plot with all analyzed errors
# - NN - obs
# - L2 - L3
# - Ice errors
# - COD errors
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#num_bins = 200
#plot_error_components_combined(direct_forcings, calc_forcings, sim_name, \
#    daily_filename, ice_filename, cod_filename, num_bins, \
#    astrofit = True, log_scale = True, save = False)
#sys.exit()


# Using 6 sfc bins, 
# ice errors: 
#   Applied ice errors: mean = 0, std dev = 15
#   Resulting forcing errors: mean = 0.011 W/m2, std dev = 3.11 W/m2
# cod errors: 
#   Applied ice errors: mean = 0, std dev = 5
#   Resulting forcing errors: mean = 0.848 W/m2, std dev = 15.9 W/m2
# L2L3 errors:
##    noland50
#       mean error (astrofit) = -2.57
#       stdv error (astrofit) = 22.2
#
# Combined error = sqrt(

# = = = = = =
#
# 2024/09/08
#
# USING 6 SURFACE BINS AND NOLAND50,
#
# Impacts of L3 ice errors on L3 forcing, from histogram
#   Applied ice errors: mean = 0, std dev = 15
#   Resulting forcing errors: mean = 0.031 W/m2, std dev = 3.24 W/m2
# Impacts of L3 COD errors on L3 forcing, from histogram
#   Applied COD errors: mean = 0, std dev = 5
#   Resulting forcing errors: mean = 1.012 W/m2, std dev = 14.374 W/m2
# L2L3 error numbers, from histogram
#   Mean L2L3 error, with astrofit = -2.57 
#   STDV L2L3 error, with astrofit = 22.2 
#
# Mean combined error = ice_error_mean + COD_error_mean + L2L3_error_mean = -1.527
# STD  combined error = sqrt(ice_error_std**2 + COD_error_std**2 + L2L3_error**2) = 26.64
# = = = = = =

# = = = = = =
#
# 2024/09/09
#
# USING 6 SURFACE BINS AND NOLAND74,
# 
# NN vs Obs errors, from histogram
#   Mean NN-ob error, with astrofit = 0.4
#   StDv NN-ob error, with astrofit = 17.3
# L2L3 error numbers, from histogram
#   Mean L2L3 error, with astrofit = -3.74 
#   STDV L2L3 error, with astrofit = 22.2 
# Impacts of L3 ice errors on L3 forcing, from histogram
#   Applied ice errors: mean = 0, std dev = 15
#   Resulting forcing errors: mean = 0.002 W/m2, std dev = 3.16 W/m2
# Impacts of L3 COD errors on L3 forcing, from histogram
#   Applied COD errors: mean = 0, std dev = 5
#   Resulting forcing errors: mean = 0.638 W/m2, std dev = 14.51 W/m2
#
# Mean combined error = ice_error_mean + COD_error_mean + L2L3_error_mean = -3.1
# STD  combined error = sqrt(ice_error_std**2 + COD_error_std**2 + L2L3_error**2) = 26.70
#
# With NN-ob errors,
# Mean combined error = nn_ob_error_mean + ice_error_mean + COD_error_mean + L2L3_error_mean = -2.7
# STD  combined error = sqrt(nn_ob_error_mean**2 + ice_error_std**2 + COD_error_std**2 + L2L3_error**2) = 31.82
#
# = = = = = =


# = = = = =
#
# TEST BULK ERROR ADDITIONS
#
# = = = = =




total_err_mean = -2.7
total_err_std = 31.82
num_sims = 300
#num_sims = 300
minlat = 65.5
maxlat = 87.5
sim_values = None

daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)

# = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = =

read_force_sim_vals = True
save_force_vals = False    
read_trend_sim_vals = False
save_trend_vals = False    

calc_region_avg_force_vals = True
calc_force_trends_from_file = False

if(read_trend_sim_vals):
   
    # NOTE: Only using 2 files here. Can change ":2" to allow it to read more files 
    # -----------------------------------------------------------------------------
    file_start = 'arctic_monthly_force_trends_count'
    force_files = glob(file_start + '*.hdf5')[:2]

    # Figure out how many simulations are in all the files
    # ----------------------------------------------------
    total_sims = 0
    for ffile in force_files:
        local_numsim = ffile.strip().split('/')[-1].split('_')[4].split('_')[0].split('.')[0][5:]   
        total_sims += int(local_numsim)

    print("Total number of sims:", total_sims)

    # Create an array to hold the values (Assume 96 months)
    # -----------------------------------------------------
    forcing_trends = np.full( (total_sims, 6, daily_dict['latitude'].shape[0], daily_dict['longitude'].shape[0]), np.nan)

    # Insert the values into the array
    # --------------------------------
    beg_idx = 0
    end_idx = 0
    for ii in range(len(force_files)):
        local_numsim = int(force_files[ii].strip().split('/')[-1].split('_')[4].split('_')[0].split('.')[0][5:])
        data = h5py.File(force_files[ii])
        end_idx = beg_idx + local_numsim
        forcing_trends[beg_idx:end_idx,:,:,:] = data['force_trends'][:,:,:,:]
        data.close()
        print(beg_idx, end_idx)
        beg_idx = end_idx
    
else:
    if(read_force_sim_vals):
        
        # Grab the filenames
        # NOTE: Only using 2 files here. Can change ":2" to allow it to read more files 
        # -----------------------------------------------------------------------------
        #file_start = 'arctic_monthly_force_values_count300'
        #force_files = glob(file_start +'*.hdf5')[:2]
        # NOTE: These files are for noland74, with old error

        force_files = ['arctic_monthly_force_values_count300.hdf5', \
                       'arctic_monthly_force_values_count300_v1.hdf5', \
                       'arctic_monthly_force_values_count300_v2.hdf5', \
                       'arctic_monthly_force_values_count300_v3.hdf5', \
                       'arctic_monthly_force_values_count300_v4.hdf5']

        force_files = force_files[:2]   
 
        # Figure out how many simulations are in all the files
        # ----------------------------------------------------
        total_sims = 0
        for ffile in force_files:
            print(ffile)
            local_numsim = ffile.strip().split('/')[-1].split('_')[4].split('_')[0].split('.')[0][5:]   
            total_sims += int(local_numsim)
    
        print("Total number of sims:", total_sims)
   
        if(calc_region_avg_force_vals): 
            # Create an array to hold the region-averaged values 
            # (Assume 96 months & 3 regions (65 - 90, 65 - 75, and 75 - 90)
            # -------------------------------------------------------------
            sim_values = np.full( (total_sims, 3, 96), np.nan)
        else:
            # Create an array to hold the values (Assume 96 months)
            # -----------------------------------------------------
            sim_values = np.full( (total_sims, 96, daily_dict['latitude'].shape[0], daily_dict['longitude'].shape[0]), np.nan)
    
        # Insert the values into the array
        # --------------------------------
        beg_idx = 0
        end_idx = 0
        for ii in range(len(force_files)):
            local_numsim = int(force_files[ii].strip().split('/')[-1].split('_')[4].split('_')[0].split('.')[0][5:])
            data = h5py.File(force_files[ii])
            end_idx = beg_idx + local_numsim
            
            if(calc_region_avg_force_vals):
                lat_mins = [65.5, 65.5, 75.5]
                lat_maxs = [89.5, 75.5, 89.5]
                for jj in range(len(lat_mins)):
                    lat_idxs = np.where( (data['latitude'][:] >= lat_mins[jj]) & (data['latitude'][:] < lat_maxs[jj])) 
                    lat_beg_idx = lat_idxs[0][0]
                    lat_end_idx = lat_idxs[0][-1] + 1
                    sim_values[beg_idx:end_idx,jj,:] = \
                        np.nanmean(data['monthly_force_vals'][:,:,lat_beg_idx:lat_end_idx,:], axis = (2,3))
            else:
                sim_values[beg_idx:end_idx,:,:,:] = data['monthly_force_vals'][:,:,:,:]

            data.close()
            print(beg_idx, end_idx)
            beg_idx = end_idx
        
    
    else:
        #for nn in range(2):
            
        return_std = False
        if(return_std):
            sim_values, std_values = calc_force_vals_bulk(daily_dict, total_err_mean, total_err_std, \
                minlat = minlat, maxlat = maxlat, num_sims = num_sims, return_std = True)
        else:
            sim_values = calc_force_vals_bulk(daily_dict, total_err_mean, total_err_std, \
                minlat = minlat, maxlat = maxlat, num_sims = num_sims, return_std = False)
            #sim_values_2 = calc_force_vals_bulk(daily_dict, total_err_mean, total_err_std, \
            #    minlat = minlat, maxlat = maxlat, num_sims = num_sims, return_std = False)
    
        if(save_force_vals):
            # Save error monthly forcing values to an output file
            # ---------------------------------------------------
            write_monthly_force_vals_sims_to_HDF5(daily_dict, sim_values, \
                total_err_mean, total_err_std, save_path = './', name_add = '')
    
    
        # Test calculating trends across all calculations and across the whole grid
        # -------------------------------------------------------------------------
        forcing_trends = np.full( (sim_values.shape[0], 6, sim_values.shape[2], sim_values.shape[3]), np.nan)
        #forcing_pvals  = np.full( (sim_values.shape[0], 6, sim_values.shape[2], sim_values.shape[3]), np.nan)
        #forcing_uncert = np.full( (sim_values.shape[0], 6, sim_values.shape[2], sim_values.shape[3]), np.nan)
        
        for ii in range(6):
            for jj in range(sim_values.shape[0]):
                print(ii,jj)
                forcing_trends[jj,ii,:,:], _, \
                    _ = calc_forcing_grid_trend(\
                    sim_values[jj,ii::6,:,:], 'standard')
                #forcing_trends[jj,ii,:,:], forcing_pvals[jj,ii,:,:], \
                #    forcing_uncert[jj,ii,:,:] = calc_forcing_grid_trend(\
                #    sim_values[jj,ii::6,:,:], 'standard')
        
        
        if(save_force_vals): 
            # Save error trends to an output file
            # -----------------------------------
            write_monthly_force_trend_sims_to_HDF5(daily_dict, forcing_trends, \
                total_err_mean, total_err_std, save_path = './', name_add = '')
        
        #del(sim_values)
        #del(forcing_trends) 

# Plot time series of the many region-averaged monthly forcing values
# for a given region index (0 = entire Arctic, 1 = low Arctic, 
# 2 = high Arctic
# -------------------------------------------------------------------
#plot_arctic_avg_region_trends(sim_values, 0)

calc_arctic_avg_region_trends(sim_values)
#trends, pvals = calc_arctic_avg_region_trends(sim_values)

sys.exit()

# Plot an 18-panel (6 row, 3 or 4 column) figure with:
# - AI trends
# - Mean of the forcing trends
# - Standard deviation of the forcing trends
# ----------------------------------------------------
#plot_bulk_force_AI_trend_v2(daily_dict, forcing_trends, shawn_file, \
#    vmax = 1.5, min_AI = 0.0, max_AI = 20.0, minlat = 65.5, \
#    maxlat = 90.5,  save = False)
#sys.exit()
   


"""
# Test calculating the standard deviations 
# ----------------------------------------
avg_std_per_sim_per_month = np.sqrt(np.sum(std_values[:,:,:,:]**2., axis = (2,3)) / std_values[:,:,:,:][0,0,:,:].flatten().shape[0])
avg_mean_per_sim_per_month = np.nanmean(sim_values, axis = (2, 3))

avg_mean_across_sims_aug = np.nanmean(avg_mean_per_sim_per_month[:,4::6], axis = 0)
avg_std_across_sims_aug = np.sqrt(np.sum(avg_std_per_sim_per_month[:,4::6]**2., axis = (0)) / avg_std_per_sim_per_month[:,4::6].shape[0])
"""

# Calculate confidence intervals for each trend
# ---------------------------------------------
# In the lib file,  from scipy.stats import norm as statnorm
# conf_intvl = statnorm.interval(alpha = 0.90, loc = np.mean(forcing_trends[:,4,2,280]), scale = st.sem(forcing_trends[:,4,2,280]))

#month_idx = 3
#lat_idx = 3
#lon_idx = 341
#plot_force_trend_mean_std_dist(daily_dict, forcing_trends, month_idx, lat_idx, lon_idx, \
#    vmin = -1.5, vmax = 1.5, conf_window = 90, save = False)
#sys.exit()


#plot_test_trends_stdevs(daily_dict, forcing_trends, meanmax = 3.0, stdmax = 1.5)

# Plot the spatial trend averages or trend standard deviations at each
# grid point
# ---------------------------------------------------------------------
plot_bulk_sim_trends(daily_dict, forcing_trends, 'mean', vmax = 1.5, \
    save = False)

sys.exit()

# Plot the distribution of trend estimates at a lat/lon idx and month
# -------------------------------------------------------------------
test_error_dist(daily_dict, forcing_trends, 4, 2, 285, 20)
sys.exit()

# Plot the distribution of mean trends across the Arctic for each month
# ---------------------------------------------------------------------
plot_bulk_trend_dist_multimonth(daily_dict, forcing_trends, bins = 50, save = False, \
    log_yscale = True, plot_single_month = None)
sys.exit()


#sys.exit()

#plot_sim_errors_bulk_arctic_avg_combined(daily_filename, total_err_mean, total_err_std, \
#    minlat = 65.5, maxlat = 90.5, num_sims = num_sims, sim_values = sim_values, \
#    plot_result_min_max_range = True, trend_type = 'linregress', \
#    flat_axs = None)

# Plot a 24-panel (6 row, 3 or 4 column) figure with:
# - AI trends
# - Mean of the forcing trends
# - Standard deviation of the forcing trends
# - (Optional) histogram of the trend values at each day
#plot_bulk_force_AI_trend(daily_dict, forcing_trends, shawn_file, \
#    vmax = 1.0, min_AI = 0.0, max_AI = 20.0, minlat = 65.5, \
#    maxlat = 90.5,  save = False)


sys.exit()

minlat = 65.5
maxlat = 90.5
plot_sim_errors_bulk_arctic_avg(daily_filename, total_err_mean, \
    total_err_std, minlat = minlat, maxlat = maxlat, num_sims = num_sims, \
    sim_values = sim_values, plot_result_min_max_range = True, flat_axs = None)



plot_sim_errors_bulk_arctic_avg_combined(daily_filename, total_err_mean, total_err_std, \
    minlat = 65.5, maxlat = 90.5, num_sims = num_sims, sim_values = sim_values, \
    plot_result_min_max_range = True, trend_type = 'linregress', \
    flat_axs = None)

fig.tight_layout()
plt.show()
sys.exit()

#


# CONTROL RUN
# -----------
cod_err_mean = None
cod_err_std = None
ice_err_mean = None
ice_err_std = None
all_month_vals_orig_alldaily = \
    calculate_type_forcing_v4_alldaily(daily_VSJ4, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
    filter_bad_vals = False, return_modis_nsidc = False, \
    mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
    use_intercept = True, debug = False)

write_daily_month_force_L2L3_error_to_HDF5(\
    all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
    minlat = minlat, maxlat = maxlat, \
    maxerr = maxerr, ai_thresh = ai_thresh, \
    L2L3_err_mean = None, L2L3_err_std = None, \
    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
    dtype = None, \
    overwrite_old_file = False, \
    write_daily_values = True, \
    OMI_daily_data = daily_VSJ4, \
    save_path = './', name_add = '_numsfcbins6')

# ICE RUN
# -----------
cod_err_mean = None
cod_err_std = None
ice_err_mean = 0.
ice_err_std = 15.
all_month_vals_orig_alldaily = \
    calculate_type_forcing_v4_alldaily(daily_VSJ4, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
    filter_bad_vals = False, return_modis_nsidc = False, \
    mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
    use_intercept = True, debug = False)

write_daily_month_force_L2L3_error_to_HDF5(\
    all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
    minlat = minlat, maxlat = maxlat, \
    maxerr = maxerr, ai_thresh = ai_thresh, \
    L2L3_err_mean = None, L2L3_err_std = None, \
    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
    dtype = None, \
    overwrite_old_file = False, \
    write_daily_values = True, \
    OMI_daily_data = daily_VSJ4, \
    save_path = './', name_add = '_numsfcbins6_iceerr')

# COD RUN
# -------
cod_err_mean = 0
cod_err_std = 5
ice_err_mean = None
ice_err_std = None
all_month_vals_orig_alldaily = \
    calculate_type_forcing_v4_alldaily(daily_VSJ4, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
    filter_bad_vals = False, return_modis_nsidc = False, \
    mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
    use_intercept = True, debug = False)

write_daily_month_force_L2L3_error_to_HDF5(\
    all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
    minlat = minlat, maxlat = maxlat, \
    maxerr = maxerr, ai_thresh = ai_thresh, \
    L2L3_err_mean = None, L2L3_err_std = None, \
    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
    dtype = None, \
    overwrite_old_file = False, \
    write_daily_values = True, \
    OMI_daily_data = daily_VSJ4, \
    save_path = './', name_add = '_numsfcbins6_coderr')


sys.exit()

all_month_vals_orig_alldaily = \
    calculate_type_forcing_v4_alldaily(daily_VSJ4, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
    filter_bad_vals = False, return_modis_nsidc = False, \
    mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
    ice_err_mean = 0., ice_err_std = 15., \
    use_intercept = True, debug = False)

write_daily_month_force_L2L3_error_to_HDF5(\
    all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
    minlat = minlat, maxlat = maxlat, \
    maxerr = maxerr, ai_thresh = ai_thresh, \
    L2L3_err_mean = None, L2L3_err_std = None, \
    ice_err_mean = 0, ice_err_std = 15, \
    dtype = None, \
    overwrite_old_file = False, \
    write_daily_values = True, \
    OMI_daily_data = daily_VSJ4, \
    save_path = './', name_add = '_iceerr')


sys.exit()

write_daily_month_force_L2L3_error_to_HDF5(\
    all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
    minlat = minlat, maxlat = maxlat, \
    maxerr = maxerr, ai_thresh = ai_thresh, \
    L2L3_err_mean = None, L2L3_err_std = None, \
    ice_err_mean = 0, ice_err_std = 15, \
    dtype = None, \
    overwrite_old_file = False, \
    write_daily_values = True, \
    OMI_daily_data = daily_VSJ4, \
    save_path = './', name_add = '')

# Control daily values
# --------------------
daily_filename = 'arctic_daily_est_forcing.hdf5'        
# Daily values with ice modifiations
# ----------------------------------
daily_filename = 'arctic_daily_est_forcing_iceerr.hdf5'
daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)
#print_L2L3_sim_info(daily_dict)

sys.exit()
# Calculate original values for reference
all_month_vals_orig = \
    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
    filter_bad_vals = False, return_modis_nsidc = False, \
    mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
    use_intercept = True, debug = False)
#all_month_vals_orig_filterbadvals = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
#    use_intercept = True, debug = False)


# Set up an array to hold the new values
# --------------------------------------
num_calc = 50
all_month_vals_err_combined = np.full( (num_calc, \
    all_month_vals_orig.shape[0], all_month_vals_orig.shape[1], \
    all_month_vals_orig.shape[2]), np.nan)

for ii in range(num_calc):
    print("\n\n# = = = = = = = = = = = = \n#\nSIMULATION NUMBER ",ii,"\n#\n# = = = = = = = ")
    all_month_vals_err_combined[ii,:,:,:] = \
        calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
        slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
        ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
        filter_bad_vals = False, return_modis_nsidc = False, \
        mod_L2_L3_error = None, L2L3_err_mean = -3.321, L2L3_err_std = 22.969, \
        use_intercept = True, debug = False)

# NEED A READER FUNCTION TOO
#count10_errfile = 'arctic_month_est_forcing_L2L3err_count10.hdf5'
write_daily_month_force_L2L3_error_to_HDF5(\
    all_month_vals_err_combined, OMI_daily_VSJ4, \
    slope_dict_lin, bin_dict, all_month_vals_orig, \
    minlat = minlat, maxlat = maxlat, \
    maxerr = maxerr, ai_thresh = ai_thresh, \
    L2L3_err_mean = -3.321, L2L3_err_std = 22.969, \
    dtype = None, \
    save_path = './', name_add = '')

current_files = ['arctic_month_est_forcing_L2L3err_count10.hdf5', \
                 'arctic_month_est_forcing_L2L3err_count40.hdf5', \
                 'arctic_month_est_forcing_L2L3err_count50_v1.hdf5']


test_plot_L2L3_err_sims(current_files, save = False)
sys.exit()

fig = plt.figure()
axs = fig.subplots(nrows = 2, ncols = 3)
flat_axs = axs.flatten()
titles = ['April','May','June','July','August','September']
for tfile in current_files:
    data_dict = read_daily_month_force_L2L3_error_from_HDF5(tfile)

    arctic_avgs = np.nanmean(data_dict['force_estimate'], axis = (2, 3))

    if('force_estimate_orig' in data_dict.keys()):
        arctic_avgs_orig = np.nanmean(data_dict['force_estimate_orig'], axis = (1, 2))

    # Plot the error results
    # ----------------------
    for ii in range(len(flat_axs)):
        flat_axs[ii].plot(arctic_avgs[:,ii::6].T)

        if('force_estimate_orig' in data_dict.keys()):
            flat_axs[ii].plot(arctic_avgs_orig[ii::6], color = 'k')
    
flat_axs[0].set_title(titles[0]) 
flat_axs[1].set_title(titles[1]) 
flat_axs[2].set_title(titles[2]) 
flat_axs[3].set_title(titles[3]) 
flat_axs[4].set_title(titles[4]) 
flat_axs[5].set_title(titles[5]) 

sys.exit()







#values = test_calculate_type_forcing_v4(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict, bin_dict, date_str, minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr,\
#    filter_bad_vals = True, \
#    reference_ice = None, \
#    reference_cld = None, \
#    mod_slopes = None, \
#    return_modis_nsidc = False)

#sys.exit()

ai_thresh = 0.7
#all_month_vals = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, \
#    reference_ice = None, reference_cld = None, mod_slopes = None, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# AS OF 2024/08/22: NEED TO GO BACK AND RECALCULATE ALL OF THESE USING
#   THE NEW ICE BINS. THE CALCULATION CODE WILL NEED TO BE TWEAKED ACCORDINGLY
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

infile_aimin1 = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4.hdf5'
infile_aimin0 = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt.hdf5'
infile_aimin0_upper = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_upperslope.hdf5'
infile_aimin0_lower = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_lowerslope.hdf5'
infile_aimin0_lin = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_lin.hdf5'
infile_aimin0_linaddslopeerror = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_linaddslopeerror.hdf5'
infile_aimin0_linsubslopeerror = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_linsubslopeerror.hdf5'
infile_aimin0_linaddintcpterror = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_linaddintcpterror.hdf5'
infile_aimin0_linsubintcpterror = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_linsubintcpterror.hdf5'
infile_aimin0_linaddbotherror = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_linaddintcpterror_addslopeerror.hdf5'
infile_aimin0_linsubbotherror = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_linsubintcpterror_subslopeerror.hdf5'
infile_aimin0_linicep5   = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_liniceplus5.hdf5'
infile_aimin0_linicem5   = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_liniceminus5.hdf5'
infile_aimin0_linicep15  = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_liniceplus15.hdf5'
infile_aimin0_linicem15  = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_liniceminus15.hdf5'
infile_aimin0_lincodp5   = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_lincodplus5.hdf5'
infile_aimin0_lincodm5   = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_lincodminus5.hdf5'
infile_aimin0_lincodp2   = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_lincodplus2.hdf5'
infile_aimin0_lincodm2   = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_lincodminus2.hdf5'
infile_aimin0_linL2L3p30 = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_linaddL2L3error30.hdf5'
infile_aimin0_linL2L3m30 = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_linsubL2L3error30.hdf5'

#infile = infile_aimin0
infile = infile_aimin0_lin
all_month_dict = read_daily_month_force_HDF5(infile)
all_month_dict_thl = read_daily_month_force_HDF5(infile_aimin0)
#all_month_dict_upper = read_daily_month_force_HDF5(infile_aimin0_linicep15)
#all_month_dict_lower = read_daily_month_force_HDF5(infile_aimin0_linicem15)
#all_month_dict_upper = read_daily_month_force_HDF5(infile_aimin0_upper)
#all_month_dict_lower = read_daily_month_force_HDF5(infile_aimin0_lower)
#all_month_dict_upper = read_daily_month_force_HDF5(infile_aimin0_linaddslopeerror)
#all_month_dict_lower = read_daily_month_force_HDF5(infile_aimin0_linsubslopeerror)
#all_month_dict_upper = read_daily_month_force_HDF5(infile_aimin0_linaddintcpterror)
#all_month_dict_lower = read_daily_month_force_HDF5(infile_aimin0_linsubintcpterror)

# Compare the results for different sea ice and cloud values
# ----------------------------------------------------------
##plot_type_forcing_v3_all_months_arctic_avg_manyrefice(all_month_dict['FORCE_EST'], \
##       OMI_daily_VSJ4, minlat = 65, trend_type = 'standard', stype = 'ice', \
##       ptype = 'forcing', vtype = 'v4', version4 = True)
#plot_type_forcing_v3_all_months_arctic_avg_manyrefice(all_month_dict_thl['FORCE_EST'], \
#       OMI_daily_VSJ4, minlat = 65, maxlat = 75., trend_type = 'standard', stype = 'ice', \
#       ptype = 'forcing', vtype = 'v4', slope_type = 'thl', version4 = True)

#plot_type_forcing_v4_all_months_arctic_avg_manyrefice_combined(all_month_dict_thl['FORCE_EST'], \
#    OMI_daily_VSJ4, slope_type = 'thl', stype = 'ice', show_trends = False, \
#    version4 = False, horiz_orient = False)
#sys.exit()

# Calculates the slopes of the refice or refcld simulations plotted
# in the "plot_type_forcing_v3_all_months_arctic_avg_manyrefice" function.
# ------------------------------------------------------------------------
#ice_slopes = calc_forcing_slopes_v3_all_months_arctic_avg_manyrefice(all_month_dict['FORCE_EST'], \
#       OMI_daily_VSJ4, minlat = 65, trend_type = 'standard', stype = 'ice', \
#       ptype = 'forcing', version4 = True)

#ice_slopes = calc_forcing_slopes_v4_all_months_arctic_avg_manyrefice(all_month_dict['FORCE_EST'], \
#    OMI_daily_VSJ4, minlat = 65., maxlat = 87., stype = 'ice', slope_type = 'lin', \
#    save = False)
#cld_slopes = calc_forcing_slopes_v4_all_months_arctic_avg_manyrefice(all_month_dict['FORCE_EST'], \
#    OMI_daily_VSJ4, minlat = 65., maxlat = 87., stype = 'cld', slope_type = 'lin', \
#    save = False)

# Print the results of the refice/cld simulation trend comparisons as 
# a table
# -------------------------------------------------------------------
#calc_print_forcing_slope_error_v3(all_month_dict['FORCE_EST'], \
#       OMI_daily_VSJ4, minlat = 65, trend_type = 'standard', \
#       vtype = 'v3', ptype = 'forcing', version4 = True)
#
#sys.exit()

# Plot the trend of daily-gridded monthly forcing values
# ------------------------------------------------------
##plot_type_forcing_v3_all_months(all_month_vals, OMI_daily_VSJ4, \
#plot_type_forcing_v3_all_months(all_month_dict['FORCE_EST'], OMI_daily_VSJ4, \
#           minlat = 65, omi_data_type = 'lin', version4 = True)
#
# Calculate the Arctic-wide average of the daily-gridded monthly forcing
# values for each month and plot them
# ----------------------------------------------------------------------
#plot_type_forcing_v3_all_months_arctic_avg(all_month_vals, OMI_daily_VSJ4, \
#plot_type_forcing_v3_all_months_arctic_avg(all_month_dict['FORCE_EST'], OMI_daily_VSJ4, \
#    minlat = 65, trend_type = 'standard', omi_data_type = 'lin', version4 = True)
#sys.exit()
#write_daily_month_force_to_HDF5(all_month_vals, OMI_daily_VSJ4, \
#    name_add = '_dayaithresh07_aipert_dataminlat70')
#write_daily_month_force_to_HDF5(all_month_vals, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt')


# Calculate daily forcing values and average the daily values into
# monthly forcing estimates, but by adding or subtracting the 
# slope error from the SZA-mean AI/SWF slopes. 
# ----------------------------------------------------------------
#all_month_vals_upper = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = 'upper', \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_upper, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_upperslope')
#all_month_vals_lower = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = 'lower', \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_lower, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_lowerslope')
#sys.exit()

# NOTE: maxerr is not used in the v4 calculations as of 2024/07/03
#all_month_vals_upper_lin = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = 'upper', \
#    mod_intercepts = None, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_upper_lin, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_linadderror')
#all_month_vals_lower_lin = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = 'lower', \
#    mod_intercepts = None, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_lower_lin, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_linsuberror')
#all_month_vals_lin = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_intercepts = None, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_lin, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_lin')
#all_month_vals_upper_lin = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_intercepts = 'upper', \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_upper_lin, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_linaddintcpterror')
#all_month_vals_lower_lin = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_intercepts = 'lower', \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_lower_lin, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_linsubintcpterror')
#all_month_vals_upper_lin = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = 'upper', \
#    mod_intercepts = 'upper', \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_upper_lin, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_linaddintcpterror_addslopeerror')
#all_month_vals_lower_lin = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = 'lower', \
#    mod_intercepts = 'lower', \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_lower_lin, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_linsubintcpterror_subslopeerror')
#sys.exit()

# Calculate monthly forcing values with L2/L3 (L2_L3) errors 
# ---------------------------------------------------------
#all_month_vals_lower_lin = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_intercepts = None, mod_L2_L3_error = 30,\
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_lower_lin, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_linaddL2L3error30')
#all_month_vals_lower_lin = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_intercepts = None, mod_L2_L3_error = -30,\
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_lower_lin, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_linsubL2L3error30')
#sys.exit()


# Test the sensitivity to ice
# ---------------------------
#all_month_vals_ice_plus5 = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_ice = 5, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_ice_plus5, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_liniceplus5')
#all_month_vals_ice_minus5 = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_ice = -5, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_ice_minus5, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_liniceminus5')
#all_month_vals_ice_plus5 = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_cod = 5, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_ice_plus5, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_lincodplus5')
#all_month_vals_ice_minus5 = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_cod = -5, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_ice_minus5, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_lincodminus5')

#all_month_vals_ice_plus5 = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_cod = 2, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_ice_plus5, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_lincodplus2')
#all_month_vals_ice_minus5 = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_cod = -2, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_ice_minus5, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_lincodminus2')
#sys.exit()


#all_month_vals_ice_plus15 = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_ice = 15, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_ice_plus15, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_liniceplus15')
#all_month_vals_ice_minus15 = \
#    calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    mod_ice = -15, \
#    filter_bad_vals = True, return_modis_nsidc = False, \
#    use_intercept = True, debug = False)
#write_daily_month_force_to_HDF5(all_month_vals_ice_minus15, OMI_daily_VSJ4, \
#    maxerr = maxerr, ai_thresh = ai_thresh, minlat = minlat, 
#    dtype = 'pert', 
#    name_add = '_dayaithresh07_v4_aimin0_useintcpt_liniceminus15')
#
#sys.exit()

# Calculate daily forcing values and average the daily values into
# monthly forcing estimates, but here using the daily ice concentration
# values from 2005 as a reference. Am using this to try to see how
# the change in sea ice affects the aerosol forcing. May need to use
# the average of the first three years of sea ice rather than one
# year...
# ----------------------------------------------------------------
##infile_aimin0 = home_dir + '/Research/Arctic_compares/arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt.hdf5'
#ref_ice_vals = np.arange(2005,2021)
#for ref_ice in ref_ice_vals:
#    #all_month_vals_ice = calculate_type_forcing_v3_monthly(daily_VSJ4, \
#    #    OMI_daily_VSJ4, lin_smth2_dict_v6, 'all', ai_thresh = ai_thresh, minlat = minlat, \
#    #    maxerr = maxerr, reference_ice = str(ref_ice))
#    all_month_vals_ice = \
#        calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#        slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#        ai_thresh = ai_thresh, maxerr = maxerr, \
#        reference_ice = str(ref_ice), reference_cld = None, mod_slopes = None, \
#        filter_bad_vals = True, return_modis_nsidc = False, \
#        use_intercept = True, debug = False)
#    write_daily_month_force_to_HDF5(all_month_vals_ice, OMI_daily_VSJ4, \
#        name_add = '_dayaithresh07_v4_aimin0_useintcpt_lin_refice' + str(ref_ice))
#
#ref_cld_vals = np.arange(2005,2021)
#for ref_cld in ref_cld_vals:
#    #all_month_vals_cld = calculate_type_forcing_v3_monthly(daily_VSJ4, \
#    #    OMI_daily_VSJ4, lin_smth2_dict_v6, 'all', ai_thresh = ai_thresh, minlat = minlat, \
#    #    maxerr = maxerr, reference_cld = str(ref_cld))
#    all_month_vals_cld = \
#        calculate_type_forcing_v4_monthly(daily_VSJ4, OMI_daily_VSJ4, \
#        slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#        ai_thresh = ai_thresh, maxerr = maxerr, \
#        reference_ice = None, reference_cld = str(ref_cld), mod_slopes = None, \
#        filter_bad_vals = True, return_modis_nsidc = False, \
#        use_intercept = True, debug = False)
#    write_daily_month_force_to_HDF5(all_month_vals_cld, OMI_daily_VSJ4, \
#        name_add = '_dayaithresh07_v4_aimin0_useintcpt_lin_refcld' + str(ref_cld))
#sys.exit()
    
# Calculate the Arctic-wide average of the daily-gridded monthly forcing
# values for each month and plot them, but also plotting the '_adderror'
# and '_suberror' results for the first look at an error analysis.
# ----------------------------------------------------------------------
#plot_type_forcing_v3_all_months_arctic_avg(all_month_dict['FORCE_EST'], 
#    OMI_daily_VSJ4, minlat = 65, trend_type = 'standard', 
#    month_values2 = all_month_dict_upper['FORCE_EST'], 
#    month_values3 = all_month_dict_lower['FORCE_EST'], \
#    omi_data_type = 'pert', labels = ['upper', 'lower'])

# Plot the 65 - 87, 65 - 75, and 75 - 87 Arctic-avg results
# ----------------------------------------------------------------------
#plot_type_forcing_v3_all_months_arctic_avg_combined(all_month_dict['FORCE_EST'], \
#    OMI_daily_VSJ4, version4 = False, max_pval = 0.05, save = False)
#sys.exit()

# Calculate Arctic-wide average for all of the provided monthly forcings.
# Meant to be a plot of the Monte Carlo-esque uncertainty analysis.
# The first file in the list is the control
# -----------------------------------------------------------------------
all_month_files = [
    infile_aimin0_lin, \
    infile_aimin0_linaddslopeerror, \
    infile_aimin0_linsubslopeerror, \
    infile_aimin0_linaddintcpterror, \
    infile_aimin0_linsubintcpterror, \
    infile_aimin0_linicep5, \
    infile_aimin0_linicem5, \
    infile_aimin0_linicep15, \
    infile_aimin0_linicem15, \
    infile_aimin0_lincodp2, \
    infile_aimin0_lincodm2, \
    infile_aimin0_lincodp5, \
    infile_aimin0_lincodm5, \
    infile_aimin0_linL2L3p30, \
    infile_aimin0_linL2L3m30
]

# This plots the 6-panel results for a specific latitude band
# -----------------------------------------------------------
#plot_type_forcing_v4_all_months_arctic_avg(all_month_files, \
#        minlat = 70., maxlat = 87., \
#        data_type = 'raw', trend_type = 'standard', \
#        save = False,debug = False, labels = None)

# This plots all the results from all 3 latitude bands
# ----------------------------------------------------
plot_type_forcing_v4_all_months_arctic_avg_combined(all_month_files, \
    version4 = False, max_pval = 0.05, \
    horiz_orient = True, save = False)
sys.exit()

#uncert_slopes =  calc_forcing_slopes_v4_all_months_arctic_avg_uncert(all_month_files, \
#        OMI_daily_VSJ4, slope_type = 'lin')

calc_print_forcing_slope_error_v4(all_month_dict['FORCE_EST'], \
    OMI_daily_VSJ4, all_month_files, \
    min_year = 2005, max_year = 2020, \
    minlat = 65., maxlat = 87., \
    trend_type = 'standard', ptype = 'forcing', \
    vtype = '', version4 = True, slope_type = 'lin', 
    save = False)

sys.exit()

"""
# = = = = = = = = = = = =
#
# This code compares how handling -999s and NaNs in the coloc data
# gets rid of missing zones in the coloc data. 
#
# = = = = = = = = = = = =

#test_file = 'test_comp_data/colocated_subset_201807052213.hdf5'
test_file = 'comp_data/testing/colocated_subset_201807052213.hdf5'

data = h5py.File(test_file)

orig_data = data['modis_cld_top_pres'][:,:]
no_999 = np.ma.masked_where(orig_data == -999., orig_data)
no_nan = np.ma.masked_where(np.isnan(orig_data) == True, orig_data)
nan_0  = np.where(np.isnan(orig_data), 1025., orig_data)
all_fix  = np.where(np.isnan(orig_data), 1025., orig_data)
all_fix  = np.ma.masked_where(  (all_fix <= 0.), all_fix)

cod_orig_data = data['modis_cod'][:,:]
cod_no_999 = np.ma.masked_where(cod_orig_data == -999., cod_orig_data)
cod_no_nan = np.ma.masked_where(np.isnan(cod_orig_data) == True, cod_orig_data)
#cod_nan_0  = np.where(np.isnan(orig_data), 1025., orig_data)

fig = plt.figure()
ax1 = fig.add_subplot(2,4,1, projection = ccrs.NorthPolarStereo())
ax2 = fig.add_subplot(2,4,2, projection = ccrs.NorthPolarStereo())
ax3 = fig.add_subplot(2,4,3, projection = ccrs.NorthPolarStereo())
ax4 = fig.add_subplot(2,4,4, projection = ccrs.NorthPolarStereo())
ax5 = fig.add_subplot(2,4,5, projection = ccrs.NorthPolarStereo())
ax6 = fig.add_subplot(2,4,6, projection = ccrs.NorthPolarStereo())
ax7 = fig.add_subplot(2,4,7, projection = ccrs.NorthPolarStereo())
ax8 = fig.add_subplot(2,4,8, projection = ccrs.NorthPolarStereo())

ax1.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], orig_data, \
    shading = 'auto', transform = ccrs.PlateCarree())
ax1.coastlines()
ax1.set_extent([-180,180,65,90], ccrs.PlateCarree())

ax2.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], no_999, \
    shading = 'auto', transform = ccrs.PlateCarree())
ax2.coastlines()
ax2.set_extent([-180,180,65,90], ccrs.PlateCarree())

ax3.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], no_nan, \
    shading = 'auto', transform = ccrs.PlateCarree())
ax3.coastlines()
ax3.set_extent([-180,180,65,90], ccrs.PlateCarree())

ax4.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], nan_0, \
    shading = 'auto', transform = ccrs.PlateCarree())
ax4.coastlines()
ax4.set_extent([-180,180,65,90], ccrs.PlateCarree())

ax5.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], cod_orig_data, \
    shading = 'auto', transform = ccrs.PlateCarree())
ax5.coastlines()
ax5.set_extent([-180,180,65,90], ccrs.PlateCarree())

ax6.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], cod_no_999, \
    shading = 'auto', transform = ccrs.PlateCarree())
ax6.coastlines()
ax6.set_extent([-180,180,65,90], ccrs.PlateCarree())

ax7.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], cod_no_nan, \
    shading = 'auto', transform = ccrs.PlateCarree())
ax7.coastlines()
ax7.set_extent([-180,180,65,90], ccrs.PlateCarree())

ax8.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], all_fix, \
    shading = 'auto', transform = ccrs.PlateCarree())
ax8.coastlines()
ax8.set_extent([-180,180,65,90], ccrs.PlateCarree())

fig.tight_layout()

plt.show()

sys.exit()

# = = = = = = = = = = = = =
#
# END OF CTP FIX STUFF
#
# = = = = = = = = = = = = =

"""

##!#date_str = '201807040819'
##!#date_str = '201807040501'
##!##date_str = '201507081837'
##!#plot_compare_colocate_spatial(date_str, minlat = 65., zoom = False, \
##!#    save = False)

data = h5py.File('comp_data/testing/colocated_subset_201908110033.hdf5')

fig = plt.figure()

ax1 = fig.add_subplot(2,3,1, projection = ccrs.NorthPolarStereo())
ax2 = fig.add_subplot(2,3,2, projection = ccrs.NorthPolarStereo())
ax3 = fig.add_subplot(2,3,3, projection = ccrs.NorthPolarStereo())
ax4 = fig.add_subplot(2,3,4, projection = ccrs.NorthPolarStereo())
ax5 = fig.add_subplot(2,3,5, projection = ccrs.NorthPolarStereo())
ax6 = fig.add_subplot(2,3,6, projection = ccrs.NorthPolarStereo())

# Plot AI Perturbed
mask_ai = np.ma.masked_invalid(data['omi_uvai_pert'])
mask_ai = np.ma.masked_where(mask_ai == -999., mask_ai)
ax1.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], mask_ai, \
    transform = ccrs.PlateCarree(), shading = 'auto') 
ax1.set_extent([-180,180,65,90], ccrs.PlateCarree())
ax1.coastlines()
ax1.set_title('AI PERT ' + str(np.max(mask_ai)))

# Plot AI Raw
mask_ai = np.ma.masked_invalid(data['omi_uvai_raw'])
mask_ai = np.ma.masked_where(mask_ai == -999., mask_ai)
ax2.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], mask_ai, \
    transform = ccrs.PlateCarree(), shading = 'auto') 
ax2.set_extent([-180,180,65,90], ccrs.PlateCarree())
ax2.coastlines()
ax2.set_title('AI RAW ' + str(np.max(mask_ai)))

# Plot MODIS CH7
mask_ai = np.ma.masked_invalid(data['modis_ch7'])
mask_ai = np.ma.masked_where(mask_ai == -999., mask_ai)
ax3.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], mask_ai, \
    transform = ccrs.PlateCarree(), shading = 'auto') 
ax3.set_extent([-180,180,65,90], ccrs.PlateCarree())
ax3.coastlines()
ax3.set_title('CH7 ' + str(np.max(mask_ai)))

# Plot NSIDC ICE
mask_ai = np.ma.masked_invalid(data['nsidc_ice'])
mask_ai = np.ma.masked_where(mask_ai == -999., mask_ai)
ax4.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], mask_ai, \
    transform = ccrs.PlateCarree(), shading = 'auto') 
ax4.set_extent([-180,180,65,90], ccrs.PlateCarree())
ax4.coastlines()
ax4.set_title('ICE ' + str(np.max(mask_ai)))

# Plot CERES ALB
mask_ai = np.ma.masked_invalid(data['ceres_alb'])
mask_ai = np.ma.masked_where(mask_ai == -999., mask_ai)
ax5.pcolormesh(data['omi_lon'][:,:], data['omi_lat'][:,:], mask_ai, \
    transform = ccrs.PlateCarree(), shading = 'auto') 
ax5.set_extent([-180,180,65,90], ccrs.PlateCarree())
ax5.coastlines()
ax5.set_title('CH7 ' + str(np.max(mask_ai)))

fig.tight_layout()
plt.show()

data.close()

sys.exit()




##!#    #'20180701','20180702','20180703','20180706', '20180707','20180708']
##!#run_list = [
##!#    '20180702','20180703','20180706', '20180707','20180708']
##!#
##!##run_list = [\
##!#
##!## RUN THIS PART FIRST, TURN ON VPN AND RUN THIS
##!#
##!##shawn_path = home_dir+ '/data/OMI/shawn_files/ltc3/'
##!##for ddate in run_list:
##!##    dt_date_str = datetime.strptime(ddate, '%Y%m%d')
##!##    shawn_files = glob(dt_date_str.strftime(shawn_path + '%Y%m%d*')) 
##!##    if(len(shawn_files) == 0):
##!##        # download shawn files
##!##        cmnd = dt_date_str.strftime(\
##!##            'scp -r bsorenson@134.129.222.68:/Research/OMI/'+\
##!##            'out_files-ltc3/new/%Y%m%d* ') + shawn_path
##!##            #'scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/'+\
##!##            #'out_files-ltc3/new/%Y%m%d* ') + shawn_path
##!##        print(cmnd)
##!##        os.system(cmnd)
##!##    
##!##        shawn_files = glob(dt_date_str.strftime(shawn_path + '%Y%m%d*')) 
##!##sys.exit()
##!#
##!##final_list = entire_wrapper(min_AI = 2.0, minlat = 70., download = True, \
##!##    images = False, process = False, run_list = run_list, copy_to_raindrop = False, \
##!##    include_tropomi = True, remove_ch2_file = True)
##!#final_list = entire_wrapper(min_AI = -2.0, minlat = 65., download = True, \
##!#    images = False, process = False, run_list = run_list, copy_to_raindrop = False, \
##!#    include_tropomi = False, remove_ch2_file = True)
##!#
##!#sys.exit()

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

calc_pcnt_aerosol_over_type(dates, 1.5, minlat = 70., dtype = 'PERT', ax = None)
sys.exit()
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

