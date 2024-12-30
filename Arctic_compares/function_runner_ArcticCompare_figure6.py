#!/usr/bin/env python

"""

  SYNTAX FOR PAPER: ./function_runner_ArcticCompare_figure6.py 
        neuralnet_output_clear_newfiles/test_calc_out_noland105_201807082244.hdf5 
        neuralnet_output/test_calc_out_noland105_201807052213.hdf5

"""

import Arctic_compare_lib
from Arctic_compare_lib import *
import random
#from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

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

#sim_name = 'noland74'
#
## Read in the desired files
## -------------------------
#if( (sim_name == 'noland103') | (sim_name == 'noland104') ):
#    files = glob('neuralnet_output_clear_newfiles/test_calc_out_' + sim_name + '*.hdf5')
#else:
#    files = glob('neuralnet_output_clear/test_calc_out_' + sim_name + '*.hdf5')
#
#if(len(files) == 0):
#    print("ERROR: NO CLEAR FILES FOUND FOR SIM " + sim_name)
#    sys.exit()
#
## Figure out the total size to insert the data
## ---------------------------------------------
##minlat = 70.
#minlat = 65.    # 2024/01/10: changed to 65.
#total_size = 0
#for ff in files:
#    data = h5py.File(ff,'r')
#    #local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
#    local_data = np.ma.masked_invalid(data['omi_uvai_pert'])
#    local_data = np.ma.masked_where((local_data < -12) | (local_data > 1.0), local_data)
#    local_data = np.ma.masked_where(data['omi_lat'][:,:] < minlat, \
#        local_data) 
#    local_data = np.ma.masked_where((data['ceres_swf'][:,:] < -200.) | \
#        (data['ceres_swf'][:,:] > 3000), \
#        local_data) 
#    local_data = np.ma.masked_where(np.isnan(data['calc_swf'][:,:]), local_data)
#    local_size = local_data.compressed().shape[0]
#    max_ai = np.max(local_data)
#    print(ff, local_size, np.round(max_ai, 3))
#    total_size += local_size
#
#    data.close()
#
#
## Set up the data structure to hold all the data
## ----------------------------------------------
#combined_data = {}
#combined_data['calc_swf']      = np.full(total_size, np.nan)
#combined_data['ceres_swf']     = np.full(total_size, np.nan)
#
#print("Loading data")
#
## Loop back over the files and insert the data into the structure
## ---------------------------------------------------------------
#total_size = 0
#beg_idx = 0
#end_idx = 0
#for ff in files:
#
#    data = h5py.File(ff,'r')
#    # NOTE: Changed the omi variable here from "pert" to "raw" on 20230623.
#    #       This move should allow for coloc data to be read after 2020
#    #local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
#    local_data = np.ma.masked_invalid(data['omi_uvai_pert'])
#    print("MAX AI FOR SWATH", ff, np.max(local_data))
#    local_data = np.ma.masked_where((local_data < -12) | (local_data > 1.0), local_data)
#    local_data = np.ma.masked_where(data['omi_lat'][:,:] < minlat, \
#        local_data) 
#    local_data = np.ma.masked_where((data['ceres_swf'][:,:] < -200.) | \
#        (data['ceres_swf'][:,:] > 3000), \
#        local_data) 
#    local_data = np.ma.masked_where(np.isnan(data['calc_swf'][:,:]), local_data)
#    local_size = local_data.compressed().shape[0]
#
#    beg_idx = end_idx
#    end_idx = beg_idx + local_size
#
#    for tkey in combined_data.keys():
#        combined_data[tkey][beg_idx:end_idx] = \
#            data[tkey][~local_data.mask]
#
#    #print(local_size)
#    total_size += local_size
#
#    data.close()
#
##errors = combined_data['calc_swf'] - combined_data['ceres_swf']
#
## Calculate RMSE for binned values
## --------------------------------
#swf_bins1 = np.arange(30, 600, 30)
#rmse_vals1 = np.full(len(swf_bins1) - 1, np.nan)
##rmse_vals2 = np.full(len(swf_bins1) - 1, np.nan)
#for ii in range(len(swf_bins1) - 1):
#    keep_idxs = np.where( (combined_data['ceres_swf'] >= swf_bins1[ii]) & \
#        (combined_data['ceres_swf'] < swf_bins1[ii + 1]))
#    if(len(keep_idxs[0]) > 2):
#        rmse1 = mean_squared_error(combined_data['ceres_swf'][keep_idxs], \
#            combined_data['calc_swf'][keep_idxs], squared = False)
#        counts1 = len(keep_idxs[0])
#        #keep_idxs = np.where( (both_orig2 >= swf_bins1[ii]) & (both_orig2 < swf_bins1[ii + 1]))
#        #rmse2 = mean_squared_error(both_orig2[keep_idxs], both_calc2[keep_idxs], squared = False)
#        #counts2 = len(keep_idxs[0])
#        #print(swf_bins1[ii], swf_bins1[ii + 1], np.round(rmse1, 1), counts1, np.round(rmse2, 1), counts2) 
#
#        rmse_vals1[ii] = rmse1
#        #rmse_vals2[ii] = rmse2
#
#swf_centers = (swf_bins1[1:] + swf_bins1[:-1])/2.
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.scatter(swf_centers, rmse_vals1)
#plt.show()
#sys.exit()





l_calc_noise_floor = False




if(l_calc_noise_floor):

    in_calc1 = h5py.File(sys.argv[1], 'r')
    mask_orig1 = np.ma.masked_where((in_calc1['ceres_swf'][:,:] == -999.) | \
                                   (in_calc1['ceres_swf'][:,:] > 3000), in_calc1['ceres_swf'][:,:])
    mask_calc1 = np.ma.masked_invalid(in_calc1['calc_swf'])
    mask_calc1 = np.ma.masked_where(mask_calc1 == -999., mask_calc1)
    
    in_calc2 = h5py.File(sys.argv[2], 'r')
    mask_orig2 = np.ma.masked_where((in_calc2['ceres_swf'][:,:] == -999.) | \
                                   (in_calc2['ceres_swf'][:,:] > 3000), in_calc2['ceres_swf'][:,:])
    mask_calc2 = np.ma.masked_invalid(in_calc2['calc_swf'])
    mask_calc2 = np.ma.masked_where(mask_calc2 == -999., mask_calc2)
    
    
    
    #mask_orig1 = np.ma.masked_where(np.isnan(in_calc1['omi_uvai_pert'][:,:]), mask_orig1)
    #mask_calc1 = np.ma.masked_where(np.isnan(in_calc1['omi_uvai_pert'][:,:]), mask_calc1)
    mask_orig1 = np.ma.masked_where(in_calc1['omi_uvai_pert'][:,:] > 1.0, mask_orig1)
    mask_calc1 = np.ma.masked_where(in_calc1['omi_uvai_pert'][:,:] > 1.0, mask_calc1)
    both_orig1 = np.ma.masked_where((mask_orig1.mask == True) | (mask_calc1.mask == True), mask_orig1)
    both_calc1 = np.ma.masked_where((mask_orig1.mask == True) | (mask_calc1.mask == True), mask_calc1)
    both_orig1 = both_orig1.compressed()
    both_calc1 = both_calc1.compressed()
    
    #mask_orig2 = np.ma.masked_where(np.isnan(in_calc2['omi_uvai_pert'][:,:]), mask_orig2)
    #mask_calc2 = np.ma.masked_where(np.isnan(in_calc2['omi_uvai_pert'][:,:]), mask_calc2)
    mask_orig2 = np.ma.masked_where(in_calc2['omi_uvai_pert'][:,:] > 1.0, mask_orig2)
    mask_calc2 = np.ma.masked_where(in_calc2['omi_uvai_pert'][:,:] > 1.0, mask_calc2)
    both_orig2 = np.ma.masked_where((mask_orig2.mask == True) | (mask_calc2.mask == True), mask_orig2)
    both_calc2 = np.ma.masked_where((mask_orig2.mask == True) | (mask_calc2.mask == True), mask_calc2)
    both_orig2 = both_orig2.compressed()
    both_calc2 = both_calc2.compressed()
    
    # Calculate RMSE for binned values
    # --------------------------------
    swf_bins1 = np.arange(30, 600, 30)
    rmse_vals1 = np.full(len(swf_bins1) - 1, np.nan)
    rmse_vals2 = np.full(len(swf_bins1) - 1, np.nan)
    for ii in range(len(swf_bins1) - 1):
        keep_idxs = np.where( (both_orig1 >= swf_bins1[ii]) & (both_orig1 < swf_bins1[ii + 1]))
        if(len(keep_idxs[0]) > 2):
            rmse1 = mean_squared_error(both_orig1[keep_idxs], both_calc1[keep_idxs], squared = False)
            counts1 = len(keep_idxs[0])
            keep_idxs = np.where( (both_orig2 >= swf_bins1[ii]) & (both_orig2 < swf_bins1[ii + 1]))
            rmse2 = mean_squared_error(both_orig2[keep_idxs], both_calc2[keep_idxs], squared = False)
            counts2 = len(keep_idxs[0])
            print(swf_bins1[ii], swf_bins1[ii + 1], np.round(rmse1, 1), counts1, np.round(rmse2, 1), counts2) 
    
            rmse_vals1[ii] = rmse1
            rmse_vals2[ii] = rmse2
    
    swf_centers = (swf_bins1[1:] + swf_bins1[:-1])/2.
    
    in_calc1.close()
    in_calc2.close()
    
    fig = plt.figure(figsize = (11, 4))
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)
    
    xy = np.vstack([both_orig1, both_calc1])
    r2 = r2_score(both_orig1, both_calc1)
    z = stats.gaussian_kde(xy)(xy)       
    ax1.scatter(both_orig1, both_calc1, c = z, s = 1)
    lims = [\
            np.min([ax1.get_xlim(), ax1.get_ylim()]),\
            np.max([ax1.get_xlim(), ax1.get_ylim()]),
    ]
    ax1.plot(lims, lims, 'r', linestyle = ':', alpha = 1.00)
    ax1.set_xlim(lims)
    ax1.set_ylim(lims)
    ax1.set_xlabel('CERES SWF [Wm$^{-2}$]')
    ax1.set_ylabel('NN SWF [Wm$^{-2}$]')
    
    xy = np.vstack([both_orig2, both_calc2])
    r2 = r2_score(both_orig2, both_calc2)
    z = stats.gaussian_kde(xy)(xy)       
    ax2.scatter(both_orig2, both_calc2, c = z, s = 1)
    lims = [\
            np.min([ax2.get_xlim(), ax2.get_ylim()]),\
            np.max([ax2.get_xlim(), ax2.get_ylim()]),
    ]
    ax2.plot(lims, lims, 'r', linestyle = ':', alpha = 1.00)
    ax2.set_xlim(lims)
    ax2.set_ylim(lims)
    ax2.set_xlabel('CERES SWF [Wm$^{-2}$]')
    ax2.set_ylabel('NN SWF [Wm$^{-2}$]')
    
    ax3.scatter(swf_centers, rmse_vals1, label = 'File1')
    ax3.scatter(swf_centers, rmse_vals2, label = 'File2')
    ax3.legend()
    ax3.set_xlabel('CERES SWF [Wm$^{-2}$]')
    ax3.set_ylabel('RMSE [Wm$^{-2}$]')
    
    fig.tight_layout()
    plt.show()
    
    sys.exit()

else:

    # Compare the OMI, CERES, NN, and NN - CERES for two swaths
    # ---------------------------------------------------------
    plot_compare_NN_output_double(sys.argv[1], \
        sys.argv[2], \
        save = True, include_scatter = False)
    #plot_compare_NN_output_double(sys.argv[1], sys.argv[2], save = False, include_scatter = False)


