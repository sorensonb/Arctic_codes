#!/usr/bin/env python

"""


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
#sim_name = 'noland74'
sim_name = 'noland105'
#sim_name = 'noland103'
#sim_name = 'noland75'
print("AS OF 2024/09/09, USING ", sim_name)

dates = [
         '200607240029',
         '200607240208',
         '200607240347',
         '200607240526',
         '200607242016',
         '200607242155',
         '200607242334',
         '200607250112',
         '200607250251',
         '200607250430',
         '200607252238',
         '200607260017',
         '200607260156',
         '200607260335',
         '200607260513',
         '200607260652',
         '200607260831',
         '200607262003',
         '200607262142',
         '200607262321',
         '200607270100',
         '200607270239',
         '200607270418',
         '200607270557',
         '200607270736',
         '200607270914',
         '200607272047',
         '200607272226',
         '200804221841',
         '200804222159',
         '201408110046',
         '201408110404',
         '201408111853',
         '201408112032',
         '201408112211',
         '201408120308',
         '201408121758',
         '201408121937',
         '201408122115',
         '201408122254',
         '201506271220',
         '201506271359',
         '201506271538',
         '201506271717',
         '201506271856',
         '201507061353',
         '201507061532',
         '201507061711',
         '201507061850',
         '201507062028',
         '201507062207',
         '201507071615',
         '201507071754',
         '201507081837',
         '201507082016',
         '201507082155',
         '201507090113',
         '201507090748',
         '201507090927',
         '201507091245',
         '201507091424',
         '201507091603',
         '201507100514',
         '201507100653',
         '201708161146',
         '201708161325',
         '201708161643',
         '201708161821',
         '201708162000',
         '201708171050',
         '201708171229',
         '201708171547',
         '201708171726',
         '201708172222',
         '201708181133',
         '201708181312',
         '201708181451',
         '201708181630',
         '201708181809',
         '201708181948',
         '201708191038',
         '201708191217',
         '201708191355',
         '201708191534',
         '201708191852',
         '201807040005',
         '201807040144',
         '201807040322',
         '201807041633',
         '201807041812',
         '201807041951',
         '201807042130',
         '201807042309',
         '201807050048',
         '201807050227',
         '201807051538',
         '201807051717',
         '201807051856',
         '201807052034',
         '201807052213',
         '201807052352',
         '201807210047',
         '201807211358',
         '201807211537',
         '201807211716',
         '201807211855',
         '201807212034',
         '201808100200',
         '201808140135',
         '201808141804',
         '201808141942',
         '201808142121',
         '201808142300',
         '201808260158',
         '201808260337',
         '201808260655',
         '201908100129',
         '201908100308',
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
         '201908112158',
         '201908112337'
    ]

#calc_pcnt_aerosol_over_type(dates, 1.5, minlat = 70., dtype = 'PERT', ax = None, local_dir = 'comp_data/new_data_20241029/')
#sys.exit()

calc_pcnt_aerosol_over_type_v2(sim_name, 1.0, minlat = 65., \
    dtype = 'PERT', ax = None, save = True)
