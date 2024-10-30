#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *
import random
#from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

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
    file_start = 'arctic_monthly_force_trends_count300_newerror'
    #file_start = 'arctic_monthly_force_trends_count'
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
        #file_start = 'arctic_monthly_force_values_count300_newerror'
        #file_start = 'arctic_monthly_force_values_count'
        #force_files = glob(file_start +'*.hdf5')[:2]

        # NOTE: These files are for noland74, with old error
        #force_files = ['arctic_monthly_force_values_count300.hdf5', \
        #               'arctic_monthly_force_values_count300_v1.hdf5', \
        #               'arctic_monthly_force_values_count300_v2.hdf5', \
        #               'arctic_monthly_force_values_count300_v3.hdf5', \
        #               'arctic_monthly_force_values_count300_v4.hdf5']

        # NOTE: These files are for noland74, with new error
        #force_files = ['arctic_monthly_force_values_count300_newerror.hdf5', \
        #               'arctic_monthly_force_values_count300_newerror_v1.hdf5']

        # NOTE: These files are for noland103
        force_files = ['arctic_monthly_force_values_count300_noland103.hdf5', \
                       'arctic_monthly_force_values_count300_noland103_v1.hdf5']

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
    
trends, pvals = calc_arctic_avg_region_trends(sim_values)

