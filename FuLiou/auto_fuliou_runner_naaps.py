#!/usr/bin/env python

"""
  NAME:
    auto_fuliou_runner.py

  PURPOSE:
    Automate the running of FuLiou for each unique point in a passed
        CCPEXCV-HALO NAAPS/lidar file.

  SYNTAX:
    python extract_naapslidar.py CCPEXCV-HALO_file

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2023/10/03:
      Written

"""

import sys
from datetime import datetime, timedelta
import h5py
import numpy as np
import os
from scipy.stats import mode

# This switch determines if an HDF5 file is generated from the
# FuLiou output files, in addition to the text file. The 
# HDF5 file is titled 'fuliou_out_YYYYMMDD.hdf5'
# ----------------------------------------------------------------------
l_GEN_HDF5_OUTPUT = False

# Check command line arguments
# ----------------------------
if(len(sys.argv) != 2):
    print("SYNTAX: python extract_naapslidar.py CCPEXCV-HALO_file")
    sys.exit(1)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Set path variables
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#base_dir = "/Research/for_blake_from_zhang/fuliou_package/
#base_dir = "/Research/fuliou_lidarflux/fuliou_package/"
base_dir = "/Research/fuliou_lidarflux/test_package/"
data_dir = "/Research/NAAPS_NOGAPS_2023/"
mhome = base_dir + "test_run"
output_dir = mhome + '/data'

homedir = mhome

########################################################
#Parameter Group 2

########################################################
#locations for executables
#naaps_fuliou
#setenv NAAPS_FULIOUCODE "/home/jzhang/fuliou_package/sourcecodes"
naaps_fulioucode = base_dir + "sourcecodes"

tabledir  = naaps_fulioucode + "/table"             ##place to store .txt files needed for the program
srcdir  =  naaps_fulioucode 

#tools
mtools = naaps_fulioucode + "/tools"
#setenv MTOOLS      "$NAAPS_FULIOUCODE/tools"


########################################################
#locations for data

#NOGAPS data
#setenv NOGAPS_SOURCE "/home/jzhang/fuliou_package/test_data/test_nogaps"
#setenv NOGAPS_SOURCE "/Research/for_blake_from_zhang/fuliou_package/test_data/test_nogaps"
#setenv NOGAPS_SOURCE "/Research/NAAPS_NOGAPS_2023/NAAPS/"

#setenv NOGAPS_SOURCE "/Research/NAAPS_NOGAPS_2023/nog_cmorph"
nogaps_source = data_dir + "nog_cmorph"

#NAAPS DATA
#setenv naapsdir  "/home/jzhang/fuliou_package/test_data/test_naaps"
#setenv naapsdir  "/Research/for_blake_from_zhang/fuliou_package/test_run/test_naaps"

#setenv naapsdir  "/Research/NAAPS_NOGAPS_2023/NAAPS"
naapsdir = data_dir + "NAAPS"
concdir = os.getcwd()
nogapsdir = os.getcwd()
predir = os.getcwd()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Write namelist file
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

namelist_file = 'namelist'
dtg = '2000010100'
#set namelist = "${runhome}/namelist"
#/bin/rm -r $namelist >& /dev/null
with open(namelist_file, 'w') as fname:
    fname.write("    &nl\n" + \
                "    tape        = '" + dtg + "',\n" + \
                "    tabledir    = '" + tabledir + "',\n" + \
                "    sulfatefile = 'anthro_farop_v3.txt',\n" + \
                "    dustfile    = 'dust_farop_v3.txt',\n" + \
                "    smokefile   = 'smoke_farop_v3.txt',\n" + \
                "    saltfile    = 'salt_farop_v3.txt',\n" + \
                "    concdir     = '" + concdir + "',\n" + \
                "    nogapsdir   = '" + nogapsdir + "',\n" + \
                "    predir      = '" + predir + "',\n" + \
                "    sceneidfile =  'global_scene_id.map',\n" + \
                "    &end\n")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Open the data file and extract the variables
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

infile = sys.argv[1]
hdf_data = h5py.File(infile)

naaps_press   = hdf_data['NAAPS_pressure'][:,:] / 100.

# If the user wants to screen out any extra NAAPS layer info,
# -----------------------------------------------------------
l_remove_doubles = True
if(l_remove_doubles):
    
    # Figure out the vertical indices
    # -------------------------------
    ##!#unique_vert_press = np.unique(naaps_press[0,:])[::-1]
    ##!#keep_vert_idxs = np.array([\
    ##!#    np.where(naaps_press[0,:] == tup)[0][0] \
    ##!#    for tup in unique_vert_press])

    # Figure out the horizontal indices
    # -------------------------------
    unique_horz_press = np.unique(naaps_press[:,0])
    res, ind = np.unique(naaps_press[:,0], return_index = True)
    keep_horz_idxs = ind[np.argsort(ind)]

    #naaps_press   = naaps_press[keep_horz_idxs,:][:,keep_vert_idxs]
    naaps_press   = naaps_press[keep_horz_idxs,:]
else:     
    keep_horz_idxs = np.arange(naaps_press.shape[0])
    keep_vert_idxs = np.arange(naaps_press.shape[1])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Convert the mass concentrations from ug/m3 to kg/m3 by  multiplying
#  each value by 1e-9.
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Extract metadata
# ----------------
split_date    = hdf_data['time'].attrs['units'].strip().split()
date_str      = ' '.join(split_date[2:])
base_date     = datetime.strptime(date_str, '%Y%m%d %H:%M:%S')
time          = hdf_data['time'][:][keep_horz_idxs]
naaps_masl    = hdf_data['NAAPS_MASL'][keep_horz_idxs,:]
naaps_dz      = hdf_data['NAAPS_dz'][keep_horz_idxs,:]
lat           = hdf_data['lat'][:][keep_horz_idxs]
lon           = hdf_data['lon'][:][keep_horz_idxs]

# Extract NAAPS metorological variables
# -------------------------------------
naaps_temp    = hdf_data['NAAPS_temperature'\
    ][keep_horz_idxs,:]
naaps_sfcpres = hdf_data['NAAPS_ps'][:][keep_horz_idxs]
naaps_spechum = hdf_data['NAAPS_q'][keep_horz_idxs,:]
naaps_relhum  = hdf_data['NAAPS_rh'][keep_horz_idxs,:]

# Extract NAAPS aerosol variables
# -------------------------------
naaps_tot_aot   = hdf_data['NAAPS_total_aot'][:][keep_horz_idxs]

naaps_dust_ext  = hdf_data['NAAPS_dust_ext'][keep_horz_idxs,:]
naaps_dust_aot  = hdf_data['NAAPS_dust_aot'][:][keep_horz_idxs]
naaps_dust_mass = \
    hdf_data['NAAPS_dust_mass_concentration'\
    ][keep_horz_idxs,:] * 1e-9

naaps_abf_ext   = hdf_data['NAAPS_ABF_ext'][keep_horz_idxs,:]
naaps_abf_aot   = hdf_data['NAAPS_ABF_aot'][:][keep_horz_idxs]
naaps_abf_mass  = \
    hdf_data['NAAPS_ABF_mass_concentration'\
    ][keep_horz_idxs,:] * 1e-9

naaps_salt_ext  = hdf_data['NAAPS_SS_ext'][keep_horz_idxs,:]
naaps_salt_aot  = hdf_data['NAAPS_SS_aot'][:][keep_horz_idxs]
naaps_salt_mass = \
    hdf_data['NAAPS_SS_mass_concentration'\
    ][keep_horz_idxs,:] * 1e-9

naaps_smoke_ext  = hdf_data['NAAPS_smoke_ext'\
    ][keep_horz_idxs,:]
naaps_smoke_aot  = hdf_data['NAAPS_smoke_aot'][:][keep_horz_idxs]
naaps_smoke_mass = \
    hdf_data['NAAPS_smoke_mass_concentration'\
    ][keep_horz_idxs,:] * 1e-9

# Extract lidar aerosol variables
# -------------------------------
ext_532       = hdf_data['532_ext'][keep_horz_idxs,:]

ext_532 = np.ma.where(np.ma.masked_invalid(ext_532).mask == True, \
    0., ext_532)

# Close the file
# --------------
hdf_data.close()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Loop over the data and print the stuff
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Variable descriptions
#  0 - napsp
#  1 - naaps_z
#  2 - naaps_dz
#  3 - napst
#  4 - napsrh
#  5 - napssh
#  6 - sulext
#  7 - dustext
#  8 - smokeext
#  9 - saltext
# 10 - sulmass
# 11 - dustmass
# 12 - smokemass
# 13 - saltmass

fmt_str = "{0:7.1f} {1:8.1f} {2:8.1f} {3:7.2f} " + \
    "{4:7.3f} {5:10.3e} " + \
    "{6:10.3e} {7:10.3e} {8:10.3e} {9:10.3e} " + \
    "{10:10.3e} {11:10.3e} {12:10.3e} {13:10.3e}\n"

# totaod = [aod1, aod2, aod3, aod4, aod5]
# totaod: [UNUSED, sulfate, dust, smoke, salt]

#foutname = base_dir + \
#            'test_data/test_naaps_new/test_naaps_file.txt' 
foutname = 'test_naaps_file.txt' 

#out_data_name = output_dir + '/test_output_file.txt'
out_data_name = 'test_output_file.txt'

# If the outputfile exists in its current name, delete the old one
# ----------------------------------------------------------------
if(os.path.exists(out_data_name)):
    print("Removing old output data file")
    cmnd = 'rm ' + out_data_name
    os.system(cmnd)

##work_file = 'test_comp_work.txt'
##if(os.path.exists(work_file)):
##    print("Removing work file")
##    cmnd = 'rm ' + work_file
##    os.system(cmnd)

# Prep the output HDF5 file, if desired
# -------------------------------------
num_layers = naaps_temp.shape[1]
if(l_GEN_HDF5_OUTPUT):
    # NOTE: In the FuLiou code, the output printed to the screen begins at 
    # index 7 of the arrays in the code. This is mirrored here
    # --------------------------------------------------------------------
    #num_pres_layers = num_layers - 6
    #num_heat_layers = num_layers - 7
    num_pres_layers = num_layers - 1
    num_heat_layers = num_layers - 2

    output_time           = np.zeros((time.shape[0]))
    output_lat            = np.zeros((time.shape[0]))
    output_lon            = np.zeros((time.shape[0]))
    output_aod            = np.zeros((time.shape[0]))
    output_sfc_dwn_sw_flx = np.zeros((time.shape[0]))
    output_toa_up_sw_flx  = np.zeros((time.shape[0]))
    output_sfc_dwn_lw_flx = np.zeros((time.shape[0]))
    output_toa_up_lw_flx  = np.zeros((time.shape[0]))

    output_alt        = np.zeros((num_pres_layers))
    output_press      = np.zeros((time.shape[0], num_pres_layers))
    output_tmp        = np.zeros((time.shape[0], num_pres_layers))
    output_relhum     = np.zeros((time.shape[0], num_pres_layers))
    output_spchum     = np.zeros((time.shape[0], num_pres_layers))
    output_abf_ext    = np.zeros((time.shape[0], num_pres_layers))
    output_dust_ext   = np.zeros((time.shape[0], num_pres_layers))
    output_smoke_ext  = np.zeros((time.shape[0], num_pres_layers))
    output_salt_ext   = np.zeros((time.shape[0], num_pres_layers))

    output_heat       = np.zeros((time.shape[0], num_heat_layers))


    output_time[:]           = np.nan
    output_lat[:]            = np.nan
    output_lon[:]            = np.nan
    output_aod[:]            = np.nan
    output_sfc_dwn_sw_flx[:] = np.nan
    output_toa_up_sw_flx[:]  = np.nan
    output_sfc_dwn_lw_flx[:] = np.nan
    output_toa_up_lw_flx[:]  = np.nan
    output_press[:,:]        = np.nan
    output_heat[:,:]         = np.nan

    output_alt[:]          = np.nan
    output_tmp[:,:]        = np.nan
    output_relhum[:,:]     = np.nan
    output_spchum[:,:]     = np.nan
    output_abf_ext[:,:]    = np.nan
    output_dust_ext[:,:]   = np.nan
    output_smoke_ext[:,:]  = np.nan
    output_salt_ext[:,:]   = np.nan


for ii in range(time.shape[0]):
#for ii in range(18,19):

    # Figure out the vertical indices
    # -------------------------------
    unique_vert_press = np.unique(naaps_press[ii,:])[::-1]
    keep_vert_idxs = np.array([\
        np.where(naaps_press[ii,:] == tup)[0][0] \
        for tup in unique_vert_press])

    local_ext_532       = ext_532[ii,:][keep_vert_idxs]
    local_naaps_press   = naaps_press[ii,:][keep_vert_idxs]
    local_naaps_masl    = naaps_masl[ii,:][keep_vert_idxs]
    local_naaps_dz      = naaps_dz[ii,:][keep_vert_idxs]
    local_naaps_temp    = naaps_temp[ii,:][keep_vert_idxs]
    local_naaps_spechum = naaps_spechum[ii,:][keep_vert_idxs]
    local_naaps_relhum  = naaps_relhum[ii,:][keep_vert_idxs]
    local_naaps_dust_ext   = naaps_dust_ext[ii,:][keep_vert_idxs]
    local_naaps_dust_mass  = naaps_dust_mass[ii,:][keep_vert_idxs] * 1e-9
    local_naaps_abf_ext    = naaps_abf_ext[ii,:][keep_vert_idxs]
    local_naaps_abf_mass   = naaps_abf_mass[ii,:][keep_vert_idxs] * 1e-9
    local_naaps_smoke_ext  = naaps_smoke_ext[ii,:][keep_vert_idxs]
    local_naaps_smoke_mass = naaps_smoke_mass[ii,:][keep_vert_idxs] * 1e-9
    local_naaps_salt_ext   = naaps_salt_ext[ii,:][keep_vert_idxs]
    local_naaps_salt_mass  = naaps_salt_mass[ii,:][keep_vert_idxs] * 1e-9

    with open(foutname, 'w') as fout:
    
        """
        
        Output format
        YYYYMMDDHHMM lat lon tau_aer
         sulfate_aod, dust_aod, smoke_aod, salt_aod
         press NAAPS_masl NAAPS_dz temp rhum spechum abf_ext dust_ext smoke_ext salt_ext abf_mas dust_mas smoke_mas salt_mas
         ...
         ...
         ...
        
         
        """
   
        # Print the metadata information
        # ------------------------------
        local_dt_date = base_date + timedelta(seconds = time[ii]) 
        local_date = local_dt_date.strftime('%Y%m%d%H%M%S')
        print(ii, local_dt_date.strftime('%Y%m%d %H:%M:%S'), lat[ii])

        fout.write("{0} {1:8.5f} {2:8.5f} {3:8.5f}\n".format(\
            local_date, lat[ii], lon[ii], naaps_tot_aot[ii]))
        
        fout.write("{0:10.3e} {1:10.3e} {2:10.3e} {3:10.3e}\n".format(\
            naaps_abf_aot[ii], naaps_dust_aot[ii], \
            naaps_smoke_aot[ii], naaps_salt_aot[ii]))
        
        # Print the profile information in reverse order, with the lowest
        # data in the profile being at the last index.
        # ----------------------------------------------------------------------
        for jj in range(local_ext_532.shape[0]-1, -1, -1):
            fout.write(fmt_str.format(\
                local_naaps_press[jj], local_naaps_masl[jj], local_naaps_dz[jj], \
                local_naaps_temp[jj], local_naaps_relhum[jj], \
                local_naaps_spechum[jj], \
                local_naaps_abf_ext[jj], local_naaps_dust_ext[jj], \
                local_naaps_smoke_ext[jj], local_naaps_salt_ext[jj], \
                local_naaps_abf_mass[jj], local_naaps_dust_mass[jj], \
                local_naaps_smoke_mass[jj], local_naaps_salt_mass[jj]
            ))
  
    
    mode_vals = mode(local_naaps_press[:])
    if(mode_vals.count[0] != 1):
        print(ii, "BAD PRESS", mode_vals.mode[0], mode_vals.count[0])
 
    ##!#os.system('cp ' + foutname + ' ' + base_dir + \
    #!#        'test_data/test_naaps_new/test_naaps_file_' + str(int(ii)) + '.txt')

    #### Copy the first couple of lines to the work file 
    #### -----------------------------------------------
    ###cmnd = 'head -n 2 ' + foutname + ' >> ' + work_file
    ###print(cmnd)
    ###os.system(cmnd)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    #  Call the FuLiou code
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    """
    #print("Calling FuLiou code") 
    print(' Execute naaps_fuliou.exe')
    cmnd = 'time ' + srcdir + '/bin/naaps_fuliou.exe'
    print(cmnd)
    os.system(cmnd)
    """
    # Concatenate the temp output file to the end of the actual
    # output file. Only do this if the code generated an output file, 
    # which only happens if the code successfully runs.
    # ---------------------------------------------------------
    if(os.path.exists('fort.22')):
        print("COPYING FILE OUTPUTS")
        cmnd = 'cat fort.22 >> ' + out_data_name
        print(cmnd)
        os.system(cmnd)

        if(l_GEN_HDF5_OUTPUT):
        
            output_time[ii]            = time[ii]

            # Read in the last output file
            # ----------------------------
            with open('fort.22','r') as fin:
                flines = fin.readlines()

            line1 = flines[0].strip().split()
            output_lat[ii]            = float(line1[0])
            output_lon[ii]            = float(line1[1])
            output_aod[ii]            = float(line1[2])
            output_sfc_dwn_sw_flx[ii] = float(line1[3])
            output_toa_up_sw_flx[ii]  = float(line1[4])
            output_sfc_dwn_lw_flx[ii] = float(line1[5])
            output_toa_up_lw_flx[ii]  = float(line1[6])

            # Figure out how to properly divide the file
            # ------------------------------------------
            num_data_lines = (len(flines) - 1) / 2
            
            output_press[ii,:] = np.array([float(tvar) for tvar in \
                ''.join(flines[1:num_data_lines+1]).strip().split()])
            output_heat[ii,:] = np.array([float(tvar) for tvar in \
                ''.join(flines[num_data_lines+1:]).strip().split()])

            # Read in the input file data, containing met profiles and aerosol
            # ----------------------------------------------------------------
            indata = np.loadtxt('test_naaps_file.txt', skiprows = 2)

            work_idx = num_layers - num_pres_layers
            output_alt[:]           = indata[work_idx:,1]
            output_tmp[ii,:]        = indata[work_idx:,3]
            output_relhum[ii,:]     = indata[work_idx:,4]
            output_spchum[ii,:]     = indata[work_idx:,5]
            output_abf_ext[ii,:]    = indata[work_idx:,6]
            output_dust_ext[ii,:]   = indata[work_idx:,7]
            output_smoke_ext[ii,:]  = indata[work_idx:,8]
            output_salt_ext[ii,:]   = indata[work_idx:,9]

    ###cmnd = 'head -n 1 fort.22 >> ' + work_file
    ###print(cmnd)
    ###os.system(cmnd)

# At the end, if desired, make the output HDF5 file
if(l_GEN_HDF5_OUTPUT):

    out_name = 'fuliou_output_naaps_' + \
        infile.strip().split('/')[-1].split('_')[2] + '.h5'

    # Here, remove the indices that contain missing values?
    # -----------------------------------------------------
    keep_idxs = np.where(~np.isnan(output_lat[:]))

    dset = h5py.File(out_name, 'w')

    cdt = dset.create_dataset('time',   data = output_time[keep_idxs])
    cdt.attrs['long_name'] = 'UTC time of measurement instant'
    cdt.attrs['units'] = 'seconds after midnight on date provided in file name'
    cdt = dset.create_dataset('latitude',   data = output_lat[keep_idxs])
    cdt.attrs['long_name'] = 'latitude'
    cdt.attrs['units'] = 'degrees north (negative for degrees south)'
    cdt = dset.create_dataset('longitude',  data = output_lon[keep_idxs])
    cdt.attrs['long_name'] = 'longitude'
    cdt.attrs['units'] = 'degrees east (negative for degrees west)'
    cdt = dset.create_dataset('AOD_550',    data = output_aod[keep_idxs])
    cdt.attrs['long_name'] = 'lidar total AOD at 532 nm'
    cdt.attrs['units'] = 'none'
    cdt = dset.create_dataset('sfc_dwn_sw_flux',  data = output_sfc_dwn_sw_flx[keep_idxs])
    cdt.attrs['long_name'] = 'Surface downward SW flux'
    cdt.attrs['units'] = 'Wm^-2'
    cdt = dset.create_dataset('toa_up_sw_flux',   data = output_toa_up_sw_flx[keep_idxs])
    cdt.attrs['long_name'] = 'TOA upward SW flux'
    cdt.attrs['units'] = 'Wm^-2'
    cdt = dset.create_dataset('sfc_dwn_lw_flux',  data = output_sfc_dwn_lw_flx[keep_idxs])
    cdt.attrs['long_name'] = 'Surface downward LW flux'
    cdt.attrs['units'] = 'Wm^-2'
    cdt = dset.create_dataset('toa_up_lw_flux',   data = output_toa_up_lw_flx[keep_idxs])
    cdt.attrs['long_name'] = 'TOA upward LW flux'
    cdt.attrs['units'] = 'Wm^-2'
    cdt = dset.create_dataset('press_lvl',        data = output_press[keep_idxs])
    cdt.attrs['long_name'] = 'NAAPS Pressure'
    cdt.attrs['units'] = 'mb'
    cdt = dset.create_dataset('altitude',         data = output_alt[:])
    cdt.attrs['long_name'] = 'lidar bin altitude (above sea level)'
    cdt.attrs['units'] = 'm'
    cdt = dset.create_dataset('temperature',      data = output_tmp[keep_idxs])
    cdt.attrs['long_name'] = 'NAAPS air temperature'
    cdt.attrs['units'] = 'K'
    cdt = dset.create_dataset('rel_hum',          data = output_relhum[keep_idxs])
    cdt.attrs['long_name'] = 'NAAPS relative humidity'
    cdt.attrs['units'] = '%'
    cdt = dset.create_dataset('spec_hum',         data = output_spchum[keep_idxs])
    cdt.attrs['long_name'] = 'NAAPS specific humidity'
    cdt.attrs['units'] = 'kg/kg'
    cdt = dset.create_dataset('abf_ext',          data = output_abf_ext[keep_idxs])
    cdt.attrs['long_name'] = 'lidar abf extinction'
    cdt.attrs['units'] = 'm^-1'
    cdt = dset.create_dataset('dust_ext',         data = output_dust_ext[keep_idxs])
    cdt.attrs['long_name'] = 'lidar dust extinction'
    cdt.attrs['units'] = 'm^-1'
    cdt = dset.create_dataset('smoke_ext',        data = output_smoke_ext[keep_idxs])
    cdt.attrs['long_name'] = 'lidar smoke extinction'
    cdt.attrs['units'] = 'm^-1'
    cdt = dset.create_dataset('salt_ext',         data = output_salt_ext[keep_idxs])
    cdt.attrs['long_name'] = 'lidar salt extinction'
    cdt.attrs['units'] = 'm^-1'
    cdt = dset.create_dataset('net_heating_rate', data = output_heat[keep_idxs])
    cdt.attrs['long_name'] = 'net daily heating rate'
    cdt.attrs['units'] = 'K/day'

    dset.close()
    print("Saved HDF5 output file",out_name)
