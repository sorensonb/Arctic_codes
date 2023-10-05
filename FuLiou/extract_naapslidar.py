#!/usr/bin/env python

"""


"""

#import importlib, FuLiouLib
#from FuLiouLib import *
import sys
from netCDF4 import Dataset
from datetime import datetime, timedelta
import h5py
import numpy as np

# Check command line arguments
# ----------------------------
if(len(sys.argv) != 2):
    print("SYNTAX: python extract_naapslidar.py CCPEXCV-HALO_file")
    sys.exit(1)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Open the file and extract the variables
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

infile = sys.argv[1]
nc_data = Dataset(infile, 'r')
#hdf_data = h5py.File(infile)

"""
  call fuliou(latt,lont,alb,ems,
    tau_aer, totaod,csza,
    dustext,smokeext,saltext,sulext,
    dustssa,smokessa,saltssa,sulssa,
    dustasy,smokeasy,saltasy,sulasy,
    wave_mex,wave_wo,wave_g,wave,kwave,
    napsp,napst,napsrh,napso3,napssh,surft )

- latt
- lont
- alb
- ems
- tau_aer
- totaod
- csza
- x dustext
- x smokeext
- x saltext
- x sulext
- dustssa
- smokessa
- saltssa
- sulssa
- dustasy
- smokeasy
- saltasy
- sulasy
- wave_mex
- wave_wo
- wave_g
- wave
- kwave
- x napsp
- x napst
- x napsrh
- napso3
- x napssh
- surft

"""

naaps_press   = nc_data['NAAPS_pressure'][:,:] / 100.

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  If the user wants to screen out any extra NAAPS layer info
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
l_remove_doubles = True
if(l_remove_doubles):
    unique_press = np.unique(naaps_press[0,:])[::-1]
    keep_idxs = np.array([\
        np.where(naaps_press[0,:] == tup)[0][0] \
        for tup in unique_press])
    naaps_press   = naaps_press[:,keep_idxs]
else:     
    keep_idxs = np.arange(naaps_press.shape[1])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  REMEMBER TO CONFIRM THAT THE UNITS MATCH THOSE USED IN THE CODE
#  Convert the mass concentrations from ug/m3 to kg/m3 by  multiplying
#  each value by 1e-9.
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Extract metadata
# ----------------
split_date    = nc_data['time'].units.strip().split()
date_str      = ' '.join(split_date[2:])
base_date     = datetime.strptime(date_str, '%Y%m%d %H:%M:%S')
time          = nc_data['time'][:]
z             = nc_data['z'][keep_idxs]
naaps_dz      = nc_data['NAAPS_dz'][:,keep_idxs]
lat           = nc_data['lat'][:]
lon           = nc_data['lon'][:]

# Extract NAAPS metorological variables
# -------------------------------------
ext_532       = nc_data['532_ext'][:,keep_idxs]
aot_532       = nc_data['532_AOT_hi'][:]
naaps_temp    = nc_data['NAAPS_temperature'][:,keep_idxs]
naaps_sfcpres = nc_data['NAAPS_ps'][:]
naaps_spechum = nc_data['NAAPS_q'][:,keep_idxs]
naaps_relhum  = nc_data['NAAPS_rh'][:,keep_idxs]

# Extract NAAPS aerosol variables
# -------------------------------
naaps_masl      = nc_data['NAAPS_MASL'][:,keep_idxs]

naaps_tot_aot   = nc_data['NAAPS_total_aot'][:]

naaps_dust_ext  = nc_data['NAAPS_dust_ext'][:,keep_idxs]
naaps_dust_aot  = nc_data['NAAPS_dust_aot'][:]
naaps_dust_mass = \
    nc_data['NAAPS_dust_mass_concentration'][:,keep_idxs] * 1e-9

naaps_abf_ext   = nc_data['NAAPS_ABF_ext'][:,keep_idxs]
naaps_abf_aot   = nc_data['NAAPS_ABF_aot'][:]
naaps_abf_mass  = \
    nc_data['NAAPS_ABF_mass_concentration'][:,keep_idxs] * 1e-9

naaps_salt_ext  = nc_data['NAAPS_SS_ext'][:,keep_idxs]
naaps_salt_aot  = nc_data['NAAPS_SS_aot'][:]
naaps_salt_mass = \
    nc_data['NAAPS_SS_mass_concentration'][:,keep_idxs] - 1e-9

naaps_smoke_ext  = nc_data['NAAPS_smoke_ext'][:,keep_idxs]
naaps_smoke_aot  = nc_data['NAAPS_smoke_aot'][:]
naaps_smoke_mass = \
    nc_data['NAAPS_smoke_mass_concentration'][:,keep_idxs] * 1e-9

# Extract lidar aerosol variables
# -------------------------------

ext_532 = np.ma.where(ext_532.mask == True, -9., ext_532)

# Close the file
# --------------
#nc_data.close()


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Loop over the data and print the stuff
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Variable descriptions
#  0 - napsp
#  X - lidar_z
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
# 14 - ext_532

#fmt_str = "{0:7.1f} {1:8.1f} {2:8.1f} {3:8.1f} {4:7.2f} " + \
#    "{5:7.3f} {6:10.3e} " + \
#    "{7:10.3e} {8:10.3e} {9:10.3e} {10:10.3e} " + \
#    "{11:10.3e} {12:10.3e} {13:10.3e} {14:10.3e} " + \
#    "{15:10.3e}"
fmt_str = "{0:7.1f} {2:8.1f} {3:8.1f} {4:7.2f} " + \
    "{5:7.3f} {6:10.3e} " + \
    "{7:10.3e} {8:10.3e} {9:10.3e} {10:10.3e} " + \
    "{11:10.3e} {12:10.3e} {13:10.3e} {14:10.3e}\n"

# totaod = [aod1, aod2, aod3, aod4, aod5]
# totaod: [UNUSED, sulfate, dust, smoke, salt]

for ii in range(time.shape[0]):
    #ii = 1567
    foutname = 'test_naaps_file.txt'
    with open(foutname, 'w') as fout:
    
        
        """
        
        Output format
        YYYYMMDDHHMM lat lon tau_aer
         sulfate_aod, dust_aod, smoke_aod, salt_aod
         press NAAPS_z dz temp rhum spechum abf_ext dust_ext smoke_ext salt_ext abf_mas dust_mas smoke_mas salt_mas
         ...
         ...
         ...
        
         
        """
     
   
        # Print the metadata information
        # ------------------------------
        local_dt_date = base_date + timedelta(seconds = time[ii]) 
        local_date = local_dt_date.strftime('%Y%m%d%H%M')
        print(ii, local_dt_date.strftime('%Y%m%d %H:%M:%S'))

        fout.write("{0} {1:8.5f} {2:8.5f} {3:8.5f}\n".format(\
            local_date, lat[ii], lon[ii], naaps_tot_aot[ii]))
        
        fout.write("{0:10.3e} {1:10.3e} {2:10.3e} {3:10.3e}\n".format(\
            naaps_abf_aot[ii], naaps_dust_aot[ii], \
            naaps_smoke_aot[ii], naaps_salt_aot[ii]))
        
        # Print the profile information. NOTE: You'll notice that there are more
        # variables in the print statement here than there are in the output
        # file. Variables can be easily removed by modifing the print format
        # statement above, so, for example, the lidar height (z) is not being
        # printed, even though it is at spot "1" here.
        # ----------------------------------------------------------------------
        for jj in range(ext_532.shape[1]):
            fout.write(fmt_str.format(\
                naaps_press[ii,jj], z[jj], naaps_masl[ii,jj], naaps_dz[ii,jj], \
                naaps_temp[ii,jj], naaps_relhum[ii,jj], \
                naaps_spechum[ii,jj], \
                naaps_abf_ext[ii,jj], naaps_dust_ext[ii,jj], \
                naaps_smoke_ext[ii,jj], naaps_salt_ext[ii,jj], \
                naaps_abf_mass[ii,jj], naaps_dust_mass[ii,jj], \
                naaps_smoke_mass[ii,jj], naaps_salt_mass[ii,jj], \
                ext_532[ii,jj]
            ))
    
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        #
        #  Call the FuLiou code
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        #print("Calling FuLiou code") 

# Close the file object
# ---------------------
print("Closing file object")
nc_data.close()
