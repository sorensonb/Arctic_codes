#!/usr/bin/env python

"""


"""

import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import sys

if(len(sys.argv) != 2):
    print("SYNTAX: ./auto_sbdart.py satellite")
    print("      satellite: modis_ch31, goes17_ch08, goees17_ch10")
    sys.exit()

sb_path = '/home/bsorenson/Research/MODIS/obs_smoke_forcing/SBDART/SBDART/'
infile = sb_path + 'src/atms_orig.dat'

satellite = sys.argv[1]
#satellite = 'modis_ch31'

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Filter function prep
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Copy the filter function to working directory
# ---------------------------------------------
os.system('cp ' + sb_path + 'data/filter/' + satellite + '_filter.dat ./filter.dat')

# Modify the filter file
# ----------------------
with open('filter.dat','r') as fin:
    indata = fin.readlines()
with open('filter.dat','w') as fout:
    fout.writelines(indata[4:])

# Read with pandas and convert to microns
# ---------------------------------------
ff_data = pd.read_csv('filter.dat', names = ['wavenum','filter'], \
    delim_whitespace = True)
ff_data['wavenum'] = (1. / (ff_data['wavenum'] * 100.)) * 1e6
ff_data = ff_data.reindex(index = ff_data.index[::-1])
ff_data.to_csv('filter.dat', \
        index = False, sep = ' ', header = None)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Input file prep
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

##isat = -1
##
##with open('INPUT','w') as fout:
##    fout.write(" &INPUT\n" + \
##        "   idatm=0\n" + \
##        "   isat={0}\n".format(isat) + \
##        "   wlinf=8.00\n" + \
##        "   wlsup=12.00\n" + \
##        "   wlinc=0.01\n" + \
##        "   iout=1\n" + \
##        "/")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# SBDART runs
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Read the atmos file
# -------------------
data = pd.read_csv(infile,skiprows = 0, names = \
    ['z','p','t','wv','o2'], delim_whitespace = True)

multipliers = np.arange(0.2,1.5,0.2)

outfile_list = []

for mtp in multipliers:
    print("processing for wv multiplier = ",np.round(mtp,2))

    # Make copies for editing
    # -----------------------
    new_data = data.copy(deep = True)
    
    # Multiply lower moisture values by multiplier
    # --------------------------------------------
    new_data['wv'][new_data['z'] <= 5] *= mtp
   
    print(len(new_data['wv']))
    header_line = [str(len(new_data['wv'])),' ',' ',' ',' ']
 
    # Save to outfile
    # ---------------
    outname = 'atms.dat'
    #outname = 'atms_' + str(int(mtp*100)) + '.dat'
    new_data.to_csv(outname, \
        index = False, sep = ' ', header = header_line)
    print("Saved file",outname)

    # Run SBDART with the new atmosphere profile
    # ------------------------------------------
    outfile = 'sbout_' + str(int(mtp*100)) + '.dat'
    cmnd = sb_path + 'bin/sbdart > ' + outfile
    print(cmnd)
    os.system(cmnd)

    # Modify the output file
    # ----------------------
    with open(outfile,'r') as fin:
        indata = fin.readlines()
    with open(outfile,'w') as fout:
        fout.writelines(indata[3:])

    outfile_list.append(outfile)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Data reading
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Read in the SBDART output 
# -------------------------
outfile_list = glob.glob(sb_path + 'src/sbout_*.dat')

names = ['wavelength','filter','down_sol_flux_toa','up_rad_flux_toa',\
    'dir_sol_flux_toa','down_rad_flux_sfc','up_rad_flux_sfc',\
    'dir_sol_flux_sfc']

data_dict = {}
for ofile in outfile_list:
    wv_pct = int(ofile.split('.')[0].split('_')[-1])
    data_dict[wv_pct] = {}
    data_dict[wv_pct]['data'] = pd.read_csv(ofile,skiprows = 0, names = names,\
            delim_whitespace = True, header = None)
    data_dict[wv_pct]['data'] = data_dict[wv_pct]['data'].set_index('wavelength')

    # Calculated weighted average wavelength
    # --------------------------------------
    data_dict[wv_pct]['avg_wavel'] = \
        np.sum(data_dict[wv_pct]['data'].index * data_dict[wv_pct]['data']['filter']) / \
        np.sum(data_dict[wv_pct]['data']['filter'])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Calculations
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
h_const = 6.626e-34 #J*s
k_const = 1.3806e-23 #J/K
c_const = 3e8 # m/s


# wavel in microns
# bght_tmp in (W/(m2*μm*Sr))
def inverse_planck(bght_tmp, wavel):
    return ((h_const * c_const) / (wavel*1e-6 * k_const)) * \
        (np.log( ((2.*h_const * (c_const**2.)) / (wavel * (wavel * 1e-6)**4. \
        * bght_tmp)) + 1.))**-1.

# radiance is a function that calculates the blackbody radiance for a given 
# wavelength (in meters) and temperature (in Kelvin)
# -------------------------------------------------------------------------
def radiance(wavelength,tmp):
    return ((2.*h_const*(c_const**2.)) / \
           ((wavelength**5.)*(np.exp((h_const*c_const)/\
            (wavelength*k_const*tmp)) - 1.))) * \
           1e-6

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Output visualization
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Plot stuff
# ----------
plt.close('all')
fig = plt.figure()
ax0 = fig.add_subplot(1,1,1)

plot_val = 'up_rad_flux_toa'

for key in sorted(data_dict.keys()):
    xvals = pd.to_numeric(data_dict[key]['data'].index).values
    yvals = pd.to_numeric(data_dict[key]['data'][plot_val]).values
    ax0.plot(xvals,yvals,label = key) 

ax0.legend()
ax0.set_title(satellite)
plt.show()
