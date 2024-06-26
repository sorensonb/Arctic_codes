#!/usr/bin/env python

"""


"""

import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import sys

h_const = 6.626e-34 #J*s
k_const = 1.3806e-23 #J/K
c_const = 3e8 # m/s

# wavel in microns
# radiance in (W/(m2*μm*Sr))
def inverse_planck(radiance, wavel):
    return ((h_const * c_const) / (wavel*1e-6 * k_const)) * \
        (np.log( ((2.*h_const * (c_const**2.)) / (wavel * (wavel * 1e-6)**4. \
        * radiance)) + 1.))**-1.

# radiance is a function that calculates the blackbody radiance for a given 
# wavelength (in meters) and temperature (in Kelvin)
# -------------------------------------------------------------------------
def radiance(wavelength,tmp):
    return ((2.*h_const*(c_const**2.)) / \
           ((wavelength**5.)*(np.exp((h_const*c_const)/\
            (wavelength*k_const*tmp)) - 1.))) * \
           1e-6


if(len(sys.argv) != 2):
    print("SYNTAX: ./auto_sbdart.py satellite")
    print("      satellite: modis_ch31, goes17_ch08, goees17_ch10")
    sys.exit()

#sb_path = '/home/bsorenson/Research/MODIS/obs_smoke_forcing/SBDART/SBDART/'
infile = sb_path + 'src/atms_orig.dat'

satellite = sys.argv[1]
calc_radiance = True 
run = False
#satellite = 'modis_ch31'

def run_sbdart(satellite, calc_radiance, run = True):
    sb_path = '/home/bsorenson/Research/SBDART/'
    infile = sb_path + 'src/atms_orig.dat'

  
    if(calc_radiance):
        file_adder = '_rad.dat'
    else:
        file_adder = '.dat'
 
    if(run): 
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
        
        isat = -1
        iout = 5
        #iout = 20
        #iout = 0
        
        ##with open('INPUT','w') as fout:
        ##    fout.write(" &INPUT\n" + \
        ##        "   idatm=0\n" + \
        ##        "   isat={0}\n".format(isat) + \
        ##        #"   wlinf=8.00\n" + \
        ##        #"   wlsup=12.00\n" + \
        ##        "   wlinc=0.01\n" + \
        ##        "   iout={0}\n".format(iout) + \
        ##        "   uzen=0\n" + \
        ##        #"   vzen=180\n" + \
        ##        "   phi=0\n" + \
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
        
        multipliers = np.array([0.2, 1.0, 1.5, 2.0, 2.5])
        #multipliers = np.arange(0.1,2.1,0.1)
        
        outfile_list = []
        
        for mtp in multipliers:
            print("processing for wv multiplier = ",np.round(mtp,2))
        
            # Make copies for editing
            # -----------------------
            new_data = data.copy(deep = True)
            
            # Multiply lower moisture values by multiplier
            # --------------------------------------------
            new_data['wv'][new_data['z'] <= 5] *= mtp
           
            header_line = [str(len(new_data['wv'])),' ',' ',' ',' ']
         
            # Save to outfile
            # ---------------
            outname = 'atms.dat'
            #outname = 'atms_' + str(int(mtp*100)) + '.dat'
            new_data.to_csv(outname, \
                index = False, sep = ' ', header = header_line)
        
            # Run SBDART with the new atmosphere profile
            # ------------------------------------------
            outfile = 'sbout_' + str(int(mtp*100)) + '_' + satellite + file_adder
            cmnd = sb_path + 'bin/sbdart > ' + outfile
            print(cmnd)
            os.system(cmnd)
        
            # Modify the output file
            # ----------------------
            if(not calc_radiance):
                with open(outfile,'r') as fin:
                    indata = fin.readlines()
                with open(outfile,'w') as fout:
                    fout.writelines(indata[3:])
        
            outfile_list.append(outfile)

    else:
        search_str = sb_path + 'src/sbout_*' + satellite + \
            file_adder
        print(search_str)
        outfile_list = glob.glob(search_str)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Data reading
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    # Read in the SBDART output 
    # -------------------------
    #outfile_list = glob.glob(sb_path + 'src/sbout_*.dat')
    
    names = ['wavelength','filter','down_sol_flux_toa','up_rad_flux_toa',\
        'dir_sol_flux_toa','down_rad_flux_sfc','up_rad_flux_sfc',\
        'dir_sol_flux_sfc']
   
    print(outfile_list)
 
    data_dict = {}
    for ofile in outfile_list:
        if(not calc_radiance):
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
        else:
            wv_pct = int(ofile.split('/')[-1].split('_')[1])
            data_dict[wv_pct] = {}
            data_dict[wv_pct]['data'] = {}
            with open(ofile,'r') as fin: 
                flines = fin.readlines()

                # Extract the number of reports
                # -----------------------------
                num_rep = int(flines[2].strip().split()[0])

                # Pull out the number of VZAs in the file
                # ---------------------------------------
                num_vzen = int(flines[4].strip().split()[1])

                # Extract the VZAs
                # ----------------
                vza = np.array([float(tval) for tval in flines[6].strip().split()])
           
                # Extract the wavelength and filter values
                # ----------------------------------------
                datalines = flines[3:][::4+num_vzen]
                wavel  = np.array([float(dstr.strip().split()[0]) for dstr in datalines])
                fvals  = np.array([float(dstr.strip().split()[1]) for dstr in datalines])
                topdwn = np.array([float(dstr.strip().split()[2]) for dstr in datalines])
                topup  = np.array([float(dstr.strip().split()[3]) for dstr in datalines])
                topdir = np.array([float(dstr.strip().split()[4]) for dstr in datalines])
                botdwn = np.array([float(dstr.strip().split()[5]) for dstr in datalines])
                botup  = np.array([float(dstr.strip().split()[6]) for dstr in datalines])
                botdir = np.array([float(dstr.strip().split()[7]) for dstr in datalines])

                # Declare an array to hold all of the radiances
                # ---------------------------------------------
                total_radiances = np.full((num_rep, num_vzen), np.nan)

                #vzen_lines = flines[6:][:
 
                # Extract the radiances
                # --------------------- 
                for ii in range(num_vzen):
                    total_radiances[:,ii] = np.array([float(tstr.strip()) for tstr in flines[(7+ii):][::4+num_vzen]])
                #radlines = flines[7:][::5]
                #radiances = np.array([float(tstr.strip()) for tstr in radlines])
    
    
                # Insert all the data into the dictionary
                # ---------------------------------------
                data_dict[wv_pct]['vza'] = vza
                data_dict[wv_pct]['data']['wavelength'] = wavel
                data_dict[wv_pct]['data']['filter']     = fvals
                data_dict[wv_pct]['data']['radiance']   = total_radiances
                data_dict[wv_pct]['data']['down_sol_flux_toa'] = topdwn
                data_dict[wv_pct]['data']['up_rad_flux_toa']   = topup
                data_dict[wv_pct]['data']['dir_sol_flux_toa']  = topdwn
                data_dict[wv_pct]['data']['down_rad_flux_sfc'] = botdwn
                data_dict[wv_pct]['data']['up_rad_flux_sfc']   = botup 
                data_dict[wv_pct]['data']['dir_sol_flux_sfc']  = botdir
    
            # Calculated weighted average wavelength
            # --------------------------------------
            data_dict[wv_pct]['avg_wavel'] = \
                np.sum(data_dict[wv_pct]['data']['wavelength'] * data_dict[wv_pct]['data']['filter']) / \
                np.sum(data_dict[wv_pct]['data']['filter'])
    
            # Calculated weighted average radiance  
            # ------------------------------------
            data_dict[wv_pct]['avg_rads'] = \
                np.array([np.sum(data_dict[wv_pct]['data']['radiance'][:,jj] * \
                data_dict[wv_pct]['data']['filter']) / np.sum(data_dict[wv_pct]['data']['filter']) \
                for jj in range(len(data_dict[wv_pct]['vza']))])

            #data_dict[wv_pct]['avg_rad'] = \
            #    np.sum(data_dict[wv_pct]['data']['radiance'] * data_dict[wv_pct]['data']['filter']) / \
            #    np.sum(data_dict[wv_pct]['data']['filter'])

            # Calculate temperatures for all wavelengths
            # ------------------------------------------
            data_dict[wv_pct]['data']['bght_tmp_all'] = \
                np.array([inverse_planck(data_dict[wv_pct]['data']['radiance'][:,jj], data_dict[wv_pct]['data']['wavelength']) \
                for jj in range(len(data_dict[wv_pct]['vza']))]).transpose()

            #data_dict[wv_pct]['data']['bght_tmp_all'] = \
            #    inverse_planck(data_dict[wv_pct]['data']['radiance'], \
            #    data_dict[wv_pct]['data']['wavelength'])

            # Calculate temperature for filter averaged values
            # ------------------------------------------------
            data_dict[wv_pct]['bght_tmp'] = \
                inverse_planck(data_dict[wv_pct]['avg_rads'], \
                data_dict[wv_pct]['avg_wavel'])

    return data_dict

# Run SBDART for the different channels
# -------------------------------------
modis31 = run_sbdart('modis_ch31', calc_radiance, run = True)
modis32 = run_sbdart('modis_ch32', calc_radiance, run = True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Calculations
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def plot_bright(modis31, modis32, vza_idx = 0, save = False):

    # Pull out the brightness temperatures for each water vapor amount
    # ----------------------------------------------------------------
    m31_bghts = np.array([modis31[key]['bght_tmp'][vza_idx] for key in modis31.keys()])
    m32_bghts = np.array([modis32[key]['bght_tmp'][vza_idx] for key in modis32.keys()])
    diffs = m31_bghts - m32_bghts
    wvs = modis31.keys()
    
    plt.close('all')
    fig1 = plt.figure(figsize = (10,5))
    ax0 = fig1.add_subplot(1,2,1)
    ax1 = fig1.add_subplot(1,2,2)
    ax0.plot(wvs, m31_bghts, label = 'MODIS Ch 31')
    ax0.plot(wvs, m32_bghts, label = 'MODIS Ch 32')
    ax1.plot(wvs, diffs, label = 'MODIS Ch 32')
    ax0.set_xlabel('Percent of original WV in lowest 5 km')
    ax0.set_ylabel('Brightness temperature [K]')
    ax1.set_xlabel('Percent of original WV in lowest 5 km')
    ax1.set_ylabel('Brightness temperature Difference [K]')

    if(save):
        outname = 'modis_ch31vch32_btmp.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_bright_vza(modis31, save = False):
    plt.close('all')
    fig1 = plt.figure(figsize = (8,4))
    ax0 = fig1.add_subplot(1,1,1)
  
    wvs = [20, 100, 150, 200, 250]
 
    for wv in wvs:
        # Plot the correct wv amounts
        ax0.plot(modis31[wv]['vza'],  modis31[wv]['bght_tmp'], label = str(wv)+'% WV')

    ax0.set_xlabel('Viewing Zenith Angle [$^{o}$]')
    ax0.set_ylabel('Brightness temperature [K]')
    ax0.set_title('MODIS Channel 31 (11.0 μm)')
    ax0.legend()
    if(save):
        outname = 'sbdart_'+satellite+'_bright_vza.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image",outname)
    else:
        plt.show() 

#plot_bright(modis31, modis32)

sys.exit()

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
