#!/usr/bin/env python

"""
  NAME:
    plot_fuliou_lidar_output.py

  PURPOSE:
    Generates a plot of the FuLiou-calculated heating rates given in the
        output file generated by the FuLiou system. The figure consists
        of two panels: the heating rates are plotted in the top panel
        while the lidar (or NAAPS) dust extinctions are plotted in
        the bottom panel. Note that the type of aerosol extinction
        being plotted may be changed by changing line 55 to be
        "abf_ext", "salt_ext", or "smoke_ext", instead of "dust_ext".

  SYNTAX:
    ./plot_fuliou_lidar_output.py fuliou_output_file [save]")
            save: optional argument, saves the image")

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2023/11/15:
      Written

"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import sys

if(len(sys.argv) < 2):
    print("SYNTAX: ./plot_fuliou_lidar_output.py fuliou_output_file [save]")
    print("         save: optional argument, saves the image")
    sys.exit()

l_SAVE = False
if("save" in sys.argv):
    l_SAVE = True 

infile = sys.argv[1]

# Parse the file name to retrieve run information
# -----------------------------------------------
split_file = infile.strip().split('/')[-1].split('.')[0].split('_')
run_type = split_file[2]
plot_title = ' '.join(split_file[3:])


# Open the file and extract variables
# -----------------------------------
data = h5py.File(infile, 'r')

time      = data['time'][:]
mask_heat = data['net_heating_rate'][:,:]
mask_ext  = data['dust_ext'][:,:]
mask_ext = np.ma.masked_where(mask_ext <= 0, mask_ext)
mask_ext = np.log10(mask_ext)

xvals = time
yvals = data['press_lvl'][0,:-1]
yvals2 = data['press_lvl'][0,:]

data.close()

# Plot the heating rate and NAAPS/lidar extinction values
# -------------------------------------------------------
fig = plt.figure(figsize = (9, 6))
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
mesh = ax1.pcolormesh(xvals, yvals, mask_heat.T, shading = 'auto', \
    vmin = -10, vmax = 10, cmap = 'bwr')
cbar = fig.colorbar(mesh, ax = ax1, label = 'Net Heating Rate [K/day]')
ax1.invert_yaxis()
ax1.set_title('Heating Rate')
ax1.set_xlabel('UTC Time [Seconds after midnight]')
ax1.set_ylabel('Pressure [mb]')

mesh = ax2.pcolormesh(xvals, yvals2, mask_ext.T, shading = 'auto', \
    vmin = -6, vmax = -3.0)
cbar = fig.colorbar(mesh, ax = ax2, label = 'Log10 of Dust Extinction')
ax2.invert_yaxis()
ax2.set_title(run_type.upper() + ' Dust extinction')
ax2.set_xlabel('UTC Time [Seconds after midnight]')
ax2.set_ylabel('Pressure [mb]')

fig.suptitle(plot_title)
fig.tight_layout()

if(l_SAVE):
    plot_title = '_'.join(split_file[3:])
    outname = 'fuliou_heating_rates_' + run_type + '_' + plot_title + '.png'
    fig.savefig(outname, dpi = 200)
    print("Saved image", outname)
else:
    plt.show()