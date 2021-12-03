#!/usr/bin/env python

"""


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

if(len(sys.argv) != 2):
    print("SYNTAX: ./plot_radiances.py sbdart_outfile")
    sys.exit()

infile = sys.argv[1]

plt.close('all')
data = pd.read_csv(infile,skiprows = 0, names = \
    ['wavelength','filter','down_sol_flux_toa','up_rad_flux_toa',\
    'dir_sol_flux_toa','down_rad_flux_sfc','up_rad_flux_sfc',\
    'dir_sol_flux_sfc'], delim_whitespace = True)
data_wave = data.set_index('wavelength')

print(data_wave)

data_wave['dir_sol_flux_toa'].plot(legend = True)
#data_wave['down_sol_flux_toa'].plot()
data_wave['up_rad_flux_toa'].plot(secondary_y = True, legend = True)
plt.show()
