#!/usr/bin/env python

"""


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

if(len(sys.argv) != 2):
    print("SYNTAX: ./plot_radiances.py atmos_file")
    sys.exit()

infile = sys.argv[1]

# Read the atmos file
# -------------------
data = pd.read_csv(infile,skiprows = 0, names = \
    ['z','p','t','wv','o2'], delim_whitespace = True)

multipliers = np.arange(0.6,1.5,0.2)

for mtp in multipliers:
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
    outname = 'atms_' + str(int(mtp*100)) + '.dat'
    new_data.to_csv(outname, \
        index = False, sep = ' ', header = header_line)

    # Append the length of the file to the beginning of the file
    # ----------------------------------------------------------
    
    print("Saved file",outname)
