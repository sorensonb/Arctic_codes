#!/usr/bin/env python
"""
  NAME:

  PURPOSE:

  SYNTAX:

  EXAMPLE:

  MODIFICATIONS:

"""

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color
import cartopy.crs as ccrs
import subprocess
import h5py

if(len(sys.argv)<2):
    print("SYNTAX: python plot_single_OMI.py filename")
    sys.exit()

filename = sys.argv[1]

# This is the path that points to the HDF5 OMI files. This must be changed
# if running on a new system.
data_path = '/home/bsorenson/data/MODIS/'


