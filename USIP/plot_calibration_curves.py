#1/usr/bin/env python


"""


"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

if(len(sys.argv) != 2):
    print("SYNTAX: ./plot_calibration_curves.py calib_file")
    sys.exit()

# Open the calibration data
# -------------------------

