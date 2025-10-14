#!/usr/bin/env python

"""


"""

import numpy as np
import sys

sys.path.append('/home/bsorenson')
from sounding_lib import *
from python_lib import *
from cn2_lib import *
from compare_cn2 import *
from datetime import datetime

date_str = '2019050403'
#plot_combined_figure(date_str, save = True)
#plot_combined_figure(date_str, save = True, show_synth_data = False)

plot_combined_figure_v2(date_str, fourpanel = True, save = True)
sys.exit()
