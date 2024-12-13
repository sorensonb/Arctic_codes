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

date_str = '2018050505'
#date_str = '2019050403'
plot_dual_vert_tempdiffs(date_str, save = True)
