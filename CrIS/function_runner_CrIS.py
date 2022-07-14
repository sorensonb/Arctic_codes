#!/usr/bin/env python

"""


"""

from CrISLib import *

date_str = '202107231030'
CrIS_data = readCrIS_granule(date_str, satellite = 'JPSS', resolution = 'NSR')
plot_CrIS_granule_test(CrIS_data, zoom = True, save = False)
