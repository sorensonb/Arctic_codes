#!/usr/bin/env python

"""


"""

from VIIRSLib import *

#
# Nighttime
#   22 July (203): 1000
#   23 July (204): 0942
# Daytime
#   22 July (203): 2124
#   23 July (204): None yet
#

plot_VIIRS_sixpanel(save = True)
#plot_VIIRS_ninepanel(save = False)

##!#filename = glob('/home/bsorenson/data/VIIRS/DNB/VNP02DNB.A2021204.0942*.nc')
##!##filename = glob('/home/bsorenson/data/VIIRS/DNB/VNP02MOD.A2021204.0942*.nc')
##!##filename = glob('/home/bsorenson/data/VIIRS/DNB/VNP02DNB.A2021204*')
##!##filename = glob('/home/bsorenson/data/VIIRS/DNB/VNP46A1.A2021204*')
##!##plot_VIIRS_DNB(filename, band = 'M05')
##!#plot_VIIRS_figure(filename, band = 'M11', zoom = True)
