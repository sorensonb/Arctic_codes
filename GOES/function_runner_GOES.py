#!/usr/bin/env python

"""


"""

from GOESLib import *

date_str = ['202107201200',\
            '202107201500',\
            '202107201800',\
            '202107202100',\
            '202107210000', \
            '202107210300',\
            '202107210600',\
            '202107210900',\
            '202107211200',\
            '202107211500',\
            '202107211800',\
            '202107212100',\
            '202107220000',\
            '202107220300'] 
channel = 2
#plot_GOES_satpy(date_str, channel, ax = None, zoom=True,save=False)
for dstr in date_str:
    plot_GOES_satpy_6panel(dstr, 2, 6, 13, 8, 9, 10, zoom = True, save = True)
