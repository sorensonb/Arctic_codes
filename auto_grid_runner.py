#!/usr/bin/env python
"""


"""

import numpy as np
import matplotlib.pyplot as plt
import gridCERESLib
import sys
from gridCERESLib import readgridCERES,calc_CERES_trend

def make_scatter(CERES_dict_alb_clr,CERES_dict_sw_clr,\
                 CERES_dict_lw_clr,CERES_dict_net_clr):
    dot_size=12
    # Generate scatter plots
    title_text = 'January 2001 - December 2018'
    fig1 = plt.figure()
    plt.scatter(CERES_dict_alb_clr['trends'].flatten(),\
        CERES_dict_sw_clr['trends'].flatten(),s=dot_size)
    plt.title(title_text)
    plt.xlabel('Terra CERES Clear-Sky Albedo Trend')
    plt.ylabel('Terra CERES Clear-Sky SWF Trend')
    outname = 'ceres_clralb_clrsw_trend_scatter.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)
    
    fig2 = plt.figure()
    plt.scatter(CERES_dict_alb_clr['trends'].flatten(),\
        CERES_dict_lw_clr['trends'].flatten(),s=dot_size)
    plt.title(title_text)
    plt.xlabel('Terra CERES Clear-Sky Albedo Trend')
    plt.ylabel('Terra CERES Clear-Sky LWF Trend')
    outname = 'ceres_clralb_clrlw_trend_scatter.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)
    
    fig3 = plt.figure()
    plt.scatter(CERES_dict_alb_clr['trends'].flatten(),\
        CERES_dict_net_clr['trends'].flatten(),s=dot_size)
    plt.title(title_text)
    plt.xlabel('Terra CERES Clear-Sky Albedo Trend')
    plt.ylabel('Terra CERES Clear-Sky Net Flux Trend')
    outname = 'ceres_clralb_clrnet_trend_scatter.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)


if(len(sys.argv)!=5):
    print("SYNTAX: ./auto_grid_runner.py start_date end_date minlat season")
    print("        Ex: start_date = '200101'")
    print("            end_date   = '201812'")
    print("            minlat     = 45.5    ")
    print("            season     = 'spring','summer','autumn','winter','all'")
    sys.exit()

start_date= sys.argv[1]
end_date  = sys.argv[2]
minlat    = float(sys.argv[3])
season    = sys.argv[4]

# Read in data for all wanted data_types

#CERES_dict_sw_all  = readgridCERES(start_date,end_date,'toa_sw_all_mon', minlat=minlat,season=season)
#CERES_dict_sw_clr  = readgridCERES(start_date,end_date,'toa_sw_clr_mon', minlat=minlat,season=season)
#CERES_dict_lw_all  = readgridCERES(start_date,end_date,'toa_lw_all_mon', minlat=minlat,season=season)
#CERES_dict_lw_clr  = readgridCERES(start_date,end_date,'toa_lw_clr_mon', minlat=minlat,season=season)
#CERES_dict_net_all = readgridCERES(start_date,end_date,'toa_net_all_mon',minlat=minlat,season=season)
#CERES_dict_net_clr = readgridCERES(start_date,end_date,'toa_net_clr_mon',minlat=minlat,season=season)
#CERES_dict_alb_all = readgridCERES(start_date,end_date,'toa_alb_all_mon',minlat=minlat,season=season)
CERES_dict_alb_clr = readgridCERES(start_date,end_date,'toa_alb_clr_mon',minlat=minlat,season=season)

# Plot all of them
#calc_CERES_trend(CERES_dict_sw_all,minlat=minlat, adjusted=True,save=True)
#calc_CERES_trend(CERES_dict_sw_clr,minlat=minlat, adjusted=True,save=True)
#calc_CERES_trend(CERES_dict_lw_all,minlat=minlat, adjusted=True,save=True)
#calc_CERES_trend(CERES_dict_lw_clr,minlat=minlat, adjusted=True,save=True)
#calc_CERES_trend(CERES_dict_net_all,minlat=minlat,adjusted=True,save=True)
#calc_CERES_trend(CERES_dict_net_clr,minlat=minlat,adjusted=True,save=True)
#calc_CERES_trend(CERES_dict_alb_all,minlat=minlat,adjusted=True,save=True)
calc_CERES_trend(CERES_dict_alb_clr,minlat=minlat,adjusted=True,save=True)

#calc_CERES_trend(CERES_dict_sw_all, minlat=minlat,save=True)
#calc_CERES_trend(CERES_dict_sw_clr, minlat=minlat,save=True)
#calc_CERES_trend(CERES_dict_lw_all, minlat=minlat,save=True)
#calc_CERES_trend(CERES_dict_lw_clr, minlat=minlat,save=True)
#calc_CERES_trend(CERES_dict_net_all,minlat=minlat,save=True)
#calc_CERES_trend(CERES_dict_net_clr,minlat=minlat,save=True)
#calc_CERES_trend(CERES_dict_alb_all,minlat=minlat,save=True)
#calc_CERES_trend(CERES_dict_alb_clr,minlat=minlat,save=True)


#make_scatter(CERES_dict_alb_clr,CERES_dict_sw_clr,\
#             CERES_dict_lw_clr,CERES_dict_net_clr)
