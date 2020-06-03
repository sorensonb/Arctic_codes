#!/usr/bin/env python
"""
  NAME:
    paper_figures.py

  PURPOSE:
    Generate figures for the Arctic positive feedback paper.

"""
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
from scipy import stats
sys.path.append('/home/bsorenson/Research/Ice_analysis/')
sys.path.append('/home/bsorenson/Research/CERES/')
from IceLib import read_ice,ice_trendCalc,grid_data,grid_data_conc,ice_gridtrendCalc
from gridCERESLib import readgridCERES,calc_CERES_trend
from comparelib import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Read in the ice and CERES data
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#inseason='sunlight'
inseason='summer'
adjusted=adj=True

# Read in NSIDC ice data
ice_data = read_ice(inseason)
ice_data = grid_data_conc(ice_data)
ice_data = ice_gridtrendCalc(ice_data)
ice_data = ice_trendCalc(ice_data)
ice_data = grid_data(ice_data)

sys.exit()

# Read in CERES data
CERES_lw_clr_dict  = readgridCERES(200012,201812, 'toa_lw_clr_mon',minlat=30.5,season=inseason)
CERES_sw_clr_dict  = readgridCERES(200012,201812, 'toa_sw_clr_mon',minlat=30.5,season=inseason)
#CERES_alb_clr_dict = readgridCERES(200012,201812,'toa_alb_clr_mon',minlat=30.5,season=inseason)
CERES_net_clr_dict = readgridCERES(200012,201812,'toa_net_clr_mon',minlat=30.5,season=inseason)
# Read in the clear data
CERES_lw_all_dict  = readgridCERES(200012,201812, 'toa_lw_all_mon',minlat=30.5,season=inseason)
CERES_sw_all_dict  = readgridCERES(200012,201812, 'toa_sw_all_mon',minlat=30.5,season=inseason)
#CERES_alb_all_dict = readgridCERES(200012,201812,'toa_alb_all_mon',minlat=30.5,season=inseason)
CERES_net_all_dict = readgridCERES(200012,201812,'toa_net_all_mon',minlat=30.5,season=inseason)

# Calculate CERES trends
calc_CERES_trend(CERES_lw_clr_dict, adjusted=adj,save=True)
calc_CERES_trend(CERES_sw_clr_dict, adjusted=adj,save=True)
#calc_CERES_trend(CERES_alb_clr_dict,adjusted=adj,save=True)
calc_CERES_trend(CERES_net_clr_dict,adjusted=adj,save=True)
calc_CERES_trend(CERES_lw_all_dict, adjusted=adj,save=True)
calc_CERES_trend(CERES_sw_all_dict, adjusted=adj,save=True)
#calc_CERES_trend(CERES_alb_all_dict,adjusted=adj,save=True)
calc_CERES_trend(CERES_net_all_dict,adjusted=adj,save=True)

# Reshape the CERES data
templon = CERES_lw_clr_dict['lon']
lon2    = np.where(templon>179.9)
goodlon = np.where(templon<=179.9)
CERES_lw_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
CERES_sw_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
#CERES_alb_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
CERES_net_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
CERES_lw_all_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
CERES_sw_all_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
#CERES_alb_all_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
CERES_net_all_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
for xi in range(len(CERES_lw_clr_dict['lat'])):
    CERES_lw_clr_dict['trends'][xi,:] = np.concatenate([CERES_lw_clr_dict['trends'][xi,:][lon2],CERES_lw_clr_dict['trends'][xi,:][goodlon]])
    CERES_sw_clr_dict['trends'][xi,:] = np.concatenate([CERES_sw_clr_dict['trends'][xi,:][lon2],CERES_sw_clr_dict['trends'][xi,:][goodlon]])
#    CERES_alb_clr_dict['trends'][xi,:] = np.concatenate([CERES_alb_clr_dict['trends'][xi,:][lon2],CERES_alb_all_dict['trends'][xi,:][goodlon]])
    CERES_net_clr_dict['trends'][xi,:] = np.concatenate([CERES_net_clr_dict['trends'][xi,:][lon2],CERES_net_clr_dict['trends'][xi,:][goodlon]])
    CERES_lw_all_dict['trends'][xi,:] = np.concatenate([CERES_lw_all_dict['trends'][xi,:][lon2],CERES_lw_all_dict['trends'][xi,:][goodlon]])
    CERES_sw_all_dict['trends'][xi,:] = np.concatenate([CERES_sw_all_dict['trends'][xi,:][lon2],CERES_sw_all_dict['trends'][xi,:][goodlon]])
#    CERES_alb_all_dict['trends'][xi,:] = np.concatenate([CERES_alb_all_dict['trends'][xi,:][lon2],CERES_alb_all_dict['trends'][xi,:][goodlon]])
    CERES_net_all_dict['trends'][xi,:] = np.concatenate([CERES_net_all_dict['trends'][xi,:][lon2],CERES_net_all_dict['trends'][xi,:][goodlon]])
##if(cloud==True):
##    CERES_lw_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
##    CERES_sw_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
##    CERES_alb_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
##    CERES_net_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
##    for xi in range(len(CERES_lw_dict['lat'])):
##        CERES_lw_clr_dict['trends'][xi,:] = np.concatenate([CERES_lw_clr_dict['trends'][xi,:][lon2],CERES_lw_clr_dict['trends'][xi,:][goodlon]])
##        CERES_sw_clr_dict['trends'][xi,:] = np.concatenate([CERES_sw_clr_dict['trends'][xi,:][lon2],CERES_sw_clr_dict['trends'][xi,:][goodlon]])
##        CERES_alb_clr_dict['trends'][xi,:] = np.concatenate([CERES_alb_clr_dict['trends'][xi,:][lon2],CERES_alb_clr_dict['trends'][xi,:][goodlon]])
##        CERES_net_clr_dict['trends'][xi,:] = np.concatenate([CERES_net_clr_dict['trends'][xi,:][lon2],CERES_net_clr_dict['trends'][xi,:][goodlon]])
##    calc_CERES_trend(CERES_lw_clr_dict,adjusted=adj,save=True)
##    calc_CERES_trend(CERES_sw_clr_dict,adjusted=adj,save=True)
##    calc_CERES_trend(CERES_alb_clr_dict,adjusted=adj,save=True)

 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Create figures for the paper
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 
# Figure 1:  Spatial trends
# 
# Three-panel (1 row, 3 column)
#     Panel 1: SW spatial trends
#     Panel 2: LW spatial trends
#     Panel 3: spatial ice trends
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#write_toASCII(ice_data,CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_net_clr_dict,\
#                  CERES_lw_all_dict,CERES_sw_all_dict,CERES_net_all_dict,\
#                  nonzero=True)
figure_1(ice_data,CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_net_clr_dict,\
         CERES_lw_all_dict,CERES_sw_all_dict,CERES_net_all_dict,adjusted=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 
# Figure 2: Trend Scatter Comparisons 
# 
# Six-panel (2 row, 3 column)
# Row 1: Clear sky trends 
#     Row 1, Panel 1: ice vs clr SW 
#     Row 1, Panel 2: ice vs clr LW
#     Row 1, Panel 3: ice vs clr net
#
# Row 2: All sky trends 
#     Row 2, Panel 1: ice vs all SW 
#     Row 2, Panel 2: ice vs all LW
#     Row 2, Panel 3: ice vs all net
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

figure_2(ice_data,CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_net_clr_dict,\
         CERES_lw_all_dict,CERES_sw_all_dict,CERES_net_all_dict,inseason,adjusted=True)



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 
# Figure 3:  Positive Feedback Model 
# 
# Five lines, showing the decrease in ice area
# for different ice thicknesses.
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#figure_3()


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 
# Figure 4:  Positive Feedback Model comparisons with observed ice area
# 
#   One line shows the model output profile for 2m thickness. Overlay
#   the observed ice area from NSIDC
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#figure_4(ice_data)
#figure_4(ice_data,model_overlay=False,zoomed=True)


