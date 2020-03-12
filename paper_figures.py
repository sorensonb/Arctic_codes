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
import glob
import subprocess
from scipy import stats
sys.path.append('/home/bsorenson/Research/Ice_analysis/')
sys.path.append('/home/bsorenson/Research/CERES/')
from IceLib import read_ice,ice_trendCalc,grid_data,grid_data_conc,ice_gridtrendCalc
from gridCERESLib import readgridCERES,calc_CERES_trend
from comparelib import plot_fourscatter, plot_cld_clr_scatter,figure_1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Read in the ice and CERES data
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

inseason='summer'
adjusted=adj=True

# Read in NSIDC ice data
ice_data = read_ice(inseason)
ice_data = grid_data_conc(ice_data)
ice_data = ice_gridtrendCalc(ice_data)
ice_data = ice_trendCalc(ice_data)
ice_data = grid_data(ice_data)

# Read in CERES data
CERES_lw_clr_dict  = readgridCERES(200012,201812, 'toa_lw_clr_mon',minlat=30.5,season=inseason)
CERES_sw_clr_dict  = readgridCERES(200012,201812, 'toa_sw_clr_mon',minlat=30.5,season=inseason)
CERES_alb_clr_dict = readgridCERES(200012,201812,'toa_alb_clr_mon',minlat=30.5,season=inseason)
CERES_net_clr_dict = readgridCERES(200012,201812,'toa_net_clr_mon',minlat=30.5,season=inseason)
# Read in the clear data
CERES_lw_all_dict  = readgridCERES(200012,201812, 'toa_lw_all_mon',minlat=30.5,season=inseason)
CERES_sw_all_dict  = readgridCERES(200012,201812, 'toa_sw_all_mon',minlat=30.5,season=inseason)
CERES_alb_all_dict = readgridCERES(200012,201812,'toa_alb_all_mon',minlat=30.5,season=inseason)
CERES_net_all_dict = readgridCERES(200012,201812,'toa_net_all_mon',minlat=30.5,season=inseason)

# Calculate CERES trends
calc_CERES_trend(CERES_lw_clr_dict, adjusted=adj,save=True)
calc_CERES_trend(CERES_sw_clr_dict, adjusted=adj,save=True)
calc_CERES_trend(CERES_alb_clr_dict,adjusted=adj,save=True)
calc_CERES_trend(CERES_net_clr_dict,adjusted=adj,save=True)
calc_CERES_trend(CERES_lw_all_dict, adjusted=adj,save=True)
calc_CERES_trend(CERES_sw_all_dict, adjusted=adj,save=True)
calc_CERES_trend(CERES_alb_all_dict,adjusted=adj,save=True)
calc_CERES_trend(CERES_net_all_dict,adjusted=adj,save=True)


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

figure_1(ice_data,CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_alb_clr_dict,\
             CERES_net_clr_dict,CERES_lw_all_dict,CERES_sw_all_dict,\
             CERES_alb_all_dict,CERES_net_all_dict,adjusted=True)
##colormap = plt.cm.bwr
##if(adjusted==True):
##    mapcrs = ccrs.NorthPolarStereo(central_longitude=45.)
##else:
##    mapcrs = ccrs.NorthPolarStereo()
##datacrs = ccrs.PlateCarree()
##fig = plt.figure(1, figsize=(22,15))
##gs = gridspec.GridSpec(nrows=1, ncols=3, hspace = 0.03, wspace = 0.03)
##sw_clr_min_val = -40.
##sw_clr_max_val = 40.
##lw_clr_min_val=-15.
##lw_clr_max_val = 15.
##title_adder = 'TOA Clear-Sky SWF'
### - - - - - - - - - - - - - 
### Plot the SW clr data
### - - - - - - - - - - - - - 
##ax1 = plt.subplot(gs[0],projection=mapcrs)
##ax1.set_extent([-180,180,45,90],datacrs)
##ax1.gridlines()
###ax = plt.axes(projection=ccrs.Miller())
##file_season = '_'+CERES_sw_clr_dict['season']
##title_season = ' '+CERES_sw_clr_dict['season']
##if(CERES_sw_clr_dict['season']=='all'):
##    file_season=''
##    title_season=''
##    
##mesh = ax1.pcolormesh(CERES_sw_clr_dict['lon'],CERES_sw_clr_dict['lat'],CERES_sw_clr_dict['trends'],\
##        transform=datacrs,vmin=sw_clr_min_val,vmax=sw_clr_max_val,cmap=colormap)
##ax1.set_title('Terra CERES '+title_adder+title_season.title()+' Trend\n'+\
##             CERES_sw_clr_dict['dates'][0]+' - '+CERES_sw_clr_dict['dates'][-1])
##ender = '.png'
##if(adjusted==True):
##    ax1.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
##    axins = inset_axes(ax1,width="5%",height="100%",loc='lower left',\
##                bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax1.transAxes,\
##                borderpad=0)
##    cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=CERES_sw_clr_dict['parm_unit'])
##    ax1.set_xlim(-5170748.535086173,5167222.438879491)
##    ax1.set_ylim(-3913488.8763307533,3943353.899053069)
##    ender = '_adjusted'+ender
##else:
##    cbar = plt.colorbar(mesh,cmap=colormap,label=CERES_sw_clr_dict['parm_unit'])
###cbar.set_label(parm_unit)
##ax1.coastlines()
##
### - - - - - - - - - - - - - 
### Plot the LW clr data
### - - - - - - - - - - - - - 
##ax2 = plt.subplot(gs[1],projection=mapcrs)
##ax2.set_extent([-180,180,45,90],datacrs)
##ax2.gridlines()
###ax = plt.axes(projection=ccrs.Miller())
##file_season = '_'+CERES_lw_clr_dict['season']
##title_season = ' '+CERES_lw_clr_dict['season']
##if(CERES_lw_clr_dict['season']=='all'):
##    file_season=''
##    title_season=''
##    
##mesh = ax2.pcolormesh(CERES_lw_clr_dict['lon'],CERES_lw_clr_dict['lat'],CERES_lw_clr_dict['trends'],\
##        transform=datacrs,vmin=lw_clr_min_val,vmax=lw_clr_max_val,cmap=colormap)
##ax2.set_title('Terra CERES '+title_adder+title_season.title()+' Trend\n'+\
##             CERES_lw_clr_dict['dates'][0]+' - '+CERES_lw_clr_dict['dates'][-1])
##ender = '.png'
##if(adjusted==True):
##    ax2.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
##    axins = inset_axes(ax2,width="5%",height="100%",loc='lower left',\
##                bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax2.transAxes,\
##                borderpad=0)
##    cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=CERES_lw_clr_dict['parm_unit'])
##    ax2.set_xlim(-5170748.535086173,5167222.438879491)
##    ax2.set_ylim(-3913488.8763307533,3943353.899053069)
##    ender = '_adjusted'+ender
##else:
##    cbar = plt.colorbar(mesh,cmap=colormap,label=CERES_lw_clr_dict['parm_unit'])
###cbar.set_label(parm_unit)
##ax2.coastlines()
##
##
##lon_ranges  = np.arange(-180.,180.,1.0)
##lat_ranges  = np.arange(30.,90.,1.0)
##
##local_grid_ice = np.copy(ice_data['grid_ice'])
##local_grid_ice_bad = np.copy(ice_data['grid_ice'])
##
##local_grid_ice[local_grid_ice==-999.] = np.nan
##local_grid_ice_bad[local_grid_ice!=-999.] = np.nan
##plot_good_data = ma.masked_invalid(local_grid_ice)
##plot_land_data = ma.masked_invalid(local_grid_ice_bad)
##
##colormap = plt.cm.bwr
###coolwarm = plt.cm.coolwarm
###ax = plt.axes(projection=ccrs.NorthPolarStereo())
##file_adder=''
##ax3 = plt.subplot(gs[2],projection=mapcrs)
##ax3.set_extent([-180,180,45,90],datacrs)
##ax3.gridlines()
##mesh = ax3.pcolormesh(lon_ranges,lat_ranges,plot_good_data,\
##        transform=datacrs,vmin=-50,vmax=50,cmap=colormap)
##if(adjusted==True):
##    ax3.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
##    axins = inset_axes(ax3,width="5%",height="100%",loc='lower left',\
##                bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax3.transAxes,\
##                borderpad=0)
##    cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label='Ice Concentration Trend [%]')
##    ax3.set_xlim(-5170748.535086173,5167222.438879491)
##    ax3.set_ylim(-3913488.8763307533,3943353.899053069)
##    #ender = '_adjusted'+ender
##else:
##    ax3.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
##    cbar = plt.colorbar(mesh,cmap=colormap,label='Ice Concentration Trend [%]')
##ax3.coastlines()
##
##
###if(ice_data['season_adder']!=''):
###    ax.set_title('NSIDC Sea Ice Concentration'+ice_data['season_adder'].title()+\
###        ' Trends\nJan 2001 to Dec 2018')
###else:
###    ax.set_title('NSIDC Sea Ice Concentration Trends\nJan 2001 to Dec 2018')
###if(save==True):
###    outname = 'ice_trend_gridded_200101_201812'+ice_data['file_adder']+file_adder+'.png'
###    plt.savefig(outname,dpi=300)
###    print("Saved image ",outname)
###else:
##plt.show()



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




# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 
# Figure 3:  Positive Feedback Model 
# 
# Five lines, showing the decrease in ice area
# for different ice thicknesses.
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Clear-sky summer values from net flux figure
dflux_dsigma = -1.3183
del_T = (86400.)*90  # number of seconds of solar heating per summer
l_f = 3.3e5  # latent heat of fusion (J/kg)
rho_i = 917  # density of ice (kg/m3)
def equation(sigma_p,thick):

    sigma_pnew = sigma_p * (1. - (dflux_dsigma * del_T) / (l_f * rho_i* thick))
    return sigma_pnew

def fill_array(start_val,thick):
    sigma_values = []
    sigma_values.append(start_val)
    # Find 1m-thickness values
    while(start_val>-100):
        sigma_new = equation(start_val,thick)
        sigma_values.append(sigma_new) 
        start_val = sigma_new
    return sigma_values

sum_init_ice  = 7820576049633.293/1e6  # km2 as of June 2001
sum_final_ice = 6481114233056.605/1e6  # km2 

# 10 years
#starting_val = -10. 
starting_val = ((sum_final_ice-sum_init_ice)/sum_init_ice)*100.
start_year = 2018

sigma_1m = np.array(fill_array(starting_val,1.))
sigma_2m = np.array(fill_array(starting_val,2.))
sigma_3m = np.array(fill_array(starting_val,3.))
sigma_4m = np.array(fill_array(starting_val,4.))
sigma_5m = np.array(fill_array(starting_val,5.))

extents_1m = (sum_init_ice*(100.+sigma_1m)/100.)/1e6 
years_1m   = np.arange(start_year,start_year+len(extents_1m))
extents_2m = (sum_init_ice*(100.+sigma_2m)/100.)/1e6 
years_2m   = np.arange(start_year,start_year+len(extents_2m))
extents_3m = (sum_init_ice*(100.+sigma_3m)/100.)/1e6 
years_3m   = np.arange(start_year,start_year+len(extents_3m))
extents_4m = (sum_init_ice*(100.+sigma_4m)/100.)/1e6 
years_4m   = np.arange(start_year,start_year+len(extents_4m))
extents_5m = (sum_init_ice*(100.+sigma_5m)/100.)/1e6 
years_5m   = np.arange(start_year,start_year+len(extents_5m))

print("1m depth gone at year:",years_1m[-1])
print("2m depth gone at year:",years_2m[-1])
print("3m depth gone at year:",years_3m[-1])
print("4m depth gone at year:",years_4m[-1])
print("5m depth gone at year:",years_5m[-1])

fig1 = plt.figure()
plt.plot(years_1m,extents_1m,label='1 meter thickness')
plt.plot(years_2m,extents_2m,label='2 meter thickness')
plt.plot(years_3m,extents_3m,label='3 meter thickness')
plt.plot(years_4m,extents_4m,label='4 meter thickness')
plt.plot(years_5m,extents_5m,label='5 meter thickness')
plt.legend()
plt.xlabel('Year')
plt.ylabel('Ice Area (millions of km2)')
plt.show()


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 
# Figure 4:  Positive Feedback Model comparisons with observed ice area
# 
#   One line shows the model output profile for 2m thickness. Overlay
#   the observed ice area from NSIDC
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



