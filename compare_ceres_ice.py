#!/usr/bin/env python
"""


"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/bsorenson/Research/Ice_analysis/')
sys.path.append('/home/bsorenson/Research/CERES/')
from IceLib import read_ice,ice_trendCalc,grid_data,grid_data_conc
from gridCERESLib import readgridCERES,calc_CERES_trend
from comparelib import plot_fourscatter, plot_cld_clr_scatter

###def correlation(x,y):
###    avg_x = np.average(x)
###    avg_y = np.average(y)
###    N = len(x)
###    if(len(x)!=len(y)):
###        print("ERROR: Arrays are not the same size.\nArray x has len=",len(x),\
###              "\nArray y has len=",len(y))
###    else:
###        cov = (np.sum((x-avg_x)*(y-avg_y)))/(N-1)
###        std_x = np.std(x)
###        std_y = np.std(y)
###        corr = cov/(std_x*std_y)
###        return(corr)
###
###def plot_fourscatter(CERES_lw_dict,CERES_sw_dict,CERES_alb_dict,CERES_net_dict,inseason):
###    # Quick comparison
###    markersize=6
###    fig, axs = plt.subplots(2,2)
###    fig.set_size_inches(15,13)
###    plt.title("CERES and NSIDC "+inseason.title()+" Trends")
###    # Plot LWF data
###    plotx = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
###    ploty = CERES_lw_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
###    plotx = plotx[ploty>-1000]
###    ploty = ploty[ploty>-1000]
###    corr_lwf = correlation(plotx,ploty)
###    print("LWF correlation: ",corr_lwf)
###    axs[0,0].scatter(plotx,ploty,s=markersize,color='black')
###    axs[0,0].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
###    axs[0,0].set_title('Clear-Sky LWF Trends')
###    axs[0,0].set_xlabel('Ice Trends')
###    axs[0,0].set_ylabel('LWF Trends')
###    x_pos = ((axs[0,0].get_xlim()[1]-axs[0,0].get_xlim()[0])/8)+axs[0,0].get_xlim()[0]
###    y_pos = ((axs[0,0].get_ylim()[1]-axs[0,0].get_ylim()[0])/8)+axs[0,0].get_ylim()[0]
###    axs[0,0].text(x_pos,y_pos,"Corr = "+str(np.round(corr_lwf,3)))
###    print("  x_pos = ",x_pos)
###    print("  y_pos = ",y_pos)
###    # Plot SWF data
###    plotx = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
###    ploty = CERES_sw_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
###    corr_swf = correlation(plotx,ploty)
###    print("SWF correlation: ",corr_swf)
###    axs[0,1].scatter(plotx,ploty,s=markersize,color='black')
###    axs[0,1].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
###    axs[0,1].set_title('Clear-Sky SWF Trends')
###    axs[0,1].set_xlabel('Ice Trends')
###    axs[0,1].set_ylabel('SWF Trends')
###    x_pos = ((axs[0,1].get_xlim()[1]-axs[0,1].get_xlim()[0])/8)+axs[0,1].get_xlim()[0]
###    y_pos = ((axs[0,1].get_ylim()[1]-axs[0,1].get_ylim()[0])/8)-axs[0,1].get_ylim()[1]
###    axs[0,1].text(x_pos,y_pos,"Corr = "+str(np.round(corr_swf,3)))
###    print("  x_pos = ",x_pos)
###    print("  y_pos = ",y_pos)
###    # Plot Albedo Data
###    plotx = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
###    ploty = CERES_alb_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
###    corr_alb = correlation(plotx,ploty)
###    print("Albedo correlation: ",corr_alb)
###    axs[1,0].scatter(plotx,ploty,s=markersize,color='black')
###    axs[1,0].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
###    axs[1,0].set_title('Clear-Sky Albedo Trends')
###    axs[1,0].set_xlabel('Ice Trends')
###    axs[1,0].set_ylabel('Albedo Trends')
###    x_pos = ((axs[1,0].get_xlim()[1]-axs[1,0].get_xlim()[0])/8)+axs[1,0].get_xlim()[0]
###    y_pos = ((axs[1,0].get_ylim()[1]-axs[1,0].get_ylim()[0])/8)-axs[1,0].get_ylim()[1]
###    axs[1,0].text(x_pos,y_pos,"Corr = "+str(np.round(corr_alb,3)))
###    print("  x_pos = ",x_pos)
###    print("  y_pos = ",y_pos)
###    # Plot Net Flux Data
###    plotx = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
###    ploty = CERES_net_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
###    corr_net = correlation(plotx,ploty)
###    print("Net flux correlation: ",corr_net)
###    axs[1,1].scatter(plotx,ploty,s=markersize,color='black')
###    axs[1,1].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
###    axs[1,1].set_title('Clear-Sky Net Flux Trends')
###    axs[1,1].set_xlabel('Ice Trends')
###    axs[1,1].set_ylabel('Net Flux Trends')
###    x_pos = ((axs[1,1].get_xlim()[1]-axs[1,1].get_xlim()[0])/8)+axs[1,1].get_xlim()[0]
###    y_pos = ((axs[1,1].get_ylim()[1]-axs[1,1].get_ylim()[0])/8)+axs[1,1].get_ylim()[0]
###    axs[1,1].text(x_pos,y_pos,"Corr = "+str(np.round(corr_net,3)))
###    print("  x_pos = ",x_pos)
###    print("  y_pos = ",y_pos)
###    outname = "ceres_nsidc_trends_four_panel_"+inseason+".png"
###    plt.savefig(outname,dpi=300)
###    print("Saved image "+outname)
###    plt.show()


if(len(sys.argv)!=3):
    print("SYNTAX: ./compare_ceres_ice.py type season")
    print("        type   = 'clr','cld','all''")
    print("        season = 'spring','summer','autumn','winter','all'")
    sys.exit()

adj=True
dtype = sys.argv[1]
inseason = sys.argv[2]

# Read in NSIDC ice data
ice_data = read_ice(inseason)
ice_data = grid_data_conc(ice_data)
ice_data = ice_trendCalc(ice_data)
ice_data = grid_data(ice_data)
sys.exit()


cloud = False
if(dtype=='cld'):
    cloud=True
    dtype='all'

# Read in CERES data
CERES_lw_dict = readgridCERES(200012,201812,'toa_lw_'+dtype+'_mon',minlat=30.5,season=inseason)
CERES_sw_dict = readgridCERES(200012,201812,'toa_sw_'+dtype+'_mon',minlat=30.5,season=inseason)
CERES_alb_dict = readgridCERES(200012,201812,'toa_alb_'+dtype+'_mon',minlat=30.5,season=inseason)
CERES_net_dict = readgridCERES(200012,201812,'toa_net_'+dtype+'_mon',minlat=30.5,season=inseason)
if(cloud==True):
    # Read in the clear data
    CERES_lw_clr_dict  = readgridCERES(200012,201812, 'toa_lw_clr_mon',minlat=30.5,season=inseason)
    CERES_sw_clr_dict  = readgridCERES(200012,201812, 'toa_sw_clr_mon',minlat=30.5,season=inseason)
    CERES_alb_clr_dict = readgridCERES(200012,201812,'toa_alb_clr_mon',minlat=30.5,season=inseason)
    CERES_net_clr_dict = readgridCERES(200012,201812,'toa_net_clr_mon',minlat=30.5,season=inseason)
    # Calculate cloudy data by subtracting the clear data from the all-sky data 
    # Longwave
    lw_cloud_data = CERES_lw_dict['data'] - CERES_lw_clr_dict['data']
    lw_cloud_data[np.where((CERES_lw_dict['data']==-999) | (CERES_lw_clr_dict['data']==-999))] = -999
    CERES_lw_dict['data'] = lw_cloud_data
    CERES_lw_dict['param'] = 'toa_lw_cld_mon' 
    CERES_lw_dict['parm_name'] = 'Observed TOA Longwave Flux - Cloudy-sky'
    # Shortwave
    sw_cloud_data = CERES_sw_dict['data'] - CERES_sw_clr_dict['data']
    sw_cloud_data[np.where((CERES_sw_dict['data']==-999) | (CERES_sw_clr_dict['data']==-999))] = -999
    CERES_sw_dict['data'] = sw_cloud_data
    CERES_sw_dict['param'] = 'toa_sw_cld_mon' 
    CERES_sw_dict['parm_name'] = 'Observed TOA Shortwave Flux - Cloudy-sky'
    # Albedo
    alb_cloud_data = CERES_alb_dict['data'] - CERES_alb_clr_dict['data']
    alb_cloud_data[np.where((CERES_alb_dict['data']==-999) | (CERES_alb_clr_dict['data']==-999))] = -999
    CERES_alb_dict['data'] = alb_cloud_data
    CERES_alb_dict['param'] = 'toa_alb_cld_mon' 
    CERES_alb_dict['parm_name'] = 'Observed TOA Albedo - Cloudy-sky'
    # Net flux
    net_cloud_data = CERES_net_dict['data'] - CERES_net_clr_dict['data']
    net_cloud_data[np.where((CERES_net_dict['data']==-999) | (CERES_net_clr_dict['data']==-999))] = -999
    CERES_net_dict['data'] = net_cloud_data
    CERES_net_dict['param'] = 'toa_net_cld_mon' 
    CERES_net_dict['parm_name'] = 'Observed TOA Net Flux - Cloudy-sky'
    dtype='cld' 

calc_CERES_trend(CERES_lw_dict,adjusted=adj,save=True)
calc_CERES_trend(CERES_sw_dict,adjusted=adj,save=True)
calc_CERES_trend(CERES_alb_dict,adjusted=adj,save=True)
calc_CERES_trend(CERES_net_dict,adjusted=adj,save=True)


# Reshape the CERES data
templon = CERES_lw_dict['lon']
lon2    = np.where(templon>179.9)
goodlon = np.where(templon<=179.9)
CERES_lw_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
CERES_sw_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
CERES_alb_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
CERES_net_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
for xi in range(len(CERES_lw_dict['lat'])):
    CERES_lw_dict['trends'][xi,:] = np.concatenate([CERES_lw_dict['trends'][xi,:][lon2],CERES_lw_dict['trends'][xi,:][goodlon]])
    CERES_sw_dict['trends'][xi,:] = np.concatenate([CERES_sw_dict['trends'][xi,:][lon2],CERES_sw_dict['trends'][xi,:][goodlon]])
    CERES_alb_dict['trends'][xi,:] = np.concatenate([CERES_alb_dict['trends'][xi,:][lon2],CERES_alb_dict['trends'][xi,:][goodlon]])
    CERES_net_dict['trends'][xi,:] = np.concatenate([CERES_net_dict['trends'][xi,:][lon2],CERES_net_dict['trends'][xi,:][goodlon]])
if(cloud==True):
    CERES_lw_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
    CERES_sw_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
    CERES_alb_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
    CERES_net_clr_dict['lon'] = np.concatenate([templon[lon2]-360.,templon[goodlon]])
    for xi in range(len(CERES_lw_dict['lat'])):
        CERES_lw_clr_dict['trends'][xi,:] = np.concatenate([CERES_lw_clr_dict['trends'][xi,:][lon2],CERES_lw_clr_dict['trends'][xi,:][goodlon]])
        CERES_sw_clr_dict['trends'][xi,:] = np.concatenate([CERES_sw_clr_dict['trends'][xi,:][lon2],CERES_sw_clr_dict['trends'][xi,:][goodlon]])
        CERES_alb_clr_dict['trends'][xi,:] = np.concatenate([CERES_alb_clr_dict['trends'][xi,:][lon2],CERES_alb_clr_dict['trends'][xi,:][goodlon]])
        CERES_net_clr_dict['trends'][xi,:] = np.concatenate([CERES_net_clr_dict['trends'][xi,:][lon2],CERES_net_clr_dict['trends'][xi,:][goodlon]])
    calc_CERES_trend(CERES_lw_clr_dict,adjusted=adj,save=True)
    calc_CERES_trend(CERES_sw_clr_dict,adjusted=adj,save=True)
    calc_CERES_trend(CERES_alb_clr_dict,adjusted=adj,save=True)
    calc_CERES_trend(CERES_net_clr_dict,adjusted=adj,save=True)



# Start comparisons
plot_fourscatter(ice_data,CERES_lw_dict,CERES_sw_dict,CERES_alb_dict,CERES_net_dict,inseason,dtype)

## Quick comparison
#plot_cld_clr_scatter(CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_alb_clr_dict,CERES_net_clr_dict,\
#                     CERES_lw_dict,CERES_sw_dict,CERES_alb_dict,CERES_net_dict,\
#                     ice_data,inseason)

### Loop over the grid, comparing the trends if the ice trends are non-missing.
### Trends over land are missing.
### Use the ice range as the loop bounds; the CERES data should have higher
### extent anyway.
##c_trend = np.full(len(CERES_dict['trends'].flatten()),-99.)
##i_trend = np.full(len(ice_data['grid_ice'].flatten()),-99.)
##count = 0
##for xi in range(ice_data['grid_ice'].shape[0]):
##    for yj in range(ice_data['grid_ice'].shape[1]):
##        if(ice_data['grid_ice'][xi,yj]!=-999.):
##            c_trend[count] = CERES_dict['trends'][xi,yj]   
##            i_trend[count] = ice_data['grid_ice'][xi,yj]   
##            count+=1
