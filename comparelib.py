#!/usr/bin/env python
"""


"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats
from IceLib import read_ice,ice_trendCalc,grid_data
from gridCERESLib import readgridCERES,calc_CERES_trend

def correlation(x,y):
    avg_x = np.average(x)
    avg_y = np.average(y)
    N = len(x)
    if(len(x)!=len(y)):
        print("ERROR: Arrays are not the same size.\nArray x has len=",len(x),\
              "\nArray y has len=",len(y))
    else:
        cov = (np.sum((x-avg_x)*(y-avg_y)))/(N-1)
        std_x = np.std(x)
        std_y = np.std(y)
        corr = cov/(std_x*std_y)
        return(corr)

def plot_fourscatter(ice_data,CERES_lw_dict,CERES_sw_dict,CERES_alb_dict,CERES_net_dict,inseason,dtype):
    if(dtype=='clr'):
        dtype_adder = "clear"
    else:
        dtype_adder = "all"
    # Quick comparison
    markersize=6
    fig, axs = plt.subplots(2,2)
    fig.set_size_inches(15,13)
    plt.title("CERES and NSIDC "+inseason.title()+" Trends")
    # Plot LWF data
    plotx = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    ploty = CERES_lw_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    plotx = plotx[ploty>-1000]
    ploty = ploty[ploty>-1000]
    corr_lwf = correlation(plotx,ploty)
    print("LWF correlation: ",corr_lwf)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    axs[0,0].scatter(plotx,ploty,s=markersize,color='black')
    axs[0,0].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
    axs[0,0].set_title(dtype_adder.title()+'-Sky LWF'+ice_data['season_adder'].title()+' Trends')
    axs[0,0].set_xlabel('Ice Trends')
    axs[0,0].set_ylabel('LWF Trends')
    x_pos = ((axs[0,0].get_xlim()[1]-axs[0,0].get_xlim()[0])/8)+axs[0,0].get_xlim()[0]
    y_pos = ((axs[0,0].get_ylim()[1]-axs[0,0].get_ylim()[0])/8)+axs[0,0].get_ylim()[0]
    axs[0,0].text(x_pos,y_pos,"Corr = "+str(np.round(corr_lwf,3)))
    ##print("  y_max = ",axs[0,0].get_ylim()[1],"  y_min = ",axs[0,0].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    # Plot SWF data
    plotx = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    ploty = CERES_sw_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    corr_swf = correlation(plotx,ploty)
    print("SWF correlation: ",corr_swf)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    axs[0,1].scatter(plotx,ploty,s=markersize,color='black')
    axs[0,1].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
    axs[0,1].set_title(dtype_adder.title()+'-Sky SWF'+ice_data['season_adder'].title()+' Trends')
    axs[0,1].set_xlabel('Ice Trends')
    axs[0,1].set_ylabel('SWF Trends')
    x_pos = ((axs[0,1].get_xlim()[1]-axs[0,1].get_xlim()[0])/8)+axs[0,1].get_xlim()[0]
    y_pos = axs[0,1].get_ylim()[1]-((axs[0,1].get_ylim()[1]-axs[0,1].get_ylim()[0])/8)
    axs[0,1].text(x_pos,y_pos,"Corr = "+str(np.round(corr_swf,3)))
    ##print("  y_max = ",axs[0,1].get_ylim()[1],"  y_min = ",axs[0,1].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    # Plot Albedo Data
    plotx = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    ploty = CERES_alb_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    corr_alb = correlation(plotx,ploty)
    print("Albedo correlation: ",corr_alb)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    axs[1,0].scatter(plotx,ploty,s=markersize,color='black')
    axs[1,0].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
    axs[1,0].set_title(dtype_adder.title()+'-Sky Albedo'+ice_data['season_adder'].title()+' Trends')
    axs[1,0].set_xlabel('Ice Trends')
    axs[1,0].set_ylabel('Albedo Trends')
    x_pos = ((axs[1,0].get_xlim()[1]-axs[1,0].get_xlim()[0])/8)+axs[1,0].get_xlim()[0]
    y_pos = axs[1,0].get_ylim()[1]-((axs[1,0].get_ylim()[1]-axs[1,0].get_ylim()[0])/8)
    axs[1,0].text(x_pos,y_pos,"Corr = "+str(np.round(corr_alb,3)))
    ##print("  y_max = ",axs[1,0].get_ylim()[1],"  y_min = ",axs[1,0].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    # Plot Net Flux Data
    plotx = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    ploty = CERES_net_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    corr_net = correlation(plotx,ploty)
    print("Net flux correlation: ",corr_net)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    axs[1,1].scatter(plotx,ploty,s=markersize,color='black')
    axs[1,1].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
    axs[1,1].set_title(dtype_adder.title()+'-Sky Net Flux'+ice_data['season_adder'].title()+' Trends')
    axs[1,1].set_xlabel('Ice Trends')
    axs[1,1].set_ylabel('Net Flux Trends')
    x_pos = ((axs[1,1].get_xlim()[1]-axs[1,1].get_xlim()[0])/8)+axs[1,1].get_xlim()[0]
    y_pos = ((axs[1,1].get_ylim()[1]-axs[1,1].get_ylim()[0])/8)+axs[1,1].get_ylim()[0]
    axs[1,1].text(x_pos,y_pos,"Corr = "+str(np.round(corr_net,3)))
    ##print("  y_max = ",axs[1,1].get_ylim()[1],"  y_min = ",axs[1,1].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    outname = "ceres_nsidc_trends_four_panel_"+inseason+"_"+dtype+"sky.png"
    plt.savefig(outname,dpi=300)
    print("Saved image "+outname)

