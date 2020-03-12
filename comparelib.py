#!/usr/bin/env python
"""


"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
from scipy import stats
from IceLib import read_ice,ice_trendCalc,grid_data
from gridCERESLib import readgridCERES,calc_CERES_trend
import cartopy
import cartopy.crs as ccrs
import glob
import subprocess

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

def plot_fourscatter(ice_data,CERES_lw_dict,CERES_sw_dict,CERES_alb_dict,CERES_net_dict,inseason,dtype,area=True,save=False):
    if(dtype=='clr'):
        dtype_adder = "clear"
    elif(dtype=='cld'):
        dtype_adder = "cloudy"
    else:
        dtype_adder = "all"

    if(area==True):
        ice_trend_key = 'grid_ice_area_trend'
        trend_checker = -3000.
        pcnt_reject = 10
        missing_trend = -9999.
    else:
        ice_trend_key = 'grid_ice'
        trend_checker = -1000.
        pcnt_reject = 0.1
        missing_trend = -999.

    starting_area = np.nansum((ice_data['grid_ice_conc'][0,:,:]/100.)*ice_data['grid_total_area'])

    # Quick comparison
    markersize=6
    fig, axs = plt.subplots(2,2)
    fig.set_size_inches(15,13)
    plt.title("CERES and NSIDC "+inseason.title()+" Trends")
    # Plot LWF data
    plotx = ice_data[ice_trend_key][(ice_data[ice_trend_key]!=missing_trend) & (np.abs(ice_data[ice_trend_key])>=pcnt_reject)]
    ploty = CERES_lw_dict['trends'][(ice_data[ice_trend_key]!=missing_trend) & (np.abs(ice_data[ice_trend_key])>=pcnt_reject)]
    plotx = plotx[ploty>-1000]
    ploty = ploty[ploty>-1000]
    plotx = plotx[ploty<1000]
    ploty = ploty[ploty<1000]
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
    if(dtype=='cld'):
        x_pos = ((axs[0,0].get_xlim()[1]-axs[0,0].get_xlim()[0])/8)+axs[0,0].get_xlim()[0]
        y_pos = axs[0,0].get_ylim()[1]-((axs[0,0].get_ylim()[1]-axs[0,0].get_ylim()[0])/8)
    axs[0,0].text(x_pos,y_pos,"Corr = "+str(np.round(corr_lwf,3)))
    # Plot SWF data
    plotx = ice_data[ice_trend_key][(ice_data[ice_trend_key]!=missing_trend) & (np.abs(ice_data[ice_trend_key])>=pcnt_reject)]
    ploty = CERES_sw_dict['trends'][(ice_data[ice_trend_key]!=missing_trend) & (np.abs(ice_data[ice_trend_key])>=pcnt_reject)]
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
    if(dtype=='cld'):
        x_pos = ((axs[0,1].get_xlim()[1]-axs[0,1].get_xlim()[0])/8)+axs[0,1].get_xlim()[0]
        y_pos = ((axs[0,1].get_ylim()[1]-axs[0,1].get_ylim()[0])/8)+axs[0,1].get_ylim()[0]
    axs[0,1].text(x_pos,y_pos,"Corr = "+str(np.round(corr_swf,3)))
    ##print("  y_max = ",axs[0,1].get_ylim()[1],"  y_min = ",axs[0,1].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    # Plot Albedo Data
    plotx = ice_data[ice_trend_key][(ice_data[ice_trend_key]!=missing_trend) & (np.abs(ice_data[ice_trend_key])>=pcnt_reject)]
    ploty = CERES_alb_dict['trends'][(ice_data[ice_trend_key]!=missing_trend) & (np.abs(ice_data[ice_trend_key])>=pcnt_reject)]
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
    if(dtype=='cld'):
        x_pos = ((axs[1,0].get_xlim()[1]-axs[1,0].get_xlim()[0])/8)+axs[1,0].get_xlim()[0]
        y_pos = ((axs[1,0].get_ylim()[1]-axs[1,0].get_ylim()[0])/8)+axs[1,0].get_ylim()[0]
    axs[1,0].text(x_pos,y_pos,"Corr = "+str(np.round(corr_alb,3)))
    ##print("  y_max = ",axs[1,0].get_ylim()[1],"  y_min = ",axs[1,0].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    # Plot Net Flux Data
    plotx = ice_data[ice_trend_key][(ice_data[ice_trend_key]!=missing_trend) & (np.abs(ice_data[ice_trend_key])>=pcnt_reject)]
    ploty = CERES_net_dict['trends'][(ice_data[ice_trend_key]!=missing_trend) & (np.abs(ice_data[ice_trend_key])>=pcnt_reject)]
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
    if(dtype=='cld'):
        x_pos = ((axs[1,1].get_xlim()[1]-axs[1,1].get_xlim()[0])/8)+axs[1,1].get_xlim()[0]
        y_pos = axs[1,1].get_ylim()[1]-((axs[1,1].get_ylim()[1]-axs[1,1].get_ylim()[0])/8)
    axs[1,1].text(x_pos,y_pos,"Corr = "+str(np.round(corr_net,3)))
    ##print("  y_max = ",axs[1,1].get_ylim()[1],"  y_min = ",axs[1,1].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    outname = "ceres_nsidc_trends_four_panel_"+inseason+"_"+dtype+"sky.png"
    if(save==True):
        plt.savefig(outname,dpi=300)
        print("Saved image "+outname)
    else:
        plt.show()


def plot_cld_clr_scatter(CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_alb_clr_dict,CERES_net_clr_dict,\
                         CERES_lw_cld_dict,CERES_sw_cld_dict,CERES_alb_cld_dict,CERES_net_cld_dict,\
                         ice_data,inseason):
    #if(dtype=='clr'):
    #    dtype_adder = "clear"
    #elif(dtype=='cld'):
    #    dtype_adder = "cloudy"
    #else:
    #    dtype_adder = "all"
    # Quick comparison
    markersize=6
    fig, axs = plt.subplots(2,2)
    fig.set_size_inches(15,13)
    plt.title("CERES Cloudy-Sky and Clear-sky "+inseason.title()+" Trends")
    # Plot LWF data
    #ice   = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    plotx = CERES_lw_cld_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    ploty = CERES_lw_clr_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    plotx = plotx[ploty>-1000]
    ploty = ploty[ploty>-1000]
    plotx = plotx[ploty<1000]
    ploty = ploty[ploty<1000]
    corr_lwf = correlation(plotx,ploty)
    print("LWF correlation: ",corr_lwf)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    axs[0,0].scatter(plotx,ploty,s=markersize,color='black')
    axs[0,0].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
    axs[0,0].set_title('CERES LWF'+ice_data['season_adder'].title()+' Trends')
    axs[0,0].set_xlabel('Cloudy-sky Trends')
    axs[0,0].set_ylabel('Clear-sky Trends')
    x_pos = ((axs[0,0].get_xlim()[1]-axs[0,0].get_xlim()[0])/8)+axs[0,0].get_xlim()[0]
    y_pos = ((axs[0,0].get_ylim()[1]-axs[0,0].get_ylim()[0])/8)+axs[0,0].get_ylim()[0]
    #if(dtype=='cld'):
    #    x_pos = ((axs[0,0].get_xlim()[1]-axs[0,0].get_xlim()[0])/8)+axs[0,0].get_xlim()[0]
    #    y_pos = axs[0,0].get_ylim()[1]-((axs[0,0].get_ylim()[1]-axs[0,0].get_ylim()[0])/8)
    axs[0,0].text(x_pos,y_pos,"Corr = "+str(np.round(corr_lwf,3)))
    # Plot SWF data
    #ice   = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    plotx = CERES_sw_cld_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    ploty = CERES_sw_clr_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    corr_swf = correlation(plotx,ploty)
    print("SWF correlation: ",corr_swf)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    axs[0,1].scatter(plotx,ploty,s=markersize,color='black')
    axs[0,1].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
    axs[0,1].set_title('CERES SWF'+ice_data['season_adder'].title()+' Trends')
    axs[0,1].set_xlabel('Cloudy-sky Trends')
    axs[0,1].set_ylabel('Clear-sky Trends')
    x_pos = ((axs[0,1].get_xlim()[1]-axs[0,1].get_xlim()[0])/8)+axs[0,1].get_xlim()[0]
    y_pos = axs[0,1].get_ylim()[1]-((axs[0,1].get_ylim()[1]-axs[0,1].get_ylim()[0])/8)
    #if(dtype=='cld'):
    #    x_pos = ((axs[0,1].get_xlim()[1]-axs[0,1].get_xlim()[0])/8)+axs[0,1].get_xlim()[0]
    #    y_pos = ((axs[0,1].get_ylim()[1]-axs[0,1].get_ylim()[0])/8)+axs[0,1].get_ylim()[0]
    axs[0,1].text(x_pos,y_pos,"Corr = "+str(np.round(corr_swf,3)))
    ##print("  y_max = ",axs[0,1].get_ylim()[1],"  y_min = ",axs[0,1].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    # Plot Albedo Data
    #plotx = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    plotx = CERES_alb_cld_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    ploty = CERES_alb_clr_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    corr_alb = correlation(plotx,ploty)
    print("Albedo correlation: ",corr_alb)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    axs[1,0].scatter(plotx,ploty,s=markersize,color='black')
    axs[1,0].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
    axs[1,0].set_title('CERES Albedo'+ice_data['season_adder'].title()+' Trends')
    axs[1,0].set_xlabel('Cloudy-sky Trends')
    axs[1,0].set_ylabel('Clear-sky Trends')
    x_pos = ((axs[1,0].get_xlim()[1]-axs[1,0].get_xlim()[0])/8)+axs[1,0].get_xlim()[0]
    y_pos = axs[1,0].get_ylim()[1]-((axs[1,0].get_ylim()[1]-axs[1,0].get_ylim()[0])/8)
    #if(dtype=='cld'):
    #    x_pos = ((axs[1,0].get_xlim()[1]-axs[1,0].get_xlim()[0])/8)+axs[1,0].get_xlim()[0]
    #    y_pos = ((axs[1,0].get_ylim()[1]-axs[1,0].get_ylim()[0])/8)+axs[1,0].get_ylim()[0]
    axs[1,0].text(x_pos,y_pos,"Corr = "+str(np.round(corr_alb,3)))
    ##print("  y_max = ",axs[1,0].get_ylim()[1],"  y_min = ",axs[1,0].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    # Plot Net Flux Data
    #plotx = ice_data['grid_ice'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    plotx = CERES_net_cld_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    ploty = CERES_net_clr_dict['trends'][(ice_data['grid_ice']!=-999.) & (ice_data['grid_ice']!=0.)]
    corr_net = correlation(plotx,ploty)
    print("Net flux correlation: ",corr_net)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    axs[1,1].scatter(plotx,ploty,s=markersize,color='black')
    axs[1,1].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
    axs[1,1].set_title('CERES Net Flux'+ice_data['season_adder'].title()+' Trends')
    axs[1,1].set_xlabel('Cloudy-sky Trends')
    axs[1,1].set_ylabel('Clear-sky Trends')
    x_pos = ((axs[1,1].get_xlim()[1]-axs[1,1].get_xlim()[0])/8)+axs[1,1].get_xlim()[0]
    y_pos = ((axs[1,1].get_ylim()[1]-axs[1,1].get_ylim()[0])/8)+axs[1,1].get_ylim()[0]
    #if(dtype=='cld'):
    #    x_pos = ((axs[1,1].get_xlim()[1]-axs[1,1].get_xlim()[0])/8)+axs[1,1].get_xlim()[0]
    #    y_pos = axs[1,1].get_ylim()[1]-((axs[1,1].get_ylim()[1]-axs[1,1].get_ylim()[0])/8)
    axs[1,1].text(x_pos,y_pos,"Corr = "+str(np.round(corr_net,3)))
    ##print("  y_max = ",axs[1,1].get_ylim()[1],"  y_min = ",axs[1,1].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    outname = "ceres_clr_cld_trends_four_panel_"+inseason+"_sky.png"
    plt.savefig(outname,dpi=300)
    print("Saved image "+outname)
    plt.show()

def figure_1(ice_data,CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_alb_clr_dict,\
             CERES_net_clr_dict,CERES_lw_all_dict,CERES_sw_all_dict,\
             CERES_alb_all_dict,CERES_net_all_dict,adjusted=True):
    colormap = plt.cm.bwr
    if(adjusted==True):
        mapcrs = ccrs.NorthPolarStereo(central_longitude=45.)
    else:
        mapcrs = ccrs.NorthPolarStereo()
    datacrs = ccrs.PlateCarree()
    fig = plt.figure(1, figsize=(20,8))
    gs = gridspec.GridSpec(nrows=1, ncols=3, hspace = 0.0, wspace = 0.06)
    sw_clr_min_val = -40.
    sw_clr_max_val = 40.
    lw_clr_min_val=-15.
    lw_clr_max_val = 15.
    title_adder = 'TOA Clear-Sky SWF'
    # - - - - - - - - - - - - - 
    # Plot the SW clr data
    # - - - - - - - - - - - - - 
    ax1 = plt.subplot(gs[0],projection=mapcrs)
    ax1.set_extent([-180,180,45,90],datacrs)
    ax1.gridlines()
    #ax = plt.axes(projection=ccrs.Miller())
    file_season = '_'+CERES_sw_clr_dict['season']
    title_season = ' '+CERES_sw_clr_dict['season']
    if(CERES_sw_clr_dict['season']=='all'):
        file_season=''
        title_season=''
        
    mesh = ax1.pcolormesh(CERES_sw_clr_dict['lon'],CERES_sw_clr_dict['lat'],CERES_sw_clr_dict['trends'],\
            transform=datacrs,vmin=sw_clr_min_val,vmax=sw_clr_max_val,cmap=colormap)
    ax1.set_title('Terra CERES '+title_adder+title_season.title()+' Trend\n'+\
                 CERES_sw_clr_dict['dates'][0]+' - '+CERES_sw_clr_dict['dates'][-1])
    ender = '.png'
    if(adjusted==True):
        ax1.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        ##axins = inset_axes(ax1,width="5%",height="100%",loc='lower left',\
        ##            bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax1.transAxes,\
        ##            borderpad=0)
        ##cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=CERES_sw_clr_dict['parm_unit'])
        cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,aspect=50,label=CERES_sw_clr_dict['parm_unit'])
        ax1.set_xlim(-5170748.535086173,5167222.438879491)
        ax1.set_ylim(-3913488.8763307533,3943353.899053069)
        ender = '_adjusted'+ender
    else:
        cbar = plt.colorbar(mesh,cmap=colormap,label=CERES_sw_clr_dict['parm_unit'])
    #cbar.set_label(parm_unit)
    ax1.coastlines()
    
    # - - - - - - - - - - - - - 
    # Plot the LW clr data
    # - - - - - - - - - - - - - 
    ax2 = plt.subplot(gs[1],projection=mapcrs)
    ax2.set_extent([-180,180,45,90],datacrs)
    ax2.gridlines()
    #ax = plt.axes(projection=ccrs.Miller())
    file_season = '_'+CERES_lw_clr_dict['season']
    title_season = ' '+CERES_lw_clr_dict['season']
    if(CERES_lw_clr_dict['season']=='all'):
        file_season=''
        title_season=''
        
    mesh = ax2.pcolormesh(CERES_lw_clr_dict['lon'],CERES_lw_clr_dict['lat'],CERES_lw_clr_dict['trends'],\
            transform=datacrs,vmin=lw_clr_min_val,vmax=lw_clr_max_val,cmap=colormap)
    ax2.set_title('Terra CERES '+title_adder+title_season.title()+' Trend\n'+\
                 CERES_lw_clr_dict['dates'][0]+' - '+CERES_lw_clr_dict['dates'][-1])
    ender = '.png'
    if(adjusted==True):
        ax2.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        ##axins = inset_axes(ax2,width="5%",height="100%",loc='lower left',\
        ##            bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax2.transAxes,\
        ##            borderpad=0)
        ##cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=CERES_lw_clr_dict['parm_unit'])
        cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,aspect=50,label=CERES_lw_clr_dict['parm_unit'])
        ax2.set_xlim(-5170748.535086173,5167222.438879491)
        ax2.set_ylim(-3913488.8763307533,3943353.899053069)
        ender = '_adjusted'+ender
    else:
        cbar = plt.colorbar(mesh,cmap=colormap,label=CERES_lw_clr_dict['parm_unit'])
    #cbar.set_label(parm_unit)
    ax2.coastlines()
    
    
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    
    local_grid_ice = np.copy(ice_data['grid_ice'])
    local_grid_ice_bad = np.copy(ice_data['grid_ice'])
    
    local_grid_ice[local_grid_ice==-999.] = np.nan
    local_grid_ice_bad[local_grid_ice!=-999.] = np.nan
    plot_good_data = ma.masked_invalid(local_grid_ice)
    plot_land_data = ma.masked_invalid(local_grid_ice_bad)
    
    colormap = plt.cm.bwr
    #coolwarm = plt.cm.coolwarm
    #ax = plt.axes(projection=ccrs.NorthPolarStereo())
    file_adder=''
    ax3 = plt.subplot(gs[2],projection=mapcrs)
    ax3.set_extent([-180,180,45,90],datacrs)
    ax3.gridlines()
    mesh = ax3.pcolormesh(lon_ranges,lat_ranges,plot_good_data,\
            transform=datacrs,vmin=-50,vmax=50,cmap=colormap)
    if(adjusted==True):
        ax3.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        ##axins = inset_axes(ax3,width="5%",height="100%",loc='lower left',\
        ##            bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax3.transAxes,\
        ##            borderpad=0)
        ##cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label='Ice Concentration Trend [%]')
        cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,aspect=50,label='Percent Ice Concentration Trend')
        ax3.set_xlim(-5170748.535086173,5167222.438879491)
        ax3.set_ylim(-3913488.8763307533,3943353.899053069)
        #ender = '_adjusted'+ender
    else:
        ax3.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
        cbar = plt.colorbar(mesh,cmap=colormap,label='Ice Concentration Trend [%]')
    ax3.coastlines()
    plt.show()
