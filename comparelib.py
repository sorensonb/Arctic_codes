#!/usr/bin/env python
"""


"""

import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats
from scipy.stats import gaussian_kde
import cartopy
import cartopy.crs as ccrs
from netCDF4 import Dataset
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
    plt.title("CERES and NSIDC " + inseason.title() + " Trends")
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
    axs[0,0].set_title(dtype_adder.title() + '-Sky LWF' + ice_data['season_adder'].title() + ' Trends')
    axs[0,0].set_xlabel('Ice Trends')
    axs[0,0].set_ylabel('LWF Trends')
    x_pos = ((axs[0,0].get_xlim()[1]-axs[0,0].get_xlim()[0])/8) + axs[0,0].get_xlim()[0]
    y_pos = ((axs[0,0].get_ylim()[1]-axs[0,0].get_ylim()[0])/8) + axs[0,0].get_ylim()[0]
    if(dtype=='cld'):
        x_pos = ((axs[0,0].get_xlim()[1]-axs[0,0].get_xlim()[0])/8) + axs[0,0].get_xlim()[0]
        y_pos = axs[0,0].get_ylim()[1]-((axs[0,0].get_ylim()[1]-axs[0,0].get_ylim()[0])/8)
    axs[0,0].text(x_pos,y_pos,"Corr = " + str(np.round(corr_lwf,3)))
    # Plot SWF data
    plotx = ice_data[ice_trend_key][(ice_data[ice_trend_key]!=missing_trend) & (np.abs(ice_data[ice_trend_key])>=pcnt_reject)]
    ploty = CERES_sw_dict['trends'][(ice_data[ice_trend_key]!=missing_trend) & (np.abs(ice_data[ice_trend_key])>=pcnt_reject)]
    corr_swf = correlation(plotx,ploty)
    print("SWF correlation: ",corr_swf)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    axs[0,1].scatter(plotx,ploty,s=markersize,color='black')
    axs[0,1].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
    axs[0,1].set_title(dtype_adder.title() + '-Sky SWF' + ice_data['season_adder'].title() + ' Trends')
    axs[0,1].set_xlabel('Ice Trends')
    axs[0,1].set_ylabel('SWF Trends')
    x_pos = ((axs[0,1].get_xlim()[1]-axs[0,1].get_xlim()[0])/8) + axs[0,1].get_xlim()[0]
    y_pos = axs[0,1].get_ylim()[1]-((axs[0,1].get_ylim()[1]-axs[0,1].get_ylim()[0])/8)
    if(dtype=='cld'):
        x_pos = ((axs[0,1].get_xlim()[1]-axs[0,1].get_xlim()[0])/8) + axs[0,1].get_xlim()[0]
        y_pos = ((axs[0,1].get_ylim()[1]-axs[0,1].get_ylim()[0])/8) + axs[0,1].get_ylim()[0]
    axs[0,1].text(x_pos,y_pos,"Corr = " + str(np.round(corr_swf,3)))
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
    axs[1,0].set_title(dtype_adder.title() + '-Sky Albedo' + ice_data['season_adder'].title() + ' Trends')
    axs[1,0].set_xlabel('Ice Trends')
    axs[1,0].set_ylabel('Albedo Trends')
    x_pos = ((axs[1,0].get_xlim()[1]-axs[1,0].get_xlim()[0])/8) + axs[1,0].get_xlim()[0]
    y_pos = axs[1,0].get_ylim()[1]-((axs[1,0].get_ylim()[1]-axs[1,0].get_ylim()[0])/8)
    if(dtype=='cld'):
        x_pos = ((axs[1,0].get_xlim()[1]-axs[1,0].get_xlim()[0])/8) + axs[1,0].get_xlim()[0]
        y_pos = ((axs[1,0].get_ylim()[1]-axs[1,0].get_ylim()[0])/8) + axs[1,0].get_ylim()[0]
    axs[1,0].text(x_pos,y_pos,"Corr = " + str(np.round(corr_alb,3)))
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
    axs[1,1].set_title(dtype_adder.title() + '-Sky Net Flux' + ice_data['season_adder'].title() + ' Trends')
    axs[1,1].set_xlabel('Ice Trends')
    axs[1,1].set_ylabel('Net Flux Trends')
    x_pos = ((axs[1,1].get_xlim()[1]-axs[1,1].get_xlim()[0])/8) + axs[1,1].get_xlim()[0]
    y_pos = ((axs[1,1].get_ylim()[1]-axs[1,1].get_ylim()[0])/8) + axs[1,1].get_ylim()[0]
    if(dtype=='cld'):
        x_pos = ((axs[1,1].get_xlim()[1]-axs[1,1].get_xlim()[0])/8) + axs[1,1].get_xlim()[0]
        y_pos = axs[1,1].get_ylim()[1]-((axs[1,1].get_ylim()[1]-axs[1,1].get_ylim()[0])/8)
    axs[1,1].text(x_pos,y_pos,"Corr = " + str(np.round(corr_net,3)))
    ##print("  y_max = ",axs[1,1].get_ylim()[1],"  y_min = ",axs[1,1].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    outname = "ceres_nsidc_trends_four_panel_" + inseason + "_" + dtype + "sky.png"
    if(save==True):
        plt.savefig(outname,dpi=300)
        print("Saved image " + outname)
    else:
        plt.show()


def plot_cld_clr_scatter(CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_net_clr_dict,\
                         CERES_lw_cld_dict,CERES_sw_cld_dict,CERES_net_cld_dict,\
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
    plt.title("CERES Cloudy-Sky and Clear-sky " + inseason.title() + " Trends")
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
    axs[0,0].set_title('CERES LWF' + ice_data['season_adder'].title() + ' Trends')
    axs[0,0].set_xlabel('Cloudy-sky Trends')
    axs[0,0].set_ylabel('Clear-sky Trends')
    x_pos = ((axs[0,0].get_xlim()[1]-axs[0,0].get_xlim()[0])/8) + axs[0,0].get_xlim()[0]
    y_pos = ((axs[0,0].get_ylim()[1]-axs[0,0].get_ylim()[0])/8) + axs[0,0].get_ylim()[0]
    #if(dtype=='cld'):
    #    x_pos = ((axs[0,0].get_xlim()[1]-axs[0,0].get_xlim()[0])/8) + axs[0,0].get_xlim()[0]
    #    y_pos = axs[0,0].get_ylim()[1]-((axs[0,0].get_ylim()[1]-axs[0,0].get_ylim()[0])/8)
    axs[0,0].text(x_pos,y_pos,"Corr = " + str(np.round(corr_lwf,3)))
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
    axs[0,1].set_title('CERES SWF' + ice_data['season_adder'].title() + ' Trends')
    axs[0,1].set_xlabel('Cloudy-sky Trends')
    axs[0,1].set_ylabel('Clear-sky Trends')
    x_pos = ((axs[0,1].get_xlim()[1]-axs[0,1].get_xlim()[0])/8) + axs[0,1].get_xlim()[0]
    y_pos = axs[0,1].get_ylim()[1]-((axs[0,1].get_ylim()[1]-axs[0,1].get_ylim()[0])/8)
    #if(dtype=='cld'):
    #    x_pos = ((axs[0,1].get_xlim()[1]-axs[0,1].get_xlim()[0])/8) + axs[0,1].get_xlim()[0]
    #    y_pos = ((axs[0,1].get_ylim()[1]-axs[0,1].get_ylim()[0])/8) + axs[0,1].get_ylim()[0]
    axs[0,1].text(x_pos,y_pos,"Corr = " + str(np.round(corr_swf,3)))
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
    axs[1,0].set_title('CERES Albedo' + ice_data['season_adder'].title() + ' Trends')
    axs[1,0].set_xlabel('Cloudy-sky Trends')
    axs[1,0].set_ylabel('Clear-sky Trends')
    x_pos = ((axs[1,0].get_xlim()[1]-axs[1,0].get_xlim()[0])/8) + axs[1,0].get_xlim()[0]
    y_pos = axs[1,0].get_ylim()[1]-((axs[1,0].get_ylim()[1]-axs[1,0].get_ylim()[0])/8)
    #if(dtype=='cld'):
    #    x_pos = ((axs[1,0].get_xlim()[1]-axs[1,0].get_xlim()[0])/8) + axs[1,0].get_xlim()[0]
    #    y_pos = ((axs[1,0].get_ylim()[1]-axs[1,0].get_ylim()[0])/8) + axs[1,0].get_ylim()[0]
    axs[1,0].text(x_pos,y_pos,"Corr = " + str(np.round(corr_alb,3)))
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
    axs[1,1].set_title('CERES Net Flux' + ice_data['season_adder'].title() + ' Trends')
    axs[1,1].set_xlabel('Cloudy-sky Trends')
    axs[1,1].set_ylabel('Clear-sky Trends')
    x_pos = ((axs[1,1].get_xlim()[1]-axs[1,1].get_xlim()[0])/8) + axs[1,1].get_xlim()[0]
    y_pos = ((axs[1,1].get_ylim()[1]-axs[1,1].get_ylim()[0])/8) + axs[1,1].get_ylim()[0]
    #if(dtype=='cld'):
    #    x_pos = ((axs[1,1].get_xlim()[1]-axs[1,1].get_xlim()[0])/8) + axs[1,1].get_xlim()[0]
    #    y_pos = axs[1,1].get_ylim()[1]-((axs[1,1].get_ylim()[1]-axs[1,1].get_ylim()[0])/8)
    axs[1,1].text(x_pos,y_pos,"Corr = " + str(np.round(corr_net,3)))
    ##print("  y_max = ",axs[1,1].get_ylim()[1],"  y_min = ",axs[1,1].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    outname = "ceres_clr_cld_trends_four_panel_" + inseason + "_sky.png"
    plt.savefig(outname,dpi=300)
    print("Saved image " + outname)
    plt.show()

def figure_1(ice_data,CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_net_clr_dict,\
             CERES_lw_all_dict,CERES_sw_all_dict,CERES_net_all_dict,adjusted=True):
    colormap = plt.cm.bwr
    if(adjusted==True):
        mapcrs = ccrs.NorthPolarStereo(central_longitude=45.)
    else:
        mapcrs = ccrs.NorthPolarStereo()
    datacrs = ccrs.PlateCarree()
    fig = plt.figure(1, figsize=(12,16))
    if(ice_data['season_adder'].strip()=='sunlight'):
        plt.suptitle('Sunlight Months (April - September  2001 - 2018)',y=0.92,fontsize=18,fontweight=4)
    else: 
        plt.suptitle('Summer Only (June 2001 - August 2018)',y=0.92,fontsize=18,fontweight=4)
    gs = gridspec.GridSpec(nrows=3, ncols=2, hspace = 0.03, wspace = 0.06)
    sw_clr_min_val = -60.
    sw_clr_max_val = 60.
    lw_clr_min_val=-10.
    lw_clr_max_val = 10.
    
    file_season = '_' + CERES_sw_clr_dict['season']
    title_season = ' ' + CERES_sw_clr_dict['season']
    if(CERES_sw_clr_dict['season']=='all'):
        file_season=''
        title_season=''
        
    # - - - - - - - - - - - - - 
    # Plot the average summer ice conc
    # - - - - - - - - - - - - - 
    # Find the average conc
    avg_ice_conc = np.zeros((len(ice_data['grid_lat']),len(ice_data['grid_lon'])))
    for xi in range(len(ice_data['grid_lat'])):
        for yj in range(len(ice_data['grid_lon'])):
            avg_ice_conc[xi,yj] = np.average(ice_data['grid_ice_conc'][:,xi,yj][ice_data['grid_ice_conc'][:,xi,yj]!=-99.])
            if(avg_ice_conc[xi,yj]<1.0):
                avg_ice_conc[xi,yj] = np.nan
       
    plot_good_data = ma.masked_invalid(avg_ice_conc)

    colormap = plt.cm.ocean
    #coolwarm = plt.cm.coolwarm
    #ax = plt.axes(projection=ccrs.NorthPolarStereo())
    file_adder=''
    ax0 = plt.subplot(gs[0,0],projection=mapcrs)
    ax0.set_extent([-180,180,45,90],ccrs.PlateCarree())
    ax0.gridlines()
    mesh = ax0.pcolormesh(ice_data['grid_lon'],ice_data['grid_lat'],plot_good_data,\
            transform=ccrs.PlateCarree(),vmin=0,vmax=100,cmap=colormap)
    ax0.set_title('NSIDC Average Percent Ice Concentration')
    #ax0.set_title('NSIDC Percent Ice Concentration' + title_season.title() + ' Trend\n' + \
    #             CERES_sw_clr_dict['dates'][0] + ' - ' + CERES_sw_clr_dict['dates'][-1])
    if(adjusted==True):
        ax0.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        cbar = plt.colorbar(mesh,ticks=[1,20,40,60,80,100],orientation='horizontal',pad=0,aspect=50,label='Percent Ice Concentration')
        cbar.ax.set_xticklabels(['1','20','40','60','80','100'])
        ax0.set_xlim(-4170748.535086173,4167222.438879491)
        ax0.set_ylim(-2913488.8763307533,2943353.899053069)
    else:
        plt.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
        cbar = plt.colorbar(mesh,cmap=colormap,label=plabel)
    ax0.coastlines()
     

    # - - - - - - - - - - - - - 
    # Plot the Ice trend data
    # - - - - - - - - - - - - - 
    
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    
    local_grid_ice = np.copy(ice_data['grid_ice'])
    local_grid_ice_bad = np.copy(ice_data['grid_ice'])
    
    local_grid_ice[local_grid_ice==-999.] = np.nan
    # Ignore the areas without ice by making any values where the ice conc
    # is less than 5%.
    local_grid_ice[np.isnan(avg_ice_conc)] = np.nan
    local_grid_ice_bad[local_grid_ice!=-999.] = np.nan
    plot_good_data = ma.masked_invalid(local_grid_ice)
    plot_land_data = ma.masked_invalid(local_grid_ice_bad)
    
    colormap = plt.cm.bwr
    #coolwarm = plt.cm.coolwarm
    #ax = plt.axes(projection=ccrs.NorthPolarStereo())
    file_adder=''
    ax1 = plt.subplot(gs[0,1],projection=mapcrs)
    ax1.set_extent([-180,180,45,90],datacrs)
    ax1.gridlines()
    mesh = ax1.pcolormesh(lon_ranges,lat_ranges,plot_good_data,\
            transform=datacrs,vmin=-50,vmax=50,cmap=colormap)
    ax1.set_title('NSIDC Percent Ice Concentration Trend')
    #ax1.set_title('NSIDC Percent Ice Concentration' + title_season.title() + ' Trend\n' + \
    #             CERES_sw_clr_dict['dates'][0] + ' - ' + CERES_sw_clr_dict['dates'][-1])
    if(adjusted==True):
        ax1.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        ##axins = inset_axes(ax3,width="5%",height="100%",loc='lower left',\
        ##            bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax3.transAxes,\
        ##            borderpad=0)
        ##cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label='Ice Concentration Trend [%]')
        cbar = plt.colorbar(mesh,extend='both',orientation='horizontal',pad=0,aspect=50,label='Percent Ice Concentration / Study Period')
        ax1.set_xlim(-4170748.535086173,4167222.438879491)
        ax1.set_ylim(-2913488.8763307533,2943353.899053069)
        #ender = '_adjusted' + ender
    else:
        ax1.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
        cbar = plt.colorbar(mesh,extend='both',cmap=colormap,label='Ice Concentration Trend [%]')
    ax1.coastlines()

    # - - - - - - - - - - - - - 
    # Plot the SW clr data
    # - - - - - - - - - - - - - 
    ax2 = plt.subplot(gs[1,0],projection=mapcrs)
    ax2.set_extent([-180,180,45,90],datacrs)
    ax2.gridlines()
    #ax = plt.axes(projection=ccrs.Miller())
    # Ignore the areas without ice by making any values where the ice conc
    # is less than 5%.
    local_SW = np.copy(CERES_sw_clr_dict['trends'])
    local_SW[np.isnan(avg_ice_conc)] = np.nan
    plot_local_SW = ma.masked_invalid(local_SW)
    mesh = ax2.pcolormesh(CERES_sw_clr_dict['lon'],CERES_sw_clr_dict['lat'],plot_local_SW,\
            transform=datacrs,vmin=sw_clr_min_val,vmax=sw_clr_max_val,cmap=colormap)
    ax2.set_title('Terra CERES TOA Clear-sky SWF Trend')
    #ax2.set_title('Terra CERES TOA Clear-sky SWF' + title_season.title() + ' Trend\n' + \
    #             CERES_sw_clr_dict['dates'][0] + ' - ' + CERES_sw_clr_dict['dates'][-1])
    ender = '.png'
    if(adjusted==True):
        ax2.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        ##axins = inset_axes(ax1,width="5%",height="100%",loc='lower left',\
        ##            bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax1.transAxes,\
        ##            borderpad=0)
        ##cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=CERES_sw_clr_dict['parm_unit'])
        cbar = plt.colorbar(mesh,extend='both',orientation='horizontal',pad=0,aspect=50,label=CERES_sw_clr_dict['parm_unit'] + ' / Study Period')
        #ax2.set_xlim(-5170748.535086173,5167222.438879491)
        #ax2.set_ylim(-3913488.8763307533,3943353.899053069)
        ax2.set_xlim(-4170748.535086173,4167222.438879491)
        ax2.set_ylim(-2913488.8763307533,2943353.899053069)
        ender = '_adjusted' + ender
    else:
        cbar = plt.colorbar(mesh,cmap=colormap,label=CERES_sw_clr_dict['parm_unit'])
    #cbar.set_label(parm_unit)
    ax2.coastlines()
    
    # - - - - - - - - - - - - - 
    # Plot the LW clr data
    # - - - - - - - - - - - - - 
    ax3 = plt.subplot(gs[1,1],projection=mapcrs)
    ax3.set_extent([-180,180,45,90],datacrs)
    ax3.gridlines()
    #ax = plt.axes(projection=ccrs.Miller())
    file_season = '_' + CERES_lw_clr_dict['season']
    title_season = ' ' + CERES_lw_clr_dict['season']
    if(CERES_lw_clr_dict['season']=='all'):
        file_season=''
        title_season=''
    # Ignore the areas without ice by making any values where the ice conc
    # is less than 5%.
    local_LW = np.copy(CERES_lw_clr_dict['trends'])
    local_LW[np.isnan(avg_ice_conc)] = np.nan
    plot_local_LW = ma.masked_invalid(local_LW)
        
    mesh = ax3.pcolormesh(CERES_lw_clr_dict['lon'],CERES_lw_clr_dict['lat'],plot_local_LW,\
            transform=datacrs,vmin=lw_clr_min_val,vmax=lw_clr_max_val,cmap=colormap)
    ax3.set_title('Terra CERES TOA Clear-sky LWF Trend')
    #ax3.set_title('Terra CERES TOA Clear-sky LWF' + title_season.title() + ' Trend\n' + \
    #             CERES_lw_clr_dict['dates'][0] + ' - ' + CERES_lw_clr_dict['dates'][-1])
    ender = '.png'
    if(adjusted==True):
        ax3.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        ##axins = inset_axes(ax2,width="5%",height="100%",loc='lower left',\
        ##            bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax2.transAxes,\
        ##            borderpad=0)
        ##cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=CERES_lw_clr_dict['parm_unit'])
        cbar = plt.colorbar(mesh,extend='both',orientation='horizontal',pad=0,aspect=50,label=CERES_lw_clr_dict['parm_unit'] + ' / Study Period')
        #ax3.set_xlim(-5170748.535086173,5167222.438879491)
        #ax3.set_ylim(-3913488.8763307533,3943353.899053069)
        ax3.set_xlim(-4170748.535086173,4167222.438879491)
        ax3.set_ylim(-2913488.8763307533,2943353.899053069)
        ender = '_adjusted' + ender
    else:
        cbar = plt.colorbar(mesh,cmap=colormap,label=CERES_lw_clr_dict['parm_unit'])
    #cbar.set_label(parm_unit)
    ax3.coastlines()

    # - - - - - - - - - - - - - 
    # Plot the SW all data
    # - - - - - - - - - - - - - 
    ax4 = plt.subplot(gs[2,0],projection=mapcrs)
    ax4.set_extent([-180,180,45,90],datacrs)
    ax4.gridlines()
    #ax = plt.axes(projection=ccrs.Miller())
    # Ignore the areas without ice by making any values where the ice conc
    # is less than 5%.
    local_SW = np.copy(CERES_sw_all_dict['trends'])
    local_SW[np.isnan(avg_ice_conc)] = np.nan
    plot_local_SW = ma.masked_invalid(local_SW)
    mesh = ax4.pcolormesh(CERES_sw_all_dict['lon'],CERES_sw_all_dict['lat'],plot_local_SW,\
            transform=datacrs,vmin=sw_clr_min_val,vmax=sw_clr_max_val,cmap=colormap)
    ax4.set_title('Terra CERES TOA All-sky SWF Trend')
    #ax2.set_title('Terra CERES TOA Clear-sky SWF' + title_season.title() + ' Trend\n' + \
    #             CERES_sw_clr_dict['dates'][0] + ' - ' + CERES_sw_clr_dict['dates'][-1])
    ender = '.png'
    if(adjusted==True):
        ax4.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        ##axins = inset_axes(ax1,width="5%",height="100%",loc='lower left',\
        ##            bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax1.transAxes,\
        ##            borderpad=0)
        ##cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=CERES_sw_clr_dict['parm_unit'])
        cbar = plt.colorbar(mesh,extend='both',orientation='horizontal',pad=0,aspect=50,label=CERES_sw_all_dict['parm_unit'] + ' / Study Period')
        #ax2.set_xlim(-5170748.535086173,5167222.438879491)
        #ax2.set_ylim(-3913488.8763307533,3943353.899053069)
        ax4.set_xlim(-4170748.535086173,4167222.438879491)
        ax4.set_ylim(-2913488.8763307533,2943353.899053069)
        ender = '_adjusted' + ender
    else:
        cbar = plt.colorbar(mesh,cmap=colormap,label=CERES_sw_all_dict['parm_unit'])
    #cbar.set_label(parm_unit)
    ax4.coastlines()
    
    # - - - - - - - - - - - - - 
    # Plot the LW all data
    # - - - - - - - - - - - - - 
    ax5 = plt.subplot(gs[2,1],projection=mapcrs)
    ax5.set_extent([-180,180,45,90],datacrs)
    ax5.gridlines()
    #ax = plt.axes(projection=ccrs.Miller())
    file_season = '_' + CERES_lw_all_dict['season']
    title_season = ' ' + CERES_lw_all_dict['season']
    if(CERES_lw_all_dict['season']=='all'):
        file_season=''
        title_season=''
    # Ignore the areas without ice by making any values where the ice conc
    # is less than 5%.
    local_LW = np.copy(CERES_lw_all_dict['trends'])
    local_LW[np.isnan(avg_ice_conc)] = np.nan
    plot_local_LW = ma.masked_invalid(local_LW)
        
    mesh = ax5.pcolormesh(CERES_lw_clr_dict['lon'],CERES_lw_clr_dict['lat'],plot_local_LW,\
            transform=datacrs,vmin=lw_clr_min_val,vmax=lw_clr_max_val,cmap=colormap)
    ax5.set_title('Terra CERES TOA All-sky LWF Trend')
    #ax3.set_title('Terra CERES TOA Clear-sky LWF' + title_season.title() + ' Trend\n' + \
    #             CERES_lw_clr_dict['dates'][0] + ' - ' + CERES_lw_clr_dict['dates'][-1])
    ender = '.png'
    if(adjusted==True):
        ax5.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        ##axins = inset_axes(ax2,width="5%",height="100%",loc='lower left',\
        ##            bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax2.transAxes,\
        ##            borderpad=0)
        ##cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=CERES_lw_clr_dict['parm_unit'])
        cbar = plt.colorbar(mesh,extend='both',orientation='horizontal',pad=0,aspect=50,label=CERES_lw_all_dict['parm_unit'] + ' / Study Period')
        #ax3.set_xlim(-5170748.535086173,5167222.438879491)
        #ax3.set_ylim(-3913488.8763307533,3943353.899053069)
        ax5.set_xlim(-4170748.535086173,4167222.438879491)
        ax5.set_ylim(-2913488.8763307533,2943353.899053069)
        ender = '_adjusted' + ender
    else:
        cbar = plt.colorbar(mesh,cmap=colormap,label=CERES_lw_all_dict['parm_unit'])
    #cbar.set_label(parm_unit)
    ax5.coastlines()
    
    if(ice_data['season_adder'].strip()=='sunlight'):
        plt.savefig('paper_figure1_sunlight.png',dpi=300)
    else:
        plt.savefig('paper_figure1.png',dpi=300)
    plt.show()

def figure_2(ice_data,CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_net_clr_dict,\
         CERES_lw_all_dict,CERES_sw_all_dict,CERES_net_all_dict,inseason,adjusted=True):
    area=False

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

    # - - - - - - - - - - - - - 
    # Calculate average summer ice conc
    # - - - - - - - - - - - - - 
    # Find the average conc
    avg_ice_conc = np.zeros((len(ice_data['grid_lat']),len(ice_data['grid_lon'])))
    for xi in range(len(ice_data['grid_lat'])):
        for yj in range(len(ice_data['grid_lon'])):
            avg_ice_conc[xi,yj] = np.average(ice_data['grid_ice_conc'][:,xi,yj][ice_data['grid_ice_conc'][:,xi,yj]!=-99.])
            if(avg_ice_conc[xi,yj]<1.0):
                avg_ice_conc[xi,yj] = np.nan
       
    #plot_good_data = ma.masked_invalid(avg_ice_conc)

    # Quick comparison
    markersize=6
    fig, axs = plt.subplots(2,3)
    fig.set_size_inches(19,12)
    plt.title("CERES and NSIDC " + inseason.title() + " Trends")

    local_ice_trend = np.copy(ice_data[ice_trend_key])
    avg_ice_conc[np.isnan(avg_ice_conc)] = -420.
    local_ice_trend[(local_ice_trend==missing_trend)] = -420.
    local_ice_trend[np.isnan(local_ice_trend)] = -420.
    # Plot Clear-sky SWF data
    plotx = local_ice_trend[(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    ploty = CERES_sw_clr_dict['trends'][(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    
    plotx = plotx[ploty>-1000]
    ploty = ploty[ploty>-1000]
    plotx = plotx[ploty<1000]
    ploty = ploty[ploty<1000]
    corr_swf = correlation(plotx,ploty)
    print("SWF correlation: ",corr_swf)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    # Modify data for the color density plot
    xy = np.vstack([plotx,ploty])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x,y,z = plotx[idx], ploty[idx], z[idx]
    axs[0,0].scatter(x,y,c=z,s=markersize,cmap=plt.cm.jet)
    #axs[0,0].scatter(plotx,ploty,s=markersize,color='black')
    axs[0,0].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='black')
    #axs[0,0].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='red')
    axs[0,0].set_title('Clear-Sky SWF' + ice_data['season_adder'].title() + ' Trends')
    axs[0,0].set_xlabel('Percent Ice Concentration Trends [%]')
    axs[0,0].set_ylabel('SWF Trends [W/m2]')
    x_pos = ((axs[0,0].get_xlim()[1]-axs[0,0].get_xlim()[0])/8) + axs[0,0].get_xlim()[0]
    y_pos = axs[0,0].get_ylim()[1]-((axs[0,0].get_ylim()[1]-axs[0,0].get_ylim()[0])/8)
    axs[0,0].text(x_pos,y_pos,"Corr = " + str(np.round(corr_swf,3)))
    # Plot LWF data
    plotx = local_ice_trend[(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    ploty = CERES_lw_clr_dict['trends'][(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    plotx = plotx[ploty>-1000]
    ploty = ploty[ploty>-1000]
    plotx = plotx[ploty<1000]
    ploty = ploty[ploty<1000]
    corr_lwf = correlation(plotx,ploty)
    print("LWF correlation: ",corr_lwf)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    xy = np.vstack([plotx,ploty])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x,y,z = plotx[idx], ploty[idx], z[idx]
    axs[0,1].scatter(x,y,c=z,s=markersize,cmap=plt.cm.jet)
    axs[0,1].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='black')
    axs[0,1].set_title('Clear-Sky LWF' + ice_data['season_adder'].title() + ' Trends')
    axs[0,1].set_xlabel('Percent Ice Concentration Trends[%]')
    axs[0,1].set_ylabel('LWF Trends [W/m2]')
    x_pos = ((axs[0,1].get_xlim()[1]-axs[0,1].get_xlim()[0])/8) + axs[0,1].get_xlim()[0]
    y_pos = axs[0,1].get_ylim()[1]-((axs[0,1].get_ylim()[1]-axs[0,1].get_ylim()[0])/8)
    axs[0,1].text(x_pos,y_pos,"Corr = " + str(np.round(corr_lwf,3)))
    ##print("  y_max = ",axs[0,1].get_ylim()[1],"  y_min = ",axs[0,1].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    # Plot Net Flux Data
    plotx = local_ice_trend[(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    ploty = CERES_net_clr_dict['trends'][(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    plotx = plotx[ploty>-1000]
    ploty = ploty[ploty>-1000]
    plotx = plotx[ploty<1000]
    ploty = ploty[ploty<1000]
    corr_net = correlation(plotx,ploty)
    print("Net flux correlation: ",corr_net)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope," + ",intercept)
    xy = np.vstack([plotx,ploty])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x,y,z = plotx[idx], ploty[idx], z[idx]
    axs[0,2].scatter(x,y,c=z,s=markersize,cmap=plt.cm.jet)
    axs[0,2].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='black')
    axs[0,2].set_title('Clear-Sky Net Flux' + ice_data['season_adder'].title() + ' Trends')
    axs[0,2].set_xlabel('Percent Ice Concentration Trends [%]')
    axs[0,2].set_ylabel('Net Flux Trends [W/m2]')
    x_pos = ((axs[0,2].get_xlim()[1]-axs[0,2].get_xlim()[0])/8) + axs[0,2].get_xlim()[0]
    y_pos = ((axs[0,2].get_ylim()[1]-axs[0,2].get_ylim()[0])/8) + axs[0,2].get_ylim()[0]
    axs[0,2].text(x_pos,y_pos,"Corr = " + str(np.round(corr_net,3)))

    # Plot the second row: all sky
    # Plot All-sky SWF data
    plotx = local_ice_trend[(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    ploty = CERES_sw_all_dict['trends'][(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    plotx = plotx[ploty>-1000]
    ploty = ploty[ploty>-1000]
    plotx = plotx[ploty<1000]
    ploty = ploty[ploty<1000]
    corr_swf = correlation(plotx,ploty)
    print("All SWF correlation: ",corr_swf)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    xy = np.vstack([plotx,ploty])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x,y,z = plotx[idx], ploty[idx], z[idx]
    axs[1,0].scatter(x,y,c=z,s=markersize,cmap=plt.cm.jet)
    axs[1,0].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='black')
    axs[1,0].set_title('All-Sky SWF' + ice_data['season_adder'].title() + ' Trends')
    axs[1,0].set_xlabel('Percent Ice Concentration Trends')
    axs[1,0].set_ylabel('SWF Trends [W/m2]')
    x_pos = ((axs[1,0].get_xlim()[1]-axs[1,0].get_xlim()[0])/8) + axs[1,0].get_xlim()[0]
    y_pos = axs[1,0].get_ylim()[1]-((axs[1,0].get_ylim()[1]-axs[1,0].get_ylim()[0])/8)
    axs[1,0].text(x_pos,y_pos,"Corr = " + str(np.round(corr_lwf,3)))
    # Plot All-sky LWF data
    plotx = local_ice_trend[(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    ploty = CERES_lw_all_dict['trends'][(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    plotx = plotx[ploty>-1000]
    ploty = ploty[ploty>-1000]
    plotx = plotx[ploty<1000]
    ploty = ploty[ploty<1000]
    corr_lwf = correlation(plotx,ploty)
    print("LWF correlation: ",corr_lwf)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*  ",intercept)
    xy = np.vstack([plotx,ploty])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x,y,z = plotx[idx], ploty[idx], z[idx]
    axs[1,1].scatter(x,y,c=z,s=markersize,cmap=plt.cm.jet)
    axs[1,1].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='black')
    axs[1,1].set_title('All-Sky LWF' + ice_data['season_adder'].title() + ' Trends')
    axs[1,1].set_xlabel('Ice Trends')
    axs[1,1].set_ylabel('LWF Trends')
    x_pos = ((axs[1,1].get_xlim()[1]-axs[1,1].get_xlim()[0])/8) + axs[1,1].get_xlim()[0]
    y_pos = axs[1,1].get_ylim()[1]-((axs[1,1].get_ylim()[1]-axs[1,1].get_ylim()[0])/8)
    axs[1,1].text(x_pos,y_pos,"Corr = " + str(np.round(corr_swf,3)))
    ##print("  y_max = ",axs[0,1].get_ylim()[1],"  y_min = ",axs[0,1].get_ylim()[0])
    ##print("  x_pos = ",x_pos)
    ##print("  y_pos = ",y_pos)
    # Plot Net Flux Data
    plotx = local_ice_trend[(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    ploty = CERES_net_all_dict['trends'][(local_ice_trend!=-420.) & (avg_ice_conc!=-420.)]
    plotx = plotx[ploty>-1000]
    ploty = ploty[ploty>-1000]
    plotx = plotx[ploty<1000]
    ploty = ploty[ploty<1000]
    corr_net = correlation(plotx,ploty)
    print("Net flux correlation: ",corr_net)
    slope,intercept,r_val,p_val,stderr = stats.linregress(plotx,ploty)
    print("  Equation: y = ",slope,"*x + ",intercept)
    xy = np.vstack([plotx,ploty])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x,y,z = plotx[idx], ploty[idx], z[idx]
    axs[1,2].scatter(x,y,c=z,s=markersize,cmap=plt.cm.jet)
    axs[1,2].plot(np.unique(plotx),np.poly1d(np.polyfit(plotx,ploty,1))(np.unique(plotx)),color='black')
    axs[1,2].set_title('All-Sky Net Flux' + ice_data['season_adder'].title() + ' Trends')
    axs[1,2].set_xlabel('Ice Trends')
    axs[1,2].set_ylabel('Net Flux Trends')
    x_pos = ((axs[1,2].get_xlim()[1]-axs[1,2].get_xlim()[0])/8) + axs[1,2].get_xlim()[0]
    y_pos = ((axs[1,2].get_ylim()[1]-axs[1,2].get_ylim()[0])/8) + axs[1,2].get_ylim()[0]
    axs[1,2].text(x_pos,y_pos,"Corr = " + str(np.round(corr_net,3)))

    if(ice_data['season_adder'].strip()=='sunlight'):
        plt.savefig('paper_figure2_sunlight.png',dpi=300)
    else:
        plt.savefig('paper_figure2.png',dpi=300)
    plt.show()
    #if(save==True):
    #    plt.savefig(outname,dpi=300)
    #    print("Saved image " + outname)
    #else:
    #    plt.show()

def figure_3():
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
   
    # These values come from total_ice_summer_area_obs, which I believe
    # comes from the ice data Excel sheets from the NSIDC website 
    sum_init_ice =  7548960.59925     # 1989
    sum_final_ice =  6102900.002686667 # 2018
    #sum_init_ice  = 7820576049633.293/1e6  # km2 as of June 2001
    #sum_final_ice = 6481114233056.605/1e6  # km2 
    
    # 10 years
    #starting_val = -10. 
    starting_val = ((sum_final_ice-sum_init_ice)/sum_init_ice)*100.
    start_year = 2018
    
    sigma_1m = np.array(fill_array(starting_val,1.))
    sigma_2m = np.array(fill_array(starting_val,2.))
    sigma_3m = np.array(fill_array(starting_val,3.))
    sigma_4m = np.array(fill_array(starting_val,4.))
    sigma_5m = np.array(fill_array(starting_val,5.))
    
    extents_1m = (sum_init_ice*(100. + sigma_1m)/100.)/1e6 
    years_1m   = np.arange(start_year,start_year + len(extents_1m))
    extents_2m = (sum_init_ice*(100. + sigma_2m)/100.)/1e6 
    years_2m   = np.arange(start_year,start_year + len(extents_2m))
    extents_3m = (sum_init_ice*(100. + sigma_3m)/100.)/1e6 
    years_3m   = np.arange(start_year,start_year + len(extents_3m))
    extents_4m = (sum_init_ice*(100. + sigma_4m)/100.)/1e6 
    years_4m   = np.arange(start_year,start_year + len(extents_4m))
    extents_5m = (sum_init_ice*(100. + sigma_5m)/100.)/1e6 
    years_5m   = np.arange(start_year,start_year + len(extents_5m))
    
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
    plt.savefig('paper_figure3.png',dpi=300)

# Plot the observed average sea ice areas over the time period.
def figure_4(ice_data,model_overlay=False,zoomed=True):
    inseason = ice_data['season_adder'].strip()
    divider=3
    num_months_pyear = 3
    if(inseason=='sunlight'):
        divider=6
        num_months_pyear = 6
    summer_averages_grid = np.zeros((int(len(ice_data['titles'])/divider)))
    summer_averages_raw  = np.zeros((int(len(ice_data['titles'])/divider)))
    #summer_averages_grid = np.zeros((int(len(ice_data['titles'])/3)))
    #summer_averages_raw  = np.zeros((int(len(ice_data['titles'])/3)))
    all_summer_avgs = []
    local_grid_ice = np.copy(ice_data['grid_ice_conc'])
    local_grid_ice[local_grid_ice==-99.] = np.nan

    ## Calculate yearly average ice fraction.
    #for ti in range(0,len(season_ice_avgs)):
    #    idx = (ti*divider)
    #    season_ice_avgs[ti] = np.average(ice_data['grid_ice_conc'][idx:idx+divider,xi,yj])

    for ti in range(0,len(summer_averages_raw)):
        # Calculate the ice area from the raw data
        idx = (ti*divider)
        print("i=",ti," idx=",idx)

        ##!#summer_averages_raw[ti] = np.sum((ice_data['data'][idx:idx + divider,:,:]\
        ##!#    [ice_data['data'][idx:idx + divider,:,:]<251]/100.)*\
        ##!#    ice_data['area'][ice_data['data'][:,:]<251]) / divider


        # Attempting to redo how the summer averages are calculated to improve
        # the calculated average values and make more efficient.
        raw_sum1 = np.sum((ice_data['data'][idx,:,:][ice_data['data'][idx,:,:]<251]/100.)*ice_data['area'][ice_data['data'][idx,:,:]<251])         
        raw_sum2 = np.sum((ice_data['data'][idx + 1,:,:][ice_data['data'][idx + 1,:,:]<251]/100.)*ice_data['area'][ice_data['data'][idx + 1,:,:]<251])         
        raw_sum3 = np.sum((ice_data['data'][idx + 2,:,:][ice_data['data'][idx + 2,:,:]<251]/100.)*ice_data['area'][ice_data['data'][idx + 2,:,:]<251])         
        all_summer_avgs.append(raw_sum1)
        all_summer_avgs.append(raw_sum2)
        all_summer_avgs.append(raw_sum3)
        if(inseason!='sunlight'):
            avg_raw = (raw_sum1 + raw_sum2 + raw_sum3)/3.
        else:
            raw_sum4 = np.sum((ice_data['data'][idx + 3,:,:][ice_data['data'][idx + 3,:,:]<251]/100.)*ice_data['area'][ice_data['data'][idx + 3,:,:]<251])         
            raw_sum5 = np.sum((ice_data['data'][idx + 4,:,:][ice_data['data'][idx + 4,:,:]<251]/100.)*ice_data['area'][ice_data['data'][idx + 4,:,:]<251])         
            raw_sum6 = np.sum((ice_data['data'][idx + 5,:,:][ice_data['data'][idx + 5,:,:]<251]/100.)*ice_data['area'][ice_data['data'][idx + 5,:,:]<251])         
            all_summer_avgs.append(raw_sum4)
            all_summer_avgs.append(raw_sum5)
            all_summer_avgs.append(raw_sum6)
            avg_raw = (raw_sum1 + raw_sum2 + raw_sum3 + raw_sum4 + raw_sum5 + raw_sum6)/6.
        summer_averages_raw[ti] = avg_raw
        # Calculate the area from the gridded data

    #all_summer_avgs = np.array(all_summer_avgs)

    ##!##print(all_summer_avgs) 
    ##!#fig1 = plt.figure()
    ##!#plt.plot(summer_averages_raw)
    ##!##plt.scatter(np.arange(len(all_summer_avgs))[::divider],summer_averages_raw)
    ##!#plt.show()
    ##!#plt.close() 
    ##!#return summer_averages_raw
 
    # Clear-sky summer values from net flux figure
    dflux_dsigma = 1.2645
    #dflux_dsigma = -1.2645
    #if(inseason=='sunlight'): dflux_dsigma = -0.3600
    del_T = (86400.)*30*6  # number of seconds of solar heating per sunlit time 
    l_f = 3.3e5  # latent heat of fusion (J/kg)
    rho_i = 917  # density of ice (kg/m3)
    def equation(sigma_p,thick):
    
        sigma_pnew = sigma_p * (1. + (dflux_dsigma * del_T) / (l_f * rho_i* thick))
        return sigma_pnew
    
    def fill_array(start_val,thick):
        sigma_values = []
        sigma_values.append(start_val)
        # Find 1m-thickness values
        while(start_val>-100.):
            sigma_new = equation(start_val,thick)
            print(sigma_new)
            sigma_values.append(sigma_new) 
            start_val = sigma_new
        return sigma_values
    
    sum_init_ice  = 7820576049633.293/1e6  # km2 as of June 2001
    sum_final_ice = 6481114233056.605/1e6  # km2 
    #sum_init_ice  = 7820576049633.293/1e6  # km2 as of June 2001
    #sum_final_ice = 6481114233056.605/1e6  # km2 
    
    # 10 years
    starting_val = -10. 
    # Use the 1991 summer average extent as the new "sum init ice" value   
    ##!## Use the 2001 summer average extent as the new "sum final ice" value
    sum_init_ice = summer_averages_raw[9] # 1989
    if(model_overlay==False):
        sum_final_ice = summer_averages_raw[38] # 2018
        #sum_final_ice = summer_averages_raw[38] # 2018
    else:
        sum_final_ice = summer_averages_raw[21] # 2001
    #starting_val = 100.
    #starting_val = ((sum_final_ice-sum_init_ice)/sum_init_ice)*100.
    #starting_val = ((summer_averages_raw[0]-sum_init_ice)/sum_init_ice)*100.
    starting_val = ((sum_final_ice-sum_init_ice)/sum_init_ice)*100.
    beginning_year = int(ice_data['titles'][0].split('_')[1][:4])
    if(model_overlay==False):
        if(inseason=='summer'):
            start_year = int(ice_data['titles'][116].split('_')[1][:4]) # 2018
        else:
            start_year = int(ice_data['titles'][232].split('_')[1][:4]) # 2018
    else:
        if(inseason=='summer'):
            start_year = int(ice_data['titles'][63].split('_')[1][:4]) # 2001
        else:
            start_year = int(ice_data['titles'][126].split('_')[1][:4]) # 2001
    
    # Use 1990 as the starting year
    # NOTE: This assumes that old data for summer are read in.
    #       old_summer_data = read_ice('summer',pre2001=True)
    #start_year = int(ice_data['titles'][30].split('_')[1][:4]) # 2001
    #start_year = 2001

    ### Write the obs to a file
    ##with open('total_ice_summer_area_obs.txt','w') as fout:
    ##    for ti in range(len(summer_averages_raw)):
    ##        out_year = beginning_year + ti
    ##        fout.write(str(out_year) + '   ' + str(summer_averages_raw[ti]) + '\n')
    
    
    sigma_1m = np.array(fill_array(starting_val,1.))
    sigma_2m = np.array(fill_array(starting_val,2.))
    sigma_3m = np.array(fill_array(starting_val,3.))
    sigma_4m = np.array(fill_array(starting_val,4.))
    sigma_5m = np.array(fill_array(starting_val,5.))
    # Use the 2001 average extent as the 
    extents_1m = (sum_init_ice*(100. + sigma_1m)/100.)/1e6 
    extents_2m = (sum_init_ice*(100. + sigma_2m)/100.)/1e6 
    extents_3m = (sum_init_ice*(100. + sigma_3m)/100.)/1e6 
    extents_4m = (sum_init_ice*(100. + sigma_4m)/100.)/1e6 
    extents_5m = (sum_init_ice*(100. + sigma_5m)/100.)/1e6 
    #extents_2m = (summer_averages_raw[0]*(100. + sigma_2m)/100.)/1e6 
    #extents_2m = (sum_init_ice*(100. + sigma_2m)/100.)/1e6 
    years_1m   = np.arange(start_year,start_year + len(extents_1m))
    years_2m   = np.arange(start_year,start_year + len(extents_2m))
    years_3m   = np.arange(start_year,start_year + len(extents_3m))
    years_4m   = np.arange(start_year,start_year + len(extents_4m))
    years_5m   = np.arange(start_year,start_year + len(extents_5m))

    # Print off the forecasted ending years
    print("1 meter end = ",years_1m[-1])
    print("2 meter end = ",years_2m[-1])
    print("3 meter end = ",years_3m[-1])
    print("4 meter end = ",years_4m[-1])
    print("5 meter end = ",years_5m[-1])

    # All summer averages x vals
    asa_x_vals = np.arange(beginning_year, 2018.9, 1 / divider)
    #print("First model extent: ",extents_2m[0])
    #print("Summer Averages: ",summer_averages_raw)
    #print("All monthly Averages: ",all_summer_avgs)
    #plt.plot(asa_x_vals,all_summer_avgs/1e6,label='All summer averages')

    fig, ax = plt.subplots()
    if(zoomed==True):
        ax.scatter(asa_x_vals[::divider],summer_averages_raw/1e6,\
            s=6,label='Observed',color='black')
        #ax.scatter(asa_x_vals[::divider],summer_averages_raw/1e6,\
        #    s=6,label='Observed',color='black')
    else:
        ax.plot(asa_x_vals[::divider],summer_averages_raw/1e6,\
            label='Observed',color='black')
    ax.plot(years_1m,extents_1m,label='1m Model extents')
    ax.plot(years_2m,extents_2m,label='2m Model extents')
    ax.plot(years_3m,extents_3m,label='3m Model extents')
    ax.plot(years_4m,extents_4m,label='4m Model extents')
    ax.plot(years_5m,extents_5m,label='5m Model extents')
    #if(model_overlay==True):
    #    ax.set_xlim(1980,2018)
    #    ax.set_ylim(5.5,8.0)
    #else:
    if(model_overlay==False):
        ax.axvline(1989,linestyle=':',ymin=0.85,color='black')
        ax.axvline(2018,linestyle=':',ymin=0.7,ymax=0.85,color='black')
        ax.axhline(sum_init_ice/1e6,xmin=0.05,xmax=0.1,linestyle=':',color='black')
        ax.axhline(sum_final_ice/1e6,xmin=0.143,xmax=0.193,linestyle=':',color='black')
        #ax.text(1985,5.5,'1989')
    if(zoomed==True):
        ax.set_xlim(1990,2018)
        ax.set_ylim(5.5,8.0)
    ax.legend()
    ax.set_ylabel('Ice Area (millions of km2)')
    ax.set_title("Average "+inseason.title()+" Arctic Ice Area")
    #if(model_overlay==True):
    #    if(zoomed==True):
    #        plt.savefig('paper_figure4_model_overlay_zoom.png',dpi=300)
    #    else:
    #        plt.savefig('paper_figure4_model_overlay.png',dpi=300)
    #else:
    #    plt.savefig('paper_figure4.png',dpi=300)
    plt.show()

# write_toASCII writes the contents of each dictionary to an ASCII file
# Columns in these files are as follows:
# lat  lon  trend
def write_toASCII(ice_data,CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_net_clr_dict,\
                  CERES_lw_all_dict,CERES_sw_all_dict,CERES_net_all_dict,\
                  nonzero=True):
    inseason=ice_data['season_adder'].strip()
    divider = 3
    if(inseason == 'sunlight'):
        divider = 6

    # Find the length of the time dimension, divide by divider to determine
    # the number of years
    num_months = int(len(ice_data['titles']))
    num_years = int(num_months / divider)
    season_ice_avgs = np.zeros(num_years)

    # Set up the format string for printing off the data
    # format_list includes the formats for all the yearly average concentrations
    # Separate each by a space
    formats = ' '.join(['{:10.6}' for item in season_ice_avgs])
    # Combine the formats string with the other format string
    total_ice_format = ' '.join(['{:.3} {:8.4} {:10.6} {:10.6}',formats,'\n'])

    # Set up the header line format string
    # Extract the individual years from the titles parameter
    string_years = [tstring[7:11] for tstring in ice_data['titles'][::divider]]
    header_formats = ' '.join(['{:10}' for item in season_ice_avgs])
    total_header_format = ' '.join(['{:7} {:8} {:11} {:11}',header_formats,'\n'])   
 
    # - - - - - - - - - - - - - 
    # Plot the average summer ice conc
    # - - - - - - - - - - - - - 
    # Find the average conc
    avg_ice_conc = np.full((len(ice_data['grid_lat']),len(ice_data['grid_lon'])),-999.)
    ##!#for xi in range(len(ice_data['grid_lat'])):
    ##!#    for yj in range(len(ice_data['grid_lon'])):
    ##!#        avg_ice_conc[xi,yj] = np.average(ice_data['grid_ice_conc'][:,xi,yj][ice_data['grid_ice_conc'][:,xi,yj]>0.])
       
    # Create the Ice data ASCII file
    ioutname = "nsidc_gridded_ice_trends_" + inseason + ".txt"
    coutname = "ceres_gridded_flux_trends_" + inseason + ".txt"
    if(nonzero==True):
        ioutname = "nsidc_gridded_ice_trends_" + inseason + "_nonzero.txt"
        coutname = "ceres_gridded_flux_trends_" + inseason + "_nonzero.txt"
    print("Generating file ",ioutname)
    with open(ioutname,'w') as ifout:
        ifout.write(total_header_format.format('Lat','Lon','Ice_trend','Avg_conc',*string_years))
        with open(coutname,'w') as cfout:
            cfout.write('{:7} {:8} {:9} {:9} {:9} {:9} {:9} {:9}\n'.format('Lat','Lon','sw_clr','sw_all','lw_clr','lw_all','net_clr','net_all'))
            for xi in range(len(ice_data['grid_lat'])):
                for yj in range(len(ice_data['grid_lon'])):

                    # Calculate total average ice fraction
                    # NOTE: Before 2020/05/28, this line calculated the average
                    #       value of all ice fractions greater than 0., so this
                    #       was screening out 0s. Likely screwed up the average
                    #       values.
                    avg_ice_conc[xi,yj] = np.average(ice_data['grid_ice_conc'][:,xi,yj])

                    # Calculate yearly average ice fraction.
                    for ti in range(0,len(season_ice_avgs)):
                        idx = (ti*divider)
                        season_ice_avgs[ti] = np.average(ice_data['grid_ice_conc'][idx:idx+divider,xi,yj])

                    # Only write to file if the ice fractions at the current
                    # location are always nonzero.
                    if(nonzero==True):
                        if(0.0 not in ice_data['grid_ice_conc'][:,xi,yj]):
                            ifout.write(total_ice_format.format(ice_data['grid_lat'][xi],\
                                    ice_data['grid_lon'][yj],ice_data['grid_ice'][xi,yj],avg_ice_conc[xi,yj],*season_ice_avgs))
                            cfout.write('{:.3} {:8.4} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f}\n'.format(CERES_sw_clr_dict['lat'][xi],\
                                        CERES_sw_clr_dict['lon'][yj],CERES_sw_clr_dict['trends'][xi,yj],CERES_sw_all_dict['trends'][xi,yj],\
                                        CERES_lw_clr_dict['trends'][xi,yj],CERES_lw_all_dict['trends'][xi,yj],CERES_net_clr_dict['trends'][xi,yj],\
                                        CERES_net_all_dict['trends'][xi,yj]))
                    else:
                        #ifout.write('{:.3} {:8.4} {:10.6} {:10.6}\n'.format(ice_data['grid_lat'][xi],\
                        ifout.write(total_ice_format.format(ice_data['grid_lat'][xi],\
                                ice_data['grid_lon'][yj],ice_data['grid_ice'][xi,yj],\
                                avg_ice_conc[xi,yj],*season_ice_avgs))
                        cfout.write('{:.3} {:8.4} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f}\n'.format(CERES_sw_clr_dict['lat'][xi],\
                                    CERES_sw_clr_dict['lon'][yj],CERES_sw_clr_dict['trends'][xi,yj],CERES_sw_all_dict['trends'][xi,yj],\
                                    CERES_lw_clr_dict['trends'][xi,yj],CERES_lw_all_dict['trends'][xi,yj],CERES_net_clr_dict['trends'][xi,yj],\
                                    CERES_net_all_dict['trends'][xi,yj]))
    

    ### Create the CERES Clear-sky data ASCII file
    ##coutname = "ceres_gridded_flux_trends_sunlight.txt"
    ##if(nonzero==True):
    ##    coutname = "ceres_gridded_flux_trends_sunlight_nonzero.txt"
    ##print("Generating file ",outname)
    ##with open(outname,'w') as cfout:
    ##    cfout.write('{:7} {:8} {:9} {:9} {:9} {:9} {:9} {:9}\n'.format('Lat','Lon','sw_clr','sw_all','lw_clr','lw_all','net_clr','net_all'))
    ##    for xi in range(len(CERES_sw_clr_dict['lat'])):
    ##        for yj in range(len(CERES_sw_clr_dict['lon'])):
    ##            # Only write to file if the ice fractions at the current
    ##            # location are always nonzero.
    ##            cfout.write('{:.3} {:8.4} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f}\n'.format(CERES_sw_clr_dict['lat'][xi],\
    ##                        CERES_sw_clr_dict['lon'][yj],CERES_sw_clr_dict['trends'][xi,yj],CERES_sw_all_dict['trends'][xi,yj],\
    ##                        CERES_lw_clr_dict['trends'][xi,yj],CERES_lw_all_dict['trends'][xi,yj],CERES_net_clr_dict['trends'][xi,yj],\
    ##                        CERES_net_all_dict['trends'][xi,yj]))

def write_Ice(ice_data):
    print("yay")

    # Print of a file for each month
    num_months = int(len(ice_data['titles']))

    # Set up the format string for printing off the data
    # format_list includes the formats for all the yearly average concentrations
    # Separate each by a space
    formats = '\t'.join(['{:6.4}' for item in ice_data['titles']])
    # Combine the formats string with the other format string
    total_ice_format = '\t'.join(['{:.3}\t{:6.4}',formats,'\n'])

    # Set up the header line format string
    # Extract the individual years from the titles parameter
    string_dates = [tstring.strip()[3:9] for tstring in ice_data['titles']]
    header_formats = '\t'.join(['{:6}' for item in ice_data['titles']])
    total_header_format = '\t'.join(['{:6}\t{:6}',header_formats,'\n'])   
   
    filename = "nsidc_grid_ice_values.txt"
    #for ti in range(num_months):
    #    str_date = ice_data['titles'][ti].strip()[3:9]
    #    filename = "nsidc_grid_ice_"+str_date+'.txt'
    with open(filename,'w') as fout:
        fout.write(total_header_format.format('Lat','Lon',*string_dates))
        for xi in range(len(ice_data['grid_lat'])):
            print("Printing latitude ",ice_data['grid_lat'][xi])
            for yj in range(len(ice_data['grid_lon'])):
                fout.write(total_ice_format.format(ice_data['grid_lat'][xi],\
                        ice_data['grid_lon'][yj],*ice_data['grid_ice_conc'][:,xi,yj]))

    print("Saved file",filename) 
    

def scatter_1a_1b(ice_data,adjusted=True):
    fig = plt.figure()
    if(ice_data['season_adder'].strip()=='sunlight'):
        plt.title('Sunlight Months (April - September  2001 - 2018)')
    else: 
        plt.title('Summer Only (June - August,  2001 - 2018)')
        
    # - - - - - - - - - - - - - 
    # Plot the average summer ice conc
    # - - - - - - - - - - - - - 
    # Find the average conc
    avg_ice_conc = np.zeros((len(ice_data['grid_lat']),len(ice_data['grid_lon'])))
    for xi in range(len(ice_data['grid_lat'])):
        for yj in range(len(ice_data['grid_lon'])):
            avg_ice_conc[xi,yj] = np.average(ice_data['grid_ice_conc'][:,xi,yj][ice_data['grid_ice_conc'][:,xi,yj]!=-99.])
            if(avg_ice_conc[xi,yj]<20.0):
                avg_ice_conc[xi,yj] = np.nan
       
    local_grid_ice = np.copy(ice_data['grid_ice'])
    local_grid_ice_bad = np.copy(ice_data['grid_ice'])
   
    markersize=8 
    local_grid_ice[local_grid_ice==-999.] = np.nan
    # Ignore the areas without ice by making any values where the ice conc
    # is less than 5%.
    local_grid_ice[np.isnan(avg_ice_conc)] = np.nan
    local_grid_ice_bad[local_grid_ice!=-999.] = np.nan
    avg_ice_conc[np.isnan(local_grid_ice)] = np.nan
    plot_good_data_avg = ma.masked_invalid(avg_ice_conc)
    plot_good_data_trend = ma.masked_invalid(local_grid_ice)
    plot_land_data = ma.masked_invalid(local_grid_ice_bad)
    plot_good_data_avg = plot_good_data_avg[np.isnan(plot_good_data_avg)==False]
    plot_good_data_trend = plot_good_data_trend[np.isnan(plot_good_data_trend)==False]
    xy = np.vstack([plot_good_data_avg,plot_good_data_trend])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x,y,z = plot_good_data_avg[idx], plot_good_data_trend[idx], z[idx]
    plt.scatter(x,y,c=z,s=markersize,cmap=plt.cm.jet)
    plt.xlabel('Average Ice Fraction [%]')
    plt.ylabel('Ice Fraction Trend [%/study period]')
    if(ice_data['season_adder'].strip()=='sunlight'):
        plt.savefig('scatter_1a_1b_sunlight.png',dpi=300)
        print("Saved image scatter_1a_1b_sunlight.png")
    else: 
        plt.savefig('scatter_1a_1b_summer.png',dpi=300)
        print("Saved image scatter_1a_1b_summer.png")
    plt.show()

def yearly_conc_diff_scatter(ice_data):
    # Make the figure
   
    sunlight=False
    summer=False
    if(ice_data['season_adder']==' sunlight'):
        sunlight=True
        title_divider = 6
    elif(ice_data['season_adder']==' summer'):
        summer=True
        title_divider = 3
    # Loop over the years in the dictionary
    num_years = int(len(ice_data['titles'])/title_divider)
    total_this_year = np.array([])
    total_next_year = np.array([])
    total_diff = np.array([])
    for ti in range(num_years-1):
        idx = ti*title_divider
        nidx = (ti + 1)*title_divider
        this_year_string = ice_data['titles'][idx][7:11]
        next_year_string = ice_data['titles'][nidx][7:11]
        print("Looking at ",next_year_string)
        if(sunlight==True):
            this_year = np.array([ice_data['grid_ice_conc'][idx,:,:], \
                                 ice_data['grid_ice_conc'][idx + 1,:,:], \
                                 ice_data['grid_ice_conc'][idx + 2,:,:], \
                                 ice_data['grid_ice_conc'][idx + 3,:,:], \
                                 ice_data['grid_ice_conc'][idx + 4,:,:], \
                                 ice_data['grid_ice_conc'][idx + 5,:,:]])
            next_year = np.array([ice_data['grid_ice_conc'][nidx,:,:], \
                                 ice_data['grid_ice_conc'][nidx + 1,:,:], \
                                 ice_data['grid_ice_conc'][nidx + 2,:,:], \
                                 ice_data['grid_ice_conc'][nidx + 3,:,:], \
                                 ice_data['grid_ice_conc'][nidx + 4,:,:], \
                                 ice_data['grid_ice_conc'][nidx + 5,:,:]])
        elif(summer==True):
            this_year = np.array([ice_data['grid_ice_conc'][idx,:,:], \
                                 ice_data['grid_ice_conc'][idx + 1,:,:], \
                                 ice_data['grid_ice_conc'][idx + 2,:,:]])
            next_year = np.array([ice_data['grid_ice_conc'][nidx,:,:], \
                                 ice_data['grid_ice_conc'][nidx + 1,:,:], \
                                 ice_data['grid_ice_conc'][nidx + 2,:,:]])
        this_year = np.average(this_year,axis=0)
        next_year = np.average(next_year,axis=0)

        diff = next_year-this_year

        #total_next_year = np.append(total_next_year[:],next_year)
        #total_diff = np.append(total_diff[:],diff)
        diff = diff[(next_year!=-99.) & (this_year!=-99.)]
        next_year = next_year[(next_year!=-99.) & (this_year!=-99.)]
        if(len(diff[diff>100])>0):
            print(diff[diff>100])

        fig1 = plt.figure()
        markersize=8
        xy = np.vstack([next_year,diff])
        z = gaussian_kde(xy)(xy)
        idx = z.argsort()
        x,y,z = next_year[idx], diff[idx], z[idx]
        plt.scatter(x,y,c=z,s=markersize,cmap=plt.cm.jet)
        #plt.scatter(next_year[next_year!=-99.],diff[next_year!=-99.])
        plt.xlabel(next_year_string + ' Ice Fraction')
        plt.ylabel(next_year_string + ' Ice Fraction - ' + this_year_string + ' Ice Fraction')
        if(ice_data['season_adder'].strip()=='sunlight'):
            plt.title('Sunlight data (April - September')
            outname = 'yearly_conc_diff_scatter_' + next_year_string + '_sunlight.png'
        elif(ice_data['season_adder'].strip()=='summer'): 
            plt.title('Summer data (June - August')
            outname = 'yearly_conc_diff_scatter_' + next_year_string + '_summer.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
####def write_toNCDF(ice_data,CERES_lw_clr_dict,CERES_sw_clr_dict,CERES_net_clr_dict,\
####                  CERES_lw_all_dict,CERES_sw_all_dict,CERES_net_all_dict):
####    #lat_ranges = np.arange(minlat,90,1.0)
####    #lon_ranges = np.arange(-180,180,1.0)
####   
####    # Create ice data nCDF file 
####    ice_nc = Dataset('./nsidc_gridded_dict.nc','w',format='NETCDF4')
####  
####    # Dimensions: For the raw data:  ntime, n25x, n25y
####    testdict = ice_data['data']
####    ntime = ice_data['data'].shape[0]
####    n25x  = ice_data['data'].shape[1]
####    n25y  = ice_data['data'].shape[2]
####    
####    dim_ntime = ice_nc.createDimension('ntime',ntime)
####    dim_n25x  = ice_nc.createDimension('n25x',n25x)
####    dim_n25y  = ice_nc.createDimension('n25y',n25y)
####
####    RDATA = ice_nc.createVariable('RDATA','f8',('ntime','n25x','n25y'))
####    RDATA.description('Raw Percent Ice Concentration Data. Converted from NSIDC files.')
####
####    MONTH=nc.createVariable('MONTH','i2',('nmth'))
####    MONTH.description='Months since October 2004'
####    LAT=nc.createVariable('Latitude','i2',('dlat','dlon'))
####    LAT.description='Latitude'
####    LAT.units='Degrees'
####    LON=nc.createVariable('Longitude','i2',('dlat','dlon'))
####    LON.description='Longitude'
####    LON.units='Degrees'
####    AI = nc.createVariable('AI','f4',('nmth','dlat','dlon'))
####    AI.description='Monthly Averaged Aerosol Index'
####    OB_COUNT=nc.createVariable('OB_COUNT','i2',('nmth','dlat','dlon'))
####    OB_COUNT.description='# of OMI AI measurements used in each monthly average'
####
####    #             For the grid data: ntime, nlat, nlon
####    # Fill in dimension variables
####
####    for i in range(num_time):
####        MONTH[i] = times[i]
####    for i in range(num_lat):
####        for j in range(num_lon):
####            LAT[i,j]=lat_ranges[i]
####            LON[i,j]=lon_ranges[j]
####
####    # Fill in actual variables
####    for i in range(num_lat):
####        print(lat_ranges[i])
####        for j in range(num_lon):
####            dictkey = (str(int(lat_ranges[i])) + 'x' + str(int(lon_ranges[j])))
####            if(dictkey not in OMI_dict):
####                # Insert missing values for AI and count
####                AI[:,i,j] = [-999.9 for m in range(num_time)]
####                OB_COUNT[:,i,j] = [-99 for m in range(num_time)]
####            else:
####                for m in range(num_time):
####                    timekey = testkeys[m]
####                    if(timekey not in OMI_dict[dictkey]):
####                        AI[m,i,j] = -999.9
####                        OB_COUNT[m,i,j] = -99
####                    else:
####                        AI[m,i,j] = OMI_dict[dictkey][timekey]['avg']
####                        OB_COUNT[m,i,j] = OMI_dict[dictkey][timekey]['#_obs']
####    nc.close()
