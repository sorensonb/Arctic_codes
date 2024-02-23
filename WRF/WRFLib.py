#!/usr/bin/env python
"""
  NAME:
    WRFLib.py

  PURPOSE:
    House functions for plotting WRF data.

  USAGE:
    >>> from WRFLib import *
    
    Then, import the WRF output file using
    >>> from netCDF4 import Dataset
    >>> data = Dataset('wrfout_d01_...')

    With the data read in, to plot the 500 mb height and wind field,
    use this. The variable 'time_index' defines which time to use from
    the file. The default is 0.
    >>> plot_500_hgt_wnd(data,0)

    To plot a model sounding for a city or town, use this syntax.
    The time index variable is the same; just use 0 to use the first
    time slot in the file. 'city_name' is the name of the city you want
    a sounding for. After a brief round of rigorous tests, basically any
    town you can think of works with this, as the geopy Python module
    pulls latitude / longitude coordinates for the city. The city 
    name must be formatted like "Grand Forks North Dakota" or 
    "Chicago Illinois". For example, to plot a sounding for Grand Forks, 
    use:
    >>> plot_WRF_sounding(data,time_index,"Grand Forks North Dakota")

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu> : 2020/02/12

"""

from wrf import getvar,interplevel,to_np,smooth2d,get_cartopy,cartopy_xlim,\
                cartopy_ylim,latlon_coords,destagger,ALL_TIMES,ll_to_xy,td
import metpy.calc as mcalc
from metpy.units import units
from metpy.plots import SkewT
from netCDF4 import Dataset
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from datetime import datetime
import numpy as np
from geopy.geocoders import Nominatim

# var: BC1, BC2, OC1, OC2
def calc_profile_total(values, action = 'sum'):

    #values = getvar(in_data, var, timeidx = timeidx)

    # Calculate the sum across the profile
    # ------------------------------------
    if(action == 'sum'):
        total_values = np.sum(values, axis = 0)
    elif(action == 'average'):
        total_values = np.mean(values, axis = 0)
    elif(action == 'max'):
        total_values = np.max(values, axis = 0)
    else:
        print("ERROR: Invalid action. Calculating sum")
        total_values = np.sum(values, axis = 0)

    return total_values

# Pass in any variable
def plot_single_data(in_data):
    ph00  = getvar(in_data,'PH',timeidx=0)
    phb00 = getvar(in_data,'PHB',timeidx=0)
    p00   = getvar(in_data,'P',timeidx=0)
    pb00  = getvar(in_data,'PB',timeidx=0)
    
    # Calculate geopotential height
    geopothgt00 = (ph00+phb00)/9.8
    ghgt00 = destagger(geopothgt00,0)
    # Calculate pressure field
    press00 = (p00+pb00)/100.

    hgt00_500 = interplevel(ghgt00,press00,500.)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Get the latitude and longitude points
    lats, lons = latlon_coords(press00)
    # Get the cartopy mapping object
    cart_proj = get_cartopy(ph00)

    fig = plt.figure(figsize=(12,9))
    ax = plt.axes(projection=cart_proj)
    states = NaturalEarthFeature(category="cultural",scale='50m',
                                 facecolor='none',
                                 name='admin_1_states_provinces_shp')
    ax.add_feature(states,linewidth=0.5,edgecolor='black')
    ax.coastlines('50m',linewidth=0.8)
    levels = np.arange(4800.,6600.,60.)
    contours00 = plt.contour(to_np(lons),to_np(lats),to_np(hgt00_500),levels=levels,
                           colors='black',
                           transform=ccrs.PlateCarree(),label='500 hPa 0-hr Forecast')
    plt.clabel(contours00,inline=1,fontsize=10,fmt="%i")
    plt.show()

#	float ACSWDNT(Time, south_north, west_east) ;
# var: 'PSFC', 'T2', 'SWDOWN', 'SWDNT'
def plot_sfc_field(in_data, var, time_index, ax = None, title = ''):

    # Extract the surface pressure and convert to mb
    psfc  = getvar(in_data, 'PSFC', timeidx = time_index)

    # Extract the desired variable
    # ----------------------------
    plt_data = getvar(in_data, var, timeidx = time_index)
    cmap = 'plasma'
    if(var == 'T2'):
        plt_data = plt_data - 273.15
    elif((var == 'BC1') | (var == 'BC2') | (var == 'OC1') | (var == 'OC2')):
        plt_data = calc_profile_total(plt_data, action = 'sum')
        cmap = 'jet'
    elif((var == 'CLDFRA')):
        plt_data = calc_profile_total(plt_data, action = 'max')
        cmap = 'jet'

    # Check the min and max values in the array
    # ------------------------------------------
    if(np.min(plt_data) < 2e-5):
        vmin = 0
    else:
        vmin = None

    if(np.max(plt_data) < 2e-5):
        vmax = 0
    else:
        vmax = None

    if(title is None):
        title = b''.join(in_data['Times'][time_index,:]).decode('utf-8') + \
            '\n' + var
 
    ## Smooth the surface field?
    ## -------------------------
    #psfc = smooth2d(psfc, 3, cenweight=4)
 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Get the latitude and longitude points
    lats, lons = latlon_coords(plt_data)
    # Get the cartopy mapping object
    cart_proj = get_cartopy(psfc)
    psfc = psfc / 100.

    in_ax = True
    if(ax is None):
        in_ax = False
        fig = plt.figure(figsize=(6,7))
        ax = fig.add_subplot(1,1,1, projection=cart_proj)

    ax.add_feature(cfeature.STATES)
    ax.coastlines('110m',linewidth=0.8)
    mesh = ax.pcolormesh(to_np(lons), to_np(lats), to_np(plt_data), \
        transform = ccrs.PlateCarree(), shading = 'auto', cmap = cmap, \
        vmin = vmin, vmax = vmax)
    cbar = plt.colorbar(mesh, ax = ax, fraction = 0.045)
    ax.set_title(title)
    if(not in_ax):
        plt.show()


def plot_500_hgt_wnd(in_data,time_index):
    ph00  = getvar(in_data,'PH',timeidx=time_index)
    phb00 = getvar(in_data,'PHB',timeidx=time_index)
    p00   = getvar(in_data,'P',timeidx=time_index)
    pb00  = getvar(in_data,'PB',timeidx=time_index)
    u00   = getvar(in_data,'U',timeidx=time_index)
    v00   = getvar(in_data,'V',timeidx=time_index)

    # Grab the forecast time
    plot_date = ''
    for xl in in_data.variables['Times'][time_index]:
        plot_date = plot_date + xl.decode('utf-8')

    # Calculate geopotential height
    geopothgt00 = (ph00+phb00)/9.8
    # destagger the height variable
    ghgt00 = destagger(geopothgt00,0)
    # destagger the u wind variable
    u00_d = destagger(u00,2)
    # destagger the u wind variable
    v00_d = destagger(v00,1)
    # Calculate pressure field
    press00 = (p00+pb00)/100.

    # Interpolate to 500 mb 
    hgt00_500 = interplevel(ghgt00,press00,500.)
    u00_500   = interplevel(u00_d,press00,500.)
    v00_500   = interplevel(v00_d,press00,500.)
   

    # Calculate wind speed from u/v
    wspd = np.sqrt(u00_500**2. + v00_500**2.)

    # Convert speeds to knots
    u00_500 = (u00_500 * 1.94384)*units('knots')
    v00_500 = (v00_500 * 1.94384)*units('knots')
    wspd = (wspd * 1.94384)*units('knots')

 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Get the latitude and longitude points
    lats, lons = latlon_coords(press00)
    # Get the cartopy mapping object
    cart_proj = get_cartopy(ph00)

    fig = plt.figure(figsize=(6,7))
    ax = plt.axes(projection=cart_proj)
    #states = NaturalEarthFeature(category="cultural",scale='50m',
    #                             facecolor='none',
    #                             name='admin_1_states_provinces_shp')
    #ax.add_feature(states,linewidth=0.5,edgecolor='black')
    ax.add_feature(cfeature.STATES)
    ax.coastlines('50m',linewidth=0.8)
    levels = np.arange(4800.,6600.,60.)
    # Plot height contours
    contours00 = plt.contour(to_np(lons),to_np(lats),to_np(hgt00_500),levels=levels,
                           colors='black',
                           transform=ccrs.PlateCarree(),label='500 hPa 0-hr Forecast')

    # Plot wind speed contours
    levels = np.arange(20,105,5)
    wspd_contours = plt.contourf(to_np(lons),to_np(lats),wspd,levels=levels,\
                                 transform = ccrs.PlateCarree(),\
                                 cmap=get_cmap('rainbow'))
    plt.colorbar(wspd_contours,ax = ax,orientation='horizontal',pad = 0.05,label=\
        'Wind speed [knots]')
    # Plot wind barbs
    skipper = 8
    barbs = plt.barbs(to_np(lons[::skipper,::skipper]),to_np(lats[::skipper,::skipper]),\
                to_np(u00_500[::skipper,::skipper]),to_np(v00_500[::skipper,::skipper]),\
                transform=ccrs.PlateCarree(),length=6)
    plt.clabel(contours00,inline=1,fontsize=10,fmt="%i")
    ax.set_title('500 hPa Geopotential Height\n'+plot_date)
    plt.show()

# NOTE: This function does not work yet.
def plot_850_hgt_tmp(in_data,time_index):
    ph00      = getvar(in_data,'PH',timeidx=time_index)
    phb00     = getvar(in_data,'PHB',timeidx=time_index)
    p00       = getvar(in_data,'P',timeidx=time_index)
    pb00      = getvar(in_data,'PB',timeidx=time_index)
    theta00   = getvar(in_data,'T',timeidx=time_index)

    # Calculate geopotential height
    geopothgt00 = (ph00+phb00)/9.8
    # destagger the height variable
    ghgt00 = destagger(geopothgt00,0)
    # Calculate pressure field
    press00 = (p00+pb00)/100.

    # Interpolate to 850 mb 
    hgt00_850     = interplevel(ghgt00,press00,850.)
    press00_850   = interplevel(press00,press00,850.)
    theta00_850   = interplevel(theta00,press00,850.)
    #u00_850   = interplevel(u00_d,press00,850.)
    #v00_850   = interplevel(v00_d,press00,850.)
  
    # Calculate temperature from potential temperature perturbation
    theta00_850 += 300 # base theta of 300 K
    print(to_np(theta00_850)[0]*units('K'))
    return(theta00_850)
    t00_850 = mcalc.temperature_from_potential_temperature(press00_850*units('hPa'),theta00_850*units('K')) 

    # Convert temp to Celsius
    t00_850 = t00_850.to('degC')
 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Get the latitude and longitude points
    lats, lons = latlon_coords(press00)
    # Get the cartopy mapping object
    cart_proj = get_cartopy(ph00)

    fig = plt.figure(figsize=(12,9))
    ax = plt.axes(projection=cart_proj)
    states = NaturalEarthFeature(category="cultural",scale='50m',
                                 facecolor='none',
                                 name='admin_1_states_provinces_shp')
    ax.add_feature(states,linewidth=0.5,edgecolor='black')
    ax.coastlines('50m',linewidth=0.8)
    levels = np.arange(1200.,2200.,60.)
    # Plot height contours
    contours00 = plt.contour(to_np(lons),to_np(lats),to_np(hgt00_850),levels=levels,
                           colors='black',
                           transform=ccrs.PlateCarree(),label='850 hPa 0-hr Forecast')

    # Plot temperature contours
    levels = np.arange(-40,40,5)
    tmp_contours = plt.contourf(to_np(lons),to_np(lats),to_np(t00_850),levels=levels,\
                                 transform = ccrs.PlateCarree(),\
                                 cmap=get_cmap('rainbow'))
    plt.colorbar(wspd_contours,ax = ax,orientation='horizontal',pad = 0.05)

    ### Plot wind barbs
    ##skipper = 8
    ##barbs = plt.barbs(to_np(lons[::skipper,::skipper]),to_np(lats[::skipper,::skipper]),\
    ##            to_np(u00_850[::skipper,::skipper]),to_np(v00_850[::skipper,::skipper]),\
    ##            transform=ccrs.PlateCarree(),length=6)
    plt.clabel(contours00,inline=1,fontsize=10,fmt="%i")
    plt.show()

def plot_WRF_combined_output(in_data, time_index, save = False, out_add = 'smoke'):

    ph00  = getvar(in_data,'PH',timeidx=time_index)
    cart_proj = get_cartopy(ph00)
    del(ph00)

    title = b''.join(in_data['Times'][time_index,:]).decode('utf-8')
    dt_date_str = datetime.strptime(title, '%Y-%m-%d_%H:%M:%S')

    fig = plt.figure(figsize = (11, 4))
    ax1 = fig.add_subplot(1,4,1, projection = cart_proj) # 2-m temp
    ax2 = fig.add_subplot(1,4,2, projection = cart_proj) # SWDOWN 
    ax3 = fig.add_subplot(1,4,3, projection = cart_proj) # accumulated aerosol profile
    ax4 = fig.add_subplot(1,4,4, projection = cart_proj) # average cloud fraction

    plot_sfc_field(in_data, 'T2', time_index, ax = ax1, title = 'T2')
    plot_sfc_field(in_data, 'SWDOWN', time_index, ax = ax2, title = 'SWDOWN')
    plot_sfc_field(in_data, 'BC1', time_index, ax = ax3, title = 'BC1')
    plot_sfc_field(in_data, 'CLDFRA', time_index, ax = ax4, title = 'CLOUDFRAC')

    plt.suptitle(title)

    fig.tight_layout()

    if(save):
        outname = dt_date_str.strftime('wrf_combined_output_%Y%m%d%H%M_' + out_add + '.png')
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_WRF_sounding(in_data,time_index,city_name):
    # Use geopy to get the city lat and lon
    geolocator = Nominatim(user_agent = 'myapplication')
    location = geolocator.geocode(city_name)
    lat = location.latitude
    lon = location.longitude

    # Grab the profile time
    plot_date = ''
    for xl in in_data.variables['Times'][time_index]:
        plot_date = plot_date + xl.decode('utf-8')

    # Extract the model data
    ph_00   = getvar(in_data,'PH',timeidx=time_index)
    phb_00  = getvar(in_data,'PHB',timeidx=time_index)
    p_00   = getvar(in_data,'P',timeidx=time_index)
    pb_00  = getvar(in_data,'PB',timeidx=time_index)
    theta_00   = getvar(in_data,'T',timeidx=time_index)
    qvapor_00   = getvar(in_data,'QVAPOR',timeidx=time_index)
    u00   = getvar(in_data,'U',timeidx=time_index)
    v00   = getvar(in_data,'V',timeidx=time_index)
    xlat  = getvar(in_data,'XLAT',timeidx=time_index)
    xlon  = getvar(in_data,'XLONG',timeidx=time_index)
    
    # Calculate geopotential height
    geopothgt_00 = (ph_00+phb_00)/9.8
    ghgt_00 = destagger(geopothgt_00,0)
    # Calculate pressure field
    press_00 = (p_00+pb_00)/100.
    theta_00 = theta_00+300.

    # destagger the u wind variable
    u00_d = destagger(u00,2)
    # destagger the u wind variable
    v00_d = destagger(v00,1)

     # Convert theta to temp
    t_00 = mcalc.temperature_from_potential_temperature(to_np(press_00)*units('hPa'),to_np(theta_00)*units('K')) 

    # Calculate dew point from press and mixing ratio
    td_00 = td(press_00,qvapor_00)

    # Convert temp to Celsius
    t_00 = t_00.to('degC')

    # Find indices closest
    aindices = ll_to_xy(in_data,lat,lon)
    a_x = aindices[1]
    a_y = aindices[0]

    # Extract profiles of data
    pres_prof = to_np(press_00[:,a_x,a_y])*units('hPa')
    t_prof    = to_np(t_00[:,a_x,a_y])
    td_prof   = to_np(td_00[:,a_x,a_y])*units('degC')
    u_prof    = to_np(u00_d[:,a_x,a_y])*units('m/s').to('knots')
    v_prof    = to_np(v00_d[:,a_x,a_y])*units('m/s').to('knots')

    # Plot the sounding
    skew = SkewT(rotation=45)
    skew.plot(pres_prof,t_prof,'r')
    skew.plot(pres_prof,td_prof,'g')
    skew.plot_barbs(pres_prof,u_prof,v_prof)
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()
    skew.ax.set_ylim(1000,100)
    skew.ax.set_xlim(-50,50)
    plt.title(city_name+'\n'+plot_date)
    plt.show()    


