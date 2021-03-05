"""
  NAME:
    CryoSat2Lib.py

  PURPOSE:
    House all the functions used for reading and working with the modis data



"""
import numpy as np
import numpy.ma as ma
import sys
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from matplotlib.lines import Line2D
from matplotlib.dates import DateFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
import subprocess
from scipy import stats

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Set up global variables
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
mapcrs = ccrs.NorthPolarStereo(central_longitude = 45.)
datacrs = ccrs.PlateCarree()

zoom_dict = {
    'Finland': [10,55,65,80]
}

proj_dict = {
    'Finland': ccrs.NorthPolarStereo(central_longitude = 35.) 
}

# Start_date and end_date must be formatted as 
# "YYYYMMDD"
def read_MODIS(start_date,latmin = 60):

    ##if(monthly == True):
    ##    dformat = '%Y%m' 
    ##    end_string_idx = -5
    ##else:
    ##    dformat = '%Y%m%d' 
    ##    end_string_idx = -3
    dformat = '%Y%m%d' 
    # Set up starting and ending datetime objects
    sdate = datetime.strptime(start_date,dformat)
    ##edate = datetime.strptime(end_date,dformat)
 
    
    # Grab all the modis files
    #cmnd = "ls /home/bsorenson/data/MODIS/CHL/A*.nc"
    cmnd = "ls /home/bsorenson/data/MODIS/CHL/A*.nc"
    
    
    file_initial = subprocess.check_output(cmnd,shell=True).decode('utf-8').strip().split('\n')
    file_names = []

    for fname in file_initial:
        fdate = datetime.strptime(fname[-34:-30],'%Y') + relativedelta(days = int(fname[-30:-27])-1)
        if(fdate == sdate):
            file_names.append(fname)
   
    lat_ranges = np.arange(latmin,90.25,0.25)
    lon_ranges = np.arange(-180.,180.25,0.25)
 
    # Read in the latitude, longitude, and area data
    num_files = len(file_names)
    modis_data = {}
    modis_data['chlor_a'] = np.full((num_files,len(lat_ranges),len(lon_ranges)),-9.)
    modis_data['lat']  = lat_ranges
    modis_data['lon']  = lon_ranges
    modis_data['titles']  = []
    
    count = 0
    for fname in file_names:
    
        total_name = fname
        #total_name = data_loc+fname
        print(total_name)
        data = Dataset(total_name,'r')
        
        for xi in range(len(lat_ranges)-1):       
            print(lat_ranges[xi])
            for yj in range(len(lon_ranges)-1):
                lat_indices = np.where((data.variables['lat'][:] >= lat_ranges[xi]) & (data.variables['lat'][:] < lat_ranges[xi+1]))
                lon_indices = np.where((data.variables['lon'][:] >= lon_ranges[yj]) & (data.variables['lon'][:] < lon_ranges[yj+1]))
                modis_data['chlor_a'][count,xi,yj] =\
                    np.average(data.variables['chlor_a'][lat_indices[0],lon_indices[0]])
        data.close() 
        modis_data['titles'].append(fname)
        #modis_data['dates'].append(fname[-11:end_string_idx])
        count+=1

    # Convert the longitude values from 0 - 360 to -180 - 180
    return modis_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Trend calculating functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def trend_calc(modis_data,x_ind,y_ind,variable,thielsen=False):
    temp_data = modis_data[variable][:,x_ind,y_ind]
    temp_data[temp_data < -999.] = np.nan
    #temp_data = np.ma.masked_where(temp_data < -999., temp_data)
    # Don't calculate trends if there are less than 2 valid values
    if(np.count_nonzero(~np.isnan(temp_data)) < 2):
        trend = pcnt_change = np.nan
    else:
        avg_modis = np.nanmean(temp_data)
        interpx = np.arange(len(temp_data))
        #interpx = years
        ##interper = np.poly1d(np.polyfit(interpx,temp_data,1)) 
        ### Normalize trend by dividing by number of years
        ##trend = (interper(interpx[-1])-interper(interpx[0]))

        slope, intercept, r_value, p_value, std_err = stats.linregress(interpx,temp_data)
        ##print(slope/len(test_dict[dictkey].keys())
        regress_y = interpx*slope+intercept
        trend = regress_y[-1] - regress_y[0]

        if(thielsen==True):
            #The slope
            S=0
            sm=0
            nx = len(temp_data)
            num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
            Sn=np.zeros(num_d)
            for si in range(0,nx-1):
                for sj in range(si+1,nx):
                    # Find the slope between the two points
                    Sn[sm] = (temp_data[si]-temp_data[sj])/(si-sj) 
                    sm=sm+1
                # Endfor
            # Endfor
            Snsorted=sorted(Sn)
            sm=int(num_d/2.)
            if(2*sm    == num_d):
                tsslope=0.5*(Snsorted[sm]+Snsorted[sm+1])
            if(2*sm+1 == num_d): 
                tsslope=Snsorted[sm+1]
            regress_ts = interpx*tsslope+intercept
            trend = regress_ts[-1]-regress_ts[0]


        pcnt_change = (trend/avg_modis)*100.
    # Find the percent change per decade
    return trend,pcnt_change


# modis_trendCalc calculates the trends over the time period at each
# grid point on the 25x25 km grid.
def modis_trendCalc(modis_data,thielSen=False):
    # Loop over the data and calculate trends
    for i in range(448):
        print(i)
        max_trend = -99.
        min_trend = 99.
        for j in range(304):
            modis_data['thick_trends'][i,j],temp_pcnt = \
                trend_calc(modis_data,i,j,'ice_thick',thielsen=thielSen)
            #modis_data['thick_land_trends'][i,j],temp_pcnt = \
            #    trend_calc(modis_data,i,j,'ice_thick',thielsen=thielSen)
            modis_data['con_trends'][i,j],temp_pcnt = \
                trend_calc(modis_data,i,j,'ice_con',thielsen=thielSen)
            #modis_data['con_land_trends'][i,j],temp_pcnt = \
            #    trend_calc(modis_data,i,j,'ice_con',thielsen=thielSen)
            ##temp_trend = modis_data['trends'][i,j]
            ##if(temp_trend>max_trend):
            ##    max_trend = temp_trend
            ##if(temp_trend<min_trend):
            ##    min_trend = temp_trend
        ##print("  max trend = ",max_trend)
        ##print("  min trend = ",min_trend)
        # Deal with land masks
        #good_indices = np.where(modis_data['data'][0,i,:]<251.)
        #land_indices = np.where(modis_data['data'][0,i,:]>=251.)
        #modis_data['trends'][i,land_indices] = np.nan
        #modis_data['land_trends'][i,good_indices] = np.nan

    return modis_data

# Calculate trends on the 1x1 degree lat/lon grid
def modis_gridtrendCalc(modis_data,area=True,thielSen=False):
    if(area==True):
        print("\nCalculating area trends\n")
    else:
        print("\nCalculating % concentration trends\n")
    modis_data['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = modis_data['grid_lat']
    lon_ranges = modis_data['grid_lon']
    ## Create array to hold monthly averages over region
    #initial_avgs = np.full(len(CERES_dict['dates']),-9999.)
    #initial_years  = np.zeros(len(initial_avgs))
    #initial_months = np.zeros(len(initial_avgs))
   
    # Loop over the lat and lon dimensions
    for xi in range(len(lat_ranges)):
        max_trend = -999
        min_trend = 999
        for yj in range(len(lon_ranges)):
            # Calculate the trend at the current box
            if(area==True):
                interp_data = (modis_data['grid_modis_conc'][:,xi,yj]/100.)*modis_data['grid_total_area'][xi,yj]
                good_indices = np.where(np.isnan(interp_data)==False)
            else:
                interp_data = modis_data['grid_modis_conc'][:,xi,yj]
                good_indices = np.where(interp_data!=-99.)
            if(len(good_indices[0])==0):
                total_trend = -9999.
            else:
                interpx = np.arange(len(interp_data))
                #print(CERES_dict['data'][:,xi,yj][good_indices])
                total_interper = np.poly1d(np.polyfit( \
                    interpx[good_indices],\
                    interp_data[good_indices],1)) 
                # Normalize trend by dividing by number of years
                total_trend = (total_interper(interpx[-1])-\
                               total_interper(interpx[0]))
                if(total_trend>max_trend):
                    max_trend = total_trend
                if(total_trend<min_trend):
                    min_trend = total_trend
            if(area==True):
                modis_data['grid_modis_area_trend'][xi,yj] = total_trend
        print(xi)
        #print("max trend = ",max_trend)
        #print("min trend = ",min_trend)

    return modis_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Gridding functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# grid_data_conc grids the 25x25 km gridded modis concentration data into
# a 1x1 degree lat/lon grid
def grid_data_values(modis_data):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_con    = np.full((len(modis_data['ice_con'][:,0,0]),len(lat_ranges),len(lon_ranges)),-99.)
    grid_thick  = np.full((len(modis_data['ice_thick'][:,0,0]),len(lat_ranges),len(lon_ranges)),-99.)
    grid_con_cc = np.full((len(modis_data['ice_con'][:,0,0]),len(lat_ranges),len(lon_ranges)),-999.)
    grid_thick_cc = np.full((len(modis_data['ice_thick'][:,0,0]),len(lat_ranges),len(lon_ranges)),-999.)
    print("Size of grid array: ",grid_con.shape)
    for nt in range(len(modis_data['ice_con'][:,0,0])):
        print(nt)
        for xi in range(448):
            # Don't grid the data if any portion of the lat/lon box is over land.
            # Don't include land data
            for yj in range(304):
                lat_index = np.where(np.floor(modis_data['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(modis_data['lon'][xi,yj])>=lon_ranges)[-1][-1]
                # Add the current pixel area into the correct grid box, no
                # matter if the current box is missing or not.
                #if((lat_index==20) & (lon_index==10)):
                #    print("Current grid area = ",grid_modis_area[lat_index,lon_index])
                #if((lat_index==20) & (lon_index==10)):
                #    print("    New grid area = ",grid_modis_area[lat_index,lon_index])
                if(modis_data['ice_con'][nt,xi,yj] != -9999.0):
                    #if(nt==0): grid_modis_area[lat_index,lon_index] += modis_data['area'][xi,yj]
                    if(grid_con_cc[nt,lat_index,lon_index]==-999.):
                        grid_con[nt,lat_index,lon_index] = modis_data['ice_con'][nt,xi,yj]
                        grid_con_cc[nt,lat_index,lon_index] = 1.
                    else:
                        grid_con[nt,lat_index,lon_index] = ((grid_con[nt,lat_index,\
                            lon_index]*grid_con_cc[nt,lat_index,lon_index])+\
                            modis_data['ice_con'][nt,xi,yj])/(grid_con_cc[nt,lat_index,\
                            lon_index]+1.)
                        grid_con_cc[nt,lat_index,lon_index]+=1
                if(modis_data['ice_thick'][nt,xi,yj] != -9999.0):
                    #if(nt==0): grid_modis_area[lat_index,lon_index] += modis_data['area'][xi,yj]
                    if(grid_thick_cc[nt,lat_index,lon_index]==-999.):
                        grid_thick[nt,lat_index,lon_index] = modis_data['ice_thick'][nt,xi,yj]
                        grid_thick_cc[nt,lat_index,lon_index] = 1.
                    else:
                        grid_thick[nt,lat_index,lon_index] = ((grid_thick[nt,lat_index,\
                            lon_index]*grid_thick_cc[nt,lat_index,lon_index])+\
                            modis_data['ice_thick'][nt,xi,yj])/(grid_thick_cc[nt,lat_index,\
                            lon_index]+1.)
                        grid_thick_cc[nt,lat_index,lon_index]+=1
                #else:
                #    if(nt==0): grid_modis_area[lat_index,lon_index] = np.nan

                    # end else
                # end if good modis check
            # end y grid loop
        # end x grid loop
    # end time loop 
    modis_data['grid_con']   = grid_con
    modis_data['grid_thick'] = grid_thick
    modis_data['grid_con_cc'] = grid_con_cc
    modis_data['grid_thick_cc'] = grid_thick_cc
    modis_data['grid_lat'] = lat_ranges
    modis_data['grid_lon'] = lon_ranges
    
    return modis_data

# grid_data averages the 25x25 km trends into a 1x1 degree grid
# The 'grid_modis' paramter in the dictionary therefore contains
# the gridded trends.
def grid_data_trends(modis_dict):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_modis    = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    grid_modis_cc = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    for xi in range(448):
        print(xi)
        for yj in range(304):
            if(not np.isnan(modis_dict['trends'][xi,yj])):
                lat_index = np.where(np.floor(modis_dict['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(modis_dict['lon'][xi,yj])>=lon_ranges)[-1][-1]
                if(grid_modis_cc[lat_index,lon_index]==-999.):
                    grid_modis[lat_index,lon_index] = modis_dict['trends'][xi,yj]
                    grid_modis_cc[lat_index,lon_index] = 1.
                else:
                    grid_modis[lat_index,lon_index] = ((grid_modis[lat_index,lon_index]*grid_modis_cc[lat_index,lon_index])+\
                                                     modis_dict['trends'][xi,yj])/(grid_modis_cc[lat_index,lon_index]+1.)
                    grid_modis_cc[lat_index,lon_index]+=1
    modis_dict['grid_modis'] = grid_modis
    modis_dict['grid_modis_cc'] = grid_modis_cc
    modis_dict['grid_lat'] = lat_ranges
    modis_dict['grid_lon'] = lon_ranges
    return modis_dict

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def plot_modis_data(modis_data,minlat=60,tind=0,zoom = None,save=False):

    data = modis_data['chlor_a'][tind,:,:]
    colormap = plt.cm.jet
    mask_data = np.ma.masked_where(data == -9., data)
    mask_data = np.ma.masked_invalid(mask_data)

    plot_date = datetime.strptime(modis_data['titles'][0][-34:-30],'%Y') + \
        relativedelta(days = int(modis_data['titles'][0][-30:-27])-1)

    if(minlat > 50):
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0)
    else:
        mapcrs = ccrs.Robinson()

    plot_lat, plot_lon = np.meshgrid(modis_data['lat'],modis_data['lon'])

    plt.close()
    fig1 = plt.figure()
    if(zoom == None):
        ax = plt.axes(projection = mapcrs)
        ax.set_extent([-180,180,minlat,90],datacrs)
        saver = ''
    else:
        ax = plt.axes(projection = proj_dict[zoom])
        ax.set_extent(zoom_dict[zoom],datacrs)
        saver = '_'+zoom 
    ax.gridlines()
    ax.coastlines()
    mesh = ax.pcolormesh(plot_lon,plot_lat,mask_data.T,\
            transform = datacrs, norm = cm.LogNorm(vmin = 0.01, vmax = 20.),\
            cmap = colormap)
    #CS = ax.contour(longitude,latitude,smooth_thick,[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],transform = datacrs)
    
    # Adjust and make it look good
    ax.add_feature(cfeature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
    ax.set_title('MODIS Chlorophyll-α')
    cbar = plt.colorbar(mesh,ticks = [0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0],orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905, label='Chlorophyll Concentration, OCI Algorithm (mg m$^{-3}$)')
    cbar.ax.set_xticklabels(['0.01','0.02','0.05','0.1','0.2','0.5','1','2','5','10','20'])
    #ax.set_xlim(-4170748.535086173,4167222.438879491)
    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
    #ax.set_title(datetime.strftime(base_dtm,'%B %Y') + ' CryoSat-2 Data')

    if(save == True):
        outname = 'modis_chlor_a_' + plot_date.strftime('%Y%m%d') + saver + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Plot a histogram of 25 x 25 km grid data
def plot_modis_hist(modis_data,tind,variable,bins = 100,save = False):
    data = modis_data[variable][tind,:,:]
    mask_data = np.ma.masked_where(data < -999., data)

    # Make a datetime object from the structure date
    # ----------------------------------------------
    if(len(modis_data['dates'][tind]) == 6):
        str_fmt = '%Y%m'
    elif(len(modis_data['dates'][tind]) == 8):
        str_fmt = '%Y%m%d'

    base_dtm = datetime.strptime(modis_data['dates'][tind],str_fmt)

    # Set up the x axis label
    # -----------------------
    if(variable == 'ice_thick'):
        xlabel = 'Derived Sea Ice Thickness [m]' 
    elif(variable == 'ice_con'):
        xlabel = 'Sea Ice Concentration [%]' 

    bin_heights,bin_borders = np.histogram(mask_data.compressed(),bins=bins)
    bin_widths = np.diff(bin_borders)
    bin_centers = bin_borders[:-1] + bin_widths / 2

    t_init = models.Gaussian1D()
    fit_t = fitting.LevMarLSQFitter()
    t = fit_t(t_init, bin_centers, bin_heights)

    print('Amplitude: ',np.round(t.amplitude.value,3))
    print('Mean:      ',np.round(t.mean.value,3))
    print('StDev:     ',np.round(t.stddev.value,3))

    x_interval_for_fit = np.linspace(bin_centers[0],bin_centers[-1],100)
    plt.close()
    plt.figure()
    plt.bar(bin_centers,bin_heights,width=bin_widths,label='histogram')
    plt.plot(x_interval_for_fit,t(x_interval_for_fit),label='fit',c='tab:red')
    plt.xlabel(xlabel)
    plt.ylabel('Counts')
    plt.title(datetime.strftime(base_dtm,'%B %Y') + ' CryoSat-2 Data')
    plt.legend()
    #plt.close()
    #plt.hist(mask_data.compressed(),bins=bins)
    #plt.title(modis_data['dates'][tind] + ' '+variable)
    #plt.xlabel(xlabel)
    #plt.ylabel('Counts')
    #plt.title(datetime.strftime(base_dtm,'%B %Y') + ' CryoSat-2 Data')

    if(save == True):
        outname = 'modissat2_' + variable + '_' + datetime.strftime(base_dtm,'%Y%m%d') + '_hist.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# plot_grid_data generates a plot of the /
def plot_grid_data(modis_data,t_ind,pvar,adjusted=False,save=False):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)

    if(pvar=='grid_con'):
        plabel = "Percent Ice Concentration"
    elif(pvar=='grid_thick'):
        plabel = "Sea Ice Thickness"

    local_grid_modis = np.copy(modis_dict['grid_modis_conc'][t_ind,:,:])
    local_grid_modis_bad = np.copy(modis_dict['grid_modis_conc'][t_ind,:,:])

    local_grid_modis[local_grid_modis==-999.] = np.nan
    local_grid_modis_bad[local_grid_modis!=-999.] = np.nan
    plot_good_data = ma.masked_invalid(local_grid_modis)
    plot_land_data = ma.masked_invalid(local_grid_modis_bad)

    colormap = plt.cm.ocean
    #coolwarm = plt.cm.coolwarm
    #ax = plt.axes(projection=ccrs.NorthPolarStereo())
    file_adder=''
    if(adjusted==True):
        file_adder='_adjusted'
        fig1 = plt.figure(figsize=(8,5))
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=45.))
    else:
        fig1 = plt.figure()
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
    ax.gridlines()
    mesh = plt.pcolormesh(lon_ranges,lat_ranges,plot_good_data,\
            transform=ccrs.PlateCarree(),vmin=0,vmax=100,cmap=colormap)
    if(adjusted==True):
        ax.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        axins = inset_axes(ax,width="5%",height="100%",loc='lower left',\
                    bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax.transAxes,\
                    borderpad=0)
        cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=plabel)
        ax.set_xlim(-5170748.535086173,5167222.438879491)
        ax.set_ylim(-3913488.8763307533,3943353.899053069)
        #ender = '_adjusted'+ender
    else:
        plt.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
        cbar = plt.colorbar(mesh,cmap=colormap,label=plabel)
    ax.coastlines()


    #fig1 = plt.figure(figsize=(9,5))
    #plt.pcolormesh(plot_good_data,cmap=plt.cm.bwr,vmin=-50,vmax=50)
    #plt.colorbar(label='Percent Ice Concentration Trend (%)')
    #plt.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
    #plt.gca().invert_xaxis()
    #plt.gca().invert_yaxis()
    ax.set_title('Gridded NSIDC Sea Ice Concentration\n'+modis_dict['titles'][t_ind])
    if(save==True):
        outname = 'modis_conc_gridded_200012_201812'+modis_dict['file_adder']+file_adder+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
   
# Plot the trends for the 25x25 km grid data 
def plot_trend(modis_data,variable):
    data = modis_data[variable][:,:]
    colormap = plt.cm.bwr
    mask_data = np.ma.masked_where(data == np.nan, data)

    plt.close()
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180,180,60,90])
    mesh = ax.pcolormesh(modis_data['lon'],modis_data['lat'],mask_data,\
            transform = datacrs, cmap = colormap,vmin=min_dict[variable],\
            vmax=max_dict[variable])
    #CS = ax.contour(modis_data['lon'],modis_data['lat'],np.ma.masked_where(modis_data['thick_trends'][:,:] == np.nan,modis_data['thick_trends'][:,:]),\
    #        np.linspace(-0.5,0.5,5),transform = datacrs)
    
    # Adjust and make it look good
    ax.add_feature(cfeature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
    cbar = plt.colorbar(mesh,ticks = tick_dict[variable],orientation='horizontal',pad=0,aspect=50,label=variable)
    cbar.ax.set_xticklabels(tick_label_dict[variable])
    ax.set_xlim(-4170748.535086173,4167222.438879491)
    ax.set_ylim(-2913488.8763307533,2943353.899053069)
    ax.set_title(variable)
    plt.show()

def plot_grid_trend(modis_dict,adjusted=False,save=False):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)

    local_grid_modis = np.copy(modis_dict['grid_modis'])
    local_grid_modis_bad = np.copy(modis_dict['grid_modis'])

    local_grid_modis[local_grid_modis==-999.] = np.nan
    local_grid_modis_bad[local_grid_modis!=-999.] = np.nan
    plot_good_data = ma.masked_invalid(local_grid_modis)
    plot_land_data = ma.masked_invalid(local_grid_modis_bad)

    colormap = plt.cm.bwr
    #coolwarm = plt.cm.coolwarm
    #ax = plt.axes(projection=ccrs.NorthPolarStereo())
    file_adder=''
    if(adjusted==True):
        file_adder='_adjusted'
        fig1 = plt.figure(figsize=(8,5))
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=45.))
    else:
        fig1 = plt.figure()
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
    ax.gridlines()
    mesh = plt.pcolormesh(lon_ranges,lat_ranges,plot_good_data,\
            transform=ccrs.PlateCarree(),vmin=-50,vmax=50,cmap=colormap)
    if(adjusted==True):
        ax.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        axins = inset_axes(ax,width="5%",height="100%",loc='lower left',\
                    bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax.transAxes,\
                    borderpad=0)
        cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label='Ice Concentration Trend [%]')
        ax.set_xlim(-5170748.535086173,5167222.438879491)
        ax.set_ylim(-3913488.8763307533,3943353.899053069)
        #ender = '_adjusted'+ender
    else:
        plt.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
        cbar = plt.colorbar(mesh,cmap=colormap,label='Ice Concentration Trend [%]')
    ax.coastlines()


    #fig1 = plt.figure(figsize=(9,5))
    #plt.pcolormesh(plot_good_data,cmap=plt.cm.bwr,vmin=-50,vmax=50)
    #plt.colorbar(label='Percent Ice Concentration Trend (%)')
    #plt.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
    #plt.gca().invert_xaxis()
    #plt.gca().invert_yaxis()
    if(modis_dict['season_adder']!=''):
        ax.set_title('NSIDC Sea Ice Concentration'+modis_dict['season_adder'].title()+\
            ' Trends\nJan 2001 to Dec 2018')
    else:
        ax.set_title('NSIDC Sea Ice Concentration Trends\nJan 2001 to Dec 2018')
    if(save==True):
        outname = 'modis_trend_gridded_200101_201812'+modis_dict['file_adder']+file_adder+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
    
def plot_grid_time_series(modis_dict,lat_ind,lon_ind,thielsen=False):
    inseason = modis_dict['season_adder'].strip()
    temp_data = modis_dict['grid_modis_conc'][:,lat_ind,lon_ind]
    interpx = np.arange(len(temp_data))
    #interpx = years
    ##interper = np.poly1d(np.polyfit(interpx,temp_data,1)) 
    ### Normalize trend by dividing by number of years
    ##trend = (interper(interpx[-1])-interper(interpx[0]))

    slope, intercept, r_value, p_value, std_err = stats.linregress(interpx,temp_data)
    ##print(slope/len(test_dict[dictkey].keys())
    regress_y = interpx*slope+intercept
    trend = regress_y[-1] - regress_y[0]

    if(thielsen==True):
        #The slope
        S=0
        sm=0
        nx = len(temp_data)
        num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
        Sn=np.zeros(num_d)
        for si in range(0,nx-1):
            for sj in range(si+1,nx):
                # Find the slope between the two points
                Sn[sm] = (temp_data[si]-temp_data[sj])/(si-sj) 
                sm=sm+1
            # Endfor
        # Endfor
        Snsorted=sorted(Sn)
        sm=int(num_d/2.)
        if(2*sm    == num_d):
            tsslope=0.5*(Snsorted[sm]+Snsorted[sm+1])
        if(2*sm+1 == num_d): 
            tsslope=Snsorted[sm+1]
        regress_ts = interpx*tsslope+intercept
        trend = regress_ts[-1]-regress_ts[0]

    fig1 = plt.figure() 
    plt.title('Gridded Ice Data: '+str(modis_dict['grid_lat'][lat_ind])+'x'+\
                str(modis_dict['grid_lon'][lon_ind])+'\n'+inseason.title()+\
                ' season of each year)')
    plt.plot(temp_data,label='observations') 
    plt.plot(regress_y,'--',label='trend')
    plt.ylabel('Ice Concentration [%]')
    plt.legend()
    outname = 'nsidc_grid_time_series_'+inseason+'_'+str(int(modis_dict['grid_lat'][lat_ind]))+'x'+\
                str(int(modis_dict['grid_lon'][lon_ind]))+'.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)   
    plt.show()

# tind indicates the reference time to check concentrations and thicknesses before
# plotting time series 
# Syntax to run:
# >>> import sys
# >>> sys.path.append('/home/bsorenson/Research/Ice_analysis')
# >>> from CryoSat2Lib import *
# >>> from IceLib import *
# >>> ice_data = read_ice('all')
# >>> modis_data = read_modis('all','201011','201911')
# >>> albedo_effect_test(modis_data,ice_data,11,1.8)
# This makes the figure "modissat2_melt_time_20120301.png"
def albedo_effect_test(modis_data,ice_data,tind,start_thick,save=False):
#def albedo_effect_test(modis_data,tind,start_thick):

    # Find the time index in the NSIDC ice structure that matches
    # with the time index in the CryoSat2 structure
    ice_ind = 0
    for xi in range(len(ice_data['titles'])):
        if(ice_data['titles'][xi][7:13] == modis_data['dates'][tind]):
            ice_ind = xi

    low_thick   = start_thick - 0.01
    high_thick = start_thick + 0.01

    # Identify regions with the thickness desired by 'start_thick'
    test_con = ice_data['data'][ice_ind,:,:][(modis_data['ice_thick'][tind,:,:] >= low_thick) & \
                    (modis_data['ice_thick'][tind,:,:] < high_thick)]
    test_lats = ice_data['lat'][:,:][(modis_data['ice_thick'][tind,:,:] >= low_thick) & \
                    (modis_data['ice_thick'][tind,:,:] < high_thick)]
    ##test_con = modis_data['ice_con'][tind,:,:][(modis_data['ice_thick'][tind,:,:] >= low_thick) & \
    ##                (modis_data['ice_thick'][tind,:,:] < high_thick)]
    ice_cons = np.zeros((12,test_con.size))

    colors = plt.cm.plasma(np.linspace(0,1,test_con.size))

    # Extract time series for the year (assuming October is the reference
    # month, the next 7 indices
    # -------------------------------------------------------------------
    count = 0
    for ti in range(ice_ind,ice_ind+12):
    #for ti in range(tind,tind+12):
        ice_cons[count,:] = ice_data['data'][ti,:,:][(modis_data['ice_thick'][tind,:,:] >= low_thick) & \
                    (modis_data['ice_thick'][tind,:,:] < high_thick)] 
        #ice_cons[count,:] = modis_data['ice_con'][ti,:,:][(modis_data['ice_thick'][tind,:,:] >= low_thick) & \
        #            (modis_data['ice_thick'][tind,:,:] < high_thick)] 
        count+=1

    # Go through and calculate the number of months to melt 50% of the
    # beginning ice concentration. For each grid point, plot starting
    # ice concentration on the x axis and the number of months to melt
    # on the y axis
    # ----------------------------------------------------------------

    # Set up colorbar to color the lines by latitude
    # ----------------------------------------------
    norm = mpl.colors.Normalize(vmin = test_lats.min(), vmax = test_lats.max())
    cmap = mpl.cm.ScalarMappable(norm = norm, cmap = mpl.cm.plasma)
    cmap.set_array([])
    c = np.arange(int(test_lats.min()),int(test_lats.max()))

    # Set up a base datetime object for the x axis
    base_dtm = datetime.strptime(modis_data['dates'][tind],'%Y%m')
    dates = [base_dtm + timedelta(days = int(xval * 31)) for xval in np.arange(12)]

    # Plot the total time series data. Also, calculate
    # the difference in ice concentration over the next four months
    # -------------------------------------------------------------
    diff = np.zeros(ice_cons.shape[1])
    melt_months = np.zeros(ice_cons.shape[1])
    melt_limit = 0  # when zero, means complete melting

    fig, axs = plt.subplots(2,dpi=100)
    fig.set_size_inches(9,5.5)
    for xi in range(ice_cons.shape[1]):
        # Plot the time series for the current grid box
        axs[0].plot(dates,ice_cons[:,xi],c = cmap.to_rgba(test_lats[xi]))

        # Calculate the difference
        diff[xi] = ice_cons[0,xi] - ice_cons[4,xi] 

        # Determine how long it takes 
        try:
            melt_months[xi] = np.argwhere(ice_cons[:,xi] <= melt_limit)[0][0]
        except:
            melt_months[xi] = -9 

    axs[0].xaxis.set_major_formatter(DateFormatter('%b %Y'))    
    
    #fig.colorbar(cmap, ticks = c)
    axs[0].set_ylabel('Ice concentration [%]')
    axs[0].set_title(datetime.strftime(base_dtm,'%B %Y') + '\n' + str(start_thick) + ' m Thickness ')

    # Plot the differences
    # ---------------------
    #fig2,ax2 = plt.subplots(dpi=100)
    p_scat = axs[1].scatter(ice_cons[0,:][melt_months > -9],melt_months[melt_months > -9],cmap = mpl.cm.plasma,c = cmap.to_rgba(test_lats[melt_months > -9]))
    axs[1].set_ylabel('Months needed to melt\n below ' + str(melt_limit) + '% ice concentration')
    #axs[1].set_ylim(bottom = 0)
    #p_scat = axs[1].scatter(ice_cons[0,:],diff[:],cmap = mpl.cm.plasma,s = 12,c = cmap.to_rgba(test_lats[:]))
    #axs[1].set_ylabel('Ice concentration decrease \nfrom '+datetime.strftime(base_dtm,'%B') + \
    #    ' to ' + datetime.strftime(base_dtm + timedelta(days = int(4 * 31)),'%B') + ' [%]')
    cbar_ax = fig.add_axes([0.92,0.11,0.02,0.77])
    fig.colorbar(cmap,cax = cbar_ax,ticks = c[::2],label = 'Latitude [$^{o}$]')
    #fig2.colorbar(cmap, ticks = c)
    axs[1].set_xlabel(datetime.strftime(base_dtm,'%B %Y') + ' Concentration')

    if(save == True):
        outname = 'modissat2_melt_time_' + datetime.strftime(base_dtm,'%Y%m%d') + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()
    return


    # Plot 
    #mask_data = np.ma.masked_where(data < -999., data)

    if(variable[:4] == 'grid'):
        lat_vals = modis_data['grid_lat']
        lon_vals = modis_data['grid_lon']
    else:
        lat_vals = modis_data['lat']
        lon_vals = modis_data['lon']

    plt.close()
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180,180,60,90])
    mesh = ax.pcolormesh(lon_vals,lat_vals,mask_data,\
            transform = datacrs, cmap = colormap,vmin=min_dict[variable],\
            vmax=max_dict[variable])
    #CS = ax.contour(longitude,latitude,smooth_thick,[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],transform = datacrs)
    
    # Adjust and make it look good
    ax.add_feature(cfeature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
    cbar = plt.colorbar(mesh,ticks = tick_dict[variable],orientation='horizontal',pad=0,aspect=50,label=variable)
    cbar.ax.set_xticklabels(tick_label_dict[variable])
    ax.set_xlim(-4170748.535086173,4167222.438879491)
    ax.set_ylim(-2913488.8763307533,2943353.899053069)
    ax.set_title(modis_data['titles'][tind])
    plt.show()



##!## This code replicates Jianglong's figure showing how the April - September
##!## modis concentration values have changed over the time period
##!## NOTE: To run, do something like
##!## >>> all_modis_data = read_modis('all')
##!## >>> all_modis_data = grid_data_conc(all_modis_data)
##!## and then
##!## >>> plot_apr_sep_changes(all_modis_data)
##!#def plot_apr_sep_changes(modis_dict):
##!#  
##!#    inseason = modis_dict['season_adder'].strip()
##!# 
##!#    upper_vals = np.array([100,80,60,40])
##!#    lower_vals = np.array([80,60,40,20])
##!#
##!#    # Generate 3-year averages for 3 time periods:
##!#    # 2001-2003
##!#    # 2008-2011
##!#    # 2016-2018
##!#   
##!#    # For now, look at April - September 
##!#    if(inseason=='sunlight'):
##!#        num_months = 6
##!#        start_idx  = 0
##!#    else:
##!#        num_months = 12
##!#        start_idx  = 4
##!#    # Dimensions of modis_avgs are
##!#    # 0 - year blocks (size = 3, see above)
##!#    # 1 - months      (size = num_months)
##!#    # 2 - modis values  (size = 4)
##!#    modis_avgs = np.zeros((3,num_months,len(upper_vals)))
##!#    indices = [0,8,15]
##!#  
##!#    # Base the locations for each modis concentration range on the average
##!#    # April concentration between 2001 and 2003
##!#
##!#    avg_Apr_old = \
##!#        np.average(modis_dict['grid_modis_conc'][start_idx::num_months,:,:][:3,:,:],axis=0)
##!# 
##!#    # Use the locations during the first April average to base everything
##!#    # on
##!#    locations_80_100 = np.where((avg_Apr_old <= 100.) & (avg_Apr_old > 80.)) 
##!#    locations_60_80  = np.where((avg_Apr_old <= 80.)  & (avg_Apr_old > 60.)) 
##!#    locations_40_60  = np.where((avg_Apr_old <= 60.)  & (avg_Apr_old > 40.)) 
##!#    locations_20_40  = np.where((avg_Apr_old <= 40.)  & (avg_Apr_old > 20.)) 
##!#
##!#    # Fill the modis_avgs array with the data for each time period 
##!#    for ri in range(len(indices)):
##!#        for mi in range(num_months):
##!#            print("Year block = ",ri,"  Month = ",mi)
##!#            # Deal with 80 to 100
##!#            ##temp_arr = np.zeros((3,len(locations_80_100)))
##!#            ##temp_arr[0] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri]][locations_80_100]
##!#            ##temp_arr[1] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri] + 1][locations_80_100]
##!#            ##temp_arr[2] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri] + 2][locations_80_100]
##!#            #modis_avgs[ri,mi,0] = np.average(temp_arr)
##!#            modis_avgs[ri,mi,0] = \
##!#                np.average(np.average(modis_dict['grid_modis_conc'][mi::num_months,:,:]
##!#                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_80_100])
##!#
##!#            ### Deal with 60 to 80
##!#            ##temp_arr = np.zeros((3,len(locations_60_80)))
##!#            ##temp_arr[0] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri]][locations_60_80]
##!#            ##temp_arr[1] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri] + 1][locations_60_80]
##!#            ##temp_arr[2] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri] + 2][locations_60_80]
##!#            ##modis_avgs[ri,mi,1] = np.average(temp_arr)
##!#            modis_avgs[ri,mi,1] = \
##!#                np.average(np.average(modis_dict['grid_modis_conc'][mi::num_months,:,:]
##!#                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_60_80])
##!#
##!#            # Deal with 40 to 60
##!#            ##temp_arr = np.zeros((3,len(locations_40_60)))
##!#            ##temp_arr[0] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri]][locations_40_60]
##!#            ##temp_arr[1] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri] + 1][locations_40_60]
##!#            ##temp_arr[2] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri] + 2][locations_40_60]
##!#            #modis_avgs[ri,mi,2] = np.average(temp_arr)
##!#            modis_avgs[ri,mi,2] = \
##!#                np.average(np.average(modis_dict['grid_modis_conc'][mi::num_months,:,:]
##!#                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_40_60])
##!#
##!#            # Deal with 40 to 60
##!#            ##temp_arr = np.zeros((3,len(locations_20_40)))
##!#            ##temp_arr[0] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri]][locations_20_40]
##!#            ##temp_arr[1] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri] + 1][locations_20_40]
##!#            ##temp_arr[2] = modis_dict['grid_modis_conc'][mi::num_months][indices[ri] + 2][locations_20_40]
##!#            ##modis_avgs[ri,mi,3] = np.average(temp_arr)
##!#            modis_avgs[ri,mi,3] = \
##!#                np.average(np.average(modis_dict['grid_modis_conc'][mi::num_months,:,:]
##!#                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_20_40])
##!#
##!#            #modis_avgs[ri,mi] = \
##!#                #np.average(modis_dict['grid_modis_conc'][0::num_months]
##!#                #    [indices[ri]:indices[ri] + 3],axis = 0)
##!#
##!#    # Generate the figure
##!#    plt.close()
##!#    fig1 = plt.figure(figsize=(8,6))
##!#    ax = plt.subplot()
##!#    if(inseason=='sunlight'):
##!#        months = ['Apr','May','June','July','Aug','Sep']
##!#    else:
##!#        months = ['Dec','Jan','Feb','Mar','Apr','May',\
##!#                  'June','July','Aug','Sep','Oct','Nov']
##!#    # Plot the 2001 - 2003 data
##!#    ax.plot(modis_avgs[0,:,0],color='black')
##!#    ax.plot(modis_avgs[0,:,1],color='tab:blue')
##!#    ax.plot(modis_avgs[0,:,2],color='tab:green')
##!#    ax.plot(modis_avgs[0,:,3],color='tab:red')
##!#    # Plot the 2009 - 2011 data
##!#    ax.plot(modis_avgs[1,:,0],'--',color='black')
##!#    ax.plot(modis_avgs[1,:,1],'--',color='tab:blue')
##!#    ax.plot(modis_avgs[1,:,2],'--',color='tab:green')
##!#    ax.plot(modis_avgs[1,:,3],'--',color='tab:red')
##!#    # Plot the 2016 - 2018 data
##!#    ax.plot(modis_avgs[2,:,0],linestyle='dotted',color='black')
##!#    ax.plot(modis_avgs[2,:,1],linestyle='dotted',color='tab:blue')
##!#    ax.plot(modis_avgs[2,:,2],linestyle='dotted',color='tab:green')
##!#    ax.plot(modis_avgs[2,:,3],linestyle='dotted',color='tab:red')
##!#    ax.set_xticks(np.arange(num_months),months)
##!#    ax.set_ylabel('Ice Concentration [%]')
##!#
##!#    # Shrink the current axis's height by 10% to make room for the legend
##!#    box = ax.get_position()
##!#    ax.set_position([box.x0, box.y0 + box.height * 0.1,\
##!#                    box.width, box.height * 0.9])
##!#
##!#    # Make the legend
##!#    custom_lines = [Line2D([0],[0],color='black'),\
##!#                    Line2D([0],[0],color='black'),\
##!#                    Line2D([0],[0],color='tab:blue'),\
##!#                    Line2D([0],[0],linestyle='dashed',color='black'),\
##!#                    Line2D([0],[0],color='tab:green'),\
##!#                    Line2D([0],[0],linestyle='dotted',color='black'),\
##!#                    Line2D([0],[0],color='tab:red')]
##!#    ax.legend(custom_lines,['80 - 100%',\
##!#                             '2001 - 2003',\
##!#                             '60 - 80%',\
##!#                             '2009 - 2011',\
##!#                             '40 - 60%',\
##!#                             '2016 - 2018',\
##!#                             '20 - 40%'],\
##!#              loc = 'upper center',bbox_to_anchor=(0.5,-0.05),\
##!#              fancybox=True,shadow=True,ncol=4)
##!#    plt.show()
##!#    plt.close()

