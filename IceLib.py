"""
  NAME:
    IceLib.py

  PURPOSE:
    House all the functions used for reading and working with the ice data



"""
import numpy as np
import numpy.ma as ma
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy
import cartopy.crs as ccrs
import glob
import subprocess
from scipy import stats

def plot_data(ice_dict,index):
    data = ice_dict['data'][index,:,:]
    plt.pcolormesh(data.T)
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.title(ice_dict['titles'][index])
    plt.show()

def trend_calc(ice_dict,x_ind,y_ind,thielsen=False):
    temp_data = ice_dict['data'][:,x_ind,y_ind]
    avg_ice = np.average(temp_data)
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


    pcnt_change = (trend/avg_ice)*100.
    # Find the percent change per decade
    return trend,pcnt_change

def read_ice(season,pre2001=False):
    spring=False
    summer=False
    autumn=False
    winter=False
    sunlight=False
    season_adder = ' '+season
    file_adder = '_'+season
    if(season=='spring'):
        spring=True
        ls_check = ['03_','04_','05_']
    elif(season=='summer'):
        summer=True
        ls_check = ['06_','07_','08_']
    elif(season=='autumn'):
        autumn=True
        ls_check = ['09_','10_','11_']
    elif(season=='winter'):
        winter=True
        ls_check = ['12_','01_','02_']
    elif(season=='sunlight'):
        sunlight=True
        ls_check = ['04_','05_','06_','07_','08_','09_']
    else:
        print("Analyzing all seasons. Option given was",season)
        season_adder = ''
        file_adder = ''
        ls_check=['01_','02_','03_','04_','05_','06_','07_','08_','09_','10_','11_','12_']
    
    
    # Grab all the ice files
    #file_names = glob.glob(data_loc+'*.bin')
    # Using this syntax, ignores 2000
    cmnd = "ls /home/bsorenson/data/NSIDC/nt_*.bin"
    #status,output = commands.getstatusoutput(cmnd)
    #file_initial = output.strip().split('\n')
    #
    #file_names = []
    #for fname in file_initial:
    #    if(fname[50:52] in ls_check):
    #        file_names.append(fname)
    
    
    file_initial = subprocess.check_output(cmnd,shell=True).decode('utf-8').strip().split('\n')
    file_names = []

    if(pre2001==True):
        cmnd = "ls /home/bsorenson/data/NSIDC/pre_2001/nt_*.bin"
        pre_file_initial = subprocess.check_output(cmnd,shell=True).decode('utf-8').strip().split('\n')
        for fname in pre_file_initial:
            if(fname[-17:-14] in ls_check):
                file_names.append(fname)
        
    for fname in file_initial:
        if(fname[-17:-14] in ls_check):
            file_names.append(fname)
    
    # Read in the latitude, longitude, and area data
    latfileo = open('/home/bsorenson/Research/Ice_analysis/psn25lats_v3.dat','r')
    lats = np.reshape(np.fromfile(latfileo,dtype=np.uint32)/100000.,(448,304))
    lonfileo = open('/home/bsorenson/Research/Ice_analysis/psn25lons_v3.dat','r')
    tlons = np.fromfile(lonfileo,dtype=np.uint32)/100000.
    tlons[tlons>180.] = tlons[tlons>180.]-42949.67296
    lons = np.reshape(tlons,(448,304))
    areafileo = open('/home/bsorenson/Research/Ice_analysis/psn25area_v3.dat','r')
    areas = np.reshape(np.fromfile(areafileo,dtype=np.uint32)/1000.,(448,304))
    
    #lons = np.reshape(np.fromfile(lonfileo,dtype=np.uint32)/100000.,(448,304))
    
    num_files = len(file_names)
    ice_data = {}
    ice_data['data'] = np.full((num_files,448,304),-9.)
    ice_data['lat']  = lats
    ice_data['lon']  = lons
    ice_data['area']  = areas
    ice_data['trends'] = np.full((448,304),-9.)
    ice_data['land_trends'] = np.full((448,304),-9.)
    ice_data['titles'] = []
    ice_data['season_adder'] = season_adder
    ice_data['file_adder'] = file_adder
    
    count = 0
    for fname in file_names:
    
        total_name = fname
        #total_name = data_loc+fname
        print(total_name)
        fileo = open(total_name,'rb')
        test = bytearray(fileo.read())
        
        ## To print the header
        #print(test[0:299])
        
        # psg = polar stereographic grid
        
        missing_val             = int(test[0:5].decode('utf-8'))
        ncol_psg                = int(test[6:11].decode('utf-8'))
        nrow_psg                = int(test[12:17].decode('utf-8'))
        lat_enclosed_psg        = float(test[24:29].decode('utf-8'))
        greenwich_orient_psg    = float(test[30:35].decode('utf-8'))
        j_coord_pole            = float(test[42:47].decode('utf-8'))
        i_coord_pole            = float(test[48:53].decode('utf-8'))
        inst_desc               = test[54:59].decode('utf-8')
        data_desc               = test[60:65].decode('utf-8')
        julian_start_day        = int(test[66:71].decode('utf-8'))
        start_hour              = int(test[72:77].decode('utf-8'))
        start_min               = int(test[78:83].decode('utf-8'))
        julian_end_day          = int(test[84:89].decode('utf-8'))
        end_hour                = int(test[90:95].decode('utf-8'))
        end_min                 = int(test[96:101].decode('utf-8'))
        year                    = int(test[102:107].decode('utf-8'))
        julian_day              = int(test[108:113].decode('utf-8'))
        channel_desc            = int(test[114:119].decode('utf-8'))
        scaling_factor          = float(test[120:125].decode('utf-8'))
        file_name               = test[126:149].decode('utf-8')
        image_title             = test[150:229].decode('utf-8')
        data_info               = test[230:299].decode('utf-8')
        
        # Move the data into a np array
        data = np.reshape(np.array(test[300:]),(448,304))
        
        #for i in range(len(data)):
        #    templine = data[i,:].astype(np.float)
        #    templine[templine<251] = (templine[templine<251]/scaling_factor)*100.
        #    templine = templine.astype(int)
        #    print(templine)
        
        #data[np.where(data>=251)]=data[np.where(data>=251)]*0.
        data[np.where(data<251)]=  \
            (data[np.where(data<251)]/scaling_factor)*100.
    
        ice_data['data'][count,:,:] = data[:,:]
        ice_data['titles'].append(file_name)
    #    total_data[:,:,count] = data[:,:]
        count+=1

    return ice_data

def ice_trendCalc(ice_data,thielSen=False):
    # Loop over the data and calculate trends
    for i in range(448):
        print(i)
        max_trend = -99.
        min_trend = 99.
        for j in range(304):
            ice_data['trends'][i,j],temp_pcnt = trend_calc(ice_data,i,j,thielsen=thielSen)
            ice_data['land_trends'][i,j],temp_pcnt = trend_calc(ice_data,i,j,thielsen=thielSen)
            temp_trend = ice_data['trends'][i,j]
            if(temp_trend>max_trend):
                max_trend = temp_trend
            if(temp_trend<min_trend):
                min_trend = temp_trend
        print("  max trend = ",max_trend)
        print("  min trend = ",min_trend)
        # Deal with land masks
        good_indices = np.where(ice_data['data'][0,i,:]<251.)
        land_indices = np.where(ice_data['data'][0,i,:]>=251.)
        ice_data['trends'][i,land_indices] = np.nan
        ice_data['land_trends'][i,good_indices] = np.nan

    return ice_data

def ice_gridtrendCalc(ice_data,area=True,thielSen=False):
    if(area==True):
        print("\nCalculating area trends\n")
    else:
        print("\nCalculating % concentration trends\n")
    ice_data['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = ice_data['grid_lat']
    lon_ranges = ice_data['grid_lon']
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
                interp_data = (ice_data['grid_ice_conc'][:,xi,yj]/100.)*ice_data['grid_total_area'][xi,yj]
                good_indices = np.where(np.isnan(interp_data)==False)
            else:
                interp_data = ice_data['grid_ice_conc'][:,xi,yj]
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
                ice_data['grid_ice_area_trend'][xi,yj] = total_trend
        print(xi)
        #print("max trend = ",max_trend)
        #print("min trend = ",min_trend)

    return ice_data

def grid_data_conc(ice_dict):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_ice_conc    = np.full((len(ice_dict['data'][:,0,0]),len(lat_ranges),len(lon_ranges)),-99.)
    grid_ice_area    = np.zeros((len(lat_ranges),len(lon_ranges)))
    grid_ice_area_trend    = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    grid_ice_conc_cc = np.full((len(ice_dict['data'][:,0,0]),len(lat_ranges),len(lon_ranges)),-999.)
    print("Size of grid array: ",grid_ice_conc.shape)
    for nt in range(len(ice_dict['data'][:,0,0])):
        print(nt)
        for xi in range(448):
            # Don't grid the data if any portion of the lat/lon box is over land.
            # Don't include land data
            for yj in range(304):
                lat_index = np.where(np.floor(ice_dict['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(ice_dict['lon'][xi,yj])>=lon_ranges)[-1][-1]
                # Add the current pixel area into the correct grid box, no
                # matter if the current box is missing or not.
                #if((lat_index==20) & (lon_index==10)):
                #    print("Current grid area = ",grid_ice_area[lat_index,lon_index])
                #if((lat_index==20) & (lon_index==10)):
                #    print("    New grid area = ",grid_ice_area[lat_index,lon_index])
                if(ice_dict['data'][nt,xi,yj] < 251):
                    if(nt==0): grid_ice_area[lat_index,lon_index] += ice_dict['area'][xi,yj]
                    if(grid_ice_conc_cc[nt,lat_index,lon_index]==-999.):
                        grid_ice_conc[nt,lat_index,lon_index] = ice_dict['data'][nt,xi,yj]
                        grid_ice_conc_cc[nt,lat_index,lon_index] = 1.
                    else:
                        grid_ice_conc[nt,lat_index,lon_index] = ((grid_ice_conc[nt,lat_index,lon_index]*grid_ice_conc_cc[nt,lat_index,lon_index])+\
                                                         ice_dict['data'][nt,xi,yj])/(grid_ice_conc_cc[nt,lat_index,lon_index]+1.)
                        grid_ice_conc_cc[nt,lat_index,lon_index]+=1
                else:
                    if(nt==0): grid_ice_area[lat_index,lon_index] = np.nan

                    # end else
                # end if good ice check
            # end y grid loop
        # end x grid loop
    # end time loop 
    ice_dict['grid_ice_conc'] = grid_ice_conc
    ice_dict['grid_ice_conc_cc'] = grid_ice_conc_cc
    ice_dict['grid_total_area'] = grid_ice_area
    ice_dict['grid_ice_area_trend'] = grid_ice_area_trend
    ice_dict['grid_lat'] = lat_ranges
    ice_dict['grid_lon'] = lon_ranges
    
    return ice_dict

def grid_data(ice_dict):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_ice    = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    grid_ice_cc = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    for xi in range(448):
        print(xi)
        for yj in range(304):
            if(not np.isnan(ice_dict['trends'][xi,yj])):
                lat_index = np.where(np.floor(ice_dict['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(ice_dict['lon'][xi,yj])>=lon_ranges)[-1][-1]
                if(grid_ice_cc[lat_index,lon_index]==-999.):
                    grid_ice[lat_index,lon_index] = ice_dict['trends'][xi,yj]
                    grid_ice_cc[lat_index,lon_index] = 1.
                else:
                    grid_ice[lat_index,lon_index] = ((grid_ice[lat_index,lon_index]*grid_ice_cc[lat_index,lon_index])+\
                                                     ice_dict['trends'][xi,yj])/(grid_ice_cc[lat_index,lon_index]+1.)
                    grid_ice_cc[lat_index,lon_index]+=1
    ice_dict['grid_ice'] = grid_ice
    ice_dict['grid_ice_cc'] = grid_ice_cc
    ice_dict['grid_lat'] = lat_ranges
    ice_dict['grid_lon'] = lon_ranges
    return ice_dict

def plot_grid_data(ice_dict,t_ind,pvar,adjusted=False,save=False):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)

    if(pvar=='ice'):
        plabel = "Percent Ice Concentration"

    local_grid_ice = np.copy(ice_dict['grid_ice_conc'][t_ind,:,:])
    local_grid_ice_bad = np.copy(ice_dict['grid_ice_conc'][t_ind,:,:])

    local_grid_ice[local_grid_ice==-999.] = np.nan
    local_grid_ice_bad[local_grid_ice!=-999.] = np.nan
    plot_good_data = ma.masked_invalid(local_grid_ice)
    plot_land_data = ma.masked_invalid(local_grid_ice_bad)

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
    ax.set_title('Gridded NSIDC Sea Ice Concentration\n'+ice_dict['titles'][t_ind])
    if(save==True):
        outname = 'ice_trend_gridded_200012_201812'+ice_dict['file_adder']+file_adder+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
    

def plot_grid_trend(ice_dict,adjusted=False,save=False):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)

    local_grid_ice = np.copy(ice_dict['grid_ice'])
    local_grid_ice_bad = np.copy(ice_dict['grid_ice'])

    local_grid_ice[local_grid_ice==-999.] = np.nan
    local_grid_ice_bad[local_grid_ice!=-999.] = np.nan
    plot_good_data = ma.masked_invalid(local_grid_ice)
    plot_land_data = ma.masked_invalid(local_grid_ice_bad)

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
    if(ice_dict['season_adder']!=''):
        ax.set_title('NSIDC Sea Ice Concentration'+ice_dict['season_adder'].title()+\
            ' Trends\nJan 2001 to Dec 2018')
    else:
        ax.set_title('NSIDC Sea Ice Concentration Trends\nJan 2001 to Dec 2018')
    if(save==True):
        outname = 'ice_trend_gridded_200101_201812'+ice_dict['file_adder']+file_adder+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
    

