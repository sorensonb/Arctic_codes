"""
  NAME:
    NSIDCLib.py

  PURPOSE:
    House all the functions used for reading and working with the ice data



"""
import numpy as np
import numpy.ma as ma
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy
import cartopy.crs as ccrs
import glob
import subprocess
import os
from scipy import stats
from datetime import datetime, timedelta
from netCDF4 import Dataset
import h5py

home_dir = os.environ['HOME']

sys.path.append(home_dir)
from python_lib import circle, plot_trend_line, nearest_gridpoint, \
    aerosol_event_dict, init_proj, plot_lat_circles, plot_figure_text, \
    plot_subplot_label

data_dir = home_dir + '/data/NSIDC/'
datacrs = ccrs.PlateCarree()
mapcrs = ccrs.NorthPolarStereo()

def plot_NSIDC_data(NSIDC_data,index):
    data = NSIDC_data['data'][index,:,:]
    plt.pcolormesh(data.T)
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.title(NSIDC_data['titles'][index])
    plt.show()

def trend_calc(NSIDC_data,x_ind,y_ind,thielsen=False):
    temp_data = NSIDC_data['data'][:,x_ind,y_ind]
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
    #cmnd = "ls /data/NSIDC/nt_*.bin"
    cmnd = "ls " + home_dir + "/data/NSIDC/nt_*.bin"
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
        cmnd = "ls " + home_dir + "/data/NSIDC/pre_2001/nt_*.bin"
        pre_file_initial = subprocess.check_output(cmnd,shell=True).decode('utf-8').strip().split('\n')
        for fname in pre_file_initial:
            if(fname[-17:-14] in ls_check):
                file_names.append(fname)
        
    for fname in file_initial:
        if(fname[-17:-14] in ls_check):
            file_names.append(fname)
    
    # Read in the latitude, longitude, and area data
    latfileo = open(home_dir + '/data/NSIDC/psn25lats_v3.dat','r')
    lats = np.reshape(np.fromfile(latfileo,dtype=np.uint32)/100000.,(448,304))
    lonfileo = open(home_dir + '/data/NSIDC/psn25lons_v3.dat','r')
    tlons = np.fromfile(lonfileo,dtype=np.uint32)/100000.
    tlons[tlons>180.] = tlons[tlons>180.]-42949.67296
    lons = np.reshape(tlons,(448,304))
    areafileo = open(home_dir + '/data/NSIDC/psn25area_v3.dat','r')
    areas = np.reshape(np.fromfile(areafileo,dtype=np.uint32)/1000.,(448,304))
    
    #lons = np.reshape(np.fromfile(lonfileo,dtype=np.uint32)/100000.,(448,304))
    
    num_files = len(file_names)
    NSIDC_data = {}
    NSIDC_data['data'] = np.full((num_files,448,304),-9.)
    NSIDC_data['lat']  = lats
    NSIDC_data['lon']  = lons
    NSIDC_data['area']  = areas
    NSIDC_data['trends'] = np.full((448,304),-9.)
    NSIDC_data['land_trends'] = np.full((448,304),-9.)
    NSIDC_data['titles'] = []
    NSIDC_data['season_adder'] = season_adder
    NSIDC_data['file_adder'] = file_adder
    
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
    
        NSIDC_data['data'][count,:,:] = data[:,:]
        NSIDC_data['titles'].append(file_name)
    #    total_data[:,:,count] = data[:,:]
        count+=1

    return NSIDC_data

# This function downloads an NSIDC netCDF file of daily sea ice concentration
# NOTE: This function may not work if not already logged in to Earthdata.
#       May need to fix this to account for URS cookies.
# ---------------------------------------------------------------------------
def download_NSIDC_daily(date_str, save_dir = data_dir):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d')
    base_link = "https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0051.002/"
    file_str = 'NSIDC0051_SEAICE_PS_N25km_%Y%m%d_v2.0.nc'

    # First check if the file is already there
    # ----------------------------------------
    if(os.path.exists(dt_date_str.strftime(save_dir + file_str))):
        print("NSIDC file already exists.")
        return
    else:

        total_cmnd = dt_date_str.strftime('wget ' + base_link + '%Y.%m.%d/' + file_str)
        print(total_cmnd)
        os.system(total_cmnd)
        total_cmnd = dt_date_str.strftime('mv ' + file_str + ' ' + save_dir)
        print(total_cmnd)
        os.system(total_cmnd)

def readNSIDC_daily(date_str, grid_data = False):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d') 
  
    # Read in the latitude, longitude, area and ice data
    # --------------------------------------------------
    latfileo = open(home_dir + '/data/NSIDC/psn25lats_v3.dat','r')
    lats = np.reshape(np.fromfile(latfileo,dtype=np.uint32)/100000.,(448,304))
    lonfileo = open(home_dir + '/data/NSIDC/psn25lons_v3.dat','r')
    tlons = np.fromfile(lonfileo,dtype=np.uint32)/100000.
    tlons[tlons>180.] = tlons[tlons>180.]-42949.67296
    lons = np.reshape(tlons,(448,304))
    areafileo = open(home_dir + '/data/NSIDC/psn25area_v3.dat','r')
    areas = np.reshape(np.fromfile(areafileo,dtype=np.uint32)/1000.,(448,304))

    # Look for both the binary file and the netCDF file.
    # --------------------------------------------------
    bin_filename = data_dir + 'nt_' + date_str + '_f17_v1.1_n.bin'
    net_filename = data_dir + 'NSIDC0051_SEAICE_PS_N25km_' + date_str + '_v2.0.nc'
    print(bin_filename)

    bin_found = os.path.exists(bin_filename)
    net_found = os.path.exists(net_filename)

    if((bin_found) and not (net_found)):
 
        print('Reading ', bin_filename)
        fileo = open(filename,'rb')
        test = bytearray(fileo.read())
        
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
        
        data[np.where(data<251)]=  \
            (data[np.where(data<251)]/scaling_factor)*100.
   
    else:
        print('Reading ', net_filename)
        # Read the netCDf file
        # --------------------
        in_data = Dataset(net_filename, 'r')
        if('F17_ICECON' in in_data.variables.keys()):
            data = in_data['F17_ICECON'][0,:,:] * 100.
        elif('F13_ICECON' in in_data.variables.keys()):
            data = in_data['F13_ICECON'][0,:,:] * 100.
        else:
            print("ERROR: incorrect NSIDC ice variable name")
            in_data.close()
            return

        in_data.close()
 
    NSIDC_data = {}
    NSIDC_data['data'] = data
    NSIDC_data['lat']  = lats
    NSIDC_data['lon']  = lons
    NSIDC_data['area']  = areas
    NSIDC_data['date'] = date_str
 
    return NSIDC_data

def plotNSIDC_daily(pax, NSIDC_data, minlat=65, \
        vmin = None, vmax = None, title = None, label = None, \
        labelsize = 11, labelticksize = 9, circle_bound = False, \
        gridlines = True, zoom = False):

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.ocean
    ##!#if(pax is None):
    ##!#    if(minlat < 45):
    ##!#        mapcrs = ccrs.Miller()
    ##!#    else:
    ##!#        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    ##!#if(grid_data):
    ##!#    plot_lon  = NSIDC_data['lon']
    ##!#    plot_lat  = NSIDC_data['lat']
    ##!#    if(param == 'total'):
    ##!#        mask_flux = NSIDC_data['swf'] + NSIDC_data['lwf']
    ##!#    else:
    ##!#        mask_flux = NSIDC_data[param.lower()]
    ##!#else:
    mask_data = np.ma.masked_where((NSIDC_data['data'] > 250) | \
        (NSIDC_data['data'] == 0), \
        NSIDC_data['data'])
 
    if(title == None):
        title =  'SSMI/S Sea Ice Concentration\n' + NSIDC_data['date']
    if(label == None):
        label = 'Concentration [%]'
    pax.set_title(title)

    #plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
    mesh = pax.pcolormesh(NSIDC_data['lon'], NSIDC_data['lat'], mask_data,\
        transform = datacrs, cmap = colormap, vmin = vmin, vmax = vmax, \
        shading = 'auto')
    cbar = plt.colorbar(mesh,ax = pax, orientation='vertical',\
        fraction = 0.046, pad = 0.04)

    cbar.set_label(label,fontsize = labelsize, weight='bold')
    cbar.ax.tick_params(labelsize=labelticksize)

    if(gridlines): 
        pax.gridlines()

    pax.coastlines(resolution = '50m')
    #pax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())

    if(circle_bound):
        pax.set_boundary(circle, transform=pax.transAxes)

# Plot just a single day 
# ----------------------
def plotNSIDC_daily_figure(date_str, minlat = 65., \
        lat_circles = None, grid_data = False, zoom = False, \
        vmin = None, vmax = None, circle_bound = True, 
        title = None, label = None, \
        ax = None, gridlines = False, save = False):

    dt_date_str = datetime.strptime(date_str, "%Y%m%d")
    
    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    ##!#if(grid_data):
    ##!#    NSIDC_data = readgridCERES_hrly_grid(date_str[:10], param, \
    ##!#        satellite = 'Aqua', minlat = minlat)

    ##!#    if(NSIDC_data is None):
    ##!#        print("ERROR: no data returned from readgridCERES_hrly_grid")
    ##!#        print("Quitting")
    ##!#        return
    ##!#else:
    NSIDC_data = readNSIDC_daily(date_str, grid_data = grid_data)

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    in_ax = True 
    if(ax is None): 
        in_ax = False
        plt.close('all')
        fig1 = plt.figure(figsize = (6,6))
        mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = init_proj(date_str)
        #mapcrs = ccrs.LambertConformal(central_longitude = -100.)
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot data 
    # ----------------------------------------------------
    if(zoom):
        circle_bound = False
    plotNSIDC_daily(ax, NSIDC_data, minlat = minlat, \
        vmin = vmin, vmax = vmax, title = title, label = label, \
        circle_bound = circle_bound, \
        gridlines = gridlines, zoom = False)

    if(zoom):
        ax.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lon'][0], \
                       aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lon'][1], \
                       aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lat'][0], \
                       aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lat'][1]],\
                       datacrs)
##!#        ax0.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
##!#                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
##!#                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
##!#                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
##!#                        datacrs)
    else:
        ax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())

    #ax.set_title(date_str)

    # ----------------------------------------------------
    # If the user wants circles along latitude lines,    
    # plot them here      
    # ----------------------------------------------------
    if(lat_circles is not None):
        plot_lat_circles(ax, lat_circles) 

    if(not in_ax):
        fig1.tight_layout()
        
        if(save == True):
            out_name = 'nsidc_daily_conc_' + NSIDC_data['date'] + '.png'
            plt.savefig(out_name,dpi=300)
            print('Saved image '+out_name)
        else:
            plt.show()

# Writes a MODIS channel dictionary to HDF5 for Fortran colocation
# ----------------------------------------------------------------
def writeNSIDC_to_HDF5(NSIDC_data, save_path = './'):

    if(isinstance(NSIDC_data, str)):
        NSIDC_data = readNSIDC_daily(NSIDC_data)

    # Convert the filename object to datetime
    # ---------------------------------------
    file_date = NSIDC_data['date']
    dt_date_str = datetime.strptime(file_date, '%Y%m%d')

    # Create a new netCDF dataset to write to the file
    # ------------------------------------------------
    outfile = save_path + 'nsidc_comp_' + file_date + '.hdf5'
    dset = h5py.File(outfile,'w')
 
    dset.create_dataset('latitude',  data = NSIDC_data['lat'][:,:])
    dset.create_dataset('longitude', data = NSIDC_data['lon'][:,:])
    dset.create_dataset('SeaIceConcentration', data = NSIDC_data['data'][:,:].data)

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)  

def writeNSIDC_toNCDF(NSIDC_data, save_path = './'):

    # Create a new netCDF dataset to write to the file
    # ------------------------------------------------
    outfile = save_path + 'nsidc_'+ NSIDC_data['date'] + '.nc'
    nc = Dataset(outfile,'w',format='NETCDF4')
  
    # Dimensions = lon, lat
    # Create the sizes of each dimension in the file. In this case,
    # the dimensions are "# of latitude" and "# of longitude"
    # -------------------------------------------------------------
    num_x = NSIDC_data['data'].shape[0]
    num_y = NSIDC_data['data'].shape[1]
    
    # Use the dimension size variables to actually create dimensions in 
    # the file.
    # ----------------------------------------------------------------- 
    n_x  = nc.createDimension('nx',num_x)
    n_y  = nc.createDimension('ny',num_y)

    # Create variables for the three dimensions. Note that since these
    # are variables, they are still given 'dimensions' using 'dlat', 
    # and 'dlon'. Latitude and longitude are each 2-d grids, so they are 
    # given 2 dimensions (dlat, dlon).
    # ------------------------------------------------------------------
    #X_DIM=nc.createVariable('x_dim','i4',('nx','ny'))
    #X_DIM.description='Across track'
    #X_DIM.units='km?'
    #Y_DIM=nc.createVariable('Latitude','i4',('nx','ny'))
    #Y_DIM.description='Along track'
    #Y_DIM.units='km?'

    # Create a variable for the AI data, dimensioned using both 
    # dimensions.
    # ---------------------------------------------------------  
    LAT = nc.createVariable('Latitude','f4',('nx','ny'))
    LAT.description='Latitude'
    LAT.units      ='Degrees'
    LON = nc.createVariable('Longitude','f4',('nx','ny'))
    LON.description='Longitude'
    LON.units      ='Degrees'

    CONC = nc.createVariable('SeaIceConcentration','f4',('nx','ny'))
    CONC.description='25 x 25 km daily sea ice concentration.'
    CONC.units = 'Percent'

    AREA = nc.createVariable('GridBoxArea','f4',('nx','ny'))
    AREA.description='Area of each grid box.'
    AREA.units = 'km2'

    # Fill in dimension variables one-by-one.
    # NOTE: not sure if you can insert the entire data array into the
    # dimension (for example, doing :
    #      LAT, LON = np.meshgrid(lat_ranges,lon_ranges) ),
    # so this could be
    # something for you to try. Might make this faster if it works
    # ---------------------------------------------------------------
    LON[:,:]      = NSIDC_data['lat'][:,:]
    LAT[:,:]      = NSIDC_data['lon'][:,:]
    CONC[:,:]     = NSIDC_data['data'][:,:]
    AREA[:,:]     = NSIDC_data['area'][:,:]

    # Save, write, and close the netCDF file
    # --------------------------------------
    nc.close()

    print("Saved file ",outfile)  


def write_Ice(NSIDC_data):
    print("yay")

    # Print of a file for each month
    num_months = int(len(NSIDC_data['titles']))

    # Set up the format string for printing off the data
    # format_list includes the formats for all the yearly average concentrations
    # Separate each by a space
    formats = '\t'.join(['{:6.4}' for item in NSIDC_data['titles']])
    # Combine the formats string with the other format string
    total_ice_format = '\t'.join(['{:.3}\t{:6.4}',formats,'\n'])

    # Set up the header line format string
    # Extract the individual years from the titles parameter
    string_dates = [tstring.strip()[3:9] for tstring in NSIDC_data['titles']]
    header_formats = '\t'.join(['{:6}' for item in NSIDC_data['titles']])
    total_header_format = '\t'.join(['{:6}\t{:6}',header_formats,'\n'])   
   
    filename = "nsidc_grid_ice_values.txt"
    #for ti in range(num_months):
    #    str_date = NSIDC_data['titles'][ti].strip()[3:9]
    #    filename = "nsidc_grid_ice_"+str_date+'.txt'
    with open(filename,'w') as fout:
        fout.write(total_header_format.format('Lat','Lon',*string_dates))
        for xi in range(len(NSIDC_data['grid_lat'])):
            print("Printing latitude ",NSIDC_data['grid_lat'][xi])
            for yj in range(len(NSIDC_data['grid_lon'])):
                fout.write(total_ice_format.format(NSIDC_data['grid_lat'][xi],\
                        NSIDC_data['grid_lon'][yj],*NSIDC_data['grid_ice_conc'][:,xi,yj]))

    print("Saved file",filename) 
    


# ice_trendCalc calculates the trends over the time period at each
# grid point on the 25x25 km grid.
def ice_trendCalc(NSIDC_data,thielSen=False):
    # Loop over the data and calculate trends
    for i in range(448):
        print(i)
        max_trend = -99.
        min_trend = 99.
        for j in range(304):
            NSIDC_data['trends'][i,j],temp_pcnt = trend_calc(NSIDC_data,i,j,thielsen=thielSen)
            NSIDC_data['land_trends'][i,j],temp_pcnt = trend_calc(NSIDC_data,i,j,thielsen=thielSen)
            temp_trend = NSIDC_data['trends'][i,j]
            if(temp_trend>max_trend):
                max_trend = temp_trend
            if(temp_trend<min_trend):
                min_trend = temp_trend
        print("  max trend = ",max_trend)
        print("  min trend = ",min_trend)
        # Deal with land masks
        good_indices = np.where(NSIDC_data['data'][0,i,:]<251.)
        land_indices = np.where(NSIDC_data['data'][0,i,:]>=251.)
        NSIDC_data['trends'][i,land_indices] = np.nan
        NSIDC_data['land_trends'][i,good_indices] = np.nan

    return NSIDC_data

def ice_gridtrendCalc(NSIDC_data,area=True,thielSen=False):
    if(area==True):
        print("\nCalculating area trends\n")
    else:
        print("\nCalculating % concentration trends\n")
    NSIDC_data['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = NSIDC_data['grid_lat']
    lon_ranges = NSIDC_data['grid_lon']
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
                interp_data = (NSIDC_data['grid_ice_conc'][:,xi,yj]/100.)*NSIDC_data['grid_total_area'][xi,yj]
                good_indices = np.where(np.isnan(interp_data)==False)
            else:
                interp_data = NSIDC_data['grid_ice_conc'][:,xi,yj]
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
                NSIDC_data['grid_ice_area_trend'][xi,yj] = total_trend
        print(xi)
        #print("max trend = ",max_trend)
        #print("min trend = ",min_trend)

    return NSIDC_data

# grid_data_conc grids the 25x25 km gridded ice concentration data into
# a 1x1 degree lat/lon grid
def grid_data_conc(NSIDC_data):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_ice_conc    = np.full((len(NSIDC_data['data'][:,0,0]),len(lat_ranges),len(lon_ranges)),-99.)
    grid_ice_area    = np.zeros((len(lat_ranges),len(lon_ranges)))
    grid_ice_area_trend    = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    grid_ice_conc_cc = np.full((len(NSIDC_data['data'][:,0,0]),len(lat_ranges),len(lon_ranges)),-999.)
    print("Size of grid array: ",grid_ice_conc.shape)
    for nt in range(len(NSIDC_data['data'][:,0,0])):
        print(nt)
        for xi in range(448):
            # Don't grid the data if any portion of the lat/lon box is over land.
            # Don't include land data
            for yj in range(304):
                lat_index = np.where(np.floor(NSIDC_data['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(NSIDC_data['lon'][xi,yj])>=lon_ranges)[-1][-1]
                # Add the current pixel area into the correct grid box, no
                # matter if the current box is missing or not.
                #if((lat_index==20) & (lon_index==10)):
                #    print("Current grid area = ",grid_ice_area[lat_index,lon_index])
                #if((lat_index==20) & (lon_index==10)):
                #    print("    New grid area = ",grid_ice_area[lat_index,lon_index])
                if(NSIDC_data['data'][nt,xi,yj] < 251):
                    if(nt==0): grid_ice_area[lat_index,lon_index] += NSIDC_data['area'][xi,yj]
                    if(grid_ice_conc_cc[nt,lat_index,lon_index]==-999.):
                        grid_ice_conc[nt,lat_index,lon_index] = NSIDC_data['data'][nt,xi,yj]
                        grid_ice_conc_cc[nt,lat_index,lon_index] = 1.
                    else:
                        grid_ice_conc[nt,lat_index,lon_index] = ((grid_ice_conc[nt,lat_index,lon_index]*grid_ice_conc_cc[nt,lat_index,lon_index])+\
                                                         NSIDC_data['data'][nt,xi,yj])/(grid_ice_conc_cc[nt,lat_index,lon_index]+1.)
                        grid_ice_conc_cc[nt,lat_index,lon_index]+=1
                else:
                    if(nt==0): grid_ice_area[lat_index,lon_index] = np.nan

                    # end else
                # end if good ice check
            # end y grid loop
        # end x grid loop
    # end time loop 
    NSIDC_data['grid_ice_conc'] = grid_ice_conc
    NSIDC_data['grid_ice_conc_cc'] = grid_ice_conc_cc
    NSIDC_data['grid_total_area'] = grid_ice_area
    NSIDC_data['grid_ice_area_trend'] = grid_ice_area_trend
    NSIDC_data['grid_lat'] = lat_ranges
    NSIDC_data['grid_lon'] = lon_ranges
    
    return NSIDC_data

# grid_data averages the 25x25 km trends into a 1x1 degree grid
# The 'grid_ice' paramter in the dictionary therefore contains
# the gridded trends.
def grid_data(NSIDC_data):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_ice    = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    grid_ice_cc = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    for xi in range(448):
        print(xi)
        for yj in range(304):
            if(not np.isnan(NSIDC_data['trends'][xi,yj])):
                lat_index = np.where(np.floor(NSIDC_data['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(NSIDC_data['lon'][xi,yj])>=lon_ranges)[-1][-1]
                if(grid_ice_cc[lat_index,lon_index]==-999.):
                    grid_ice[lat_index,lon_index] = NSIDC_data['trends'][xi,yj]
                    grid_ice_cc[lat_index,lon_index] = 1.
                else:
                    grid_ice[lat_index,lon_index] = ((grid_ice[lat_index,lon_index]*grid_ice_cc[lat_index,lon_index])+\
                                                     NSIDC_data['trends'][xi,yj])/(grid_ice_cc[lat_index,lon_index]+1.)
                    grid_ice_cc[lat_index,lon_index]+=1
    NSIDC_data['grid_ice'] = grid_ice
    NSIDC_data['grid_ice_cc'] = grid_ice_cc
    NSIDC_data['grid_lat'] = lat_ranges
    NSIDC_data['grid_lon'] = lon_ranges
    return NSIDC_data

# plot_grid_data generates a plot of the /
def plot_grid_data(NSIDC_data,t_ind,pvar,adjusted=False,save=False):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)

    if(pvar=='ice'):
        plabel = "Percent Ice Concentration"

    local_grid_ice = np.copy(NSIDC_data['grid_ice_conc'][t_ind,:,:])
    local_grid_ice_bad = np.copy(NSIDC_data['grid_ice_conc'][t_ind,:,:])

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
    ax.set_title('Gridded NSIDC Sea Ice Concentration\n'+NSIDC_data['titles'][t_ind])
    if(save==True):
        outname = 'ice_conc_gridded_200012_201812'+NSIDC_data['file_adder']+file_adder+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
    

def plot_grid_trend(NSIDC_data,adjusted=False,save=False):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)

    local_grid_ice = np.copy(NSIDC_data['grid_ice'])
    local_grid_ice_bad = np.copy(NSIDC_data['grid_ice'])

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
    if(NSIDC_data['season_adder']!=''):
        ax.set_title('NSIDC Sea Ice Concentration'+NSIDC_data['season_adder'].title()+\
            ' Trends\nJan 2001 to Dec 2018')
    else:
        ax.set_title('NSIDC Sea Ice Concentration Trends\nJan 2001 to Dec 2018')
    if(save==True):
        outname = 'ice_trend_gridded_200101_201812'+NSIDC_data['file_adder']+file_adder+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
    
def plot_grid_time_series(NSIDC_data,lat_ind,lon_ind,thielsen=False):
    inseason = NSIDC_data['season_adder'].strip()
    temp_data = NSIDC_data['grid_ice_conc'][:,lat_ind,lon_ind]
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
    plt.title('Gridded Ice Data: '+str(NSIDC_data['grid_lat'][lat_ind])+'x'+\
                str(NSIDC_data['grid_lon'][lon_ind])+'\n'+inseason.title()+\
                ' season of each year)')
    plt.plot(temp_data,label='observations') 
    plt.plot(regress_y,'--',label='trend')
    plt.ylabel('Ice Concentration [%]')
    plt.legend()
    outname = 'nsidc_grid_time_series_'+inseason+'_'+str(int(NSIDC_data['grid_lat'][lat_ind]))+'x'+\
                str(int(NSIDC_data['grid_lon'][lon_ind]))+'.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)   
    plt.show()

# This code replicates Jianglong's figure showing how the April - September
# ice concentration values have changed over the time period
# NOTE: To run, do something like
# >>> all_NSIDC_data = read_ice('all')
# >>> all_NSIDC_data = grid_data_conc(all_NSIDC_data)
# and then
# >>> plot_apr_sep_changes(all_NSIDC_data)
def plot_apr_sep_changes(NSIDC_data):
  
    inseason = NSIDC_data['season_adder'].strip()
 
    upper_vals = np.array([100,80,60,40])
    lower_vals = np.array([80,60,40,20])

    # Generate 3-year averages for 3 time periods:
    # 2001-2003
    # 2008-2011
    # 2016-2018
   
    # For now, look at April - September 
    if(inseason=='sunlight'):
        num_months = 6
        start_idx  = 0
    else:
        num_months = 12
        start_idx  = 4
    # Dimensions of ice_avgs are
    # 0 - year blocks (size = 3, see above)
    # 1 - months      (size = num_months)
    # 2 - ice values  (size = 4)
    ice_avgs = np.zeros((3,num_months,len(upper_vals)))
    indices = [0,8,15]
  
    # Base the locations for each ice concentration range on the average
    # April concentration between 2001 and 2003

    avg_Apr_old = \
        np.average(NSIDC_data['grid_ice_conc'][start_idx::num_months,:,:][:3,:,:],axis=0)
 
    # Use the locations during the first April average to base everything
    # on
    locations_80_100 = np.where((avg_Apr_old <= 100.) & (avg_Apr_old > 80.)) 
    locations_60_80  = np.where((avg_Apr_old <= 80.)  & (avg_Apr_old > 60.)) 
    locations_40_60  = np.where((avg_Apr_old <= 60.)  & (avg_Apr_old > 40.)) 
    locations_20_40  = np.where((avg_Apr_old <= 40.)  & (avg_Apr_old > 20.)) 

    # Fill the ice_avgs array with the data for each time period 
    for ri in range(len(indices)):
        for mi in range(num_months):
            print("Year block = ",ri,"  Month = ",mi)
            # Deal with 80 to 100
            ##temp_arr = np.zeros((3,len(locations_80_100)))
            ##temp_arr[0] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri]][locations_80_100]
            ##temp_arr[1] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 1][locations_80_100]
            ##temp_arr[2] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 2][locations_80_100]
            #ice_avgs[ri,mi,0] = np.average(temp_arr)
            ice_avgs[ri,mi,0] = \
                np.average(np.average(NSIDC_data['grid_ice_conc'][mi::num_months,:,:]
                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_80_100])

            ### Deal with 60 to 80
            ##temp_arr = np.zeros((3,len(locations_60_80)))
            ##temp_arr[0] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri]][locations_60_80]
            ##temp_arr[1] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 1][locations_60_80]
            ##temp_arr[2] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 2][locations_60_80]
            ##ice_avgs[ri,mi,1] = np.average(temp_arr)
            ice_avgs[ri,mi,1] = \
                np.average(np.average(NSIDC_data['grid_ice_conc'][mi::num_months,:,:]
                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_60_80])

            # Deal with 40 to 60
            ##temp_arr = np.zeros((3,len(locations_40_60)))
            ##temp_arr[0] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri]][locations_40_60]
            ##temp_arr[1] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 1][locations_40_60]
            ##temp_arr[2] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 2][locations_40_60]
            #ice_avgs[ri,mi,2] = np.average(temp_arr)
            ice_avgs[ri,mi,2] = \
                np.average(np.average(NSIDC_data['grid_ice_conc'][mi::num_months,:,:]
                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_40_60])

            # Deal with 40 to 60
            ##temp_arr = np.zeros((3,len(locations_20_40)))
            ##temp_arr[0] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri]][locations_20_40]
            ##temp_arr[1] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 1][locations_20_40]
            ##temp_arr[2] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 2][locations_20_40]
            ##ice_avgs[ri,mi,3] = np.average(temp_arr)
            ice_avgs[ri,mi,3] = \
                np.average(np.average(NSIDC_data['grid_ice_conc'][mi::num_months,:,:]
                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_20_40])

            #ice_avgs[ri,mi] = \
                #np.average(NSIDC_data['grid_ice_conc'][0::num_months]
                #    [indices[ri]:indices[ri] + 3],axis = 0)

    # Generate the figure
    plt.close()
    fig1 = plt.figure(figsize=(8,6))
    ax = plt.subplot()
    if(inseason=='sunlight'):
        months = ['Apr','May','June','July','Aug','Sep']
    else:
        months = ['Dec','Jan','Feb','Mar','Apr','May',\
                  'June','July','Aug','Sep','Oct','Nov']
    # Plot the 2001 - 2003 data
    ax.plot(ice_avgs[0,:,0],color='black')
    ax.plot(ice_avgs[0,:,1],color='tab:blue')
    ax.plot(ice_avgs[0,:,2],color='tab:green')
    ax.plot(ice_avgs[0,:,3],color='tab:red')
    # Plot the 2009 - 2011 data
    ax.plot(ice_avgs[1,:,0],'--',color='black')
    ax.plot(ice_avgs[1,:,1],'--',color='tab:blue')
    ax.plot(ice_avgs[1,:,2],'--',color='tab:green')
    ax.plot(ice_avgs[1,:,3],'--',color='tab:red')
    # Plot the 2016 - 2018 data
    ax.plot(ice_avgs[2,:,0],linestyle='dotted',color='black')
    ax.plot(ice_avgs[2,:,1],linestyle='dotted',color='tab:blue')
    ax.plot(ice_avgs[2,:,2],linestyle='dotted',color='tab:green')
    ax.plot(ice_avgs[2,:,3],linestyle='dotted',color='tab:red')
    ax.set_xticks(np.arange(num_months),months)
    ax.set_ylabel('Ice Concentration [%]')

    # Shrink the current axis's height by 10% to make room for the legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,\
                    box.width, box.height * 0.9])

    # Make the legend
    custom_lines = [Line2D([0],[0],color='black'),\
                    Line2D([0],[0],color='black'),\
                    Line2D([0],[0],color='tab:blue'),\
                    Line2D([0],[0],linestyle='dashed',color='black'),\
                    Line2D([0],[0],color='tab:green'),\
                    Line2D([0],[0],linestyle='dotted',color='black'),\
                    Line2D([0],[0],color='tab:red')]
    ax.legend(custom_lines,['80 - 100%',\
                             '2001 - 2003',\
                             '60 - 80%',\
                             '2009 - 2011',\
                             '40 - 60%',\
                             '2016 - 2018',\
                             '20 - 40%'],\
              loc = 'upper center',bbox_to_anchor=(0.5,-0.05),\
              fancybox=True,shadow=True,ncol=4)
    plt.show()
    plt.close()





