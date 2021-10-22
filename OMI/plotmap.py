import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color
import cartopy.crs as ccrs

name_dict = {
    'SurfaceAlbedo':           'alb',\
    'NormRadiance':            'nrad',\
    'Reflectivity':            'refl',\
    'CloudFraction':           'cldfrc',\
    'UVAerosolIndex':          'uvai',\
    'PixelQualityFlags':       'pxlqf',\
    'MeasurementQualityFlags': 'msmtqf',\
    'ViewingZenithAngle':      'vza',\
    'SolarZenithAngle':        'sza',\
    'RelativeZenithAngle':     'rza',\
    'GroundPixelQualityFlags': 'grndqf'
}

var_dict = {
    'SurfaceAlbedo':           {'min': 0.0,  'max': 1.0},\
    'Reflectivity':            {'min': 0.0,  'max': 1.0},\
    'NormRadiance':            {'min': 0.0,  'max': 0.1},\
    'CloudFraction':           {'min': None, 'max': None},\
    'UVAerosolIndex':          {'min': -2.0, 'max': 3.0 },\
    'PixelQualityFlags':       {'min': None, 'max': None},\
    'MeasurementQualityFlags': {'min': None, 'max': None},\
    'ViewingZenithAngle':      {'min': 0.0,  'max': 180.},\
    'SolarZenithAngle':        {'min': None, 'max': None},\
    'RelativeZenithAngle':     {'min': None, 'max': None},\
    'GroundPixelQualityFlags': {'min': 0, 'max': 7}\
}

def plotmap(variable,data,masked_data,str_type = '', plot_time = '', save = False):
    
    # definitions for the axes
    left, width = 0.1, 0.8
    bottom, height = 0.1, 0.7
    spacing = 0.005

    rect_map = [left, bottom + 0.15, width, height]
    rect_his = [left, bottom, width, 0.15]
    
    # Set up values for gridding the AI data
    latmin = 60 
    lat_gridder = latmin * 4.

    lat_ranges = np.arange(latmin,90.1,0.25)
    lon_ranges = np.arange(-180,180.1,0.25)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #
    # Plot the gridded data
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #print("Time to plot")
    mapcrs = ccrs.NorthPolarStereo()
    datacrs = ccrs.PlateCarree()
    colormap = plt.cm.jet

    # Set up the polar stereographic projection map
    plt.clf()
    plt.figure(figsize=(8, 8))

    if(latmin<45):
        axs_map = plt.axes(rect_map, projection = ccrs.Miller())
    else:
        axs_map = plt.axes(rect_map, projection = ccrs.NorthPolarStereo(central_longitude = 0.))
    axs_map.gridlines()
    axs_map.coastlines(resolution = '50m')
    
    # Use meshgrid to convert the 1-d lat/lon arrays into 2-d, which is needed
    # for pcolormesh.
    plot_lat, plot_lon = np.meshgrid(lat_ranges,lon_ranges)

    plt.title('OMI ' + variable + str_type + ' ' + plot_time)
    mesh = axs_map.pcolormesh(plot_lon, plot_lat,masked_data,transform = datacrs,cmap = colormap,\
            vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])

    # Center the figure over the Arctic
    axs_map.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
        
    axs_his = plt.axes(rect_his)
    
    axs_his.hist(masked_data.flatten(), 1000, log=True, range=(var_dict[variable]['min'], var_dict[variable]['max']))
    axs_his.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    # Depending on the desired variable, set the appropriate colorbar ticks
    if(variable == 'NormRadiance'):
        tickvals = np.arange(0.0,0.101,0.01)
    elif(variable == 'Reflectivity'):
        tickvals = np.arange(0.0,1.1,0.10)
    else:
        tickvals = np.arange(-2.0,4.1,0.5)

    cbar = plt.colorbar(mesh,ticks = tickvals,orientation='horizontal',pad=.01,\
        aspect=50,shrink = 0.905,label=variable)

    if(save == True):
        out_name = 'figs/omi_climo_' + name_dict[variable] + str_type + '_' + plot_time + '.png'
        plt.savefig(out_name)
        print('Saved image '+out_name)
        plt.close()
    else:
        plt.show()
        plt.close()
    
