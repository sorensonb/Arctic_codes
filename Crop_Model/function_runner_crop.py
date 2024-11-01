#!/usr/bin/env python

"""


"""

from Crop_Model_Lib import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Read the integer MUKEY & lat/lon pairs from the netCDF files
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

l_AUTO_SELECT_FILES = False

if(l_AUTO_SELECT_FILES):

    # Grab all the netCDF files with integer MUKeys and lat/lons
    # ----------------------------------------------------------
    base_path = '/home/bsorenson/Research/Crop_Model/soil_nd*_1000.nc'
    #base_path = '/home/bsorenson/Research/Crop_Model/soil_nd*_100.nc'
    files = glob(base_path)

else:

    files = ['/home/bsorenson/Research/Crop_Model/soil_nd001_mukey_int_ncdf_100.nc',\
             '/home/bsorenson/Research/Crop_Model/soil_nd003_mukey_int_ncdf_100.nc',\
             '/home/bsorenson/Research/Crop_Model/soil_nd005_mukey_int_ncdf_100.nc',\
             '/home/bsorenson/Research/Crop_Model/soil_nd007_mukey_int_ncdf_100.nc',\
             '/home/bsorenson/Research/Crop_Model/soil_nd009_mukey_int_ncdf_100.nc',\
             '/home/bsorenson/Research/Crop_Model/soil_nd035_mukey_int_ncdf_1000.nc']

total_len = 0

beg_idx = 0
for tfile in files: 

    # Open the file
    # -------------
    data = Dataset(tfile, 'r')

    filename = tfile.strip().split('/')[-1]

    # Extract the lats, lons, and mukeys
    # ----------------------------------
    mukey = data['Band1'][:,:]
    lats = data['lat'][:]
    lons = data['lon'][:]

    # Close the file
    # --------------
    data.close()

    meshlon, meshlat = np.meshgrid(lons, lats)

    # Check distances
    # ---------------
    idx = 50
    dist1 = np.round(find_distance_between_points(lats[idx], lons[idx], lats[idx + 1], lons[idx]), 4)
    dist2 = np.round(find_distance_between_points(lats[idx], lons[idx], lats[idx], lons[idx + 1]), 4)
    dist3 = np.round(find_distance_between_points(lats[idx], lons[idx], lats[idx + 1], lons[idx + 1]), 4)
    

    # Remove any locations with missing mukeys
    # ----------------------------------------
    keep_idxs = ~mukey.flatten().mask
    good_mukey = mukey.flatten()[~mukey.flatten().mask]
    good_lat   = meshlat.flatten()[~mukey.flatten().mask]
    good_lon   = meshlon.flatten()[~mukey.flatten().mask]
    
    total_len += good_mukey.shape[0] 
    end_idx = beg_idx + good_mukey.shape[0]
    print(filename, good_mukey.shape, total_len, dist1, dist2, beg_idx, end_idx) 

    beg_idx = end_idx

# Set up arrays to hold all the data
# ----------------------------------
total_lat = np.full( (total_len), np.nan)
total_lon = np.full( (total_len), np.nan)
total_muk = np.full( (total_len), np.nan)

# Loop back over the files and insert into the arrays
# ---------------------------------------------------
beg_idx = 0
end_idx = 0
total_len = 0
for tfile in files: 

    # Open the file
    # -------------
    data = Dataset(tfile, 'r')

    filename = tfile.strip().split('/')[-1]

    # Extract the lats, lons, and mukeys
    # ----------------------------------
    mukey = data['Band1'][:,:]
    lats = data['lat'][:]
    lons = data['lon'][:]

    # Close the file
    # --------------
    data.close()

    meshlon, meshlat = np.meshgrid(lons, lats)

    # Remove any locations with missing mukeys
    # ----------------------------------------
    keep_idxs = ~mukey.flatten().mask
    good_mukey = mukey.flatten()[~mukey.flatten().mask]
    good_lat   = np.round(meshlat.flatten()[~mukey.flatten().mask], 5)
    good_lon   = np.round(meshlon.flatten()[~mukey.flatten().mask], 5)
    
    total_len += good_mukey.shape[0] 
    end_idx = beg_idx + good_mukey.shape[0]
    print(filename, good_mukey.shape, total_len) 
    total_lat[beg_idx:end_idx] = good_lat
    total_lon[beg_idx:end_idx] = good_lon
    total_muk[beg_idx:end_idx] = good_mukey
    beg_idx = end_idx

total_muk = total_muk.astype(int)

# Save to a .csv file
# -------------------
df = pd.DataFrame({'Lat': total_lat, 'Lon': total_lon, 'Mukey': total_muk})
df.to_csv('test_mukey.csv', index = False)


matchfile = 'test_mukey.csv'
inlat = 47.945
inlon = -97.55231
#matchfile = 'test_mukey_1000.csv'
#find_mukey_from_latlon(inlat, inlon, matchfile = matchfile)

sys.exit()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Work with the GDB file
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

print("hi")

gdb_file = 'gSSURGO_ND.gdb'

layers = read_layer_names_from_gdb(gdb_file)

print(layers)

sa_data = read_layer_from_gdb(gdb_file, 'SAPOLYGON')
#yld_data = read_layer_from_gdb(gdb_file, 'mucropyld')

# Get projection information
# --------------------------
print(sa_data['geometry'].crs)

# Convert the projection to EPSG 4326?
sa_data['geometry'].to_crs(epsg = 4326)

# First XY pair from MULTIPOLYGON: -146413.000 2589662.200 (NAD83 projection)
# To convert to lat/lon: 
#lon_coord = -146413.000
#lat_coord = 2589662.200
#convert_projections('EPSG:5070','EPSG:4326', lon_coord, lat_coord)

sys.exit()

# EXAMPLES FROM JZ
# ----------------
# Grand Forks:
#   mukey: 2642530
#   lat:  47.869052461271
#   lon: -97.1139702381136

#   mukey: 2642602
#   lat:   47.9844126247407
#   lon:  -97.201850380271

# Set up lat/lon pair to test
# ---------------------------
lat_val = 43.8469
lon_val = -91.2522
print("INITIAL COORDS:", lat_val, lon_val)

# Convert from EPSG 4326 (WSG-84) to EPSG 5070
# --------------------------------------------
nad83_coords = convert_projections('EPSG:4326', 'EPSG:5070', lat_val, lon_val)
print("\nNAD83 COORDS:", nad83_coords)

# Convert from EPSG 5070 to EPSG 4326 (WSG-84) 
# --------------------------------------------
wsg84_coords = convert_projections('EPSG:5070', 'EPSG:4326', nad83_coords[0], nad83_coords[1])
print("\nWSG84 COORDS:", wsg84_coords)
