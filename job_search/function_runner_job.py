#!/usr/bin/env python

"""


"""

from job_search_lib import *

atsci_file = 'atsci_programs.csv'
nws_file   = 'nws_locations.csv'
lab_file   = 'labs_centers.csv'
priv_file  = 'private_sector.csv'

atsci = pd.read_csv(atsci_file)
nws   = pd.read_csv(nws_file)
labs  = pd.read_csv(lab_file)
priv  = pd.read_csv(priv_file)

lat_lon_dict = read_lat_lon_from_csv()

radius_miles = 50
radius_km = radius_miles / 0.611


def check_point_distance(data, threshold, plat, plon, lat_lon_dict):

    if('University' in data.keys()):
        loop_var = 'University'
    elif('WFO' in data.keys()):
        loop_var = 'WFO'
    elif('Lab/Center' in data.keys()):
        loop_var = 'Lab/Center'
    elif('Company' in data.keys()):
        loop_var = 'Company'

    loop_list = list(data[loop_var])
    # Grab all the lats and lons for the loop data
    # --------------------------------------------
    loop_lats = np.full(len(loop_list), np.nan)
    loop_lons = np.full(len(loop_list), np.nan)

    for ii in range(len(loop_list)):

        xx = data['City'][ii] 
        yy = data['State'][ii]
        if(len(yy) > 2):
            ll_yy = us_state_to_abbrev[yy]
        else:
            ll_yy = yy
        city_name = xx + ' ' + ll_yy

        if(city_name in lat_lon_dict.keys()):
            llat = lat_lon_dict[city_name]['Lat']
            llon = lat_lon_dict[city_name]['Lon']

            dist = np.round(find_distance_between_points(\
                plat, plon, llat, llon), 1)

            if(dist < threshold):
                print(loop_var,':', loop_list[ii])

def mouse_event(event):
    ix, iy = event.xdata, event.ydata
    event_list = [ix, iy]
    # convert from display coordinates to data coordinates
    p_a_cart = datacrs.transform_point(event_list[0], event_list[1], src_crs=mapcrs)
    #p_a_cart[1] = lat, p_a_cart[0] = lon
    plat = p_a_cart[1]
    plon = p_a_cart[0]

    
    print("\n")
    check_point_distance(atsci, 20, plat, plon, lat_lon_dict)
    check_point_distance(nws,   20, plat, plon, lat_lon_dict)
    check_point_distance(labs,  20, plat, plon, lat_lon_dict)
    check_point_distance(priv,  20, plat, plon, lat_lon_dict)



#count_dict = calc_distances(nws, atsci, labs, priv, threshold = radius_km, lat_lon_dict = lat_lon_dict)

fig1 = plt.figure(figsize = (9, 9))
ax2 = fig1.add_subplot(2,2,1, projection = mapcrs)
ax3 = fig1.add_subplot(2,2,2, projection = mapcrs)
ax4 = fig1.add_subplot(2,2,3, projection = mapcrs)
ax5 = fig1.add_subplot(2,2,4, projection = mapcrs)

plot_job_sites_on_map(ax2, atsci, size = 7, color = None, lat_lon_dict = lat_lon_dict)
plot_job_sites_on_map(ax3, nws, size = 7, color = None, lat_lon_dict = lat_lon_dict)
plot_job_sites_on_map(ax4, labs, size = 7, color = None, lat_lon_dict = lat_lon_dict)
plot_job_sites_on_map(ax5, priv, size = 7, color = None, lat_lon_dict = lat_lon_dict)
ax2.set_title('Universities')
ax3.set_title('NWS WFOs')
ax4.set_title('Labs/Centers')
ax5.set_title('Companies')

fig1.tight_layout()

fig2 = plt.figure(figsize = (10,6))
ax1 = fig2.add_subplot(1,1,1, projection = mapcrs)
plot_all_sites(atsci, nws, labs, priv, ax = ax1, lat_lon_dict = lat_lon_dict, radius = radius_km)
ax1.set_title('Combined')
fig2.tight_layout()

fig4 = plt.figure()
ax1 = fig4.add_subplot(1,1,1)
#plot_all_sites(atsci, nws, labs, priv, ax = ax1, lat_lon_dict = lat_lon_dict, radius = radius_km)
plot_change_in_cwa_with_dist(ax1, nws, atsci, labs, priv, \
    lat_lon_dict = lat_lon_dict, min_threshold = 2, max_threshold = 100, passed_thresh = radius_miles)
ax1.grid(color = 'grey', alpha = 0.25, linestyle = '--')
#ax1.set_title('Colocated')
fig4.tight_layout()

fig3 = plt.figure(figsize = (10,6))
ax1 = fig3.add_subplot(1,1,1, projection = mapcrs)
#plot_all_sites(atsci, nws, labs, priv, ax = ax1, lat_lon_dict = lat_lon_dict, radius = radius_km)
plot_coloc_nws_sites(ax1, nws, atsci, labs, priv, lat_lon_dict = lat_lon_dict, draw_circles = False, threshold = radius_km)

cid1 = fig1.canvas.mpl_connect('button_press_event', mouse_event)
cid2 = fig2.canvas.mpl_connect('button_press_event', mouse_event)
cid3 = fig3.canvas.mpl_connect('button_press_event', mouse_event)

ax1.set_title('Colocated')
fig3.tight_layout()

plt.show()



#calc_distances(nws, atsci, threshold = 60)
#city_name = 'Quad Cities, IL'
#geolocator = Nominatim(user_agent = 'myapplication')
#location = geolocator.geocode(city_name)
#llat = location.latitude
#llon = location.longitude
#print(llat, llon)

# Loop over each of the NWS sites / universities and find
#     1. the closest university / NWS site
#     2. the number of universities / NWS sites, national labs, and companies within X kms


