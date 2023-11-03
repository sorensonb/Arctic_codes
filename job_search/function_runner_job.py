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

radius_miles = 40
radius_km = radius_miles / 0.611

calc_distances(nws, atsci, labs, priv, threshold = radius_km, lat_lon_dict = lat_lon_dict)

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


