#!/usr/bin/env python

"""


"""

from MODISLib import *


base_date = datetime(2020,9,22)
end_date  = datetime(2020,9,30)
local_date = base_date
#tester = identify_MODIS_MYD08(date_str, dest_dir = myd08_dir, dtype = 'daily')
#write_MODIS_MYD08(date_str, minlat = 65.5)

while(local_date <= end_date):
    
    date_str = local_date.strftime('%Y%m%d') 
    print(date_str)
    
    write_MODIS_MYD08(date_str, minlat = 65.5)
    ##write_MODIS_MYD08_monthly(date_str, minlat = 65.5)

    # Remove the old single-day files 
    cmnd = local_date.strftime('rm ' + myd08_dir + 'daily/MYD08_D3*.A%Y%j.*.hdf')
    #cmnd = local_date.strftime('rm ' + myd08_dir + 'MYD08_M3*.A%Y%j.*.hdf')
    print(cmnd)
    os.system(cmnd)

    # Increment the month and continue
    local_date = local_date + timedelta(days = 1)

    if(local_date.month == 10):
        local_date = local_date + relativedelta(months = 6)


sys.exit()





date_str = '202107222110'

MODIS_ch31 = read_MODIS_channel(date_str, 31)

filename = '/home/bsorenson/data/MODIS/Aqua/MYD/MYD021KM.A2021203.2110.061.2021204155922.hdf'
data1 = Dataset(filename)
data2 = SD.SD(filename)

sys.exit()

#dt_date_str = '201807052305'
#dt_date_str = '201807052310'
dt_date_str = '201507082105'
MODIS_ch1  = read_MODIS_channel(dt_date_str, 1,  zoom = False, swath = True, include_cloud_mask = True, include_myd06 = True)
MODIS_ch7  = read_MODIS_channel(dt_date_str, 7,  zoom = False, swath = True, include_cloud_mask = True, include_myd06 = True)

fig = plt.figure(figsize = (7,6))
ax1 = fig.add_subplot(2,2,1, projection = ccrs.NorthPolarStereo())
ax2 = fig.add_subplot(2,2,2, projection = ccrs.NorthPolarStereo())
ax3 = fig.add_subplot(2,2,3, projection = ccrs.NorthPolarStereo())
ax4 = fig.add_subplot(2,2,4, projection = ccrs.NorthPolarStereo())

mesh = ax1.pcolormesh(MODIS_ch1['lon'], MODIS_ch1['lat'], MODIS_ch1['data'], \
    transform = datacrs, shading = 'auto', cmap = MODIS_ch1['colors'], vmax = 0.7)
cbar = plt.colorbar(mesh, ax = ax1, label = 'Reflectance')
mesh = ax2.pcolormesh(MODIS_ch1['lon'], MODIS_ch1['lat'], MODIS_ch1['cloud_mask'], \
    transform = datacrs, \
    shading = 'auto', cmap = 'jet')
cbar = plt.colorbar(mesh, ax = ax2, label = 'Cloud Mask')
#ax2.pcolormesh(cldmsk_lon, cldmsk_lat, cldmsk_mask, transform = datacrs, \
#    shading = 'auto', cmap = 'jet')
mesh = ax3.pcolormesh(MODIS_ch1['lon'], MODIS_ch1['lat'], \
    MODIS_ch1['cloud_optical_depth'], transform = datacrs, \
    shading = 'auto', vmax = 40)
cbar = plt.colorbar(mesh, ax = ax3, label = 'Cloud Optical Depth')
mesh = ax4.pcolormesh(MODIS_ch7['lon'], MODIS_ch7['lat'], MODIS_ch7['data'], \
    transform = datacrs, vmax = 0.4, shading = 'auto', cmap = MODIS_ch1['colors'])
cbar = plt.colorbar(mesh, ax = ax4, label = 'Reflectance')
#ax4.pcolormesh(myd06_lon, myd06_lat, myd06_dat, transform = datacrs, \
#    shading = 'auto')

ax1.coastlines()
ax2.coastlines()
ax3.coastlines()
ax4.coastlines()

ax1.set_title('0.64 μm')
ax2.set_title('Cloud Mask\nBlue = Cloud, Red = Clear')
ax3.set_title('MODIS COD')
ax4.set_title('MODIS 2.1 μm Refl.')

ax1.set_boundary(circle, transform=ax1.transAxes)
ax2.set_boundary(circle, transform=ax2.transAxes)
ax3.set_boundary(circle, transform=ax3.transAxes)
ax4.set_boundary(circle, transform=ax4.transAxes)

ax1.set_extent([-180,180,65,90], datacrs)
ax2.set_extent([-180,180,65,90], datacrs)
ax3.set_extent([-180,180,65,90], datacrs)
ax4.set_extent([-180,180,65,90], datacrs)

plt.suptitle(dt_date_str)

fig.tight_layout()

#fig.savefig('modis_cod_compare_' + dt_date_str + '.png', dpi = 200)

plt.show()

sys.exit()






date_strs = [
    "201708172135",
    "201708172140",
    "201708172145",
    "201708172150"]
for dstr in date_strs:
    identify_MODIS_MYD06(dstr)

sys.exit()


#
#MODIS_data = '201507082105'
#minlat = 70
#write_MODIS_to_HDF5(MODIS_data, channel = 1, swath = True, \
#    save_path = './', minlat = minlat, remove_empty_scans = True, \
#    include_cloud_mask = True, include_myd06 = True)
#
#sys.exit()



CERES_date_str =  "2015070822"
download = True
download_dict = download_MODIS_swath(CERES_date_str,download = download)

sys.exit()

#myd06_files = [
#    'MYD06_L2.A2018186.2305.061.2018187160426.hdf',\
#    'MYD06_L2.A2018186.2310.061.2018187160720.hdf',\
#    'MYD06_L2.A2018186.2315.061.2018187155509.hdf']
##!#
##!#date_strs = [
##!#    '201507082105',\
##!#    '201507082110',\
##!#    '201507082115',\
##!#    '201507082120',\
##!#    ]
##!#
##!#for dstr in date_strs:
##!#    identify_MODIS_MYD06(dstr)
##!#
##!#sys.exit()
#start_date = '200207'
#end_date   = '201909'
start_date = '200504'
end_date   = '202009'
myd08_data = read_MODIS_MYD08_monthrange(start_date,end_date,\
    minlat=65.5, calc_month = False)
cloud_data = read_MODIS_CLDL3_monthrange(start_date,end_date,\
    minlat=65.5, calc_month = True)
plot_compare_CLDL3_MYD08_trend(cloud_data, myd08_data, month_idx = None, \
    minlat = 65.5, trend_type = 'standard', season = '', myd08_var = 'all', ax = None, colorbar = True, \
    norm_to_decade = True, save = True, show_pval = True, comp_xidx = 10, comp_yidx = 185)
sys.exit()

plotMODIS_MYD08_MonthTrend(cloud_data,month_idx=None,save=False,\
    trend_type='theil-sen',pvar = 'day', season='',minlat=65.5,return_trend=False, \
    colorbar = True, colorbar_label_size = None,title = None, \
    ax = None, show_pval = False, uncert_ax = None, norm_to_decade = True, vmin = -0.10, vmax = 0.10)

#start_date = '200504'
#end_date   = '202009'
#cloud_data = read_MODIS_CLDL3_monthrange(start_date,end_date,\
#    minlat=65.5, calc_month = True)
#plotMODIS_CLDL3_MonthTrend(cloud_data,month_idx=None,save=False,\
#    trend_type='theil-sen',season='',minlat=65.5,return_trend=False, \
#    colorbar = True, colorbar_label_size = None,title = None, \
#    ax = None, show_pval = False, uncert_ax = None, norm_to_decade = True)
plt.show()
sys.exit()

fig = plt.figure(figsize = (9, 6))
ax1 = fig.add_subplot(2,3,1, projection = ccrs.NorthPolarStereo())
ax2 = fig.add_subplot(2,3,2, projection = ccrs.NorthPolarStereo())
ax3 = fig.add_subplot(2,3,3, projection = ccrs.NorthPolarStereo())
ax4 = fig.add_subplot(2,3,4, projection = ccrs.NorthPolarStereo())
ax5 = fig.add_subplot(2,3,5, projection = ccrs.NorthPolarStereo())
ax6 = fig.add_subplot(2,3,6, projection = ccrs.NorthPolarStereo())
plotMODIS_MYD08_MonthTrend(cloud_data,month_idx=0,save=False,\
    trend_type='standard',pvar = 'day', season='',minlat=65.5,return_trend=False, \
    colorbar = True, colorbar_label_size = None,title = None, \
    ax = ax1, show_pval = False, uncert_ax = None, norm_to_decade = True)
plotMODIS_MYD08_MonthTrend(cloud_data,month_idx=1,save=False,\
    trend_type='standard',pvar = 'day', season='',minlat=65.5,return_trend=False, \
    colorbar = True, colorbar_label_size = None,title = None, \
    ax = ax2, show_pval = False, uncert_ax = None, norm_to_decade = True)
plotMODIS_MYD08_MonthTrend(cloud_data,month_idx=2,save=False,\
    trend_type='standard',pvar = 'day', season='',minlat=65.5,return_trend=False, \
    colorbar = True, colorbar_label_size = None,title = None, \
    ax = ax3, show_pval = False, uncert_ax = None, norm_to_decade = True)
plotMODIS_MYD08_MonthTrend(cloud_data,month_idx=3,save=False,\
    trend_type='standard',pvar = 'day', season='',minlat=65.5,return_trend=False, \
    colorbar = True, colorbar_label_size = None,title = None, \
    ax = ax4, show_pval = False, uncert_ax = None, norm_to_decade = True)
plotMODIS_MYD08_MonthTrend(cloud_data,month_idx=4,save=False,\
    trend_type='standard',pvar = 'day', season='',minlat=65.5,return_trend=False, \
    colorbar = True, colorbar_label_size = None,title = None, \
    ax = ax5, show_pval = False, uncert_ax = None, norm_to_decade = True)
plotMODIS_MYD08_MonthTrend(cloud_data,month_idx=5,save=False,\
    trend_type='standard',pvar = 'day', season='',minlat=65.5,return_trend=False, \
    colorbar = True, colorbar_label_size = None,title = None, \
    ax = ax6, show_pval = False, uncert_ax = None, norm_to_decade = True)
plt.suptitle('MODIS MYD08 Terra/Aqua Cloud Frac Trend\n' +
    'Percent per decade\n2003 - 2019')
ax1.set_title('April')
ax2.set_title('May')
ax3.set_title('June')
ax4.set_title('July')
ax5.set_title('August')
ax6.set_title('September')
fig.tight_layout()
plt.show()
sys.exit()


base_date = datetime(2003,6,1)
end_date  = datetime(2005,3,1)
local_date = base_date

while(local_date <= end_date):
    
    date_str = local_date.strftime('%Y%m') 
    print(date_str)
    
    # Download the data for this month
    download_MODIS_CLDL3_monthly(date_str)

    # Create the monthly average and write to an out file
    write_MODIS_CLDL3_monthly(date_str, minlat = 65.5)

    # Remove the old single-day files 
    cmnd = local_date.strftime('rm ' + cloudL3_daily_dir + '*.A%Y*.nc')
    print(cmnd)
    os.system(cmnd)

    # Increment the month and continue
    local_date = local_date + relativedelta(months = 1)

sys.exit()




date_str = '201807'
#identify_MODIS_MYD08(date_str)
write_MODIS_MYD08_monthly(date_str, minlat = 65.5)
sys.exit()
start_date = '200504'
end_date   = '202009'
cloud_data = read_MODIS_CLDL3_monthrange(start_date,end_date,\
    minlat=65.5, calc_month = True)

##!#fig = plt.figure()
##!#ax = fig.add_subplot(1,1,1, projection = ccrs.NorthPolarStereo())
##!#plotMODIS_CLDL3_MonthTrend(cloud_data,month_idx=None,save=False,\
##!#    trend_type='theil-sen',season='',minlat=65.5,return_trend=False, \
##!#    colorbar = True, colorbar_label_size = None,title = None, \
##!#    ax = ax, show_pval = False, uncert_ax = None, norm_to_decade = True)
##!#ax.coastlines()
##!#ax.set_extent([-180,180,65.5, 90], datacrs)
##!#ax.set_title('MODIS CLDL3 Terra/Aqua Cloud Frac Trend\n' +
##!#    'Percent per decade\nApril - September, 2005 - 2020')
##!#fig.tight_layout()
##!#plt.show()
##!#sys.exit()



date_str = '201807'
#download_MODIS_CLDL3_monthly(date_str)
#cloud_data = read_MODIS_CLDL3_single_month(date_str, minlat = 65.5)
#found_file = identify_MODIS_CLDL3(date_str)
#cloud_data = read_MODIS_CLDL3_daily(date_str)
#date_str = '201807'
#cloud_data = read_MODIS_CLDL3_daily_allmonth(date_str)
#write_MODIS_CLDL3_monthly(cloud_data, minlat = 65.5)
plot_MODIS_CLDL3_2panel(date_str, var1 = 'cld_frac_mean', \
    var2 = 'cld_frac_std',minlat = 65.5, ax = None, colorbar = True,\
    save = True)
sys.exit()


sys.exit()

#date_str = '201807'
#tester = read_MODIS_CLDL3_daily_allmonth(date_str, minlat = 65.5)
#sys.exit()


dt_date_str = datetime.strptime(date_str, '%Y%m%d')

base_url = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/62/MCD06COSP_D3_MODIS/'

total_url = dt_date_str.strftime(base_url + '%Y/%j')

files = listFD(dt_date_str.strftime(total_url), ext = '.nc')

print(files[0])

final_url = 'wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 \"' + \
    total_url + '\" --header \"Authorization: Bearer ' + laads_daac_key + \
    '\" -P .'

print(final_url)
#os.system(final_url)

sys.exit()

testfile = '/home/bsorenson/data/MODIS/combined/MCD06COSP_D3_MODIS.A2018186.062.2022124045334.nc'
data = Dataset(testfile,'r')

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection = ccrs.NorthPolarStereo())

ax.pcolormesh(data['longitude'][:], data['latitude'][:], data['Cloud_Mask_Fraction/Mean'][:,:].T,\
    transform = datacrs, shading = 'auto')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
ax.set_extent([-180,180,65,90], datacrs)
data.close()
plt.show()

sys.exit()

times = [
        #"201507081750",
        "201507081925",
        "201507082105",
        "201507082245",
        "201807050140",
        "201807050315",
        "201807051630",
        "201807051805",
        "201807051945",
        "201807052125",
        "201807052305",
        ]

for time in times:
    plot_compare_MODIS_cloud(time, swath = True, save = True)

sys.exit()

#def cloud_check(cloud_list):
#    count_0 = cloud_list[0]
#    count_1 = cloud_list[1]
#    count_2 = cloud_list[2]
#    count_3 = cloud_list[3]
#    if(    ((count_0 > count_1) & (count_0 > count_2) & (count_0 > count_3))):
#            out_cloud = 0.
#    elif(  ((count_1 > count_0) & (count_1 > count_2) & (count_1 > count_3)) |
#           ((count_0 > count_3) & (count_1 >= count_0) & (count_1 >= count_2)) |
#           ((count_0 == count_3) & (count_1 > count_0) & (count_1 == count_2))
#           ):
#            out_cloud = 1.
#    elif(  ((count_2 > count_0) & (count_2 > count_1) & (count_2 > count_3)) |
#           ((count_3 > count_0) & (count_2 >= count_3) & (count_1 <= count_2)) |
#           ((count_0 == count_2) & (count_0 > count_3) & (count_2 > count_1))
#           ):
#            out_cloud = 2.
#    elif(  (count_3 >= count_0) & (count_3 >= count_1) & (count_3 >= count_2)):
#            out_cloud = 3.
#    return out_cloud
#
#
#sys.exit()

download = True
#CERES_date_str = '2018070523'
CERES_dates = [
                "2015070811",
                "2015070813",
                "2015070814",
                "2015070816",
                "2015070817",
             ]

for CERES_date_str in CERES_dates:
    download_dict = download_MODIS_swath(CERES_date_str, \
            dest_dir = modis_dir, download = download)

sys.exit()

#MODIS_data = "201807050140"
#MODIS_data = "201807050315"
#MODIS_data = "201807051630"
MODIS_data = "201807051805"
#MODIS_data = "201807051945"
#MODIS_data = "201807052125"

sys.exit()

#date_str = '201506011300'
#date_str = '201807052305'
#download_MODIS_file(date_str, dest_dir = modis_dir, download_cloud_mask = True)


       
#MODIS_data = '201908110125'
MODIS_data = '201807052305'
minlat = 70.
swath = True
write_MODIS_to_HDF5(MODIS_data, channel = 1, swath = True, \
    save_path = './', minlat = minlat, remove_empty_scans = True, \
    include_cloud_mask = False)
write_MODIS_to_HDF5(MODIS_data, channel = 2, swath = True, \
    save_path = './', minlat = minlat, remove_empty_scans = True, \
    include_cloud_mask = False)


sys.exit()

#swath = True
##!#channel = 2 
##!#minlat = 65.
##!#MODIS_ch1 = read_MODIS_channel(MODIS_data, 1, swath = swath, include_cloud_mask = True)
##!#MODIS_ch7 = read_MODIS_channel(MODIS_data, 7, swath = swath, include_cloud_mask = True)
##!#
##!##cldmsk = Dataset('CLDMSK_L2_MODIS_Aqua.A2019224.0035.001.2019224173701.nc')
##!##cldmsk = Dataset('CLDMSK_L2_MODIS_Aqua.A2019224.0030.001.2019224173707.nc')
##!#cldmsk = Dataset('/home/bsorenson/data/MODIS/Aqua/CLDMSK/CLDMSK_L2_MODIS_Aqua.A2018186.2305.001.2019064032332.nc')
##!#testmask = cldmsk['geophysical_data/Integer_Cloud_Mask'][::5,::5]
##!#testlat  = cldmsk['geolocation_data/latitude'][::5,::5]
##!#testlon  = cldmsk['geolocation_data/longitude'][::5,::5]
##!##mydmsk = Dataset('MYD35_L2.A2018186.2305.061.2018187153528.hdf')
##!#
##!#cldmsk.close()
##!##cldmsk['geophysical_data/Integer_Cloud_Mask']
##!##


plot_figure2(save=False, composite = True)
#plot_combined_figure1_v6(save = True)

sys.exit()

date_str = '202107222110'
plot_ceres_scatter(date_str, zoom=True,save=True,composite=True,\
    avg_pixel=True,plume_only=False)

sys.exit()


date1 = '202107222110'
channel = 31
dt_date_str1 = datetime.strptime(date1,"%Y%m%d%H%M")
#dt_date_str2 = datetime.strptime(date2,"%Y%m%d%H%M")
filename1 = aerosol_event_dict[dt_date_str1.strftime('%Y-%m-%d')][dt_date_str1.strftime('%H%M')]['modis']
#filename2 = aerosol_event_dict[dt_date_str2.strftime('%Y-%m-%d')][dt_date_str2.strftime('%H%M')]['modis']

MODIS_data1 = read_MODIS_channel(dt_date_str1.strftime('%Y%m%d%H%M'), channel, zoom = True)

compare_data1 = nearest_grid_values(MODIS_data1)

sys.exit()

#colocate_comparison(date1, date2, channel = 31)

sys.exit()

plot_CERES_swaths(date_str = '202107222110', save = False)
sys.exit()

#date_str = '202107202125'
#date_str = '202107222110'
#date_str = '202108062025'
#date_str = '202109012105'
channel1 = 1
channel2 = 5
channel3 = 31

#date_str = '201807051950'
date_str = '201807052127'

CERES_date_str = '2008042219'
download_MODIS_swath(CERES_date_str, \
    dest_dir = '/home/bsorenson/data/MODIS/Aqua/')

#date_str = '201807052125'
MODIS_data = date_str
#MODIS_data = '201908110125'
write_MODIS_to_HDF5(MODIS_data, channel = 2, swath = True, \
    save_path = '/home/bsorenson/Research/Arctic_compares/comp_data/' + date_str[:8] + '/')
write_MODIS_to_HDF5(MODIS_data, channel = 7, swath = True, \
    save_path = '/home/bsorenson/Research/Arctic_compares/comp_data/' + date_str[:8] + '/')
sys.exit()


sys.exit()

#plot_MODIS_asos_sites(date_str, sites = None, save = False)

#plot_combined_figure1_v4(date_str = '202107222110', zoom = True, show_smoke = False, composite = True, \
#        double_fig = True, save = True)
plot_spatial_scatter(date_str, zoom = True, composite = True,\
    avg_pixel = True, plume_only = False, save = False)
sys.exit()

#plot_viewing_geometry(date_str = '202107222110', zoom = True, show_smoke = False, composite = True, \
#        save=False)
#plot_MODIS_VIIRS_SBDART(save=True, composite = True, calc_radiance = True)
#plot_spatial_scatter_wAI(date_str, zoom=True,save=True,composite=True,\
#    avg_pixel=True,plume_only=False)
#plot_combined_scatter(date_str,channel0 = 31, channel1 = 1, channel2 = 5,\
#        zoom=True,save=False,composite=True,avg_pixel=True,plume_only=False)


#plot_combined_figure1_v3(date_str = '202107202125', zoom = True, show_smoke = True, composite = True, \
#        save=False)
plot_MODIS_GOES_SBDART(save=True, composite = True, calc_radiance = True)
#plot_combined_figure1_v2(date_str = '202107202125', zoom = True, show_smoke = False, composite = True, \
#    save=True)
#plot_true_color_satpy(date_str, ax = None, zoom = True, save = False, composite = False)
#plot_combined_figure1(date_str = date_str, zoom = True, show_smoke = True, composite = True, \
#        save=True)

#plot_scatter_OMI_CERES_figure(zoom = True, show_smoke = False, composite = True, \
#        plume_only = False, avg_pixel = True, save=True)
#plot_MODIS_detection(date_str, zoom = True, save = False)
#plot_MODIS_CERES_3panel(zoom = True, show_smoke = False, composite = True, \
#        save=False)

#compare_MODIS_3panel(date_str,channel1,channel2,channel3,zoom=True,save=True,\
#        plot_ASOS_loc = False, show_smoke = True, compare_OMI = False, \
#        compare_CERES = False, return_MODIS = False)
#plot_MODIS_temporary_4panel('202107222110', \
#    zoom = True, composite = True, show_smoke = True, save = True)
#plot_MODIS_temporary('202107222110', zoom = True, save = True)
#compare_MODIS_3panel('202107222110',31,1,5,zoom=True,save=True,\
#        plot_ASOS_loc = False, show_smoke = True)
#plot_combined_figure3(save = True)
#plot_figure2(save=True, composite = True)
#plot_figureS1(save=True, composite = True)
#plot_combined_figureS2(save=True)
