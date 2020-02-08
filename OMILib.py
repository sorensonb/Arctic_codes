#!/usr/bin/env python
"""
  NAME:

  PURPOSE:
    Plot trends in OMI-measured aerosol index across the northern hemisphere.
    The average AI values are read from files produced by 
    Average_AerosolIndexCalculator.pro (included in ai_trend_codes_YYMMDD.tar.gz)

  PYTHON VERSION:
    2.6.6

  MODULES:
    - Custom module AILib
    - Matplotlib
    - mpl_toolkits.basemap
    - datetime
    - sys
    - numpy
    
"""

import numpy as np
import sys
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.colors as color
##import matplotlib.colors as colors
from matplotlib.colors import rgb2hex,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase
from scipy import stats
from netCDF4 import Dataset
import gzip

# Class MidpointNormalize is used to center the colorbar on zero
class MidpointNormalize(Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

# Function onclick performs an operation at the location of a click on the
# map. 
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    #print ix, iy  
    mapLon, mapLat = m(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = [OMI_dict[dictkey][date]['avg'] for date in                           \
        sorted(OMI_dict[dictkey].keys())]
    dates = sorted(OMI_dict[dictkey].keys())
    x_vals = np.arange(0,len(OMI_dict[dictkey].keys()))

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals,avgs)
    print(dictkey, slope*len(x_vals))
    #print slope/len(OMI_dict[dictkey].keys())
    regress_y = x_vals*slope+intercept

    #The slope
    S=0
    sm=0
    nx = len(avgs)
    num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
    Sn=np.zeros(num_d)
    for si in range(0,nx-1):
        for sj in range(si+1,nx):
            # Find the slope between the two points
            Sn[sm] = (avgs[si]-avgs[sj])/(si-sj) 
            sm=sm+1
        # Endfor
    # Endfor
    Snsorted=sorted(Sn)
    sm=int(num_d/2.)
    if(2*sm    == num_d):
        tsslope=0.5*(Snsorted[sm]+Snsorted[sm+1])
    if(2*sm+1 == num_d): 
        tsslope=Snsorted[sm+1]
    regress_ts = x_vals*tsslope+intercept

    # Convert the dates to datetime objects
    ddates = [datetime.strptime(date,'%Y%m') for date in dates] 

    label_dates = [ddate.strftime('%b %Y') for ddate in ddates]

    tickNumber = 12

    fig1 = plt.figure()
    plt.plot(avgs)
    plt.plot(x_vals,regress_y,linestyle='--',color='red',label='Control')
    plt.plot(x_vals,regress_ts,linestyle='--',color='blue',label='Thiel-Sen')
    plt.title(dictkey)
    plt.xticks(np.arange(0,len(avgs))[::-int(len(avgs)/tickNumber)],label_dates[::\
        -int(len(avgs)/tickNumber)],rotation=45)
    plt.ylabel('Average Aerosol Index')
    plt.legend()
    plt.show()

def onclick_climo(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = m(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([OMI_dict[dictkey][date] for date in sorted(OMI_dict[dictkey].keys())])
    ai_avg = np.average(avgs)
    print(dictkey, ai_avg)

def readOMI(inputfile,start_date,end_date,key=None):
    global OMI_dict
    OMI_dict = {}

    if(key is not None):
        OMI_dict[key]={}

    if(inputfile.strip().split('/')[-1].split('.')[-1]=='gz'):
        f = gzip.open(inputfile,'rb')
    else:
        f = open(inputfile,'r')
    #with open(inputfile,'r') as f:
    # Skip the header line
    for line in f:
        templine = line.strip().split()
        if(len(templine)>1):
            if((int(templine[0])>=start_date) & (int(templine[0])<=end_date)):
                if(key is not None):
                    if(templine[1]==key):
                        OMI_dict[key][templine[0].decode('utf-8')] = {}
                        OMI_dict[key][templine[0].decode('utf-8')]['avg']=float(templine[2])
                        OMI_dict[key][templine[0].decode('utf-8')]['#_obs']=int(templine[3])
                else:
                    # If the current lat/lon pair are not found in the dictionary's
                    # keys, then make a new subdictionary for it.
                    if(templine[1].decode('utf-8') not in OMI_dict.keys()):
                        OMI_dict[templine[1].decode('utf-8')] = {}
                    # If the current lat/lon pair are already in the dictionary's
                    # keys, then add the new data to the subdictionary
                    OMI_dict[templine[1].decode('utf-8')][templine[0].decode('utf-8')]={}
                    OMI_dict[templine[1].decode('utf-8')][templine[0].decode('utf-8')]['avg']=float(templine[2])
                    OMI_dict[templine[1].decode('utf-8')][templine[0].decode('utf-8')]['#_obs']=int(templine[3])
    f.close()    

    return OMI_dict

def writeOMI_toNCDF(OMI_dict,minlat=30):
    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)
    
    nc = Dataset('./omi_gridded_ai.nc','w',format='NETCDF4')
  
    # Dimensions = lat, lon, time
    testdict = OMI_dict['40x5']
    testkeys = list(testdict.keys())
    num_lat = len(lat_ranges)
    num_lon = len(lon_ranges)
    num_time = len(testdict.keys())
    times = np.arange(num_time)
    
    n_time = nc.createDimension('nmth',num_time)
    n_lat  = nc.createDimension('dlat',num_lat)
    n_lon  = nc.createDimension('dlon',num_lon)

    MONTH=nc.createVariable('MONTH','i2',('nmth'))
    MONTH.description='Months since October 2004'
    LAT=nc.createVariable('Latitude','i2',('dlat','dlon'))
    LAT.description='Latitude'
    LAT.units='Degrees'
    LON=nc.createVariable('Longitude','i2',('dlat','dlon'))
    LON.description='Longitude'
    LON.units='Degrees'
    AI = nc.createVariable('AI','f4',('nmth','dlat','dlon'))
    AI.description='Monthly Averaged Aerosol Index'
    OB_COUNT=nc.createVariable('OB_COUNT','i2',('nmth','dlat','dlon'))
    OB_COUNT.description='# of OMI AI measurements used in each monthly average'

    # Fill in dimension variables
    for i in range(num_time):
        MONTH[i] = times[i]
    for i in range(num_lat):
        for j in range(num_lon):
            LAT[i,j]=lat_ranges[i]
            LON[i,j]=lon_ranges[j]

    # Fill in actual variables
    for i in range(num_lat):
        print(lat_ranges[i])
        for j in range(num_lon):
            dictkey = (str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j])))
            if(dictkey not in OMI_dict):
                # Insert missing values for AI and count
                AI[:,i,j] = [-999.9 for m in range(num_time)]
                OB_COUNT[:,i,j] = [-99 for m in range(num_time)]
            else:
                for m in range(num_time):
                    timekey = testkeys[m]
                    if(timekey not in OMI_dict[dictkey]):
                        AI[m,i,j] = -999.9
                        OB_COUNT[m,i,j] = -99
                    else:
                        AI[m,i,j] = OMI_dict[dictkey][timekey]['avg']
                        OB_COUNT[m,i,j] = OMI_dict[dictkey][timekey]['#_obs']
    nc.close()
    

def plotOMI_MK(OMI_dict,start_date,end_date,save=False,file_type='XR123',season='',minlat=30.):
    if(file_type=='NXAR'):
        title_flabel = 'No XTrack and All Rows'
        outname_flabel = 'noX_allRows'
    elif(file_type=='XAR'):
        title_flabel = 'XTrack and All Rows'
        outname_flabel = 'xtrack_allRows'
    elif(file_type=='XR123'):
        title_flabel = 'XTrack and Rows 1-23'
        outname_flabel = 'xtrack_rows1to23'
    lat_ranges = np.arange(minlat,90,1.0)
    #lat_ranges = np.arange(lowest_lat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False
    if(season=='spring'):
        spring = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    OMI_dict[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                    OMI_dict[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    OMI_dict[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    OMI_dict[lkey].pop(tkey)

    # Find the lowest lat in the file
    #lowest_lat = float(sorted(OMI_dict.keys())[0].split('x')[0])
     
    # Set up the polar stereographic map
    fig1 = plt.figure()
    ax = plt.subplot(111)
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
   
    # Define the colors for each range of P-Value
    # Ranges: >=0.4,0.4-0.3,0.3-0.2,0.2-0.15,0.15-0.1,0.1-0.05,<0.05
    # Seven colors
    # - White,yellow,orange-yellow,orange,red-orange,red,dark red
    color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
    pRED   = np.array([255.,255.,255.,255.,255.,150.,150.])
    pGREEN = np.array([255.,255.,150., 75.,  0.,  0.,  0.])
    pBLUE  = np.array([255.,  0.,  0.,  0.,  0.,  0.,150.])
    pRED   = pRED/256.
    pGREEN = pGREEN/256.
    pBLUE  = pBLUE/256.

    noRED = 0.
    noGREEN = 0.
    noBLUE = 0.
    
    # - White,blue-green1,orange-yellow,orange,red-orange,red,dark red
    nRED   = np.array([255.,100.,  0.,  0.,  0.,  0.,  0.])
    nGREEN = np.array([255.,255.,255.,255.,200.,100.,  0.])
    nBLUE  = np.array([255.,  0.,  0.,255.,255.,255.,150.])
    nRED   = nRED/256.
    nGREEN = nGREEN/256.
    nBLUE  = nBLUE/256.
 
    ##cmap = plt.get_cmap('bwr')
    ### Center the colorbar on zero
    ##norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
    ###norm = Normalize(vmin=vmin,vmax=vmax)
    ##mapper = ScalarMappable(norm=norm,cmap=cmap)

    # Set up initial values for the analysis regions
    canada_avg = 0.
    canada_counts = 0
    siberia_avg = 0.
    siberia_counts = 0 
    eus_avg = 0.
    eus_counts = 0
    wus_avg = 0.
    wus_counts = 0
    ea_avg = 0.
    ea_counts = 0
    wa_avg = 0.
    wa_counts = 0
    europe_avg = 0.
    europe_counts = 0
    
    # Loop over all the keys and print(the regression slopes 
    # Grab the averages for the key
    max_pval = -10.
    min_pval = 10.
    #for i in range(15,20):
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = (str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j])))
            #print(dictkey)
            if(dictkey=='48x-97'):
                #print("sorted(OMI_dict[dictkey].keys())=",sorted(OMI_dict[dictkey].keys()))
                min_date = sorted(OMI_dict[dictkey].keys())[0]
                max_date = sorted(OMI_dict[dictkey].keys())[-1]
    
            # If no data are present for the curent lat/lon box, fill it with
            # black
            if(dictkey not in OMI_dict.keys()):
                colorRED = 0.
                colorGREEN = 0.
                colorBLUE = 0.
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(OMI_dict[dictkey].keys())==0):
                print(dictkey,pval,color_index_value,'NO DATA')
                #print('Here 2')
                colorRED = 0.
                colorGREEN = 0.
                colorBLUE = 0.
            else:
                avgs = [OMI_dict[dictkey][date]['avg'] for date in \
                    sorted(OMI_dict[dictkey].keys())]
                if(len(avgs)<2):
                    #print('Here 3')
                    colorRED = 0.
                    colorGREEN = 0.
                    colorBLUE = 0.
                else:
                    # Check the current max and min
                    x_vals = np.arange(0,len(OMI_dict[dictkey].keys()))
                    avgs = np.ma.masked_array([OMI_dict[dictkey][date]['avg'] for date in sorted(OMI_dict[dictkey].keys())])
                    #avgs = avgs[np.where(avgs.mask==False)[0]]
                    temp_dates = np.array(sorted(OMI_dict[dictkey].keys()))
                    dates = temp_dates
                    #dates = temp_dates[np.where(avgs.mask==False)[0]]
                    x_vals = np.arange(0,len(dates))
    
                    nx=len(avgs)
                    num_d=int(nx*(nx-1)/2)  # The number of elements in d
                    Sn=np.zeros(num_d)
                    S=0
                    for ti in range(nx-1):
                        for tj in range(ti+1,nx):
                            S+=np.sign(avgs[tj]-avgs[ti])
                        # Endfor
                    # Endfor

                    # Find the unique values in the data
                    uniq = np.unique(avgs)
                    g = len(uniq)
                    if(nx==g):
                        # No ties
                        Vs = (nx*(nx-1.)*(2.*nx-5.))/18.
                    else:
                        tp = np.zeros(uniq.shape)
                        for si in range(g):
                            tp[si] = sum(avgs==uniq[si])
                        Vs = (nx*(nx-1.)*(2.*nx-5.)-np.sum(tp*(tp-1.)*(2.*tp-5.)))/18.
                    if (S > 0.): 
                        z=(S-1.)/np.sqrt(Vs)
                    elif (S < 0.): 
                        z=(S+1.)/np.sqrt(Vs)
                    else: 
                        z=0.
                    # Calculate the p value of the trend
                    pval=2*(1.-stats.norm.cdf(abs(z)))  # (two-side)
                    alpha=0.05
                    ##h = abs(z) > stats.norm.ppf(1.-alpha/2.)
                    ### Determine the trend type 
                    ##if(z>0) and h:
                    ##    trend='increasing'
                    ##elif(z<0) and h:
                    ##    trend='decreasing'
                    ##else:
                    ##    trend='no trend'

                    if(pval<min_pval):
                        min_pval=pval
                    elif(pval>max_pval):
                        max_pval=pval

                    #color = mapper.to_rgba(slope)
    #color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
                    color_index_value = 0
                    color_index = np.where(pval<color_range)[0]
                    if(len(color_index)>0):
                        color_index_value = color_index[-1]+1
                    if(S>0):
                        colorRED = pRED[color_index_value]
                        colorGREEN = pGREEN[color_index_value]
                        colorBLUE = pBLUE[color_index_value]
                    else:
                        colorRED = nRED[color_index_value]
                        colorGREEN = nGREEN[color_index_value]
                        colorBLUE = nBLUE[color_index_value]
    
                    if(np.isnan(pval) == False):
                        # Add value to analysis regions (if possible)
                        if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                            canada_avg+=pval
                            canada_counts+=1
                        elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                            siberia_avg+=pval
                            siberia_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                            eus_avg+=pval
                            eus_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                            wus_avg+=pval
                            wus_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                            ea_avg+=pval
                            ea_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                            wa_avg+=pval
                            wa_counts+=1
                        elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                            europe_avg+=pval
                            europe_counts+=1
    
                        print(dictkey,pval,color_index_value)
            # Find the y coordinates of the current LatxLon box
            y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
            # Find the x coordinates of the current LatxLon box
            x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
            # Convert x and y into map coordinates
            mx, my = m(x,y)
            mxy = zip(mx,my)
            pair1 = (mx[0],my[0])
            pair2 = (mx[1],my[1])
            pair3 = (mx[2],my[2])
            pair4 = (mx[3],my[3])
            #mxy = zip(mx,my)
            # Plot the box on the map using color
            colors = [colorRED,colorGREEN,colorBLUE]
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors)
            #poly = Polygon(mxy,facecolor=color2,edgecolor=color2)
            plt.gca().add_patch(poly)
    
    y = [0.,0.,0.,0.]
    # Find the x coordinates of the current LatxLon box
    x = [0.,0.,0.,0.]
    color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
    # Convert x and y into map coordinates
    mx, my = m(x,y)
    pair1 = (mx[0],my[0])
    pair2 = (mx[1],my[1])
    pair3 = (mx[2],my[2])
    pair4 = (mx[3],my[3])
    # Positive
    # Add legend for 0.3-0.2
    colors = [pRED[1],pGREEN[1],pBLUE[1]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.4-0.3')
    plt.gca().add_patch(poly)
    # Add legend for 0.3-0.2
    colors = [nRED[1],nGREEN[1],nBLUE[1]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.4-0.3')
    plt.gca().add_patch(poly)
    # Add legend for 0.2-0.15
    colors = [pRED[2],pGREEN[2],pBLUE[2]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.3-0.2')
    plt.gca().add_patch(poly)
    # Add legend for 0.15-0.1
    colors = [nRED[2],nGREEN[2],nBLUE[2]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.3-0.2')
    plt.gca().add_patch(poly)
    # Add legend for 0.1-0.05
    colors = [pRED[3],pGREEN[3],pBLUE[3]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.2-0.15')
    plt.gca().add_patch(poly)
    # Add legend for <0.05
    colors = [nRED[3],nGREEN[3],nBLUE[3]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.2-0.15')
    plt.gca().add_patch(poly)
    colors = [pRED[4],pGREEN[4],pBLUE[4]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.15-0.1')
    plt.gca().add_patch(poly)
    # Negative
    # Add legend for 0.2-0.15
    colors = [nRED[4],nGREEN[4],nBLUE[4]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.15-0.1')
    plt.gca().add_patch(poly)
    # Add legend for 0.15-0.1
    colors = [pRED[5],pGREEN[5],pBLUE[5]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.1-0.05')
    plt.gca().add_patch(poly)
    # Add legend for 0.1-0.05
    colors = [nRED[5],nGREEN[5],nBLUE[5]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.1-0.05')
    plt.gca().add_patch(poly)
    # Add legend for <0.05
    colors = [pRED[6],pGREEN[6],pBLUE[6]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ <0.05')
    plt.gca().add_patch(poly)
    colors = [nRED[6],nGREEN[6],nBLUE[6]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- <0.05')
    plt.gca().add_patch(poly)

    box = ax.get_position()
    ax.set_position([box.x0,box.y0+box.height*0.1,box.width,box.height*0.9])
    ax.legend(loc='lower center',prop={'size':5},bbox_to_anchor=(0.48,-0.1),ncol=6)

    minmax_diff = (max_pval-min_pval) 
    minmax_range = np.arange(min_pval,max_pval,minmax_diff/6.0)
    print("min pval = ",min_pval)
    print("max pval = ",max_pval)
    print("Range = ",minmax_range)

    if(canada_counts>0):
        canada_avg = canada_avg/canada_counts
        print("Canada       avg = ",canada_avg,"  counts = ",canada_counts)
    if(siberia_counts>0):
        siberia_avg = siberia_avg/siberia_counts
        print("Siberia      avg = ",siberia_avg,"  counts = ",siberia_counts)
    if(eus_counts>0):
        eus_avg = eus_avg/eus_counts
        print("Eastern US   avg = ",eus_avg,"  counts = ",eus_counts)
    if(wus_counts>0):
        wus_avg = wus_avg/wus_counts
        print("Western US   avg = ",wus_avg,"  counts = ",wus_counts)
    if(ea_counts>0):
        ea_avg = ea_avg/ea_counts
        print("Eastern Asia avg = ",ea_avg,"  counts = ",ea_counts)
    if(wa_counts>0):
        wa_avg = wa_avg/wa_counts
        print("Western Asia avg = ",wa_avg,"  counts = ",wa_counts)
    if(europe_counts>0):
        europe_avg = europe_avg/europe_counts
        print("Europe       avg = ",europe_avg,"  counts = ",europe_counts)

    #start_date = min_date.decode("utf-8")
    #end_date = max_date.decode("utf-8")
    title_string = 'OMI Average Ultraviolet Aerosol Index Trend P-Values\n'+\
                   start_date+' to '+end_date
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    elif(summer is True):
        title_string = title_string+'\nJune, July, August'
    elif(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    elif(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)
    # Set up the colorbar
    #cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    #cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    #color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
    
    plt.xticks(rotation=45,fontsize=6)
    #plt.legend(prop={'size':10})
    fig1.canvas.mpl_connect('button_press_event',onclick)
    if(save is True):
        if(spring is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_spring_newData.png'
        elif(summer is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_summer_newData.png'
        elif(autumn is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_autumn_newData.png'
        elif(winter is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_winter_newData.png'
        else:
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_whole_year_newData.png'
        
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

def plotOMI(OMI_dict,start_date,end_date,save=False,trend_type='standard',file_type='XR123',season='',minlat=30.):
    if(file_type=='NXAR'):
        title_flabel = 'No XTrack and All Rows'
        outname_flabel = 'noX_allRows'
    elif(file_type=='XAR'):
        title_flabel = 'XTrack and All Rows'
        outname_flabel = 'xtrack_allRows'
    elif(file_type=='XR123'):
        title_flabel = 'XTrack and Rows 1-23'
        outname_flabel = 'xtrack_rows1to23'
    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='_thielSen'
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False

    # If only summer months are being analyzed, remove all data except 
    # in summer
    if(season=='spring'):
        spring = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    OMI_dict[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                    OMI_dict[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    OMI_dict[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    OMI_dict[lkey].pop(tkey)


    # Find the lowest lat in the file
    #lowest_lat = 50.

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    # Set up the polar stereographic map
    fig1 = plt.figure()
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,         \
                resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    

    cmap = plt.get_cmap('bwr')
    #bcmap = plt.cm.set_cmap(cmap)
    if(summer is True):
        v_min = -0.350 # Summer values
        mid_val = 0
        v_max = 0.900
    else:
        v_min = -0.350  # Whole-year values
        mid_val = 0
        v_max = 0.900

    if(minlat>30):
        v_max = 0.6
        mid_val = 0
        v_min = -0.6

    # Center the colorbar on zero
    norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
    #norm = Normalize(vmin=vmin,vmax=vmax)
    mapper = ScalarMappable(norm=norm,cmap=cmap)

    # Set up initial values for the analysis regions
    canada_avg = 0.
    canada_counts = 0
    siberia_avg = 0.
    siberia_counts = 0 
    eus_avg = 0.
    eus_counts = 0
    wus_avg = 0.
    wus_counts = 0
    ea_avg = 0.
    ea_counts = 0
    wa_avg = 0.
    wa_counts = 0
    europe_avg = 0.
    europe_counts = 0


    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    max_slope = -10.
    min_slope = 10.
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = (str(int(lat_ranges[i]))+'x'+                            \
                       str(int(lon_ranges[j])))

            # If no data are present for the curent lat/lon box, fill it with
            # black
            keylist = [ky for ky in OMI_dict.keys()]
            if(dictkey not in OMI_dict.keys()):
                color=(0,0,0,0)
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(OMI_dict[dictkey].keys())==0):
                color=(0,0,0,0)
            else:
                avgs = [OMI_dict[dictkey][date]['avg'] for date in \
                    sorted(OMI_dict[dictkey].keys())]
                # Check the current max and min
                x_vals = np.arange(0,len(OMI_dict[dictkey].keys()))
                temp_dates = np.array(sorted(OMI_dict[dictkey].keys()))
                # Find the slope of the line of best fit for the time series of
                # average data
                if(trend_type=='standard'): 
                    slope, intercept, r_value, p_value, std_err = \
                        stats.linregress(x_vals,avgs)
                    slope *= len(x_vals)
                else:
                    #The slope
                    S=0
                    sm=0
                    nx = len(avgs)
                    num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
                    Sn=np.zeros(num_d)
                    for si in range(0,nx-1):
                        for sj in range(si+1,nx):
                            # Find the slope between the two points
                            Sn[sm] = (avgs[si]-avgs[sj])/(si-sj) 
                            sm=sm+1
                        # Endfor
                    # Endfor
                    Snsorted=sorted(Sn)
                    sm=int(num_d/2.)
                    print(dictkey,len(Snsorted))
                    if(len(Snsorted)==1):
                        color=(0,0,0,0) 
                    else:
                        if(2*sm    == num_d):
                            slope=0.5*(Snsorted[sm]+Snsorted[sm+1])
                        if((2*sm)+1 == num_d): 
                            slope=Snsorted[sm+1]
                        slope = slope*len(avgs)

                if(slope>max_slope):
                    max_slope=slope
                elif(slope<min_slope):
                    min_slope=slope
        
                color = mapper.to_rgba(slope)

                # Add value to analysis regions (if possible)
                if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                    canada_avg+=slope
                    canada_counts+=1
                elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                    siberia_avg+=slope
                    siberia_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                    eus_avg+=slope
                    eus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                    wus_avg+=slope
                    wus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                    ea_avg+=slope
                    ea_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                    wa_avg+=slope
                    wa_counts+=1
                elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                    europe_avg+=slope
                    europe_counts+=1
                


            # Find the y coordinates of the current LatxLon box
            y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
            # Find the x coordinates of the current LatxLon box
            x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
            # Convert x and y into map coordinates
            mx, my = m(x,y)
            mxy = zip(mx,my)
            pair1 = (mx[0],my[0])
            pair2 = (mx[1],my[1])
            pair3 = (mx[2],my[2])
            pair4 = (mx[3],my[3])
            color2 = rgb2hex(color)
            #color2 = rgb2hex(color[:2]+color[3:])
            # Plot the box on the map using color
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=color2,         \
                edgecolor=color2)
            plt.gca().add_patch(poly)
    
    #start_date = str(temp_dates[0].decode("utf-8"))
    #end_date = str(temp_dates[-1].decode("utf-8"))
    minmax_diff = (max_slope-min_slope) 
    minmax_range = np.arange(min_slope,max_slope,minmax_diff/8.0)
    print("min slope = ",min_slope)
    print("max slope = ",max_slope)
    print("Range = ",minmax_range)

    if(canada_counts>0):
        canada_avg = canada_avg/canada_counts
        print("Canada       avg = ",canada_avg,"  counts = ",canada_counts)
    if(siberia_counts>0):
        siberia_avg = siberia_avg/siberia_counts
        print("Siberia      avg = ",siberia_avg,"  counts = ",siberia_counts)
    if(eus_counts>0):
        eus_avg = eus_avg/eus_counts
        print("Eastern US   avg = ",eus_avg,"  counts = ",eus_counts)
    if(wus_counts>0):
        wus_avg = wus_avg/wus_counts
        print("Western US   avg = ",wus_avg,"  counts = ",wus_counts)
    if(ea_counts>0):
        ea_avg = ea_avg/ea_counts
        print("Eastern Asia avg = ",ea_avg,"  counts = ",ea_counts)
    if(wa_counts>0):
        wa_avg = wa_avg/wa_counts
        print("Western Asia avg = ",wa_avg,"  counts = ",wa_counts)
    if(europe_counts>0):
        europe_avg = europe_avg/europe_counts
        print("Europe       avg = ",europe_avg,"  counts = ",europe_counts)

    title_string = 'Change in OMI Ultraviolet Average Aerosol Index\n'+        \
        title_flabel+'\n'+str(start_date)+' to '+str(end_date)
        ##start_date+' to '+end_date+'\nXTrack and all rows'
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    if(summer is True):
        title_string = title_string+'\nJune, July, August'
    if(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    if(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)
    # Set up the colorbar
    cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    plt.xticks(rotation=45,fontsize=6)
    fig1.canvas.mpl_connect('button_press_event',onclick)
    if(save is True):
        if(spring is True):
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_spring_'+outname_flabel+'_newData.png'
        elif(summer is True):     
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_summer_'+outname_flabel+'_newData.png'
        elif(autumn is True):     
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_autumn_'+outname_flabel+'_newData.png'
        elif(winter is True):    
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_winter_'+outname_flabel+'_newData.png'
        else:                   
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_whole_year_'+outname_flabel+'_newData.png'
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

def plotOMI_Climo(OMI_dict,start_date,end_date,save=False,trend_type='standard',file_type='XR123',season='',minlat=30.):
    if(file_type=='NXAR'):
        title_flabel = 'No XTrack and All Rows'
        outname_flabel = 'noX_allRows'
    elif(file_type=='XAR'):
        title_flabel = 'XTrack and All Rows'
        outname_flabel = 'xtrack_allRows'
    elif(file_type=='XR123'):
        title_flabel = 'XTrack and Rows 1-23'
        outname_flabel = 'xtrack_rows1to23'
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False

    # If only summer months are being analyzed, remove all data except 
    # in summer
    if(season=='spring'):
        spring = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    OMI_dict[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                    OMI_dict[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    OMI_dict[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(OMI_dict.keys()):
            for tkey in sorted(OMI_dict[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    OMI_dict[lkey].pop(tkey)


    # Find the lowest lat in the file

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    # Set up the polar stereographic map
    fig1 = plt.figure()
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,         \
                resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    

    cmap = plt.get_cmap('jet')
    #bcmap = plt.cm.set_cmap(cmap)
    if(summer is True):
        v_min = -0.350 # Summer values
        mid_val = 0
        v_max = 0.900
    else:
        v_min = -0.350  # Whole-year values
        mid_val = 0
        v_max = 0.900
    v_min = -1.000  # Whole-year values
    mid_val = 0
    v_max = 1.500
    if(minlat>30.):
        v_min = 0.0
        mid_val = 0
        v_max = 1.00

    # Center the colorbar on zero
    #norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
    norm = Normalize(vmin=v_min,vmax=v_max)
    mapper = ScalarMappable(norm=norm,cmap=cmap)

    # Set up initial values for the analysis regions
    canada_avg = 0.
    canada_counts = 0
    siberia_avg = 0.
    siberia_counts = 0 
    eus_avg = 0.
    eus_counts = 0
    wus_avg = 0.
    wus_counts = 0
    ea_avg = 0.
    ea_counts = 0
    wa_avg = 0.
    wa_counts = 0
    europe_avg = 0.
    europe_counts = 0

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    max_avg_uvai = -10.
    min_avg_uvai = 10.
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = (str(int(lat_ranges[i]))+'x'+                            \
                       str(int(lon_ranges[j])))

            # If no data are present for the curent lat/lon box, fill it with
            # black
            keylist = [ky for ky in OMI_dict.keys()]
            if(dictkey not in OMI_dict.keys()):
                color=(0,0,0,0)
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(OMI_dict[dictkey].keys())==0):
                color=(0,0,0,0)
            else:
                avgs = np.array([OMI_dict[dictkey][date]['avg'] for date in \
                    sorted(OMI_dict[dictkey].keys())])
                counts = np.array([OMI_dict[dictkey][date]['#_obs'] for date in \
                    sorted(OMI_dict[dictkey].keys())])
                avg_uvai = sum(avgs*counts)/sum(counts)
                #avg_uvai = np.average(avgs)
                temp_dates = np.array(sorted(OMI_dict[dictkey].keys()))
                # Check the current max and min
                #x_vals = np.arange(0,len(OMI_dict[dictkey].keys()))
                # Find the slope of the line of best fit for the time series of
                # average data
                #slope, intercept, r_value, p_value, std_err = \
                #    stats.linregress(x_vals,avgs)
                #slope *= len(x_vals)

                if(avg_uvai>max_avg_uvai):
                    max_avg_uvai=avg_uvai
                elif(avg_uvai<min_avg_uvai):
                    min_avg_uvai=avg_uvai
        
                color = mapper.to_rgba(avg_uvai)

                # Add value to analysis regions (if possible)
                if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                    canada_avg+=avg_uvai
                    canada_counts+=1
                elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                    siberia_avg+=avg_uvai
                    siberia_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                    eus_avg+=avg_uvai
                    eus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                    wus_avg+=avg_uvai
                    wus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                    ea_avg+=avg_uvai
                    ea_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                    wa_avg+=avg_uvai
                    wa_counts+=1
                elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                    europe_avg+=avg_uvai
                    europe_counts+=1

            # Find the y coordinates of the current LatxLon box
            y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
            # Find the x coordinates of the current LatxLon box
            x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
            # Convert x and y into map coordinates
            mx, my = m(x,y)
            mxy = zip(mx,my)
            pair1 = (mx[0],my[0])
            pair2 = (mx[1],my[1])
            pair3 = (mx[2],my[2])
            pair4 = (mx[3],my[3])
            color2 = rgb2hex(color)
            #color2 = rgb2hex(color[:2]+color[3:])
            # Plot the box on the map using color
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=color2,         \
                edgecolor=color2)
            plt.gca().add_patch(poly)
    
    #start_date = str(temp_dates[0].decode("utf-8"))
    #end_date = str(temp_dates[-1].decode("utf-8"))
    minmax_diff = (max_avg_uvai-min_avg_uvai) 
    minmax_range = np.arange(min_avg_uvai,max_avg_uvai,minmax_diff/8.0)
    print("min avg_uvai = ",min_avg_uvai)
    print("max avg_uvai = ",max_avg_uvai)
    print("Range = ",minmax_range)

    if(canada_counts>0):
        canada_avg = canada_avg/canada_counts
        print("Canada       avg = ",canada_avg,"  counts = ",canada_counts)
    if(siberia_counts>0):
        siberia_avg = siberia_avg/siberia_counts
        print("Siberia      avg = ",siberia_avg,"  counts = ",siberia_counts)
    if(eus_counts>0):
        eus_avg = eus_avg/eus_counts
        print("Eastern US   avg = ",eus_avg,"  counts = ",eus_counts)
    if(wus_counts>0):
        wus_avg = wus_avg/wus_counts
        print("Western US   avg = ",wus_avg,"  counts = ",wus_counts)
    if(ea_counts>0):
        ea_avg = ea_avg/ea_counts
        print("Eastern Asia avg = ",ea_avg,"  counts = ",ea_counts)
    if(wa_counts>0):
        wa_avg = wa_avg/wa_counts
        print("Western Asia avg = ",wa_avg,"  counts = ",wa_counts)
    if(europe_counts>0):
        europe_avg = europe_avg/europe_counts
        print("Europe       avg = ",europe_avg,"  counts = ",europe_counts)

    start_date = str(start_date)
    end_date = str(end_date)
    title_string = 'OMI Ultraviolet Average Aerosol Index Climatology\n'+      \
        title_flabel+'\n'+start_date+' to '+end_date
        ##start_date+' to '+end_date+'\nXTrack and All Rows'
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    if(summer is True):
        title_string = title_string+'\nJune, July, August'
    if(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    if(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)
    # Set up the colorbar
    cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    plt.xticks(rotation=45,fontsize=6)
    plt.xlabel('Ultraviolet Aerosol Index',fontsize=6)
    fig1.canvas.mpl_connect('button_press_event',onclick_climo)
    if(save is True):
        if(spring is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_spring_'+outname_flabel+'_newData.png'
        elif(summer is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_summer_'+outname_flabel+'_newData.png'
        elif(autumn is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_autumn_'+outname_flabel+'_newData.png'
        elif(winter is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_winter_'+outname_flabel+'_newData.png'
        else:
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_whole_year_'+outname_flabel+'_newData.png'
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

