#Program to average CATS vertical profiles
#Written by Logan Lee May 2018 Modified May 2020
#Modified by J.Zhang to output aerosol profiles

import numpy
import os
#import matplotlib.pylab as plt
#import matplotlib
import math
from numpy import genfromtxt
import h5py
import fnmatch
import datetime
import itertools
from numpy import meshgrid
#from mpl_toolkits.basemap import Basemap
import glob


pass_count=numpy.zeros((200,4))
all_ext=numpy.zeros((200,4))
all_ext2=numpy.zeros((200,4))
all_extpoints=numpy.zeros((200,4))

#Box Boundaries
upboundlat=0+90
lowboundlat=-20+90
eastboundlon=80+180
westboundlon=60+180
def myround(x, base=5):
    return int(base * round(float(x)/base))
def myround2(x,base=2):
    return int(base * round(float(x)/base))
def myround3(x,base=6):
    return int(base*round(float(x)/base))

#Read in data files
cats_init=glob.glob('/data/CATS_V3/L2_Profiles/*/*.hdf5')
cats_list=[f for f in cats_init if (fnmatch.fnmatch(f,'*/CATS-*2015-12*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2015-01*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2015-02*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2016-12*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2016-01*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2016-02*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2017-12*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2017-01*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2017-02*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2015-03*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2015-04*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2015-05*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2016-03*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2016-04*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2016-05*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2017-03*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2017-04*.hdf5') | fnmatch.fnmatch(f,'*/CATS-*2017-05*.hdf5'))]


for r in range(0,len(cats_list)):   #Loop through list of CATS files
                all_sum=numpy.zeros((200,4))
                all_count=numpy.zeros((200,4))
		#Get variable data
                filename2=cats_list[r]
                print(filename2)
                file=h5py.File(filename2,'r')
                dataset=file['/geolocation/CATS_Fore_FOV_Latitude']
                catlat=dataset[:]
                dataset2=file['/geolocation/CATS_Fore_FOV_Longitude']
                catlon=dataset2[:]
                dataset6=file['/profile/Aerosol_Optical_Depth_1064_Fore_FOV']
                aod=dataset6[:]
                dataset4=file['/profile/Profile_UTC_Time']
                time=dataset4[:]
                dataset133=file['/profile/Profile_UTC_Date']
                date=dataset133[:]
                dataset8=file['/profile/Feature_Type_Fore_FOV']
                feattype=dataset8[:]
                dataset5=file['/profile/Sky_Condition_Fore_FOV']
                sky=dataset5[:]
                dataset7=file['/profile/Extinction_QC_Flag_1064_Fore_FOV']
                qc=dataset7[:]
                dataset10=file['/profile/Extinction_Coefficient_Uncertainty_1064_Fore_FOV']
                extuncat=dataset10[:]
                dataset9=file['/profile/Feature_Type_Score_Fore_FOV']
                featscore=dataset9[:]
                dataset11=file['/profile/Day_Night_Flag']
                dnfcat=dataset11[:]
                dataset13=file['/profile/Percent_Opacity_Fore_FOV']
                opac=dataset13[:]
                dataset17=file['/profile/DEM_Surface_Altitude_Fore_FOV']
                demsfc=dataset17[:]
                dataset18=file['/profile/Extinction_Coefficient_1064_Fore_FOV']
                ext=dataset18[:]
                dataset32=file['/metadata_parameters/Bin_Altitude_Array']
                bin_alt=dataset32[:]
                file.close()
		#lon,lat, and time have , we want mean value
                catlon=catlon[:,1]
                catlat=catlat[:,1]
                time=time[:,1]

		#Construct julian day
                date2=date.astype('S8')
                fmt='%Y%m%d'
                dt=[datetime.datetime.strptime(datey,fmt).strftime('%j') for datey in date2]
                day=numpy.array(dt).astype(float)
                tyme=numpy.array(time).astype(float)
                juldayc=day+time
                yearcats=numpy.array([datetime.datetime.strptime(datey,fmt).strftime('%Y') for datey in date2]).astype(float)
	

		#Variables for QA procedure
                featscore2=numpy.copy(featscore).astype(float)  #Valid range is from -10 to 10; -10 is confidently aerosol, 10 is confidently cloud, 0 could be either.
                featscore2[featscore2 < -10]=numpy.nan   #Don't want undefined values
                catcadmax=numpy.nanmax(featscore2,axis=0)  #Max featscore
                catcadmin=numpy.nanmin(featscore2,axis=0)  #Min featscore
                catextunmax=numpy.nanmax(extuncat,axis=0)  #Max Extinction Uncertainty
                catextunmin=numpy.nanmin(extuncat,axis=0)  #Min Extintion Uncertainty
                catextmax=numpy.nanmax(ext,axis=0)  #Max Extinction
                featscoremax=numpy.nanmax(featscore,axis=0)


		#CATS QC
		#Feature Type Values 0=, 1=, 2=, 3=, 4=
                clouds=numpy.unique(numpy.where(feattype == 1)[1])  #Get the x-values of cloudy pixels (not the heights)
                cloud_param=numpy.zeros(aod.shape)
                cloud_param[clouds]=1   #clouds present = 1, no clouds = 0
                        #Excluding profiles with underfined feature types...	
                unds=numpy.unique(numpy.where(feattype == 2)[1])
                unds_param=numpy.zeros(aod.shape)
                unds_param[unds]=1  #undefined types present = 1, not present = 0
                qcmax=numpy.nanmax(qc,axis=0)  #QC Flag values, 
                latp=numpy.copy(catlat)+90     #Making lat/lon all positive for simplicity
                lonp=numpy.copy(catlon)+180
                lonp[lonp == 360]=0  #Where lon is 360, it is also 0, so choose one
		#Valid Values based on CATS QC
                valid5=(catextmax < 1.25) & (qcmax == 0) & (sky == 1) & (cloud_param < 1) & (unds_param < 1) & (catextunmax < 10) & (catcadmax < -1) & (catcadmin > -11) & (dnfcat != 1) & ((latp < upboundlat) & (latp > lowboundlat) & (lonp < eastboundlon) & (lonp > westboundlon))	
                aod[~valid5]=-999.0
		
		#AGL CORRECT VERTICAL EXTINCTION PROFILES
                if (aod[valid5].shape[0] > 0): 
                        print(aod[valid5].shape)
                        catalt2=numpy.matrix(numpy.copy(bin_alt)).transpose()
                        cat3D=numpy.repeat(catalt2,len(demsfc),axis=1)  #3D grid of bin altitude values
                        indcat=((cat3D-demsfc) >= 0.0) #Subtract surface elevation from bin altitude, where is this 0 or greater (i.e. above ground level)
                        temp_hgt=numpy.matrix(range(0,533))
                        yvals2=numpy.repeat(numpy.transpose(temp_hgt),len(demsfc),axis=1).astype(numpy.float32)  #Grid of y values
                        yvals2[~indcat]=-1000 #Set all y-values that are below ground level to -1000
                        sort_indices2=numpy.argsort(yvals2,axis=0)
                        static_indices2=numpy.indices(indcat.shape)
                        ext_temp=numpy.copy(ext)
                        ext_temp[~indcat]=-999.		
                        catext2=numpy.copy(ext_temp[sort_indices2,static_indices2[1]]) #Move all those invalid below ground y-values to the top of the grid so that the bottom starts at ground level


			#Now, average the valid vertical profiles for this CATS data file	
                        for k in range(0,aod.shape[0]):  #Loop through vertical profiles in the data file; You could probably do this faster in some sort of array version of this than looping through each individual valid profile
                            if ((aod[k] >= 0)):  #i.e. if this profile is valid
                                timepp=myround3(tyme[k]*24.)/6.  #Round time to nearest 6 hours (0, 6, 12, 18, 24)
                                if timepp==(4.):
                                    timepp=0  #Time = 24 is the same as time = 0 (24 hours in a day)
                                    temp_ext=numpy.copy(catext2[332:532,k])  #Get valid vertical extinction profile (Note, if you wanted the full profile, just use catext2[:,k])
                                    good_count=(temp_ext > -900.).astype(int)  #All valid points in that profile
                                    temp_ext[temp_ext < -900.]=0  #Set the invalid points to 0
                                    all_sum[:,timepp]=all_sum[:,timepp]+temp_ext  #Running sum of the profile values in data fil
                                    all_count[:,timepp]=all_count[:,timepp]+good_count  #Running sum of the number of valid points in data file

                        #OK, so we now have the sum and valid count of all profiles in this data file (i.e. in this "pass" of CATS). Compute the pass average.
                        temp_count=numpy.zeros(all_count.shape)   
                        temp_count[all_count > 0]=1   #Pass Count (1, because this is 1 data file)
                        pass_count=pass_count+temp_count  #Add to previous number of pass counts to sum how many passes we've had so far
                        all_mean=numpy.zeros(all_sum.shape)
                        goodext=(all_count > 0)  
                        all_mean[goodext]=all_sum[goodext]/all_count[goodext]  #Mean for this data file( i.e. for this pass) 
                        all_ext=all_ext+all_mean   #Running Sum of the Pass Means
                        all_ext2=all_ext2+all_sum  #Running Sum of all the valid values
                        all_extpoints=all_extpoints+all_count  #Running total number of valid points
bestext=(pass_count > 0)
all_meanext=numpy.zeros(pass_count.shape)
all_meanext[bestext]=all_ext[bestext]/pass_count[bestext]   #Mean of the Pass Means after all data files have been looped through

numpy.savetxt('may2020testfixed_loc_profsix_spring.txt',all_meanext.reshape(200,4))   #Pass Mean
numpy.savetxt('may2020testfixed_loc_profsixpasscount_spring.txt',pass_count.reshape(200,4))  #Number of Passes
numpy.savetxt('may2020testfixed_loc_profsixsum_spring.txt',all_ext2.reshape(200,4))  #Sum of all values
numpy.savetxt('may2020testfixed_loc_profsixpointcount_spring.txt',all_extpoints.reshape(200,4))  #Total number of points
numpy.savetxt('may2020testfixed_loc_profsix_testzero_spring.txt',all_meanext[:,0])  #This should be timepp = 0, so you can double check which order the arrays are in, just compare to the all_meanext.reshape(200,4) one to find out 

##So, to summarize: Pass Mean (all_meanext) = sum of means that were calculated at each pass/total number of passes
##You can also calculate a point mean as: Total sum of all values (all_ext2)/total number of points (all_extpoints)


