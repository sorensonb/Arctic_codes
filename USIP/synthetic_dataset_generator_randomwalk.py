#!/usr/bin/env python

"""
  NAME:
    synthetic_dataset_generator_randomwalk.py
  
  PURPOSE:
  
  MODIFICATIONS:
  
"""

# Import functions
import numpy as np
import sys
import os
import random
import warnings
import matplotlib.pyplot as plt
from commands import getstatusoutput
from scipy import interpolate
from time import strptime
from calendar import month_name
from readsounding import *
from compare_cn2 import *

if(len(sys.argv)<3):
    print "SYNTAX: ./synthetic_dataset_generator.py time monthly_sounding_archive_file [meter]"
    print "        time - 00 or 12"
    print ""
    print "        OPTIONAL ARGUMENTS:"
    print "        meter - If meter is added as the last argument, the code"
    print "                will calculate average and standard deviations over"
    print "                500 meter windows instead of pressure levels"
    sys.exit()

METER=False
if(len(sys.argv)==4):
    if(sys.argv[3]=='meter'):
        METER=True

work_path='/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/'
# Add check for compare_cn2.py, plot_cn2.py, and readsounding.py
#if (os.path.isfile(work_path+'sounding_splitter.py') is False):
#    print "ERROR: "+work_path+"sounding_splitter.py not found"
#    print "       Please ensure the necessary program is in the correct location"
#    sys.exit(1)

# Extract the individual soundings from the archive file
infile = sys.argv[2]
# Grab the month information from the monthly archive file name
date = infile.split('/')[-1].split('_')
month = date[0]
month_no = strptime(month[:3],'%b').tm_mon
str_monthno = str(month_no)
if(month_no<10):
    str_monthno = '0'+str_monthno
tmonth = month_name[month_no]
year = date[1]
location = date[2]
file_path = '/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/sounding_archive/'+year+'/'+location+'/soundings/'
archive_path = '/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/sounding_archive/'+year+'/'+location+'/'
# If the files are already extracted and located in 
#   /nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/2017/***/soundings/,
# ignore this step
status,files_there = getstatusoutput('ls '+file_path+year[2:]+str_monthno+'*.txt')
files_there = files_there.split('\n')
if(len(files_there)==1):
    print "Extracting soundings from archive file"
    os.system(work_path+'sounding_splitter.py '+archive_path+infile.strip().split('/')[-1])
else:
    print "Soundings already extracted from archive file. Proceeding"
stime = sys.argv[1]

# Get the file names from the recently split archive file
status,output = getstatusoutput('ls '+file_path+year[2:]+str_monthno+'*_'+stime+'0000_'+location+'_UWYO.txt') 
files = output.split('\n')
list.sort(files)

# Calculate Cn2 for each sounding and calculate the average and standard deviation
# for each sounding
avgs = []
stdevs = []
cn2_sfcto925 = []
cn2_925to850 = []
cn2_850to700 = []
cn2_700to500 = []
cn2_500to400 = []
cn2_400to300 = []
cn2_300to200 = []
cn2_200to150 = []
cn2_150to100 = []
cn2_100to50 = []
cn2_50totop = []

temp_sfcto925 = []
temp_925to850 = []
temp_850to700 = []
temp_700to500 = []
temp_500to400 = []
temp_400to300 = []
temp_300to200 = []
temp_200to150 = []
temp_150to100 = []
temp_100to50 = []
temp_50totop = []

height_sfc = []
height_925 = []
height_850 = []
height_700 = []
height_500 = []
height_400 = []
height_300 = []
height_200 = []
height_150 = []
height_100 = []
height_50 = []
height_top = []

sfc_press = []
top_press = []
j=0
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for f in files:
        data = readsounding(f)
        data = cn2_calc(data) 
    #    log_cn2 = np.log10(data['CN2'])
        cn2_sfcto925.append(np.average(data['CN2'][np.where(data['PRESS']>925)]))
        cn2_925to850.append(np.average(data['CN2'][np.where((data['PRESS']<=925) & (data['PRESS']>850))]))
        cn2_850to700.append(np.average(data['CN2'][np.where((data['PRESS']<=850) & (data['PRESS']>700))]))
        cn2_700to500.append(np.average(data['CN2'][np.where((data['PRESS']<=700) & (data['PRESS']>500))]))
        cn2_500to400.append(np.average(data['CN2'][np.where((data['PRESS']<=500) & (data['PRESS']>400))]))
        cn2_400to300.append(np.average(data['CN2'][np.where((data['PRESS']<=400) & (data['PRESS']>300))]))
        cn2_300to200.append(np.average(data['CN2'][np.where((data['PRESS']<=300) & (data['PRESS']>200))]))
        cn2_200to150.append(np.average(data['CN2'][np.where((data['PRESS']<=200) & (data['PRESS']>150))]))
        cn2_150to100.append(np.average(data['CN2'][np.where((data['PRESS']<=150) & (data['PRESS']>100))]))
        cn2_100to50.append(np.average(data['CN2'][np.where((data['PRESS']<=100) & (data['PRESS']>50))]))
        cn2_50totop.append(np.average(data['CN2'][np.where(data['PRESS']<=50)]))

        temp_sfcto925.append(np.average(data['TEMP'][np.where(data['PRESS']>925)]))
        temp_925to850.append(np.average(data['TEMP'][np.where((data['PRESS']<=925) & (data['PRESS']>850))]))
        temp_850to700.append(np.average(data['TEMP'][np.where((data['PRESS']<=850) & (data['PRESS']>700))]))
        temp_700to500.append(np.average(data['TEMP'][np.where((data['PRESS']<=700) & (data['PRESS']>500))]))
        temp_500to400.append(np.average(data['TEMP'][np.where((data['PRESS']<=500) & (data['PRESS']>400))]))
        temp_400to300.append(np.average(data['TEMP'][np.where((data['PRESS']<=400) & (data['PRESS']>300))]))
        temp_300to200.append(np.average(data['TEMP'][np.where((data['PRESS']<=300) & (data['PRESS']>200))]))
        temp_200to150.append(np.average(data['TEMP'][np.where((data['PRESS']<=200) & (data['PRESS']>150))]))
        temp_150to100.append(np.average(data['TEMP'][np.where((data['PRESS']<=150) & (data['PRESS']>100))]))
        temp_100to50.append(np.average(data['TEMP'][np.where((data['PRESS']<=100) & (data['PRESS']>50))]))
        temp_50totop.append(np.average(data['TEMP'][np.where(data['PRESS']<=50)]))

        height_sfc.append(data['ALT'][0]) 
        if(len(np.where(data['PRESS']==925)[0]) != 0):
            height_925.append(data['ALT'][np.where(data['PRESS']==925)])
        if(len(np.where(data['PRESS']==850)[0]) != 0):
            height_850.append(data['ALT'][np.where(data['PRESS']==850)])
        if(len(np.where(data['PRESS']==700)[0]) != 0):
            height_700.append(data['ALT'][np.where(data['PRESS']==700)])
        if(len(np.where(data['PRESS']==500)[0]) != 0):
            height_500.append(data['ALT'][np.where(data['PRESS']==500)])
        if(len(np.where(data['PRESS']==400)[0]) != 0):
            height_400.append(data['ALT'][np.where(data['PRESS']==400)])
        if(len(np.where(data['PRESS']==300)[0]) != 0):
            height_300.append(data['ALT'][np.where(data['PRESS']==300)])
        if(len(np.where(data['PRESS']==200)[0]) != 0):
            height_200.append(data['ALT'][np.where(data['PRESS']==200)])
        if(len(np.where(data['PRESS']==150)[0]) != 0):
            height_150.append(data['ALT'][np.where(data['PRESS']==150)])
        if(len(np.where(data['PRESS']==100)[0]) != 0):
            height_100.append(data['ALT'][np.where(data['PRESS']==100)])
        if(len(np.where(data['PRESS']==50)[0]) != 0):
            height_50.append(data['ALT'][np.where(data['PRESS']==50)])
        height_top.append(data['ALT'][-1])

        sfc_press.append(data['PRESS'][0])
        top_press.append(data['PRESS'][-1])

log_data_list = [np.log10(cn2_sfcto925),np.log10(cn2_925to850),\
    np.log10(cn2_850to700),np.log10(cn2_700to500),\
    np.log10(cn2_500to400),np.log10(cn2_400to300),\
    np.log10(cn2_300to200),np.log10(cn2_200to150),\
    np.log10(cn2_150to100),np.log10(cn2_100to50),\
    np.log10(cn2_50totop)]

temp_list = [np.array(temp_sfcto925),np.array(temp_925to850),\
    np.array(temp_850to700),np.array(temp_700to500),\
    np.array(temp_500to400),np.array(temp_400to300),\
    np.array(temp_300to200),np.array(temp_200to150),\
    np.array(temp_150to100),np.array(temp_100to50),\
    np.array(temp_50totop)]

avg_heights = [np.average(height_sfc),np.average(height_925),\
    np.average(height_850),np.average(height_700),np.average(height_500),\
    np.average(height_400),np.average(height_300),np.average(height_200),\
    np.average(height_150),np.average(height_100),np.average(height_50),\
    np.average(height_top)]
#avg_heights = [np.average(height_sfc),np.average(height_top)]

avg_sfc_press = np.average(sfc_press)
avg_top_press = np.average(top_press)

# Get rid of missing values
for i,s in enumerate(log_data_list):
    log_data_list[i] = log_data_list[i][np.logical_not(np.isnan(log_data_list[i]))]
for i,s in enumerate(temp_list):
    temp_list[i] = temp_list[i][np.logical_not(np.isnan(temp_list[i]))]

avg_cn2 = [np.average(x) for x in log_data_list]
std_cn2 = [np.std(y) for y in log_data_list]

avg_temp = [np.average(z) for z in temp_list]
std_temp = [np.std(a) for a in temp_list]

press_list = [avg_sfc_press,925.,850.,700.,500.,400.,300.,200.,150.,100.,50.,avg_top_press]

# Create the synthetic pressure and altitude profiles
interp_alt = [0]*12
interp_press = [0]*11
avg_heights = np.array(avg_heights)
for i in range(1,len(avg_heights)):
    talt = np.arange(avg_heights[i-1],avg_heights[i],5.0)
    tpress = interpolate.interp1d([talt[0],talt[-1]],[press_list[i-1],press_list[i]])
    ipress = tpress(talt)
    interp_press[i-1]=ipress
    interp_alt[i-1]=talt

# Concatenate the pressure and altitude data into single arrays
interp_alt = interp_alt[:-1]
for i in range(1,len(interp_alt)):
    interp_alt[i] = np.delete(interp_alt[i],0)
    interp_press[i] = np.delete(interp_press[i],0)
all_press = np.concatenate(interp_press)
all_alt = np.concatenate(interp_alt)

# Create synthetic Cn2 and temp datasets
# s_temp and s_cn2 are lists of lists that will be concatenated later
s_temp = [0]*11
s_cn2 = [0]*11
# Loop over each average and standard deviation and create synthetic cn2 and 
# temperature data for each layer. After the loop, s_temp and s_cn2 will
# contain eleven lists containing synthetic data for each layer
        # # # # # # # # # # # # # # 
        # RANDOM WALK STUFF
        # # # # # # # # # # # # # # 
###st1 = avg_temp[0]
###test = avg_temp[0]
###sc1 = avg_cn2[0]
        # # # # # # # # # # # # # # 
        # END RANDOM WALK STUFF
        # # # # # # # # # # # # # # 
for i in range(0,len(avg_temp)):
    ###s_temp[i] = [0]*len(interp_press[i])
    ###s_cn2[i]  = [0]*len(interp_press[i])
    ###s_temp[i][0] = st1
    ###s_cn2[i][0]  = sc1

    # Interpolate the cn2 average values above and below
    cn2_bottom = avg_cn2[i]
    if(i==(len(avg_temp)-1)):
        cn2_top = avg_cn2[i]
    else:
        cn2_top = avg_cn2[i+1]
    interp_cn2 = interpolate.interp1d([interp_press[i][0],interp_press[i][-1]],[cn2_bottom,cn2_top])
    cn2_interp = interp_cn2(interp_press[i])

    # Interpolate the cn2 average values above and below
    temp_bottom = avg_temp[i]
    if(i==(len(avg_temp)-1)):
        temp_top = avg_temp[i]
    else:
        temp_top = avg_temp[i+1]
    interp_temp = interpolate.interp1d([interp_press[i][0],interp_press[i][-1]],[temp_bottom,temp_top])
    temp_interp = interp_temp(interp_press[i])

    #for j in range(0,len(interp_press[i])):
        # Use the distance from the mean to determine how
        # much to change
        # # # # # # # # # # # # # # 
        # RANDOM WALK STUFF
        # # # # # # # # # # # # # # 
        ###t_distance_factor = -((s_temp[i][j]-avg_temp[i])/std_temp[i])*(1.0/len(interp_press[i]))
        ###c_distance_factor = -((s_cn2[i][j]-avg_cn2[i])/std_cn2[i])*(1.0/len(interp_press[i]))
        ###if(s_temp[i][j]==avg_temp[i]):
        ###    t_distance_factor = 1.
        ###    c_distance_factor = 1.

        #### Create the new random value
        ###r_temp = (2*std_temp[i]*random.random())-std_temp[i]
        ###r_cn2 = (2*std_cn2[i]*random.random())-std_cn2[i]

        ###if(test==avg_temp[i]): test += r_temp
        ###else: test += r_temp*((test-avg_temp[i])/std_temp[i])
        ###print avg_temp[i], std_temp[i], ((s_temp[i][j]-avg_temp[i])/std_temp[i]), test

        ###s_temp[i][j] = s_temp[i][j]+r_temp*t_distance_factor
        ###s_cn2[i][j] = s_cn2[i][j]+r_cn2*c_distance_factor

        #print s_temp[i][j],avg_temp[i],std_temp[i],t_distance_factor

        # # # # # # # # # # # # # # 
        # END RANDOM WALK STUFF
        # # # # # # # # # # # # # # 

    s_temp[i] = [temp_interp[j]+(2*std_temp[i]*random.random())-std_temp[i] for j in range(0,len(interp_press[i]))]
    s_cn2[i] = [cn2_interp[j]+(2*std_cn2[i]*random.random())-std_cn2[i] for j in range(0,len(interp_press[i]))]
        #s_temp[i] = [avg_temp[i]+(2*std_temp[i]*random.random())-std_temp[i] for j in range(0,len(interp_press[i]))]
        #s_cn2[i] = [avg_cn2[i]+(2*std_cn2[i]*random.random())-std_cn2[i] for j in range(0,len(interp_press[i]))]
        ###st1 = s_temp[i][j]
        ###sc1 = s_cn2[i][j]

# Concatenate the lists in s_temp and s_cn2 into single arrays
synthetic_temp = np.concatenate(s_temp)
synthetic_cn2 = np.concatenate(s_cn2)

# Convert the synthetic cn2 from logarithmic to linear units
linear_cn2 = 10.**synthetic_cn2

# Convert pressure to Pascals and temperature to Kelvin
t_all_press=all_press*100.
t_synthetic_temp=synthetic_temp+273.15
# Use the pressure, temp, and cn2 data to find a synthetic CT2 dataset
synthetic_ct2 = (linear_cn2/(((79.*10.**-8.)*(t_all_press/(t_synthetic_temp**2.)))**2.))**0.5

recalculated_cn2 = cn2_calc_thermo(synthetic_ct2,all_press,synthetic_temp)
rcn2 = np.log10(recalculated_cn2['CN2T'])


# Smooth the scn2, rcn2, and tempdiff data
smooth_scn2 = np.array([])
smooth_rcn2 = np.array([])
smooth_tempdiff = np.array([])
total_s = synthetic_cn2[0]
total_r = rcn2[0]
total_t = synthetic_ct2[0]
mid_alt = np.array([])
mid_press = np.array([])
for i in range(1,len(synthetic_cn2)):
    if(i%11==0.0):
        smooth_scn2 = np.append(smooth_scn2[:],total_s/11.0)
        smooth_rcn2 = np.append(smooth_rcn2[:],total_r/11.0)
        smooth_tempdiff = np.append(smooth_tempdiff[:],total_t/11.0)
        mid_alt = np.append(mid_alt[:],all_alt[i-6])
        mid_press = np.append(mid_press[:],all_press[i-6])
        total_s=0
        total_r=0
        total_t=0

    # Add the current data to the totals
    total_s+=synthetic_cn2[i]
    total_r+=rcn2[i]
    total_t+=synthetic_ct2[i]

# Find the year and month from the archive file name
date = infile.split('/')[-1].split('.')[0].split('_')
month = date[0]
month_no = strptime(month[:3],'%b').tm_mon
tmonth = month_name[month_no]
year = date[1]
location = date[2]

# Make plots
# Make plot of all synthetic and smoothed Cn2
fig1 = plt.figure()
plt.plot(synthetic_cn2,all_alt/1000.,color='black',label='synthetic $C_n^2$')
plt.plot(smooth_scn2,mid_alt/1000.,color='red',label='smoothed synthetic $C_n^2$')
plt.legend()
plt.xlabel('$log_{10}$[$C_n^2$ [$m^{2/3}$]]',fontsize=11)
#plt.xlim(-20,-13)
plt.ylabel('Altitude [km]',fontsize=11)
#plt.ylim(0,35)
plt.title('Synthetic $C_n^2$ Profile for '+location+' '+stime+'Z '+tmonth+' '+year)
newname = 'synthetic_cn2_'+stime+'Z_'+month+'_'+year+'_interp.png'
plt.savefig(newname,dpi=300)
print "Saved image: "+newname


# Make plot of subset of synthetic and smoothed Cn2
fig2 = plt.figure()
plt.plot(synthetic_cn2,all_alt/1000.,color='black',label='synthetic $C_n^2$')
plt.plot(smooth_scn2,mid_alt/1000.,color='red',label='smoothed synthetic $C_n^2$')
plt.legend()
plt.xlabel('$log_{10}$[$C_n^2$]')
plt.ylabel('Altitude [km]')
plt.xlim(-16.5,-13.5)
plt.ylim(0.5,4.0)
plt.title('Synthetic $C_n^2$ Profile for '+location+' '+stime+'Z '+tmonth+' '+year)
newname = 'synthetic_cn2_'+stime+'Z_'+month+'_'+year+'_zoomed_interp.png'
plt.savefig(newname,dpi=300)
print "Saved image: "+newname

# Make plot of all synthetic and smoothed Cn2
fig3 = plt.figure()
plt.plot(synthetic_ct2,all_alt/1000.,color='black',label='synthetic $\Delta$T')
plt.plot(smooth_tempdiff,mid_alt/1000.,color='red',label='smoothed synthetic $\Delta$T')
plt.legend()
plt.xlabel('Temperature Difference [K]',fontsize=11)
#plt.xlim(0,.18)
plt.ylabel('Altitude [km]',fontsize=11)
#plt.ylim(0,35)
plt.title('Synthetic $\Delta$T Profile for '+location+' '+stime+'Z '+tmonth+' '+year)
newname = 'synthetic_tempdiff_'+stime+'Z_'+month+'_'+year+'_interp.png'
plt.savefig(newname,dpi=300)
print "Saved image: "+newname


# Make plot of subset of synthetic and smoothed Cn2
fig4 = plt.figure()
plt.plot(synthetic_ct2,all_alt/1000.,color='black',label='synthetic $\Delta$T')
plt.plot(smooth_tempdiff,mid_alt/1000.,color='red',label='smoothed synthetic $\Delta$T')
plt.legend()
plt.xlabel('$\Delta$T [K]')
plt.ylabel('Altitude [km]')
plt.xlim(0,0.1)
plt.ylim(0.5,4.0)
plt.title('Synthetic $\Delta$T Profile for '+location+' '+stime+'Z '+tmonth+' '+year)
newname = 'synthetic_tempdiff_'+stime+'Z_'+month+'_'+year+'_zoomed_interp.png'
plt.savefig(newname,dpi=300)
print "Saved image: "+newname

# Want to get a simulated profile of CT2, so need simulated CN2 and T as well as P values
# Interpolate P values between each pressure level

"""
# Generate a plot title from the input file name
date = infile.split('/')[-1].split('.')[0].split('_')
month = date[0]
month_no = strptime(month[:3],'%b').tm_mon
tmonth = month_name[month_no]
year = date[1]
location = date[2]
plot_title='Daily '+stime+'Z Averaged $C_n^2$ Profiles: '+location+' '+tmonth+' '+year

fig1 = plt.figure()
ax = fig1.gca()
bx = ax.boxplot(log_data_list)

ax.set_xticklabels(['Sfc-925 mb','925-850 mb','850-700 mb','700-500 mb','500-400 mb',\
    '400-300 mb','300-200 mb','200-150 mb','150-100 mb','100-50 mb','50 mb-top'])
plt.xticks(rotation=45,fontsize=6)
plt.title(plot_title)
plt.ylabel('$log_{10}$ $C_n^2$')
plt.xlabel('Pressure levels [mb]',fontsize=8)
plt.ylim(-20,-12)

if(METER is True):
    image_name = month+'_'+year+'_'+location+'_'+stime+'Z_cn2_meter_boxplots.png'
else:
    image_name = month+'_'+year+'_'+location+'_'+stime+'Z_cn2_boxplots.png'
fig1.savefig(image_name,dpi=200)
print "Saved image: "+image_name
"""
