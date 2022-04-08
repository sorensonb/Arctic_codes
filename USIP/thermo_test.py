#!/usr/bin/env python

"""
NAME:
  plot_cn2.py

PURPOSE:
  Calculate the refractive index structure parameter for a maximum of 

MODIFICATIONS:
  Blake Sorenson  - 170628: Written 

"""

# Import functions
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter
from readsounding import *
from compare_cn2 import *
from adpaa import ADPAA

if(len(sys.argv)<3):
    print "SYNTAX: ./thermo_test.py GRAW_file thermo_file"
    sys.exit()

pre2019 = False

radio = readsounding(sys.argv[1],allGraw=True,keepMissing=True)
thermo = ADPAA()
thermo.ReadFile(sys.argv[2])
# Mask any missing data in the thermosonde data
thermo.mask_MVC()
#k = kts2ms(a['UWIND'])
#v = kts2ms(a['VWIND'])

# Delete any constant altitude data or descent data from the radiosonde
# data

###############################################################################
#
# TRANSPLANTED FROM comparison_checker.py
#
# Stretch the thermosonde times to match the radiosonde's times. 
#
###############################################################################

# Delete radiosonde data where the altitudes are missing (only at the end of the file)
radio['Time'] = np.delete(radio['Time'],np.where(radio['ALT']==999999.9999))
radio['UTCTime'] = np.delete(radio['UTCTime'],np.where(radio['ALT']==999999.9999))
radio['TEMP'] = np.delete(radio['TEMP'],np.where(radio['ALT']==999999.9999))
radio['PRESS'] = np.delete(radio['PRESS'],np.where(radio['ALT']==999999.9999))
radio['UWIND'] = np.delete(radio['UWIND'],np.where(radio['ALT']==999999.9999))
radio['VWIND'] = np.delete(radio['VWIND'],np.where(radio['ALT']==999999.9999))
radio['SPD'] = np.delete(radio['SPD'],np.where(radio['ALT']==999999.9999))
radio['DIR'] = np.delete(radio['DIR'],np.where(radio['ALT']==999999.9999))
radio['ALT'] = np.delete(radio['ALT'],np.where(radio['ALT']==999999.9999))

if(pre2019 == True):
    # Thermo data
    start_time = 23331.0
    stop_time = 27380.0 
    # Radio data
    radio_ascent_start = 20154.5
    radio_ascent_stop  = 24628.5
    radio_last_contact = 26398.5
     
    radio_last_index = 5964 # burst

    thermo_start_index = 4626
    thermo_last_index = 8675
    thermo_final_last_index = 10232  # last index in the file
else:
    # Thermo data
    start_time = 13152.0
    stop_time = 14220.0 
    # After 2019/05/03
    radio_ascent_start = 13334.7
    radio_ascent_stop  = 14410.7
    radio_last_contact = 14410.7

    # After 2019/05/03
    radio_last_index = 1076 # burst

    thermo_start_index= 2622
    thermo_last_index = 3690
    thermo_final_last_index = 3918  # last index in the file

# After 2019/05/03
#start_time = 13152.0
#stop_time = 14220.0 
# Before 2019/05/03
#start_time = 23331.0
#stop_time = 27380.0 
# Radio data
# After 2019/05/03
#radio_ascent_start = 13334.7
#radio_ascent_stop  = 14410.7
#radio_last_contact = 14410.7
# Before 2019/05/03
#radio_ascent_start = 20154.5
#radio_ascent_stop  = 24628.5
#radio_last_contact = 26398.5
#radio_last_index = 7734 # last contact
# After 2019/05/03
#radio_last_index = 1076 # burst
# Before 2019/05/03
#radio_last_index = 5964 # burst

#thermo_last_index = 10226 # last contact 
# After 2019/05/03
#thermo_start_index= 2622
#thermo_last_index = 3690
#thermo_final_last_index = 3918  # last index in the file
# Before 2019/05/03
#thermo_last_index = 8675
#thermo_final_last_index = 10232  # last index in the file



#radio_last_contact = 26383.0
radio_start_time  = radio_ascent_start 
radio_last_time  = radio_last_contact 
#radio_diff = (radio_last_time-radio_start_time)/(7054.0-1490.0)
radio_diff = (radio_last_time-radio_start_time)/(np.where(radio['UTCTime']==radio_last_time)[0][0]-np.where(radio['UTCTime']==radio_start_time)[0][0])

thermo_diff = 1.0/radio_diff

diff = start_time-radio_start_time # first launch 2018 05 05 
#plotTime = radio['UTCTime']+diff

# Match the launch times between the thermosonde and radiosonde
plotTime = thermo.data['Time']-diff
thermo.data['Time'] = thermo.data['Time']-diff

if(pre2019 == True):
    # Before 2019/05/03
    ascent_rtime = radio['UTCTime'][1490:radio_last_index]
    ascent_ralt = radio['ALT'][1490:radio_last_index]
    ascent_ttime = plotTime[4626:thermo_last_index]
    ascent_talt = thermo.data['Alt'][4626:thermo_last_index]
else:
    # After 2019/05/03
    ascent_rtime = radio['UTCTime'][0:radio_last_index]
    ascent_ralt = radio['ALT'][0:radio_last_index]
    ascent_ttime = plotTime[thermo_start_index:thermo_last_index]
    ascent_talt = thermo.data['Alt'][thermo_start_index:thermo_last_index]

descent_rtime = radio['UTCTime'][radio_last_index:len(radio['UTCTime'])]
descent_ralt = radio['ALT'][radio_last_index:len(radio['UTCTime'])]
descent_ttime = plotTime[thermo_last_index:thermo_final_last_index]

# Find the ratio of the number of radiosonde ascent times to the number of
# thermosonde ascent times
ascent_diff_ratio = float(len(ascent_rtime))/float(len(ascent_ttime))
# Use that ratio to create new thermosonde ascent times that match with the
# radiosonde's ascent times. Now, the radiosonde and thermosonde ascent times
# have launches at the same time and bursts at the same time, although the 
# burst altitudes differ by a few dozen meters.
ascent_ttime = np.arange(ascent_ttime[0],ascent_rtime[-1]+1,ascent_diff_ratio)

# Unneeded for the 2019/05/03 launch

rplot = np.arange(0,len(ascent_rtime))
tplot = np.arange(0,len(ascent_ttime))

if(pre2019 == True):
    # Repeat the process for the descent
    descent_rtime = radio['UTCTime'][radio_last_index:len(radio['UTCTime'])]
    descent_ralt = radio['ALT'][radio_last_index:len(radio['UTCTime'])]
    #descent_ttime = plotTime[thermo_last_index:len(thermo.data['Time'])]
    descent_talt = thermo.data['Alt'][thermo_last_index:thermo_final_last_index]
    
    descent_diff_ratio = float(len(descent_rtime))/float(len(descent_ttime)-1.0)
    descent_ttime = np.arange(ascent_ttime[-1]+1,descent_rtime[-1]+1,descent_diff_ratio)
    # Combine the new thermosonde times into a single array
    ad_ttime = np.concatenate([ascent_ttime,descent_ttime])
    ad_talt  = np.concatenate([ascent_talt,descent_talt])

    # Before 2019/05/03
    thermo.data['Time'][4626:10232] = ad_ttime
    thermo.data['Alt'][4626:10232] = ad_talt
else:
    # After 2019/05/03
    thermo.data['Time'][thermo_start_index:thermo_last_index] = ascent_ttime
    thermo.data['Alt'][thermo_start_index:thermo_last_index] = ascent_talt

#thermo.data['Time'][2622:thermo_final_last_index] = ascent_ttime
#thermo.data['Alt'][2622:thermo_final_last_index] = ascent_talt

thermo.mask_MVC()

###############################################################################
#
# END OF comparison_checker.py TRANSPLANT
#
###############################################################################

# Account for the fact that the 
###time_offset = thermo.data['Time'][-1]-radio['UTCTime'][-1]
#thermo_start_time = 86332.0
#thermo_stop_time = 86606.0   # these three are for the tethered test
#radio_start_time  = 83093.0
#thermo_start_time = 23331.0
#thermo_stop_time = 27380.0   # these three are for the tethered test
#radio_start_time  = 20154.0

# These three are for the second full launch
thermo_start_time = 13152.0
thermo_stop_time = 14220.0
radio_start_time = 13334.7

if(pre2019 == True):
    # Before 2019/05/03
    combined_start_time = 20154.5
    combined_stop_time = 24629.5
else:
    # After 2019/05/03
    combined_start_time = 13334.7
    combined_stop_time = 14410.7

# This part is made obselete by the comparison_checker.py transplant
####time_offset = thermo_start_time-radio_start_time-23.0 # tethered test
###time_offset = thermo_start_time-radio_start_time
####radio['UTCTime']+=time_offset
###thermo.data['Time']-=time_offset

# Grab the "matching" times to plot a time series of the data later
closetime_thermo = np.array([])
matchalt_thermo = np.array([])
closetime_radio = np.array([])

# Get rid of data that is outside the ascent time
#closeindex_radio = np.where((radio['UTCTime']>=86322.0) & (radio['UTCTime']<=86606.0))[0] # tethered test
closeindex_radio = np.where((radio['UTCTime']>=combined_start_time) & (radio['UTCTime']<=combined_stop_time))[0]
for key in radio.keys():
    if((key is not 'UNITS') & (type(radio[key]) is not str)):
        radio[key] = radio[key][closeindex_radio]

tempdiff = np.array([])
for time in radio['UTCTime']:
    closetime=thermo.data['Time'][:].flat[np.abs(thermo.data['Time'][:]-time).argmin()]
    closetime_thermo = np.append(closetime_thermo[:],closetime)
    close_index = np.where(thermo.data['Time']==closetime)[0]
    matchalt_thermo = np.append(matchalt_thermo[:],thermo.data['Alt'][close_index])
#    print radio['ALT'][np.where(radio['UTCTime']==time)],thermo.data['Alt'][close_index]
#    print thermo.data['TempDiff'][np.where(thermo.data['Time']==closetime)]
    tempdiff = np.append(tempdiff[:],thermo.data['TempDiff'][np.where(thermo.data['Time']==closetime)])


thermo_cn2 = dict(radio)
thermo_cn2 = cn2_calc_thermo(tempdiff,thermo_cn2)
# Adding the method='thermo' flag causes the radiosonde Cn2 to be several
# of magnitude higher than without
sradio = dict(radio)
sradio = smooth(sradio)
sradio = cn2_calc(sradio,method='thermo')
radio = cn2_calc(radio)

# Ignore the missing thermosonde_data
masked_tempdiff = ma.masked_values(tempdiff,1e6)
masked_indices = np.where(masked_tempdiff!=ma.masked)
for key in thermo_cn2:
    if type(thermo_cn2[key]) is np.ndarray:
        thermo_cn2[key] = thermo_cn2[key][masked_indices]


# SMOOTHER
# Put an 11-point smoother on the thermosonde data
# Declare arrays to hold smoothed data
avg_t = np.array([])
avg_u = np.array([])
avg_v = np.array([])
mid_alt = np.array([])
mid_press = np.array([])
avg_cn2t = np.array([])
temp_cn2t = np.zeros(11)
j=0
# Set the totals equal to the first elements of the t, u, and v. When
# the for loop started at 0, an extra 0 was showing up at the beginning
# of the averaged data arrays. Starting the loop at 1 removes the 0, 
# but requires the first elements of the data arrays to be added
# to the total before the start.

# Convert thermo_cn2['CN2T'] to logarithmic for smoothing
thermo_scn2 = dict(thermo_cn2)
thermo_scn2['LABEL'] = 'Smoothed Thermosonde'
thermo_scn2['CN2T'] = np.log10(thermo_scn2['CN2T'])
total_t=thermo_scn2['TEMP'][0]
total_cn2t=thermo_scn2['CN2T'][0]
total_u=thermo_scn2['UWIND'][0]
total_v=thermo_scn2['VWIND'][0]
# Loop through the t, u, and v data
for i in range(1, len(thermo_cn2['CN2T'])):
    # If 11 elements have been summed, average the current total and 
    # append the averages to arrays. 
    if(i%11==0):
        avg_t = np.append(avg_t[:],total_t/11)
        avg_u = np.append(avg_u[:],total_u/11)
        avg_v = np.append(avg_v[:],total_v/11)
        mid_alt = np.append(mid_alt[:],thermo_cn2['ALT'][i-6])
        mid_press = np.append(mid_press[:],thermo_cn2['PRESS'][i-6])
        #avg_cn2t = np.append(avg_cn2t[:],total_cn2t/11)
        avg_cn2t = np.append(avg_cn2t[:],np.average(temp_cn2t))
        j=0
        total_t=0
        total_u=0
        total_v=0
    # Add the current data to the totals
    total_t+=thermo_scn2['TEMP'][i]
    #total_cn2t+=thermo_cn2['CN2T'][i]
    temp_cn2t[j] = thermo_scn2['CN2T'][i]
    j+=1
    total_u+=thermo_scn2['UWIND'][i]
    total_v+=thermo_scn2['VWIND'][i]

# REMOVE to prevent resetting the data with the smoothed data
thermo_scn2['CN2T']=avg_cn2t
thermo_scn2['ALT']=mid_alt

# Convert thermo_cn2['CN2T'] back to linear
thermo_scn2['CN2T'] = 10.**(thermo_scn2['CN2T'])

## Plot the thermosonde Cn2 and radiosonde Cn2 on a graph
#plot_cn2(thermo_scn2,sradio,'save')
##plot_cn2()
##plt.plot(np.log10(radio['CN2']),radio['ALT']/1000.)
##plt.plot(np.log10(sradio['CN2']),sradio['ALT']/1000.)
##plt.plot(np.log10(thermo_cn2['CN2T']),thermo_cn2['ALT']/1000.)
#plt.show()

# Make a plot of the temperature differences versus height
fig4 = plt.figure()
plt.plot(tempdiff[masked_indices],matchalt_thermo[masked_indices]/1000.,color='black')
if(pre2019 == True):
    plt.xlim(0,0.18)
else:
    plt.xlim(0,0.70)
if(pre2019 == True):
    plt.title('Free-Flight Launch 2018/05/05 05:00 UTC',fontsize=12)
else:
    plt.title('Free-Flight Launch 2019/05/04 03:00 UTC',fontsize=12)
plt.xlabel('$\Delta$T [K]')
plt.ylabel('Altitude [km]')
if(pre2019 == True):
    plt.ylim(0,35)
else:
    plt.ylim(0,8)
filename = radio['NAME']+'_tempdiff_vs_alt_newRange.png'
plt.savefig(filename,dpi=300)
print "Saved image",filename

calc_abs = False

# Plot the radiosonde temperature differences versus the thermosonde
# temperature differences.
fig5 = plt.figure(figsize = (8,6))
good_indices = np.where(np.diff(radio['TEMP'])!=0.)
radio_diffs = np.diff(radio['TEMP'])[good_indices]/np.diff(radio['ALT'])[good_indices]
if(calc_abs == True):
    radio_diffs = abs(radio_diffs)
#radio_diffs = abs(np.diff(radio['TEMP'])/np.diff(radio['ALT']))
plt.plot(closetime_thermo[masked_indices],tempdiff[masked_indices],label='Thermosonde [horizontal]')
plt.plot(radio['UTCTime'][:-1][good_indices],radio_diffs,label='Radiosonde [vertical]')
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.xlabel('Seconds [UTC]',fontsize=12)
plt.ylabel('Temperature Difference [K]',fontsize=12)
if(pre2019 == True):
    plt.title('Free-Flight Launch 2018/05/05 05:00 UTC',fontsize=12)
else:
    plt.title('Free-Flight Launch 2019/05/04 03:00 UTC',fontsize=12)
plt.legend()
filename = radio['NAME']+'_tempdiff_radio_vs_thermo.png'
if(calc_abs == True):
    filename = radio['NAME']+'_tempdiff_radio_vs_thermo_abs.png'
plt.savefig(filename,dpi=300)
print "Saved image",filename
sys.exit()

# Plot a scatter plot of thermosonde temp diffs to radiosonde temp diffs
#fig6 = plt.figure()
#plt.scatter(tempdiff[masked_indices][:-3],radio_diffs[masked_indices[0][:-3]])
#plt.xlabel('Thermosonde Temperature Differences [K]')
#plt.xlabel('Radiosonde Temperature Differences [K]')
#plt.show()

if(len(sys.argv)==4):
    model = readsounding(sys.argv[3])
    model = cn2_calc(model,method='thermo')
    plot_cn2(sradio,thermo_scn2,model)

sys.exit(1)

plot_cn2(thermo_cn2,radio)

fig2 = plt.figure()
plt.plot(radio['TEMP'],radio['ALT'])
plt.show()


press, temp, alt, u, v = smooth(a['PRESS'],a['TEMP'],a['ALT'],u,v) 
cn2 = cn2_calc(press,temp,alt,u,v,"thermo")
#cn2 = cn2_calc(a['PRESS'],a['TEMP'],a['ALT'],u,v,"thermo")
#alt = np.delete(a['ALT'],-1)
alt = np.delete(alt,-1)
alt = np.delete(alt,0)

# Calculate test thermosonde cn2
at = np.array([0.004,0.05,0.001,0.003,0.001,0.005,0.001,0.005,0.015,0.003,0.002,0.0015,0.005,0.0078,0.003,0.01,0.08,0.0375,0.01567,0.01009791])
t = np.tile(at,int(len(alt)/len(at)))
diff = len(cn2)-len(t)
i=0
#while(len(t)!=len(a['PRESS'])):
while(len(t)!=len(press)):
    t = np.append(t[:],at[i%len(at)])
    i+=1
#cn2t = cn2_calc_thermo(a['PRESS'],a['TEMP'],t)
cn2t = cn2_calc_thermo(press,temp,t)
cn2t = np.delete(cn2t,-1)
cn2t = np.delete(cn2t,0)
cn2 = np.delete(cn2, np.argwhere(alt<3000))
cn2t = np.delete(cn2t,np.argwhere(alt<3000))
alt = np.delete(alt, np.argwhere(alt<3000))
generate_plot()
plt.plot(cn2,alt,color='black')
plt.plot(cn2t,alt,color='red')
plt.show()
#plt.savefig("17_11_17_23_02_35_CKN_corrected_compared.png",dpi=300)
