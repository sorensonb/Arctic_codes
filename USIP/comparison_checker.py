#!/usr/bin/env python

"""
NAME:
  comparison_checker.py

PURPOSE:
  Compare the altitudes between the raw data file and the RTS table to 
  match times.

MODIFICATIONS:
  Blake Sorenson  - 170628: Written 

"""

# Import functions
import numpy as np
import sys
import matplotlib.pyplot as plt
from readsounding import *
from adpaa import ADPAA

if(len(sys.argv)!=3):
    print("SYNTAX: ./thermo_test.py Graw_RTS_file thermo_file")
    sys.exit()

radio = readsounding(sys.argv[1],allGraw=True,keepMissing=True)

flight2019 = False

if(radio['NAME'][:2] == '19'):
    flight2019 = True
else:
    flight2019 = False

thermo = ADPAA()
thermo.ReadFile(sys.argv[2])
# Mask any missing data in the thermosonde data
thermo.mask_MVC()

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


#start_time = 86332.0
#stop_time  = 86606.0
if(flight2019):
    # Thermo data
    # After 2019/05/03
    start_time = 13152.0
    stop_time = 14220.0 
    # Radio data
    # After 2019/05/03
    radio_ascent_start = 13334.7
    radio_ascent_stop  = 14410.7
    radio_last_contact = 14410.7
    # After 2019/05/03
    radio_last_index = 1076 # burst
    # After 2019/05/03
    thermo_start_index= 2622
    thermo_last_index = 3690
    thermo_final_last_index = 3918  # last index in the file

else:
    # Before 2019/05/03
    start_time = 23331.0
    stop_time = 27380.0 
    # Before 2019/05/03
    radio_ascent_start = 20154.5
    radio_ascent_stop  = 24628.5
    radio_last_contact = 26398.5
    #radio_last_index = 7734 # last contact
    # Before 2019/05/03
    radio_last_index = 5964 # burst
    
    #thermo_last_index = 10226 # last contact 
    # Before 2019/05/03
    thermo_start_index = 4626
    thermo_last_index = 8675
    thermo_final_last_index = 10232  # last index in the file

#radio_last_contact = 26383.0
radio_start_time  = radio_ascent_start 
radio_last_time  = radio_last_contact 
#radio_diff = (radio_last_time-radio_start_time)/(7054.0-1490.0)
radio_diff = (radio_last_time-radio_start_time)/(np.where(radio['UTCTime']==radio_last_time)[0][0]-np.where(radio['UTCTime']==radio_start_time)[0][0])

thermo_diff = 1.0/radio_diff

#thermo.data['Time'] = np.arange(radio_start_time,radio_last_time,thermo_diff)
    
#ascent_start = np.where(thermo.data['Time'][:]==thermo.data['Time'][:].flat[np.abs(thermo.data['Time'][:]-radio_start_time).argmin()])
#ascent_stop  = np.where(thermo.data['Time'][:]==thermo.data['Time'][:].flat[np.abs(thermo.data['Time'][:]-radio_last_time).argmin()])

#ascent_start = np.where(thermo.data['Time']==radio_start_time)[0][0]
#ascent_stop  = np.where(thermo.data['Time']==radio_last_time)[0][0]

# Get times to cooperate
#diff = 86332.0-83093.0-23.0 # tethered test

diff = start_time-radio_start_time # first launch 2018 05 05 
#plotTime = radio['UTCTime']+diff

# Match the launch times between the thermosonde and radiosonde
plotTime = thermo.data['Time']-diff
thermo.data['Time'] = thermo.data['Time']-diff

if(flight2019):
    # After 2019/05/03
    ascent_rtime = radio['UTCTime'][0:radio_last_index]
    ascent_ralt = radio['ALT'][0:radio_last_index]
    ascent_ttime = plotTime[thermo_start_index:thermo_last_index]
    ascent_talt = thermo.data['Alt'][thermo_start_index:thermo_last_index]
else:
    # Before 2019/05/03
    ascent_rtime = radio['UTCTime'][1490:radio_last_index]
    ascent_ralt = radio['ALT'][1490:radio_last_index]
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

# Repeat the process for the descent
descent_rtime = radio['UTCTime'][radio_last_index:len(radio['UTCTime'])]
descent_ralt = radio['ALT'][radio_last_index:len(radio['UTCTime'])]
#descent_ttime = plotTime[thermo_last_index:len(thermo.data['Time'])]
descent_talt = thermo.data['Alt'][thermo_last_index:thermo_final_last_index]

descent_diff_ratio = float(len(descent_rtime))/float(len(descent_ttime)-1.0)
descent_ttime = np.arange(ascent_ttime[-1]+1,descent_rtime[-1]+1,descent_diff_ratio)

rplot = np.arange(0,len(ascent_rtime))
tplot = np.arange(0,len(ascent_ttime))

##### Combine the new thermosonde times into a single array
ad_ttime = np.concatenate([ascent_ttime,descent_ttime])
ad_talt  = np.concatenate([ascent_talt,descent_talt])

"""
# Try to colocate the data based off of altitude and not time
for time in radio['UTCTime']:
    closetime=thermo.data['Time'][:].flat[np.abs(thermo.data['Time'][:]-time).argmin()]
    closetime_thermo = np.append(closetime_thermo[:],closetime)
    close_index = np.where(thermo.data['Time']==closetime)[0]
    matchalt_thermo = np.append(matchalt_thermo[:],thermo.data['Alt'][close_index])
#    print thermo.data['TempDiff'][np.where(thermo.data['Time']==closetime)]
    tempdiff = np.append(tempdiff[:],thermo.data['TempDiff'][np.where(thermo.data['Time']==closetime)])
"""
### Combine the new thermosonde times into a single array
##ad_ttime = np.concatenate([ascent_ttime,descent_ttime])
##ad_talt  = np.concatenate([ascent_talt,descent_talt])

if(flight2019):
    # After 2019/05/03
    thermo.data['Time'][thermo_start_index:thermo_last_index] = ascent_ttime
    thermo.data['Alt'][thermo_start_index:thermo_last_index] = ascent_talt
    #thermo.data['Time'][2622:thermo_final_last_index] = ascent_ttime
    #thermo.data['Alt'][2622:thermo_final_last_index] = ascent_talt
else:
    # Before 2019/05/03
    thermo.data['Time'][thermo_start_index:thermo_final_last_index] = ad_ttime
    thermo.data['Alt'][thermo_start_index:thermo_final_last_index] = ad_talt

thermo.mask_MVC()

fig1 = plt.figure()
#plt.plot(plotTime,radio['ALT'],color='orange',label='radio')
plt.plot(rplot,ascent_rtime,color='orange',label='radio')
#plt.plot(thermo.data['Time'],thermo.data['Alt'],color='blue',label='thermo')
plt.plot(tplot,ascent_ttime,color='blue',label='thermo')
plt.legend()


fig2 = plt.figure()
#plt.plot(plotTime,radio['ALT'],color='orange',label='radio')
plt.plot(ascent_rtime,ascent_ralt,color='orange',label='radio')
#plt.plot(thermo.data['Time'],thermo.data['Alt'],color='blue',label='thermo')
plt.plot(ascent_ttime,ascent_talt,color='blue',label='thermo')
#plt.plot(plotTime,radio['ALT'],color='orange',label='radio')
###plt.plot(descent_rtime,descent_ralt,color='orange',label='radio')
####plt.plot(thermo.data['Time'],thermo.data['Alt'],color='blue',label='thermo')
###plt.plot(descent_ttime,descent_talt,color='blue',label='thermo')
plt.legend()

plt.show()
sys.exit()

fig2 = plt.figure(figsize=(8,6))
plt.xlabel('Time [seconds]',fontsize=11)
plt.xticks(fontsize=10)
ax1 = plt.gca() 
ax1.plot(thermo.data['Time'][thermo_start_index:thermo_last_index],thermo.data['TempDiff'][thermo_start_index:thermo_last_index],color='blue',label='Temperature Difference')
ax1.tick_params(axis='y',labelsize=10)
#ax1.plot(thermo.data['Time'][4626:8675],thermo.data['TempDiff'][4626:8675],color='blue',label='Temperature Difference')
ax1.set_ylabel('Temperature Difference [K]',color='blue',fontsize=11)
ax2 = ax1.twinx() 
ax2.plot(thermo.data['Time'][thermo_start_index:thermo_last_index],thermo.data['Alt'][thermo_start_index:thermo_last_index],color='black',label='Altitude')
ax2.tick_params(axis='y',labelsize=10)
#ax2.plot(thermo.data['Time'][4626:8675],thermo.data['Alt'][4626:8675],color='black',label='Altitude')
ax2.set_ylabel('Altitude [m]',color='black',fontsize=11)
plt.title('2018-05-05 05:00 UTC Thermosonde Flight')
plt.savefig('18_05_05_05_11_04_thermo_altTempDiff_ascent.png')
print("Saved image 18_05_05_05_11_04_thermo_altTempDiff_ascent.png")
plt.show()

sys.exit(1)

thermo.data['Alt']-=thermo.data['Alt'][ascent_start]

fig1 = plt.figure()
plt.xlabel('Time [UTC]')
ax1 = plt.gca()
ax1.plot(thermo.data['Time'][ascent_start:ascent_stop+1],thermo.data['Alt'][ascent_start:ascent_stop+1],color='black')
ax1.set_ylabel('Altitude AGL [m]',color='black')

ax2 = ax1.twinx()
ax2.set_ylim(0.001,0.035)
ax2.plot(thermo.data['Time'][ascent_start:ascent_stop+1],thermo.data['TempDiff'][ascent_start:ascent_stop+1],color='blue')
ax2.set_ylabel('Temperature Difference [degC]',color='blue')
ax2.spines['right'].set_color('blue')

plt.title('2018/05/05 05 UTC')
#fig1.savefig('17_11_17_22_51_30_Graw_ascent_altTemp.png',dpi=300)
#fig1.savefig('18_05_05_04_04_29_Graw_ascent_altTemp.png',dpi=300)
plt.show()

