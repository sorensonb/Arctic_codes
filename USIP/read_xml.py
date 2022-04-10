#!/usr/bin/env python

"""
  NAME:
    read_xml.py
  
  PURPOSE:
    Read 5 Hz temperature difference measurements from GRAW raw data files
    (.gsf), convert them from resistance differnces to temperature differences,  
    and write the data to a UND NASA formatted file.
  
  PYTHON VERSION:
    2.7.5
  
  SYNTAX:
    ./read_xml.py <GRAW .gsf file>
  
    Note:  If the GRAW sounding session saved multiple files (_1.gsf1, _2.gsf2..)
           read_xml.py automatically looks for the next file to read from.
  
  EXAMPLE
    To read XDATA from a GRAW .gsf file titled 2017080821513208.gsf into a UND
    NASA formatted file:
      ./read_xml.py 2017080821513208.gsf
  
    This creates a file titled 17_08_08_21_11_34.GRAW.tempdiff.1Hz
  
  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>  - 2017/08/09:
      Written
  
  Data order:
      Temp (SensorId>69
      RH (SensorId>68
      Ti Vi
      Other thing
      Press (SensorID>74
  
"""

"""
xmldoc = minidom.parse(sys.argv[1])
xdatalist = xmldoc.getElementsByTagName('XDATA')
for item in xdatalist:
    print item.attributes['byte05'].value

(move through frames)
for frame in root[1]:
    for field in child:
        print field.tag, field.attrib
"""

################################################################################
#
# findData searches through the xml file for every XDATA data message. A 
# FrameTime for each XDATA is saved as well. 
#
# NOTE: Can grab altitude data from the file because XDATA is reported with 
#       every POSITION. Cannot grab temp and press data because XDATA and
#       POSITION are reported more often than press and temp, making it
#       nearly impossible to colocate the temp and press with the XDATA.
#
################################################################################
def findData(root,start_time):
    # Read data from xml file

    # Find the number of XDATA occurences in the file
    length = sum(1 for xdata in root.iter('XDATA'))

    # Declare arrays with the number of XDATA occurences in file
    voltdiff  = np.full(length, 999999.9999)
    alt       = np.full(length, 999999.9999)
    ttime     = np.full(length, 999999.9999)
    i=0

    no_time=False
    atime=start_time
    if(start_time==0):
        no_time=True

    for frame in root.find('SENSORDATA'):
        """    
        # Grab temperature
        if(int(frame.findtext('SensorId'))==69):
            ###temp[i] = frame.find('MeasValue').attrib['Value']
            tempcount+=1
        # Grab pressure 
        elif(int(frame.findtext('SensorId'))==74):
            ###press[i] = frame.find('MeasValue').attrib['Value']
            presscount+=1
        """
        # Grab altitude
        if(frame.find('POSITION') is not None):
            alt[i] = frame.find('POSITION').attrib['Altit_MSL']
        # Grab the XDATA from the current frame if it's there

        # Check if the current frame has XDATA in it
        if(frame.find('XDATA') is not None):
            # If GPS data is needed, it only appears when XDATA comes through,
            # so add code to extract GPS data
        
            # Declare emtpy list to hold the data
            bytelist = [0]*26
            j=0
            # Sort dictionary keys to get bytes in order
            for key in sorted(frame.find('XDATA').attrib):
                # Grab the current byte and insert it into the array
                bytelist[j]=int(frame.find('XDATA').attrib[key])
                j+=1

            # Pull the five byte strings from the complete byte string, reverse
            # each string, and use struct.unpack to convert byte strings to 
            # an integer value. This integer represents the voltage difference
            # between the two probes in microvolts
            value1 = struct.unpack("I",bytearray(bytelist[2:6][::-1]))[0]
            value2 = struct.unpack("I",bytearray(bytelist[6:10][::-1]))[0]
            value3 = struct.unpack("I",bytearray(bytelist[10:14][::-1]))[0]
            value4 = struct.unpack("I",bytearray(bytelist[14:18][::-1]))[0]
            value5 = struct.unpack("I",bytearray(bytelist[18:22][::-1]))[0]
            # Average the 5 measurements for the current second
            avg= (value1+value2+value3+value4+value5)/5 
            # Convert the average value from microvolts to volts
            avg = avg*1.e-6

            if(no_time==True):
                # Convert time from milliseconds to seconds
                atime=int(float(frame.findtext('FrameTime'))/1000.)
                no_time=False
            else:
                atime+=1

            # Calculate checksum to compare with the checksum provided by the
            # thermosonde. 
            ck_a = 0
            ck_b = 0
            # Only use the bytes of good data, ignoring the device ID, chain
            # index and checksum provided by the thermosonde
            for byte in bytelist[2:22]:
                ck_a += byte
                ck_b += ck_a

            final_a = ck_a % 256
            final_b = ck_b % 256
            # If the calculated checksum is the same as the checksum at the end
            # of the data, then the data is good. If not, something went wrong
            ttime[i]    = atime
            if((final_a == bytelist[22]) & (final_b==bytelist[23])):
                # Insert time and measurement into arrays
                voltdiff[i] = avg
            else:
                voltdiff[i] = 999999.9999
                print "Checksum disagreement" 
            i+=1

    outdict = dict()
    outdict['VOLT'] = voltdiff
    outdict['ALT']   = alt
    outdict['TTIME'] = ttime
    outdict['ATIME'] = atime 
    # Return temperature measurements and time
    ###return voltdiff, ttime, atime    
    return outdict

################################################################################
#
# Main
#
################################################################################

import xml.etree.ElementTree as ET
import datetime
import numpy as np
import sys
import re
import struct
import math
from adpaa import ADPAA

if(len(sys.argv)!=2):
    print "SYNTAX: ./read_xml.py <GRAW .gsf file containing XDATA>"
    sys.exit(1)

infile = sys.argv[1]
filename=infile.strip().split('/')[-1]
inPath='/'.join(infile.strip().split('/')[:-1])+'/'
# Parse the XML data from the current file
root = ET.parse(infile).getroot()

# Get UTC start time from the header data in the XML file
times = root.find('Header').find('AscentDate').text.split(' ')[1].split(':')
starttime = datetime.timedelta(hours=int(times[0]),minutes=int(times[1]),\
            seconds=int(times[2])).seconds
# Get the year, month, and date from the file name
year  = filename[2:4] 
month = filename[4:6]
day   = filename[6:8]

voltdiff = np.array([])
press    = np.array([])
alt      = np.array([])
temp     = np.array([])
time     = np.array([])
stop = False
start_time = 0
print "Reading data from " + sys.argv[1] 
# Read data as long as there is another file to read data from. If a GRAW
# data file fills up, it generates another one and inserts a section
# to the first file with the name of the next file. If the current file
# has a new filename at the end, read that data
while stop is False:
    # Grab the temperature difference and time data from the current file
    ###tvoltdiff, ttime, start_time = findData(root,start_time)     
    newdict = findData(root,start_time)
    start_time = newdict['ATIME']
    # Add the new data to the total data arrays
    voltdiff = np.concatenate([voltdiff, newdict['VOLT']])
    alt      = np.concatenate([alt,newdict['ALT']])
    time = np.concatenate([time, newdict['TTIME']])
    # Check if the current .gsf file has a 'NextFile' attribute
    # If it does, use that as the root for the next run of the while loop
    if root.find('NextFile') is not None:
        nextfile = root.find('NextFile')[0].text
        print "Reading data from " + nextfile
        root = ET.parse(inPath+nextfile).getroot()
    # If not, stop reading data
    else:
        stop = True 

# Add the UTC start time to the seconds from the start of the file to get the 
# UTC time of each temperature difference measurement
utctime = time+starttime

###############################################################################
#
# Convert voltage measurements to temperature differences
#
###############################################################################

# From Murphy paper, max deltaR = 0.052234 Ohms, results in voltage of 1.95 V
# Universal calibration formula: R0*delT = 3.572*delV
# With R0=27, delT = 0.258 K

# Convert voltage measurements to temperature differences

tempdiff = np.zeros(len(voltdiff))
good_indices = np.where(voltdiff!=999999.9999)
bad_indices = np.where((voltdiff==999999.9999) | (alt==0.0))

"""
#voltdiff2 = np.copy(voltdiff)
#tempdiff2 = np.copy(voltdiff2)
# Convert Vrms to Vcorrected using parabolic fit
voltdiff2[good_indices] = -0.8002*(voltdiff2[good_indices])**6+\
    5.736*(voltdiff2[good_indices])**5-16.295*(voltdiff2[good_indices])**4+\
    23.346*(voltdiff2[good_indices])**3-17.698*(voltdiff2[good_indices])**2+\
    7.7444*voltdiff2[good_indices]-1.0586
"""

voltdiff3 = np.copy(voltdiff)
tempdiff3 = np.copy(voltdiff3)
# Convert Vrms to Vcorrected using logarithmic fit for data below 0.66 Volts
# and linear fit for data above 0.66 Volts
# Use different equation for different launches and noise floors
# Use the new equations for the full launch on May 4th, 2018.
# For anything else, use the equations for the tethered test by default
if(infile==inPath+'2018050504512973.gsf'):
    print "new stuff"
    # After 2019/05/03
    log_indices = np.where((voltdiff3<0.467) & (voltdiff3>0))
    linear_indices = np.where(voltdiff3>=0.467)
    # Before 2019/05/03
    #log_indices = np.where((voltdiff3<0.435913982) & (voltdiff3>0))
    #linear_indices = np.where(voltdiff3>=0.435913982)
   
    # After 2019/05/03 
    voltdiff3[log_indices] = [0.6117*math.log(x)+0.8698 for x in voltdiff3[log_indices]]
    voltdiff3[linear_indices] = 1.0323*voltdiff3[linear_indices]-0.0731
    # Before 2019/05/03
    #voltdiff3[log_indices] = [0.4696*math.log(x)+0.7847 for x in voltdiff3[log_indices]]
    #voltdiff3[linear_indices] = 1.0206*voltdiff3[linear_indices]-0.046
else:

    # After 2019/05/03
    log_indices = np.where((voltdiff3<0.467) & (voltdiff3>0))
    linear_indices = np.where(voltdiff3>=0.467)
    
    voltdiff3[log_indices] = [0.6117*math.log(x)+0.8698 for x in voltdiff3[log_indices]]
    voltdiff3[linear_indices] = 1.0323*voltdiff3[linear_indices]-0.0731

    # Before 2019/05/03
    #log_indices = np.where((voltdiff3<0.66) & (voltdiff3>0))
    #linear_indices = np.where(voltdiff3>=0.66)
    #
    #voltdiff3[log_indices] = [0.5531*math.log(x)+0.8442 for x in voltdiff3[log_indices]]
    #voltdiff3[linear_indices] = 1.0189*voltdiff3[linear_indices]-0.048

voltdiff3[bad_indices]=999999.9999
voltdiff[bad_indices]=999999.9999
alt[bad_indices]=999999.9999

# Use a linear relationship to convert the corrected voltages into temperature
# differences
tempdiff3[good_indices] = 0.129*voltdiff3[good_indices]-0.001
tempdiff3[bad_indices] = voltdiff[bad_indices]
#tempdiff3[good_indices] = 0.1288*voltdiff3[good_indices]-0.0005
#tempdiff3[bad_indices] = voltdiff[bad_indices]

#tempdiff[good_indices] = (voltdiff[good_indices])/(alpha*R0)
#tempdiff[good_indices] = (3.572*voltdiff[good_indices])/R0


# Print to file
mvc_arr = ['999999.9999','999999.9999','999999.9999','999999.9999']
vscal   = ['1.0000','1.0000','1.0000','1.0000']
vmiss = mvc_arr 
vname   = ['Temperature difference between the two platinum wire temperature'\
    ' probes (2 um diameter)','Voltage measurement from the thermosonde'\
    ' Wheatstone bridge','Raw voltage measurements before the noise floor'+\
    ' correction','Altitude [m]']
vdesc   = ['Time','TempDiff','VoltDiff','VoltUncrt','Alt']
vunits  = ['second','degC','Volt','Volt','m']
outdict = dict({'Time':utctime, 'TempDiff':tempdiff3, 'VoltDiff':voltdiff3,\
                'VoltUncrt':voltdiff,'Alt':alt})

a = ADPAA()
a.NLHEAD  = '22'
a.FFI     = 1001
a.ONAME   = 'Delene, David'
a.ORG     = 'University of North Dakota'
a.SNAME   = 'Thermosonde'
a.MNAME   = 'NASA USIP 2016-2018'
a.IVOL    = 1
a.NVOL    = 1
a.DATE    = '20'+str(year)+' '+str(month)+' '+str(day)
a.RDATE   = re.sub("-"," ",str(datetime.date.today()))+""
a.DX      = 1.0
a.XNAME   = 'Time [second]; Time in UT seconds from midnight on day of the '\
    'ballon launch'
a.NV      = 4
a.VSCAL   = vscal
a.VMISS   = vmiss
a.VNAME   = vname
a.NSCOML  = 0
a.NNCOML  = 4
a.DTYPE   = 'Preliminary Data'
a.VFREQ   = '1 Hz Data'
a.VDESC   = vdesc
a.VUNITS  = vunits
a.data    = outdict
# Create format array
#formats = a.create_format_arr(mvc_arr)
#a.format = formats

# Write to file
outfile = year+'_'+month+'_'+day+'_'+times[0]+'_'+times[1]+'_'+times[2]+\
    '.GRAW.tempdiff.1Hz'
a.WriteFile(outfile)
print "Saved file "+outfile
 
"""

data = [list of 4 bytes backwards]
values = struct.unpack("I",bytearray(data))
print(values[0])
"""
