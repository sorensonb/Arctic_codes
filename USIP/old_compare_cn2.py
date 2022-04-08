#!/usr/bin/env python

"""
  NAME:
    compare_cn2.py


"""

import numpy
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from adpaa import ADPAA
#import readsounding as rs

if(len(sys.argv)!=3):
    print "SYNTAX: ./compare_cn2.py <cn2 file 1> <cn2 file 2>"
    sys.exit(1)

profile1 = ADPAA()
profile1.ReadFile(sys.argv[1])
profile2 = ADPAA()
profile2.ReadFile(sys.argv[2])

cn2_diff    = numpy.array([])
alt         = numpy.array([])
press       = numpy.array([])
percent_dif = numpy.array([])

length1=len(profile1.data['Refrct_prm'])
length2=len(profile2.data['Refrct_prm'])

print length1, length2
if(length1>length2):
    outermax=length2
    outeralt=profile2.data['Altitude']  
    outercn2=profile2.data['Refrct_prm']
    innermax=length1
    inneralt=profile1.data['Altitude']
    innercn2=profile1.data['Refrct_prm']
else:
    outermax=length1
    outeralt=profile1.data['Altitude']
    outercn2=profile1.data['Refrct_prm']
    innermax=length2
    inneralt=profile2.data['Altitude']
    innercn2=profile2.data['Refrct_prm']

# Find altitudes that are comparible between the two profiles
for i in range (0, outermax):
    temp_outeralt = 0
    temp_outercn2 = 0
    temp_inneralt = inneralt[0] 
    temp_innercn2 = innercn2[0] 
    temp_cn2diff  = 0
    for j in range(0, innermax):
        if(abs(outeralt[i] - inneralt[j]) < abs(outeralt[i] - temp_inneralt)):
 
            temp_outeralt = outeralt[i]
            temp_outercn2 = outercn2[i]
            temp_inneralt = inneralt[j]
            temp_innercn2 = innercn2[j]
            
        else:
            continue

    temp_cn2diff = temp_outercn2 - temp_innercn2
    percnt_df    = (temp_cn2diff/temp_outercn2)*100.
    cn2_diff = numpy.append(cn2_diff[:], temp_cn2diff)
    alt = numpy.append(alt[:], temp_outeralt)
    percent_dif = numpy.append(percent_dif[:], percnt_df)
    print "Outer Alt: " + str(temp_outeralt) + "\tInner Alt: " + str(temp_inneralt) + \
          "\tInside Cn2: " + str(temp_innercn2) + "\tOuter Cn2: " + str(temp_outercn2) + "\t% diff: " + str(percnt_df)

fig1 = plt.figure()
plt.plot(alt, percent_dif)
#plt.yscale('log')
#plt.ylim(10**-24, 10**-10)
plt.xlim(1, 35000)

plt.show()
