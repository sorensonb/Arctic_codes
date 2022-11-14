#!/usr/bin/env python

"""


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Read in river crests
infile = sys.argv[1]

data = pd.read_csv(infile, delim_whitespace = True, names = ['idx','hgt','l1','l2','date'])
data['date'] = pd.to_datetime(data['date'])
data.sort_values(by = 'date', inplace = True)

# Plot them
fig = plt.figure(figsize = (9,4))
ax = fig.add_subplot(1,1,1)

ax.plot(data['date'],data['hgt'])
ax.axhline(y = 28., color = 'orange')
ax.axhline(y = 40., color = 'red')
ax.axhline(y = 46., color = 'purple')
ax.set_ylabel('Stage height [ft]')
ax.set_title('Grand Forks Red River Historical Crests')
plt.show()
