#!/usr/bin/env python

"""


"""

import numpy as np
import sys

sys.path.append('/home/bsorenson')
from sounding_lib import *
from python_lib import *
from cn2_lib import *
from compare_cn2 import *
from datetime import datetime

date_str = '2018050505'
#date_str = '2019050403'
###plot_combined_figure(date_str, save = True)
#plot_combined_figure(date_str, save = True, show_synth_data = False)
#sys.exit()

plot_combined_figure_v2(date_str, fourpanel = True, save = True)
sys.exit()

thermo_scn2 = read_temp_diffs(file_dict[date_str]['radio_file_orig'], \
    file_dict[date_str]['thermo_file'])

# Make a plot of the temperature differences versus height
ax = None
if(int(date_str[:4]) < 2019):
    pre2019 = True
else:
    pre2019 = False
in_ax = True
if(ax is None):
    in_ax = False
    fig = plt.figure(figsize = (11, 4))
    axs = fig.subplots(1, 4, sharey = True)
    ax1 = axs[0]
    ax2 = axs[1]
    ax3 = axs[2]
    ax4 = axs[3]
    ax5 = ax2.twiny()

# Plot wind data
# --------------
ax1.plot(thermo_scn2['TEMP'], \
    thermo_scn2['matchalt_thermo'][thermo_scn2['masked_indices']] / 1000., \
    color = 'tab:blue')
ax1.set_xlabel('Temperature [degC]')
ax1.set_ylabel('Altitude [km]')
ax2.plot(thermo_scn2['SPD'], \
    thermo_scn2['matchalt_thermo'][thermo_scn2['masked_indices']] / 1000., \
    color = 'tab:blue')
ax2.set_xlabel('Wind Speed [kts]', color = 'tab:blue')
ax5.scatter(thermo_scn2['DIR'], \
    thermo_scn2['matchalt_thermo'][thermo_scn2['masked_indices']] / 1000., \
    color = 'tab:orange', s = 3)
ax5.set_xlabel('Wind Direction [deg]', color = 'tab:orange')

if(date_str == '2018050505'):
    synth_date = '201805 00'
    ylim = [0, 30]
else:
    synth_date = '201905 00'
    ylim = [0, 6.5]

# Plot thermosonde temperature differences
# ----------------------------------------
ax3.plot(thermo_scn2['tempdiff'][thermo_scn2['masked_indices']],\
    thermo_scn2['matchalt_thermo'][thermo_scn2['masked_indices']]/1000.,color='tab:blue', \
    label = 'Thermosonde')
if(pre2019 == True):
    ax3.set_xlim(0,0.18)
else:
    ax3.set_xlim(0,0.65)
##!#if(pre2019):
##!#    ax.set_title('Free-Flight Launch 2018/05/05 05:00 UTC',fontsize=9)
##!#else:
##!#    ax.set_title('Free-Flight Launch 2019/05/04 03:00 UTC',fontsize=9)
#ax3.set_title('Thermosonde Temperature Differences', fontsize = 8)
ax3.set_xlabel('Thermosonde $\Delta$T [K]', fontsize = 9)
#ax3.set_ylabel('Altitude [km]', fontsize = 9)
ax3.tick_params(axis = 'both', labelsize=9)


plot_cn2_figure(date_str, ax = ax4, raw = False, old = False, \
        save = False, no_3km = False, ymax = ylim[1])

ax1.grid(alpha = 0.50)
ax2.grid(alpha = 0.50)
ax3.grid(alpha = 0.50)
ax4.grid(alpha = 0.50)


ax1.set_ylim(ylim)
ax2.set_ylim(ylim)
ax3.set_ylim(ylim)
ax4.set_ylim(ylim)

plt.suptitle(date_str)

fig.tight_layout()

plt.show()
