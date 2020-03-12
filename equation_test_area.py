#!/usr/bin/env python
"""


"""
import sys
import matplotlib.pyplot as plt
import numpy as np

# Clear-sky summer values from net flux figure
dflux_darea = -0.039128 # Use the value for ice area change here. To find with
                        # respect to water area, take the negative of this
                        # value. Unis are (W/m2)/km2
del_T = (86400.)*90  # number of seconds of solar heating per summer
l_f = 3.3e5  # latent heat of fusion (J/kg)
rho_i = 917  # density of ice (kg/m3)

def equation(area,thick):

    area_new = area - (Pm(area) * del_T) / (l_f * rho_i * thick)
#    area_pnew = area_p * (1. - (dflux_darea * del_T) / (l_f * rho_i* thick))
    return area_new

def fill_array(start_val,thick):
    area_values = []
    area_values.append(start_val)
    # Find 1m-thickness values
    while(start_val>0):
        area_new = equation(start_val,thick)
        area_values.append(area_new) 
        start_val = area_new
    return area_values

###relation = (40./1000.)
###deltas = np.array([1000,200,400,850,1600,1400,600,900])
###fluxes = deltas*relation
###areas = np.array([3000,2900,2800,2700,2600,2500,2400,2300])
###
###flux_input = sum(fluxes*areas)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# This code section is written to calculate the power_water_ratio using 
# the ice area trends. If, in the future, a way is figured out to run this
# code with the ice dictionary (whether this code is moved to a function 
# in IceLib or comparelib or something else), this code can be switched on.
# Until then, the initial and final ice areas for the June 2001 to 
# August 2018 time period are hard coded in.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
auto_script = True
if(auto_script == False):

    # Find initial ice areas (ice areas at 0 time index)
    # Use this area as the base area for all flux calculations.
    # Units are km2
    ice_areas_initial = (ice_data['grid_ice_conc'][0,:,:]/100.)*ice_data['grid_total_area']
    
    # Add all trends (mostly negative) to the initial ice areas to find the final
    # ice areas.
    # Use the last month that corresponds to the initial month.
    # For example, compare June 2001 with June 2018, do not compare
    # June 2001 with August 2018.
    # On second thought, this doesn't matter. A larger area decrease will
    # mean a larger melting flux increase. A smaller area increase will
    # mean a smaller melting flux increase, so comparing June with August
    # should be okay.
    # Units are km2
    ice_areas_final = ice_areas_initial+ice_data['grid_ice_area_trend']
    water_areas_final = ice_areas_initial-ice_areas_final
    
    # Take the sums of both the initial ice areas and final ice areas.
    # Units are converted to m2
    sum_init_ice  = np.nansum(ice_areas_initial)*1e6
    sum_final_ice = np.nansum(ice_areas_final)*1e6
    sum_final_water = np.nansum(water_areas_final)*1e6
    
    
    # Find the fluxes associated with each ice area decrease in each box
    # In this case the trends act as a 'darea' value
    # Units are ((W/m2)/km2)*km2 = W/m2
    fluxes = dflux_darea*ice_data['grid_ice_area_trend']
    
    # Multiply the fluxes by the base area to find the power input in each
    # box. Convert the ice areas to square meters
    # Units are (W/m2) * m2 = W
    power_input_per_box = fluxes*(ice_areas_initial)
    
    # Taking the sum of the power input per box yields the total power input into
    # the whole base area.
    # Units are W with respect to the base area
    total_power_input = np.nansum(power_input_per_box)
    
    # Take the ratio of the total power input to the total water area 
    power_ice_ratio   = total_power_input/sum_final_ice
    power_water_ratio = total_power_input/sum_final_water
    
    #### To find the amount of melting power at a current ice extent, use:
    ###curr_ice_extent = sum_final_ice*0.95
    #### NOTE: This one doesn't work because a decrease in ice extent relates to
    ####       a decrease in melting flux. Must use the wwater term because 
    ####       a decrease in ice extent leads to an increase in melting flux
    ####melting_power_wice   = flux_ice_ratio*curr_ice_extent
    ###melting_power_wwater = power_water_ratio*(sum_init_ice-curr_ice_extent)

else:

    # Final total ice extent for August 2018, using the trends and the June 2001 
    # extent.
    sum_init_ice  = 7820576049633.293  # m2
    sum_final_ice = 6481114233056.605  # m2
    power_water_ratio = 9.584361746105995e-05  # W/m2

# Set up lambda function to give the melting power for a given ice
# extent (alpha).
Pm = lambda alpha : power_water_ratio*(sum_init_ice-alpha)


starting_val = sum_final_ice 
area_1m = np.array(fill_array(starting_val,1.))
area_2m = np.array(fill_array(starting_val,2.))
area_3m = np.array(fill_array(starting_val,3.))
area_4m = np.array(fill_array(starting_val,4.))
area_5m = np.array(fill_array(starting_val,5.))

fig1 = plt.figure()
plt.plot(area_1m/1e9,label='1 meter thickness')
plt.plot(area_2m/1e9,label='2 meter thickness')
plt.plot(area_3m/1e9,label='3 meter thickness')
plt.plot(area_4m/1e9,label='4 meter thickness')
plt.plot(area_5m/1e9,label='5 meter thickness')
plt.legend()
plt.xlabel('Years')
plt.ylabel('Ice Area (1000s of km2)')
plt.show()
