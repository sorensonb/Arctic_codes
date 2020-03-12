#!/usr/bin/env python
"""


"""
import sys
import matplotlib.pyplot as plt
import numpy as np

# Clear-sky summer values from net flux figure
dflux_dsigma = -1.3183
del_T = (86400.)*90  # number of seconds of solar heating per summer
l_f = 3.3e5  # latent heat of fusion (J/kg)
rho_i = 917  # density of ice (kg/m3)
def equation(sigma_p,thick):

    sigma_pnew = sigma_p * (1. - (dflux_dsigma * del_T) / (l_f * rho_i* thick))
    return sigma_pnew

def fill_array(start_val,thick):
    sigma_values = []
    sigma_values.append(start_val)
    # Find 1m-thickness values
    while(start_val>-100):
        sigma_new = equation(start_val,thick)
        sigma_values.append(sigma_new) 
        start_val = sigma_new
    return sigma_values

sum_init_ice  = 7820576049633.293/1e6  # km2 as of June 2001
sum_final_ice = 6481114233056.605/1e6  # km2 

# 10 years
#starting_val = -10. 
starting_val = ((sum_final_ice-sum_init_ice)/sum_init_ice)*100.
start_year = 2018

sigma_1m = np.array(fill_array(starting_val,1.))
sigma_2m = np.array(fill_array(starting_val,2.))
sigma_3m = np.array(fill_array(starting_val,3.))
sigma_4m = np.array(fill_array(starting_val,4.))
sigma_5m = np.array(fill_array(starting_val,5.))

extents_1m = (sum_init_ice*(100.+sigma_1m)/100.)/1e6 
years_1m   = np.arange(start_year,start_year+len(extents_1m))
extents_2m = (sum_init_ice*(100.+sigma_2m)/100.)/1e6 
years_2m   = np.arange(start_year,start_year+len(extents_2m))
extents_3m = (sum_init_ice*(100.+sigma_3m)/100.)/1e6 
years_3m   = np.arange(start_year,start_year+len(extents_3m))
extents_4m = (sum_init_ice*(100.+sigma_4m)/100.)/1e6 
years_4m   = np.arange(start_year,start_year+len(extents_4m))
extents_5m = (sum_init_ice*(100.+sigma_5m)/100.)/1e6 
years_5m   = np.arange(start_year,start_year+len(extents_5m))

print("1m depth gone at year:",years_1m[-1])
print("2m depth gone at year:",years_2m[-1])
print("3m depth gone at year:",years_3m[-1])
print("4m depth gone at year:",years_4m[-1])
print("5m depth gone at year:",years_5m[-1])

fig1 = plt.figure()
plt.plot(years_1m,extents_1m,label='1 meter thickness')
plt.plot(years_2m,extents_2m,label='2 meter thickness')
plt.plot(years_3m,extents_3m,label='3 meter thickness')
plt.plot(years_4m,extents_4m,label='4 meter thickness')
plt.plot(years_5m,extents_5m,label='5 meter thickness')
plt.legend()
plt.xlabel('Year')
plt.ylabel('Ice Area (millions of km2)')
plt.show()
