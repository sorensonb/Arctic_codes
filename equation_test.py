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

starting_val = -10. 
sigma_1m = fill_array(starting_val,1.)
sigma_2m = fill_array(starting_val,2.)
sigma_3m = fill_array(starting_val,3.)
sigma_4m = fill_array(starting_val,4.)
sigma_5m = fill_array(starting_val,5.)

fig1 = plt.figure()
plt.plot(sigma_1m,label='1 meter thickness')
plt.plot(sigma_2m,label='2 meter thickness')
plt.plot(sigma_3m,label='3 meter thickness')
plt.plot(sigma_4m,label='4 meter thickness')
plt.plot(sigma_5m,label='5 meter thickness')
plt.legend()
plt.xlabel('Years')
plt.ylabel('Percent Ice Concentration')
plt.show()
