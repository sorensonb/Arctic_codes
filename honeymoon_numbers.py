#!/usr/bin/env python

"""


"""

import numpy as np
import matplotlib.pyplot as plt

extra_car_mileage = {
    'Kalispell': 0.,
    'Missoula': 250.
}

car_dict = {
    'Kalispell': {
        'Honda CR-V': {
            'Cost': 795.,
            'Mileage': {
                'Highway': 31,
                'Park': 20,
            },
        },
        'Jeep Wrangler': {
            'Cost': 850.,
            'Mileage': {
                'Highway': 18,
                'Park': 10,
            },
        },
        ##!#'Honda Passport': {
        ##!#    'Cost': 723.,
        ##!#    'Mileage': {
        ##!#        'Highway': 22,
        ##!#        'Park': 15,
        ##!#    },
        ##!#},
    },
    'Missoula': {
        # Slade
        'Honda CR-V': {
            'Cost': 691.,
            'Mileage': {
                'Highway': 30,
                'Park': 20,
            },
        },
        'Honda Passport': {
            'Cost': 865.,
            'Mileage': {
                'Highway': 22,
                'Park': 15,
            },
        },
        'Nissan Rogue': {
            'Cost': 414.,
            'Mileage': {
                'Highway': 26,
                'Park': 15,
            },
        },
        'Jeep Cherokee 2017': {
            'Cost': 618.,
            'Mileage': {
                'Highway': 22,
                'Park': 15,
            },
        },
        'Jeep Cherokee 2015': {
            'Cost': 550.,
            'Mileage': {
                'Highway': 26,
                'Park': 16,
            },
        },

    }
}

outbound_flight_dict = {
    'Kalispell': {
        'Non': 365.,
        'Refund':  435.
    },
    'Missoula': {
        'Non': 338.,
        'Refund':  388.
    }
    #'Missoula': {
    #    'Non': 339.,
    #    'Refund':  389.
    #}
}

return_flight_dict = {
    'Kalispell': {
        # Delta, 6 am depart, 2:15 arrive
        'La Crosse': {
            'Non': 279.,
            'Refund':    344.,
        }, 
        ## Delta, 3:15 am depart, 8:59 arrive
        #'La Crosse': {
        #    'Non': 329.,
        #    'Refund': 379.,
        #}, 
        'MSP - Delta': {
            'Non': 249.,
            'Refund':    279.,
        }, 
        'MSP - SC': {
            'Non':    154.,
            'Refund': 154.,
        }, 
    },
    'Missoula': {
        'La Crosse': {
            'Non': 269.,
            'Refund':  319.,
        }, 
        'MSP - Delta': {
            'Non': 219.,
            'Refund': 249.,
        }, 
    }
}

fuel_cost = 4.50
extra_airport_car_mileage = 35.
extra_airport_mileage = {
    'La Crosse': {
        'miles': 0
    },
    'MSP - Delta': {
        'miles': 158.
    },
    'MSP - SC': {
        'miles': 158.
    },
}

card_points = 329.54
airbnb = 1554.14


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

total_flight_dict = {}

# outbound and return have the same first keys
for of_key in outbound_flight_dict.keys():

    tf_key = of_key
    total_flight_dict[tf_key] = {}

    for rf_key in return_flight_dict[of_key].keys():
        #print(of_key, ' - ', rf_key)

        #tf_key = of_key+ ' - '+ rf_key

        total_flight_dict[of_key][rf_key] = {}

        # Calculate the extra driving cost
        # --------------------------------
        extra_cost =  (extra_airport_mileage[rf_key]['miles'] / \
            extra_airport_car_mileage) * fuel_cost

        for ft_key in return_flight_dict[of_key][rf_key].keys():

            cost = (outbound_flight_dict[of_key][ft_key] + \
                return_flight_dict[of_key][rf_key][ft_key]) * 2. \
                - card_points + extra_cost

            total_flight_dict[of_key][rf_key][ft_key] = cost
        
            #print('\t',ft_key,np.round(cost,2))


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Calculate total Turo car costs
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

total_car_dict = {}
driving_miles = 430
for of_key in car_dict.keys():

    total_car_dict[of_key] = {}
    for ct_key in car_dict[of_key].keys():
    
        #print(of_key, ct_key)

        # Calculate the gas money needed to drive to/from the airport
        # -----------------------------------------------------------
        extra_money = (extra_car_mileage[of_key] / \
            car_dict[of_key][ct_key]['Mileage']['Highway']) * fuel_cost
   
        #print('\tExtra money = ', extra_money) 
        # Calculate the total gas money needed, including driving and airport
        # -------------------------------------------------------------------
        total_gas_money = (driving_miles / \
            car_dict[of_key][ct_key]['Mileage']['Park']) * fuel_cost + \
            extra_money
        #print('\ttotal money = ', total_gas_money) 
 
        # Add this to the rental price
        # ----------------------------
        #print('\trental cost = ', car_dict[of_key][ct_key]['Cost'])

        total_car_money = total_gas_money + car_dict[of_key][ct_key]['Cost']

        #print('\ttotal cost = ',total_car_money)

        total_car_dict[of_key][ct_key] = total_car_money

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Calculate total trip cost, including airbnb
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

max_cost = -99999
min_cost = 99999
total_trip_dict = {}
for of_key in outbound_flight_dict.keys():

    total_trip_dict[of_key] = {}

    for rf_key in return_flight_dict[of_key].keys():
        total_trip_dict[of_key][rf_key] = {}
        for ft_key in return_flight_dict[of_key][rf_key].keys():
            total_trip_dict[of_key][rf_key][ft_key] = {}
            for ct_key in car_dict[of_key].keys():
                #print(of_key, ' - ', rf_key, ' - ', ft_key, ' - ', ct_key)

                local_cost = total_flight_dict[of_key][rf_key][ft_key] + \
                    total_car_dict[of_key][ct_key] + \
                    airbnb
                total_trip_dict[of_key][rf_key][ft_key][ct_key] =  local_cost

                if(local_cost < min_cost):
                    min_cost = local_cost
                if(local_cost > max_cost):
                    max_cost = local_cost

num_subplots = len(total_trip_dict.keys())
fig = plt.figure(figsize = (5 * num_subplots,5))
#ax = fig.add_subplot(1,1,1)
# Convert the dictionary to an array
# ----------------------------------
for ii, of_key in enumerate(total_trip_dict.keys()):
    ax = fig.add_subplot(1,num_subplots, ii + 1)
#of_key = 'Kalispell'
    data_arr = np.array([[[total_trip_dict[of_key][rf_key][ft_key][ct_key] for \
        ct_key in total_trip_dict[of_key][rf_key][ft_key].keys()] for \
        ft_key in total_trip_dict[of_key][rf_key].keys()] for rf_key in \
         total_trip_dict[of_key].keys()])
   
    x_vals = np.arange(len(data_arr[:,0,0]) * len(data_arr[ii,0,:]))
 
    for ii in range(len(data_arr[:,0,0])):
        ax.bar(x_vals[ii * len(data_arr[ii,0,:]) : ii * len(data_arr[ii,0,:]) + \
            len(data_arr[ii,0,:])] + 0.0, data_arr[ii,0,:], width = 0.25, \
            label = list(total_trip_dict[of_key].keys())[ii] + ': Non')
        ax.bar(x_vals[ii * len(data_arr[ii,0,:]) : ii * len(data_arr[ii,0,:]) + \
            len(data_arr[ii,0,:])] + 0.25, data_arr[ii,1,:], width = 0.25, \
            label = list(total_trip_dict[of_key].keys())[ii] + ': Refund')

    ax.set_xticks(x_vals)
    label_list = []
    for kk in range(len(data_arr[:,0,0])):
        label_list = label_list + list(car_dict[of_key].keys())
    ax.set_xticklabels(label_list, rotation = 45, fontsize = 8)   
    ax.legend(fontsize = 8)
    ax.set_title(of_key)
    ax.grid(axis='y')
    if(ii + 1 == num_subplots):
        ax.set_yticklabels([])
        #ax.get_yaxis().set_visible(False)
    ax.set_ylim(min_cost - 100, max_cost + 100)

fig.tight_layout() 
#ax.bar(x_vals[:len(data_arr[0,0,:])] + 0.0, data_arr[0,0,:], width = 0.25, label = 'Non') 
#ax.bar(x_vals[:len(data_arr[0,0,:])] + 0.25, data_arr[0,1,:], width = 0.25, label = 'Refund')
#ax.bar(x_vals[len(data_arr[0,0,:]):len(data_arr[1,0,:])] + 0.0, data_arr[1,0,:], width = 0.25, label = 'Non') 
#ax.bar(x_vals[len(data_arr[0,0,:]):len(data_arr[1,0,:])] + 0.25, data_arr[1,1,:], width = 0.25, label = 'Refund')
#ax.bar(x_vals[:len(data_arr[2,0,:])] + 0.0, data_arr[2,0,:], width = 0.25, label = 'Non') 
#ax.bar(x_vals[:len(data_arr[2,0,:])] + 0.25, data_arr[2,1,:], width = 0.25, label = 'Refund')
plt.show() 
