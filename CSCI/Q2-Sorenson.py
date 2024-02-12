#!/usr/bin/env python

"""
  NAME:
    Q2-Sorenson.py

  PURPOSE:
    Calculate the weights and thresholds for the figure shown in question 2.

  SYNTAX:
    $ python Q2-Sorenson.py

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2024/01/26:
      Written

"""

import numpy as np

# Set initial theta and weight values for the input values
# to each of the middle neurons
# --------------------------------------------------------
x_θ_arr  = np.zeros(3)
x_w1_arr = np.zeros(3)
x_w2_arr = np.zeros(3)

# Set initial threshold and weight values from the middle
# neurons to the output neuron
# -------------------------------------------------------
y_w_arr  = np.zeros(3)
y_θ      = 0

# Set the learning rate to 1
# --------------------------
η     = 1

# Set up the x1 and x2 values, as well as the correct x1 --> x2 values
# --------------------------------------------------------------------
x1_arr = np.array( [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3])
x2_arr = np.array( [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3])
o_arr  = np.array([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], \
                   [0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1], \
                   [1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,1]])

# Set up each possible combinations of outputs from each middle neuron
# and set up the expected, correct system output for these inputs
# --------------------------------------------------------------------
y1_arr = np.array([0,0,0,0,1,1,1,1])
y2_arr = np.array([0,0,1,1,0,0,1,1])
y3_arr = np.array([0,1,0,1,0,1,0,1])
yo_arr = np.array([0,0,0,0,0,0,0,1])

# Initialize arrays to hold the errors for each epoch
# ---------------------------------------------------
errors = np.full((x1_arr.shape[0]), 2)
y_errors = np.full((y1_arr.shape[0]), 2)
epoch  = 1

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# train_neuron calculates the output from a neuron with input values
# x1 and x2 (and an optional x3), weights w1 and w2 (with an optional w3),
# and a threshold θ for a given coordinate pair (x1, x2), and corrects
# the neuron's weights and threshold using the Delta Rule. The function 
# takes as input:
#
#   ii:       the index of the input values to use
#   nx:       denotes the index of the neuron to calculate
#   x1_arr:   the array of x-axis coordinates
#   x2_arr:   the array of y-axis coordinate
#   o_arr:    the array of correct outputs to use in training
#   x_w1_arr: the array of weights from the x1 value to each neuron
#   x_w2_arr: the array of weights from the x2 value to each neuron
#   x_θ_arr:  the array of threshold values for each neuron
#   errors:   the error between the calculated neuron output and
#                 the true output.
#
# After calculating the neuron output using the current weights and
# thresholds, the error between the true value and the output is calculated,
# and then the error is used in Delta Rule to update the weights and 
# threshold values.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def train_neuron(ii, nx, x1_arr, x2_arr, o_arr, x_w1_arr, x_w2_arr, \
                    x_θ_arr, errors):

    x1 = x1_arr[ii]
    x2 = x2_arr[ii]
    w1 = x_w1_arr[nx]
    w2 = x_w2_arr[nx]
    θ  = x_θ_arr[nx]
    oo = o_arr[nx,ii]

    y_val = calc_neuron(x1, x2, θ, w1, w2)

    errors[ii] = oo - y_val

    Δθ  = -η * errors[ii]
    Δw1 = η * errors[ii] * x1
    Δw2 = η * errors[ii] * x2
    
    x_θ_arr[nx]  += Δθ
    x_w1_arr[nx] += Δw1
    x_w2_arr[nx] += Δw2

    return x_w1_arr, x_w2_arr, x_θ_arr, errors, y_val

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# calc_neuron calculates the output from a neuron with input values
# x1 and x2 (and an optional x3), weights w1 and w2 (with an optional w3),
# and a threshold θ for a given coordinate pair (x1, x2). The function 
# takes as input:
#   x1:       the x-axis coordinate
#   x2:       the y-axis coordinate
#   θ:        the threshold value for the neuron
#   w1:       the weight from the x1 value to the neuron
#   w2:       the weight from the x2 value to the neuron
#
#   OPTIONAL INPUTS:
#   x3:       the 3rd input value to the neuron
#   w3:       the weight from the 3rd input value to the neuron
#
# The function outputs either 0 or 1: 0 if the weighted sum of the input
# values is less than the threshold value, and 1 if the weighted sum of
# the input values is greater than or equal to the input value.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def calc_neuron(x1, x2, θ, w1, w2, x3 = None, w3 = None):

    if(x3 is None):
        weighted_avg = x1 * w1 + x2 * w2
    else:
        weighted_avg = x1 * w1 + x2 * w2 + x3 * w3

    if(weighted_avg >= θ): 
        y_val = 1
    else:
        y_val = 0

    return y_val

print("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
print("\n                   TRAINING THE INPUT NEURONS")
print("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
print('epoch  x1    x2     o     e   x1w1  x1w2  x1w3  x2w1  x2w2  x2w3   '+\
        'θ1    θ2    θ3')

while(len(np.where(errors == 0)[0]) != len(errors)):
    for ii in range(len(x1_arr)):
        x1 = x1_arr[ii]
        x2 = x2_arr[ii]
        oo = o_arr[0,ii]
  
        # = = = = = = = = = = = = = = = = = = = = = = = = = 
        #
        # Train neuron 1
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = 
        x_w1_arr, x_w2_arr, x_θ_arr, errors, y_val1 = \
            train_neuron(ii, 0, x1_arr, x2_arr, o_arr, \
                            x_w1_arr, x_w2_arr, x_θ_arr, errors)

        # = = = = = = = = = = = = = = = = = = = = = = = = = 
        #
        # Train neuron 2
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = 
        x_w1_arr, x_w2_arr, x_θ_arr, errors, y_val2 = \
            train_neuron(ii, 1, x1_arr, x2_arr, o_arr, \
                            x_w1_arr, x_w2_arr, x_θ_arr, errors)


        # = = = = = = = = = = = = = = = = = = = = = = = = = 
        #
        # Train neuron 3
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = 
        x_w1_arr, x_w2_arr, x_θ_arr, errors, y_val3 = \
            train_neuron(ii, 2, x1_arr, x2_arr, o_arr, \
                            x_w1_arr, x_w2_arr, x_θ_arr, errors)


        print('{0:3d} {1:5d} {2:5d} {3:5d} {4:5d} {5:5d} {6:5d} {7:5d} {8:5d} {9:5d} {10:5d} {11:5d} {12:5d} {13:5d}'.format(\
            epoch, x1, x2, oo, errors[ii], \
            int(x_w1_arr[0]), int(x_w1_arr[1]), int(x_w1_arr[2]), \
            int(x_w2_arr[0]), int(x_w2_arr[1]), int(x_w2_arr[2]), \
            int(x_θ_arr[0]), int(x_θ_arr[1]), int(x_θ_arr[2])))

    epoch += 1
    print('')

    if(epoch > 100):
        break  

# Now, train the links from the 3 middle neurons to the ending
# neurons
# ------------------------------------------------------------
print("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ")
print("\n                   TRAINING THE OUTPUT NEURON")
print("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ")
print('epoch  y1    y2    y3     o     y     e    w1    w2    w3    yθ')
epoch = 1
while(len(np.where(y_errors == 0)[0]) != len(y_errors)):
    for ii in range(len(y1_arr)):

        y1  = y1_arr[ii]
        y2  = y2_arr[ii]
        y3  = y3_arr[ii]
        yw1 = y_w_arr[0]
        yw2 = y_w_arr[1]
        yw3 = y_w_arr[2]
        yo  = yo_arr[ii]

        # = = = = = = = = = = = = = = = = = = = = = = = = = 
        #
        # Calculate weights and errors for the y neuron
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = 
        y_val = calc_neuron(y1, y2, y_θ, yw1, yw2, x3 = y3, w3 = yw3)

        y_errors[ii] = yo - y_val

        # Correct the threshold and weights using the error
        # -------------------------------------------------
        Δθ  = -η * y_errors[ii]
        Δw1 = η * y_errors[ii] * y1
        Δw2 = η * y_errors[ii] * y2
        Δw3 = η * y_errors[ii] * y3
        
        y_θ  += Δθ
        y_w_arr[0] += Δw1
        y_w_arr[1] += Δw2
        y_w_arr[2] += Δw3

        print('{0:3d} {1:5d} {2:5d} {3:5d} {4:5d} {5:5d} {6:5d} {7:5d} {8:5d} {9:5d} {10:5d}'.format(\
            epoch, y1, y2, y3, yo, y_val, y_errors[ii], \
            int(y_w_arr[0]), int(y_w_arr[1]), int(y_w_arr[2]), \
            int(y_θ)))

    epoch += 1
    print('')

    if(epoch > 100):
        break  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# calc_final_network calculates the output from the trained neural network 
# for a given coordinate pair (x1, x2). The function takes as input:
#   x1:       the x-axis coordinate
#   x2:       the y-axis coordinate
#   x_w1_arr: the trained weights from the x1 value to each of the three 
#                 middle neurons.
#   x_w2_arr: the trained weights from the x2 value to each of the three 
#                 middle neurons.
#   x_θ_arr:  the trained thresholds for each of the middle neurons
#   y_w_arr:  the trained weights from each of the middle neurons
#                 to the output neuron.
#   y_θ:      the trained threshold for the output neuron
#
# The function outputs either 0 or 1: 0 if the given coordinate pair does not
# fall within the desired region given in the question, and 1 if the pair
# does fall within that region. 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def calc_final_network(x1, x2, x_w1_arr, x_w2_arr, x_θ_arr, y_w_arr, y_θ):

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Calculate the outputs for the first three neurons
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Neuron 1
    # --------
    y_val1 = calc_neuron(x1, x2, x_θ_arr[0], x_w1_arr[0], x_w2_arr[0])

    # Neuron 2
    # --------
    y_val2 = calc_neuron(x1, x2, x_θ_arr[1], x_w1_arr[1], x_w2_arr[1])

    # Neuron 3
    # --------
    y_val3 = calc_neuron(x1, x2, x_θ_arr[2], x_w1_arr[2], x_w2_arr[2])

    # Use the previous outputs to calculate the last neuron
    # -----------------------------------------------------
    y_final = calc_neuron(y_val1, y_val2, y_θ, y_w_arr[0], y_w_arr[1], \
        x3 = y_val3, w3 = y_w_arr[2])

    return y_final


print("Final functions:")
print(str(x_w1_arr[0])+'x1 + '+str(x_w2_arr[0])+'x2 >= '+str(x_θ_arr[0]))
print(str(x_w1_arr[1])+'x1 + '+str(x_w2_arr[1])+'x2 >= '+str(x_θ_arr[1]))
print(str(x_w1_arr[2])+'x1 + '+str(x_w2_arr[2])+'x2 >= '+str(x_θ_arr[2]))
print(str(y_w_arr[0])+'y1 + '+str(y_w_arr[1])+'y2 + '+str(y_w_arr[2])+'y3 >= '+str(y_θ))

# Test the system here
# --------------------
print("\nFinal system output")
print("   x1    x2     y")
for ii in range(4):
    for jj in range(4):
        y_val = calc_final_network(ii, jj, x_w1_arr, x_w2_arr, x_θ_arr, \
            y_w_arr, y_θ)
        
        print("{0:5d} {1:5d} {2:5d}".format(ii,jj,y_val))
