#!/usr/bin/env python

"""
  NAME:
    Q1-Sorenson.py

  PURPOSE:
    Plot the TLU for the x1 --> x2 implication.

  SYNTAX:
    $ python Q1-Sorenson.py

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2024/01/24:
      Written

"""

import numpy as np
import matplotlib.pyplot as plt

# Set initial theta and weight values
# -----------------------------------
θ  = 0 
w1 = 0
w2 = 0

# Set the learning rate to 1
# --------------------------
η     = 1

# Set up the x1 and x2 values, as well as the correct x1 --> x2 values
# --------------------------------------------------------------------
x1_arr = np.array([0,0,1,1])
x2_arr = np.array([0,1,0,1])
o_arr  = np.array([1,1,0,1])

errors = np.full((4), 2)
epoch  = 1

print('epoch  x1    x2     o    Σxw    y     e    Δθ    Δw1   Δw2    '+\
    'θ    w1    w2')

fmt_str = '{0:3d} {1:5d} {2:5d} {3:5d} {4:5d} {5:5d} {6:5d} {7:5d} {8:5d}'+\
    ' {9:5d} {10:5d} {11:5d} {12:5d}'
print('{0:63d} {1:5d} {2:5d}'.format(θ,w1,w2))
while(len(np.where(errors == 0)[0]) != len(errors)):
    #for x1, x2, oo in zip(x1_arr, x2_arr, o_arr):
    for ii in range(len(x1_arr)):
        x1 = x1_arr[ii]
        x2 = x2_arr[ii]
        oo = o_arr[ii]
    
        # Calculate the weighted average of the neuron input values, using 
        # the weights 
        # ------------------------------------------------------------------- 
        weighted_avg = x1 * w1 + x2 * w2
    
        # Check how the weighted average compares to the threshold (θ). If
        # the weighted average equals or exceeds the threshold, the neuron
        # output value is 1. Otherwise, the output is 0
        # ----------------------------------------------------------------
        if(weighted_avg >= θ): 
            y_val = 1
        else:
            y_val = 0
      
        # Calculate the error between the neuron output and the corret answer
        # ------------------------------------------------------------------- 
        errors[ii] = oo - y_val
       
        # Use the Delta Rule to correct the threshold and weights
        # -------------------------------------------------------
        Δθ  = -η * errors[ii]
        Δw1 = η * errors[ii] * x1
        Δw2 = η * errors[ii] * x2
    
        θ  += Δθ
        w1 += Δw1
        w2 += Δw2
   
        # Print the output 
        print(fmt_str.format(epoch,x1,x2,oo,weighted_avg,y_val,errors[ii],\
            Δθ,Δw1,Δw2,θ,w1,w2))

    #print(errors, not (0 not in errors))
    #print(len(np.where(errors == 0)) != len(errors))
    epoch += 1
    print('')

    if(epoch > 10):
        break  

print("Final function:")
print(str(w1)+'x1 + '+str(w2)+'x2 >= '+str(θ))

# Calculate an explicit equation for the division line
# ----------------------------------------------------
line_slope = (-1 * w1) / w2
intercept  = θ / w2
print("Explicit equation:")
print("x2 = "+str(line_slope)+"x1 + ",intercept)

# Prepare values for plotting the line dividing the decision regions
# ------------------------------------------------------------------
xvals = np.arange(-1,2.5,0.5)
calc_values = xvals * line_slope + intercept

# Calculate the 'y' values for each x1-x2 pair
# --------------------------------------------
y_vals = np.array([(x1 * w1 + x2 * w2) for x1, x2 in zip(x1_arr, x2_arr)])
beyond_threshold = np.where(y_vals >= θ)
within_threshold = np.where(y_vals < θ)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Generate a figure to show the neural network output
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.axhline(0, color = 'k')
ax.axvline(0, color = 'k')
ax.scatter(x1_arr[beyond_threshold], x2_arr[beyond_threshold], c = 'red')
ax.scatter(x1_arr[within_threshold], x2_arr[within_threshold], c = 'green')
ax.fill_between(xvals, calc_values, y2 = 2, alpha = 0.25, color = 'grey')
ax.plot([xvals[0], xvals[-1]], [calc_values[0], calc_values[-1]], ':', \
    color = 'black')
ax.grid()
ax.text(1.5,0.50,'0', backgroundcolor = 'white', weight = 'bold')
ax.text(1.5,1.00,'1', backgroundcolor = 'white', weight = 'bold')
ax.text(-0.5, 0.5, str(int(w1)) + 'x$_{1}$ + ' + str(int(w2)) + \
    'x$_{2}$ = ' + str(int(θ)), backgroundcolor = 'white', weight = 'bold')
ax.set_title('TLU for implication x$_{1}$ --> x$_{2}$')
ax.set_xlim(-1, 2)
ax.set_ylim(-1, 2)
ax.set_xlabel('x$_{1}$')
ax.set_ylabel('x$_{2}$')

# Save the figure
# ---------------
fig.tight_layout()
outname = 'Q1-Sorenson.png'
fig.savefig(outname,dpi = 200)
print('Saved image',outname)
