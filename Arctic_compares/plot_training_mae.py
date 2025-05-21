#!/usr/bin/env python

"""
  NAME:
    plot_training_mae.py

  PURPOSE:
    Plots the mean absolute error (MAE) of each simulation as a function
    of the training epoch. 

  SYNTAX:
    ./plot_training_mae.py any_number_of_training_out_files

   ./plot_training_mae.py training_log_files/force_ai_test.109791_noland108.txt \
        training_log_files/force_ai_test.167022_noland112.txt \
        training_log_files/force_ai_test.167265_noland115.txt \
        training_log_files/force_ai_test.167486_noland117.txt \
        training_log_files/force_ai_test.167530_noland118.txt

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

l_save_fig = False

if(len(sys.argv) < 2):
    print("SYNTAX: plot_training_mae.py any_number_of_training_files [save]")
    print("        [save]: if added, saves the figure") 
    print("        Can also just use 'default' as an argument") 
    sys.exit()
else:
    
    # Check if the user wants to save the figure
    # ------------------------------------------
    arg_list = sys.argv[1:] 
    if('save' in arg_list):
        arg_list.remove('save')
        l_save_fig = True


    # Check if the user wants the default file list
    # ---------------------------------------------
    if(arg_list[0] == 'default'):
        file_list = ['training_log_files/force_ai_test.109791_noland108.txt',\
            'training_log_files/force_ai_test.166958_noland110.txt', \
            'training_log_files/force_ai_test.167022_noland112.txt', \
            #'training_log_files/force_ai_test.167230_noland113.txt', \
            'training_log_files/force_ai_test.167265_noland115.txt', \
            'training_log_files/force_ai_test.167486_noland117.txt']  
            #'training_log_files/force_ai_test.167530_noland118.txt']
    else:
        file_list = arg_list
    
def train_vals(infile):
    # Parse the data
    # --------------
    with open(infile, 'r') as fin:
        flines = np.array(fin.readlines())
    
    checkers = np.array([tline[:5] for tline in flines])
    good_idxs = np.where(checkers == 'Epoch')[0] + 1
    good_idxs[0] = good_idxs[1] - 2
    mae_vals  = np.array([float(tline.strip().split()[8]) for tline in flines[good_idxs]])
    return mae_vals

type_dict = {
    'noland108': 'LeakyReLU',
    'noland109': 'Linear', 
    'noland110': 'elu', 
    'noland112': 'ReLU', 
    'noland113': 'SeLU', 
    'noland115': 'Softplus', 
    'noland116': 'Tanh', 
    'noland117': 'Softsign', 
    'noland118': 'LeakyReLU 2',
}

# Set up the figure
# -----------------
fig = plt.figure()
ax = fig.add_subplot(1,1,1)


for infile in file_list:

    dtype = infile.strip().split('/')[-1].split('_')[-1].split('.')[0]

    mae_lrelu = train_vals(infile)
    
    xvals = np.arange(mae_lrelu.shape[0])
    
    ax.plot(xvals, mae_lrelu, label = type_dict[dtype])

ax.set_ylim(2.8, 3.7)
ax.set_xlabel('Training epochs')
ax.set_ylabel('Mean Absolute Error [Wm$^{-2}$]')
ax.grid(alpha = 0.5)
ax.set_title('Neural Network Error During Training')

ax.legend()

if(l_save_fig):
    outname = 'training_mae.png'
    fig.savefig(outname, dpi = 300)
    print("Saved image ", outname)
else:
    plt.show()

