#!/usr/bin/env python

"""


"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

if(len(sys.argv) < 2):
    print("SYNTAX: plot_training_mae.py any_number_of_training_files")
    sys.exit()


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
    'noland112': 'ReLU', 
    'noland115': 'Softplus', 
}

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

for infile in sys.argv[1:]:

    dtype = infile.strip().split('/')[-1].split('_')[-1].split('.')[0]

    mae_lrelu = train_vals(infile)
    
    xvals = np.arange(mae_lrelu.shape[0])
    
    ax.plot(xvals, mae_lrelu, label = type_dict[dtype])

ax.legend()

plt.show()

