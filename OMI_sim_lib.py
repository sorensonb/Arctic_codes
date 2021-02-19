#!/home/bsorenson/anaconda3/envs/py37_1/bin/python

"""



"""

import numpy as np
import sys
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons


def read_sim_text(infile):

    OMI_sim_dict = {}

    # Open up the text file
    with open(infile,'r') as fin:
        flines = fin.readlines()

        # Determine the number of lines in the file to set up data arrays
        num_lines = len(flines)-1

        AI = np.zeros(num_lines)
        SZA = np.zeros(num_lines)
        VZA = np.zeros(num_lines)
        AZM = np.zeros(num_lines)
        
        for ni, line in enumerate(flines[1:]):
            templine = line.strip().split()
            AI[ni]  = float(templine[0])
            SZA[ni] = float(templine[1])
            VZA[ni] = float(templine[2])
            AZM[ni] = float(templine[3])

    OMI_sim_dict['AI'] = AI 
    OMI_sim_dict['SZA'] = SZA
    OMI_sim_dict['VZA'] = VZA
    OMI_sim_dict['AZM'] = AZM
    OMI_sim_dict['source'] = 'text'

    return OMI_sim_dict 

def convert_to_NCDF(OMI_sim_dict):
    file_name = 'omi_sim_out.nc'
    
    # Open the dataset
    ds = Dataset(file_name,'w',format='NETCDF4')

    # Figure out sizes for each dimension
    sza_range = int(np.round(np.max(OMI_sim_dict['SZA']),1) - \
                    np.round(np.min(OMI_sim_dict['SZA']),1))+1
    vza_range = int(np.round(np.max(OMI_sim_dict['VZA']),1) - \
                    np.round(np.min(OMI_sim_dict['VZA']),1))+1
    azm_range = int(np.round(np.max(OMI_sim_dict['AZM']),1) - \
                    np.round(np.min(OMI_sim_dict['AZM']),1))+1

    # Create dimensions for the dataset. One dimension each for 
    # solar zenith angle, viewing zenith angle, and azimuth angle
    SZA = ds.createDimension('SZA', sza_range)
    VZA = ds.createDimension('VZA', vza_range)
    AZM = ds.createDimension('AZM', azm_range)

    sza_vals = ds.createVariable('SZA', 'f4', ('SZA',))
    #sza_vals.description('Solar Zenith Angle')
    #sza_vals.units = 'Degrees'
    vza_vals = ds.createVariable('VZA', 'f4', ('VZA',))
    #vza_vals.description('Viewing Zenith Angle')
    #vza_vals.units = 'Degrees'
    azm_vals = ds.createVariable('AZM', 'f4', ('AZM',))
    #azm_vals.description('Sensor Azimuth Angle')
    #azm_vals.units = 'Degrees'
    ai_vals  = ds.createVariable('AI',  'f4', ('SZA','VZA','AZM',))
    #ai_vals.units = 'Unitless'


    # Insert values into each dimension
    sza_vals[:] = np.arange(np.round(np.min(OMI_sim_dict['SZA']),1), \
                  np.round(np.max(OMI_sim_dict['SZA']),1) + 1)
    vza_vals[:] = np.arange(np.round(np.min(OMI_sim_dict['VZA']),1), \
                  np.round(np.max(OMI_sim_dict['VZA']),1) + 1)
    azm_vals[:] = np.arange(np.round(np.min(OMI_sim_dict['AZM']),1), \
                  np.round(np.max(OMI_sim_dict['AZM']),1) + 1)

    # Insert the AI data
    ai_vals[:,:,:] = \
        OMI_sim_dict['AI'].reshape((sza_range,vza_range,azm_range))
    
    ds.close() 

def read_sim_NCDF(infile):
    OMI_sim_dict = {}

    data = Dataset(infile,'r')
 
    OMI_sim_dict['SZA'] = data.variables['SZA'][:]
    OMI_sim_dict['VZA'] = data.variables['VZA'][:]
    OMI_sim_dict['AZM'] = data.variables['AZM'][:]
    OMI_sim_dict['AI'] = data.variables['AI'][:,:,:]
    OMI_sim_dict['source'] = 'NCDF'

    return OMI_sim_dict

def plot_interactive(OMI_sim_dict):
    
    # Find maximum values of values
    max_SZA = np.max(OMI_sim_dict['SZA'][:])
    max_VZA = np.max(OMI_sim_dict['VZA'][:])
    max_AZM = np.max(OMI_sim_dict['AZM'][:])

    fig, ax = plt.subplots()
    plt.subplots_adjust(left = 0.25, bottom = 0.35)
   
    l, = plt.plot(OMI_sim_dict['SZA'][:],OMI_sim_dict['AI'][:,0,0])
    plt.xlabel('SZA')
    plt.ylabel('AI')
 
    axcolor = 'lightgoldenrodyellow'
    ax_sza = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    ax_vza = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    ax_azm = plt.axes([0.25, 0.20, 0.65, 0.03], facecolor=axcolor)

    s_sza = Slider(ax_sza, 'SZA', 0.0, max_SZA, valinit = 0, valstep = 1.0)
    s_vza = Slider(ax_vza, 'VZA', 0.0, max_VZA, valinit = 0, valstep = 1.0)
    s_azm = Slider(ax_azm, 'AZM', 0.0, max_AZM, valinit = 0, valstep = 1.0)

    #def update_SZA(val):
        

    def plotSZA(event):
        print('SZA', s_sza.val, s_vza.val, s_azm.val)
        l.set_xdata(OMI_sim_dict['SZA'][:])
        l.set_ydata(OMI_sim_dict['AI'][:,int(s_vza.val), int(s_azm.val)])
        ax.set_xlabel('SZA')
        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw_idle()
    sza_button_loc = plt.axes([0.5,0.025,0.1,0.04])
    sza_button = Button(sza_button_loc, 'SZA', color=axcolor, hovercolor = '0.975')
    sza_button.on_clicked(plotSZA)

    def plotVZA(event):
        print('VZA', s_sza.val, s_vza.val, s_azm.val)
        l.set_xdata(OMI_sim_dict['VZA'][:])
        l.set_ydata(OMI_sim_dict['AI'][int(s_sza.val), :, int(s_azm.val)])
        ax.set_xlabel('VZA')
        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw_idle()
    vza_button_loc = plt.axes([0.60,0.025,0.1,0.04])
    vza_button = Button(vza_button_loc, 'VZA', color=axcolor, hovercolor = '0.975')
    vza_button.on_clicked(plotVZA)

    def plotAZM(event):
        print('AZM', s_sza.val, s_vza.val, s_azm.val)
        l.set_xdata(OMI_sim_dict['AZM'][:])
        l.set_ydata(OMI_sim_dict['AI'][int(s_sza.val), int(s_vza.val), :])
        ax.set_xlabel('AZM')
        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw_idle()
    azm_button_loc = plt.axes([0.70,0.025,0.1,0.04])
    azm_button = Button(azm_button_loc, 'AZM', color=axcolor, hovercolor = '0.975')
    azm_button.on_clicked(plotAZM)

    plt.show()
