# this file is modified from Supp/T2/T2HahnPersist/T2star_persistmode_all.py
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigSupp10and11\\T2HahnPersist','')
import sys
sys.path.insert(1, path_notebook)
import helper_functions
from helper_functions import data1d_attr, data2d_attr, data_setvaribles, hist_auto_su, hist_auto_1d


import numpy as np
import matplotlib.pyplot as plt
from projects.notebook_tools.notebook_tools import get_data_from, fit_data
from qcodes.data.data_set import load_data
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
from matplotlib import pyplot as plt
import qtt
import sys
import glob
from datetime import datetime
#%%
def expdecay(x, A, x0, alpha, y0):
    return A*np.exp(-(x/x0)**alpha) + y0

plt.rcParams.update({'font.size':15})
MRSIZE = 3
LW = 1
#%%
datakey = [
['qubit A  5mT',  138, '2023-09-02\\16-29-15'],
['qubit B  5mT',  138, '2023-09-02\\16-20-28'],
['qubit A 10mT',  163, '2023-08-02\\15-27-35'],
['qubit B 10mT',  163, '2023-08-02\\15-44-33'],
['qubit A 20mT',  163, '2023-07-25\\11-08-21'],
['qubit B 20mT',  163, '2023-07-25\\10-57-28'],
['qubit A 25mT',  152, '2023-08-30\\18-25-06'],
['qubit B 25mT',  152, '2023-08-30\\18-28-23'],
 ]

for description, vth, start_time   in datakey:

    
    end_time = start_time
    
    datadir = os.getcwd() 
    datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
    datfile = datfiles[0]
    x = 2*datfile.waittime_set.ndarray/1e3
    
    histlist =    [datfile.sensor_val_S_North_set.ndarray, datfile.hist_S_North.ndarray]  
    ydata_autothresh = 1- hist_auto_1d(histlist[0], histlist[1], x, plot_thresholding=True )[0]
    
    
    if vth == 'raw':
        y = datfile.su_S_North.ndarray
    elif vth == 'auto':
        y = ydata_autothresh
    elif isinstance(vth, int) or isinstance(vth, float):
        mask = histlist[0][0,:] < vth
        ydata_rethresh = np.zeros((histlist[0].shape[0],))
        for i in range(histlist[0].shape[0]):
            temp_var = histlist[1][i,:]
            ydata_rethresh[i] = 1-np.sum(temp_var[mask])
        y = ydata_rethresh
    else:
        raise ValueError
        
    p0 = [0.45, 50, 1.5, 0.5]
    _, p1 = fit_data(x, y, func=expdecay,p0=p0, plot=False, return_cov=True)
    
    plt.figure(figsize=(8,6))
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top = 0.9)
    # plt.subplot(121)
    # plt.scatter(x, y, s=2)
    plt.plot(x, datfile.su_S_North.ndarray, 'o-', markersize=MRSIZE, linewidth=LW)
    plt.plot(x, ydata_autothresh, 'o-', markersize=MRSIZE, linewidth=LW)
    if vth not in ['raw', 'auto']:
        plt.plot(x, ydata_rethresh,  marker='o', markersize=MRSIZE, linewidth=LW)
    xx = np.linspace(min(x), max(x), 10*len(x))
    plt.plot(xx, expdecay(xx, *p1), linewidth=LW)
    plt.xlabel('time '+ r'$2\tau$'+'(us)')
    plt.ylabel('Hahn echo amplitude ')
    
    # plt.subplot(122)
    
    plt.title(start_time + ' ' + description + '\n' + f'T2H = {p1[1]:.3g} us, alpha={p1[2]:.3g}')
    print(description + f' T2Hahn = {p1[1]:.4g} us \n')








