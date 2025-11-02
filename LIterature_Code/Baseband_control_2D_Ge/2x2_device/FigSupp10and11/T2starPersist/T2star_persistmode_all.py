# this file is modified from Supp/T2/T2starPersist/T2star_persistmode_all.py
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigSupp10and11\\T2starPersist','')
import sys
sys.path.insert(1, path_notebook)
import helper_functions
from helper_functions import data1d_attr, data2d_attr, data_setvaribles, hist_auto_su, hist_auto_1d


import numpy as np
import matplotlib.pyplot as plt
from projects.notebook_tools.notebook_tools import get_data_from, fit_data
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import datetime
import qtt
import sys
import glob
from qcodes.data.data_set import load_data


#%%
plt.rcParams.update({'font.size':15})


datakey = [
['qubit A  5mT',  8.63, '2023-09-02\\15-13-32'],
['qubit B  5mT', 18.11, '2023-09-02\\15-26-02'],
['qubit A 10mT',  16.9, '2023-08-01\\14-22-49'],
['qubit B 10mT',  35.4, '2023-08-01\\14-39-10'],    
['qubit A 20mT', 33.84, '2023-07-27\\10-18-37'],
['qubit B 20mT', 70.84, '2023-07-27\\10-42-29'],
['qubit A 25mT',  42.6, '2023-08-30\\15-31-53'],
['qubit B 25mT',  89.4, '2023-08-30\\15-58-03'],
['qubit A 40mT',  67.6, '2023-07-31\\10-48-39'],
['qubit B 40mT', 141.7, '2023-07-31\\11-09-11'],
 ]

for description, f_expect, start_time   in datakey:

    end_time = start_time
    
    datadir = os.getcwd() 
    datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
    
    data = datfiles[0]
    
    
    xdata = data.waittime_set.ndarray[0,:]*1e-3
    zdatas = data.su_S_North.ndarray#[:60]
    
    
    from scipy.fft import fft, fftfreq
    
    def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
        return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0
    
    def osc_with_decay(t, A, f, y0, phi, T2):
        return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0
    
    def osc_without_decay(t, A, f, y0, phi):
        return A*np.cos(2*np.pi*t*f+phi)+y0

    zdata = np.average(zdatas, axis=0)
        
    zf = fft(zdata)
    xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()
    plt.subplots_adjust(top=0.8)
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 5, 1.5]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay_alpha,p0=p0, plot=True, return_cov=True)

    plt.ylim(0,1)
    plt.yticks([0, 0.25, 0.5, 0.75, 1])

    
    f_guess = np.zeros((10,2))
    for n in range(f_guess.shape[0]):
        f_guess[n,:] = n/abs(xdata[1]-xdata[0])- p1[1], n/abs(xdata[1]-xdata[0])+ p1[1]
    temp_var = f_guess.ravel()    
    f_guess_f =  temp_var[np.argsort(np.abs(temp_var - f_expect))[0:2]]

    
    print(description + f' T2* = {p1[4]:.4g} us')
    
    plt.title(data.metadata['location'][:19] + ' ' + description + '\n' +\
                      f'f={f_guess_f[0]:.3g}MHz,  T2*={p1[4]:.3g}us,  alpha={p1[5]:.3g}'  )
    
    t_exp_start = datetime.datetime.strptime(data.metadata['loop']['ts_start'],'%Y-%m-%d %H:%M:%S')
    t_exp_end = datetime.datetime.strptime(data.metadata['loop']['ts_end'],'%Y-%m-%d %H:%M:%S')
    t_exp = (t_exp_end-t_exp_start).total_seconds()
    print('experimental time = ', str(t_exp//60) + ' min ' + str(t_exp % 60) + ' sec \n'  )    

