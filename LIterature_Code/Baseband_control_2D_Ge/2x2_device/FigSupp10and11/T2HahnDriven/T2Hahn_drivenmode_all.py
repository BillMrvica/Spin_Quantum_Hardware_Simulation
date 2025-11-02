# this file is modified from Supp/T2/T2HahnDriven/T2Hahn_drivenmode_5mT_to_40mT.py and 
# 2023-09-04_q1q2_Hahn_1mT.py 
# 2023-09-04_q1q2_Hahn_2mT.py
# 2023-09-03_q1q2_Hahn_3mT.py
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigSupp10and11\\T2HahnDriven','')
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
['qubit A  5mT', 'auto', '2023-09-04\\04-26-26'],
['qubit B  5mT', 'auto', '2023-09-04\\04-33-42'],
['qubit A 10mT', 143, '2023-09-03\\20-02-33'],
['qubit B 10mT', 'auto', '2023-09-03\\20-07-24'],
['qubit A 25mT', 141, '2023-09-03\\18-55-32'],
['qubit B 25mT', 'auto', '2023-09-03\\18-58-00'],
['qubit A 40mT', 147, '2023-09-03\\21-59-43'],
['qubit B 40mT', 143, '2023-09-03\\22-01-39'],
 ]

for description, vth, start_time   in datakey:

    end_time = start_time
    
    datadir = os.getcwd()  
    datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
    datfile = datfiles[0]
    x = 2*datfile.waittime_set.ndarray/1e3
    
    histlist =    [datfile.sensor_val_S_North_set.ndarray, datfile.hist_S_North.ndarray]  
    ydata_autothresh = 1- hist_auto_1d(histlist[0], histlist[1], x, plot_thresholding=False )[0]

    
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

    plt.plot(x, datfile.su_S_North.ndarray, 'o-', markersize=MRSIZE, linewidth=LW)
    plt.plot(x, ydata_autothresh, 'o-', markersize=MRSIZE, linewidth=LW)
    if vth not in ['raw', 'auto']:
        plt.plot(x, ydata_rethresh,  marker='o', markersize=MRSIZE, linewidth=LW)
    xx = np.linspace(min(x), max(x), 10*len(x))
    plt.plot(xx, expdecay(xx, *p1), linewidth=LW)
    plt.xlabel('time '+ r'$2\tau$'+'(us)')
    plt.ylabel('Hahn echo amplitude ')
    

    
    plt.title(start_time + ' ' +  description + '\n' + f'T2H = {p1[1]:.3g} us, alpha={p1[2]:.3g}')
    print(description + f' T2Hahn = {p1[1]:.4g} us \n')



#%%
def expdecay(x, A, x0, alpha, y0):
    return A*np.exp(-(x/x0)**alpha) + y0
MRSIZE = 3
LW = 1
plt.rcParams.update({'font.size':15})
#%% 1mT , modified from 2023-09-04_q1q2_Hahn_1mT.py 
description = 'qubit A 1mT' 
start_time = '2023-09-04\\02-30-39'  # very noisy
vth = 'auto'
end_time = start_time
datadir = os.getcwd()  
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


start_time = '2023-09-04\\02-31-49'  # very noisy
vth = 'auto'
end_time = start_time
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile2 = datfiles[0]

x = 2*np.hstack([datfile.waittime_set.ndarray,datfile2.waittime_set.ndarray])/1e3

histlist =    [np.vstack([datfile.sensor_val_S_North_set.ndarray,datfile2.sensor_val_S_North_set.ndarray,]), 
               np.vstack([datfile.hist_S_North.ndarray, datfile2.hist_S_North.ndarray])    ]
ydata_autothresh = 1- hist_auto_1d(histlist[0], histlist[1], x, plot_thresholding=False )[0]

if vth == 'raw':
    y = np.hstack([datfile.su_S_North.ndarray,datfile2.su_S_North.ndarray])
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

plt.figure(figsize=(8,5))
plt.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top = 0.9)

plt.plot(x, np.hstack([datfile.su_S_North.ndarray,datfile2.su_S_North.ndarray]), 'o-', markersize=MRSIZE, linewidth=LW)
plt.plot(x, ydata_autothresh, 'o-', markersize=MRSIZE, linewidth=LW)
if vth not in ['raw', 'auto']:
    plt.plot(x, ydata_rethresh,  marker='o', markersize=MRSIZE, linewidth=LW)
xx = np.linspace(min(x), max(x), 10*len(x))
plt.plot(xx, expdecay(xx, *p1), linewidth=LW)
plt.xlabel('time '+ r'$2\tau$'+'(us)')
plt.ylabel('Hahn echo amplitude ')
plt.title(start_time + ' ' + description + '\n' + f'T2H = {p1[1]:.3g} us, alpha={p1[2]:.3g}',fontsize=12)
print(description + f' T2Hahn = {p1[1]:.4g} us \n')





description = 'qubit B 1mT' 
start_time = '2023-09-04\\03-29-15'
vth = 'auto'
end_time = start_time

datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]
x = 2*datfile.waittime_set.ndarray/1e3

histlist =    [datfile.sensor_val_S_North_set.ndarray, datfile.hist_S_North.ndarray]  
ydata_autothresh = 1- hist_auto_1d(histlist[0], histlist[1], x, plot_thresholding=False )[0]

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

plt.figure(figsize=(8,5))
plt.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top = 0.9)
plt.plot(x, datfile.su_S_North.ndarray, 'o-', markersize=MRSIZE, linewidth=LW)
plt.plot(x, ydata_autothresh, 'o-', markersize=MRSIZE, linewidth=LW)
if vth not in ['raw', 'auto']:
    plt.plot(x, ydata_rethresh,  marker='o', markersize=MRSIZE, linewidth=LW)
xx = np.linspace(min(x), max(x), 10*len(x))
plt.plot(xx, expdecay(xx, *p1), linewidth=LW)
plt.xlabel('time '+ r'$2\tau$'+'(us)')
plt.ylabel('Hahn echo amplitude ')
plt.title(start_time + ' ' + description + '\n' + f'T2H = {p1[1]:.3g} us, alpha={p1[2]:.3g}',fontsize=12)
print(description + f' T2Hahn = {p1[1]:.4g} us \n')


#%% 2mT, modified from 2023-09-04_q1q2_Hahn_2mT.py
description = 'qubit A 2mT' 
start_time = '2023-09-04\\00-35-16'
vth = 'auto'
end_time = start_time

datadir = os.getcwd()  
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]
x = 2*datfile.waittime_set.ndarray/1e3

histlist =    [datfile.sensor_val_S_North_set.ndarray, datfile.hist_S_North.ndarray]  
ydata_autothresh = 1- hist_auto_1d(histlist[0], histlist[1], x, plot_thresholding=False )[0]

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

plt.figure(figsize=(8,5))
plt.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top = 0.9)
plt.plot(x, datfile.su_S_North.ndarray, 'o-', markersize=MRSIZE, linewidth=LW)
plt.plot(x, ydata_autothresh, 'o-', markersize=MRSIZE, linewidth=LW)
if vth not in ['raw', 'auto']:
    plt.plot(x, ydata_rethresh,  marker='o', markersize=MRSIZE, linewidth=LW)
xx = np.linspace(min(x), max(x), 10*len(x))
plt.plot(xx, expdecay(xx, *p1), linewidth=LW)
plt.xlabel('time '+ r'$2\tau$'+'(us)')
plt.ylabel('Hahn echo amplitude ')
plt.title(start_time + ' ' + description + '\n' + f'T2H = {p1[1]:.3g} us, alpha={p1[2]:.3g}')
print(description + f' T2Hahn = {p1[1]:.4g} us \n')




description = 'qubit B 2mT' 
start_time = '2023-09-04\\00-44-16'
vth = 'auto'
end_time = start_time

datadir = os.getcwd()  
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]
x = 2*datfile.waittime_set.ndarray/1e3

histlist =    [datfile.sensor_val_S_North_set.ndarray, datfile.hist_S_North.ndarray]  
ydata_autothresh = 1- hist_auto_1d(histlist[0], histlist[1], x, plot_thresholding=False )[0]

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

plt.figure(figsize=(8,5))
plt.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top = 0.9)
plt.plot(x, datfile.su_S_North.ndarray, 'o-', markersize=MRSIZE, linewidth=LW)
plt.plot(x, ydata_autothresh, 'o-', markersize=MRSIZE, linewidth=LW)
if vth not in ['raw', 'auto']:
    plt.plot(x, ydata_rethresh,  marker='o', markersize=MRSIZE, linewidth=LW)
xx = np.linspace(min(x), max(x), 10*len(x))
plt.plot(xx, expdecay(xx, *p1), linewidth=LW)
plt.xlabel('time '+ r'$2\tau$'+'(us)')
plt.ylabel('Hahn echo amplitude ')
plt.title(start_time + ' ' + description + '\n' + f'T2H = {p1[1]:.3g} us, alpha={p1[2]:.3g}')
print(description + f' T2Hahn = {p1[1]:.4g} us \n')

#%% 3mT, modified from 2023-09-03_q1q2_Hahn_3mT.py
description = 'qubit A 3mT'
start_time = '2023-09-03\\17-29-46'
vth = 'auto'

end_time = start_time

datadir = os.getcwd()
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]
x = 2*datfile.waittime_set.ndarray/1e3

histlist =    [datfile.sensor_val_S_North_set.ndarray, datfile.hist_S_North.ndarray]  
ydata_autothresh = 1- hist_auto_1d(histlist[0], histlist[1], x, plot_thresholding=False )[0]

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

bounds = ((0.25,  10, 0.5, 0.585), ( 0.5, 200, 4, 0.61)   )
p0 = [0.3, 50, 1.5, 0.605]
_, p1 = fit_data(x, y, func=expdecay,p0=p0, plot=False, return_cov=True, bounds=bounds)

plt.figure(figsize=(8,6))
plt.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top = 0.9)
plt.plot(x, datfile.su_S_North.ndarray, 'o-', markersize=MRSIZE, linewidth=LW)
plt.plot(x, ydata_autothresh, 'o-', markersize=MRSIZE, linewidth=LW)
if vth not in ['raw', 'auto']:
    plt.plot(x, ydata_rethresh,  marker='o', markersize=MRSIZE, linewidth=LW)
xx = np.linspace(min(x), max(x), 10*len(x))
plt.plot(xx, expdecay(xx, *p1), linewidth=LW)
plt.xlabel('time '+ r'$2\tau$'+'(us)')
plt.ylabel('Hahn echo amplitude ')
plt.title(start_time + ' ' + description + '\n' + f'T2H = {p1[1]:.3g} us, alpha={p1[2]:.3g}',fontsize=12)
print(description + f' T2Hahn = {p1[1]:.4g} us \n')





description = 'qubit B 3mT'
start_time = '2023-09-03\\17-44-11'
vth = 'auto'
end_time = start_time

datadir = os.getcwd()
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]
x = 2*datfile.waittime_set.ndarray/1e3

histlist =    [datfile.sensor_val_S_North_set.ndarray, datfile.hist_S_North.ndarray]  
ydata_autothresh = 1- hist_auto_1d(histlist[0], histlist[1], x, plot_thresholding=False )[0]

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
plt.plot(x, datfile.su_S_North.ndarray, 'o-', markersize=MRSIZE, linewidth=LW)
plt.plot(x, ydata_autothresh, 'o-', markersize=MRSIZE, linewidth=LW)
if vth not in ['raw', 'auto']:
    plt.plot(x, ydata_rethresh,  marker='o', markersize=MRSIZE, linewidth=LW)
xx = np.linspace(min(x), max(x), 10*len(x))
plt.plot(xx, expdecay(xx, *p1), linewidth=LW)
plt.xlabel('time '+ r'$2\tau$'+'(us)')
plt.ylabel('Hahn echo amplitude ')
plt.title(start_time + ' ' + description + '\n' + f'T2H = {p1[1]:.3g} us, alpha={p1[2]:.3g}',fontsize=12)
print(description + f' T2Hahn = {p1[1]:.4g} us \n')

