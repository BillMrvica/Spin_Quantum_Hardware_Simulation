# this file is modified from Supp/T2/CPMG/2023-09-02_q1q2_CPMG_5mT.py
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigSupp10and11\\T2CPMG','')
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
import datetime
import glob
        
#%%
plt.rcParams.update({'font.size':15})

datadir = os.getcwd() 



start_time = [
'\\2023-09-02\\18-00-20_sweep1D_waittime_q1_CPMG_1_q2_down', 
'\\2023-09-02\\19-06-48_sweep1D_waittime_q1_CPMG_2_q2_down',
# '\\2023-09-02\\18-01-52_sweep1D_waittime_q1_CPMG_4_q2_down', 
'\\2023-09-02\\18-40-41_sweep1D_waittime_q1_CPMG_4_q2_down',
# '\\2023-09-02\\19-10-46_sweep1D_waittime_q1_CPMG_8_q2_down',
'\\2023-09-02\\20-35-12_sweep1D_waittime_q1_CPMG_8_q2_down',
# '\\2023-09-02\\18-04-48_sweep1D_waittime_q1_CPMG_16_q2_down',
'\\2023-09-02\\20-37-10_sweep1D_waittime_q1_CPMG_16_q2_down',
'\\2023-09-02\\19-16-09_sweep1D_waittime_q1_CPMG_32_q2_down', 
'\\2023-09-02\\18-08-41_sweep1D_waittime_q1_CPMG_64_q2_down', 
'\\2023-09-02\\19-25-27_sweep1D_waittime_q1_CPMG_128_q2_down',
'\\2023-09-02\\18-19-19_sweep1D_waittime_q1_CPMG_256_q2_down',
'\\2023-09-02\\19-49-41_sweep1D_waittime_q1_CPMG_512_q2_down'
               ]

factors = np.array([ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512  ])



names = [glob.glob(datadir + x)[0]  for x in start_time]
datfiles = []
for i in range(len(names)):
    datfiles += [load_data(names[i])]

xdatalist = []
zdatalist = []
histlist = []
for i in range(len(names)):
    
    xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su_trig0', ['ndarray'] )
    zdatalist += [  zdata ]
    xdatalist += [ xdata * factors[i] * 2 *1e-3 ]
    histlist += [ [datfiles[i].sensorval_S_North_trig0_set.ndarray, datfiles[i].hist_trig0.ndarray] ]
    print(max(xdata))
def expdecay(x, A, x0, alpha, y0):
    return A*np.exp(-(x/x0)**alpha) + y0

#%%

zdata_newlist = []
for line_id in range(len(names)):
    zdata_new = 1- hist_auto_1d(histlist[line_id][0], histlist[line_id][1], xdatalist[line_id], plot_thresholding=False)[0]
    zdata_newlist += [zdata_new]



#%%

p1s = []
T2s = []
alphas = []
plt.figure(figsize=(10,8))
for i in range(len(factors)):

    zoffset = i*0.0
    xdata, zdata = xdatalist[i], zdata_newlist[i]
    

    plt.scatter(xdata, zdata+zoffset, s=40, label = f'CPMG {factors[i]}')

    
    if factors[i] == 512:
        remove_pts =  list(range(40, 51))
        mask = list( set(range(51)) - set(remove_pts) )
        xdata = xdata[mask]
        zdata = zdata[mask]
        

    if factors[i] == 512:
        from functools import partial
        bounds = ((0.2,  1, 0.5, 0.45), ( 0.5, np.inf, 4, 0.55)   )
        p0 = [0.35, 50*np.sqrt(factors[i]), 1.5, 0.47]
        _, p1 = fit_data(xdata, zdata, func=expdecay, p0=p0, plot=False, return_cov=True, bounds=bounds)
        xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
        plt.plot(xx, expdecay(xx, *p1)+ zoffset )        
    else:
        p0 = [0.25, 50*np.sqrt(factors[i]), 1.5, 0.38]
        _, p1 = fit_data(xdata, zdata, func=expdecay,p0=p0, plot=False, return_cov=True)
        xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
        plt.plot(xx, expdecay(xx, *p1)+ zoffset )        
    
    
    p1s += [p1]
    T2s += [p1[1]]
    alphas += [p1[2]]
# plt.gca().axes.yaxis.set_ticklabels([])
plt.yticks([0.5, 1])
# plt.xscale('log')    
plt.xlim(1,)    
print(  ''.join(['{:.3g}, ']*len(T2s)).format(*T2s))
print(  ''.join(['{:.3g}, ']*len(alphas)).format(*alphas))




#%%
plt.figure(figsize=(4,3))
plt.scatter(factors, T2s)
plt.xscale('log')
plt.yscale('log')
plt.plot(factors, 10**np.poly1d(np.polyfit( np.log10(factors), np.log10(T2s), deg=1))(np.log10(factors)))
beta = np.polyfit( np.log10(factors), np.log10(T2s), deg=1)[0]
plt.title(r'$T_2^{CPMG}$ ~ ' + r'$n_{\pi}^{\beta}$, ' + r'$\beta = $' + f'{beta:.3g}')

plt.plot( (3,500), ( 20* 3**0.5, 20* 500**0.5), 'k--')
plt.plot( (8,128), ( 25* 8**(3/4), 25* 128**(3/4)  ), 'k--')

#%%
plt.figure(figsize=(4,3))
plt.scatter(factors, alphas)
plt.xscale('log')
# plt.ylim(1,)
plt.grid()


#%%

plt.rcParams.update({'font.size':15})

datadir = os.getcwd() 

start_time = [    
# '\\2023-09-02\\21-00-28_sweep1D_waittime_q2_CPMG_1_q1_down', 
# '\\2023-09-03\\13-00-10_sweep1D_waittime_test',  # CPMG 1
'\\2023-09-03\\13-19-28_sweep1D_waittime_test',   # CPMG 1
# '\\2023-09-02\\21-02-20_sweep1D_waittime_q2_CPMG_2_q1_down', 
# '\\2023-09-03\\13-02-14_sweep1D_waittime_test',  # CPMG 2
'\\2023-09-03\\13-25-55_sweep1D_waittime_test',    # CPMG 2
# '\\2023-09-02\\21-05-07_sweep1D_waittime_q2_CPMG_4_q1_down', 
'\\2023-09-03\\13-05-15_sweep1D_waittime_test',   # CPMG 4
# '\\2023-09-02\\21-09-18_sweep1D_waittime_q2_CPMG_8_q1_down', 
'\\2023-09-03\\13-09-38_sweep1D_waittime_test',   # CPMG 8
# '\\2023-09-02\\21-11-48_sweep1D_waittime_q2_CPMG_16_q1_down', 
'\\2023-09-03\\13-12-21_sweep1D_waittime_test',   # CPMG 16
'\\2023-09-02\\23-08-45_sweep1D_waittime_q2_CPMG_32_q1_down',
'\\2023-09-02\\23-17-19_sweep1D_waittime_q2_CPMG_64_q1_down',
'\\2023-09-02\\23-30-24_sweep1D_waittime_q2_CPMG_128_q1_down', 
'\\2023-09-02\\23-51-23_sweep1D_waittime_q2_CPMG_256_q1_down',
'\\2023-09-03\\00-23-58_sweep1D_waittime_q2_CPMG_512_q1_down'
              ]


factors = np.array([1, 2,  4, 8, 16, 32, 64, 128, 256, 512 ])



names = [glob.glob(datadir + x)[0]  for x in start_time]
datfiles = []
for i in range(len(names)):
    datfiles += [load_data(names[i])]

xdatalist = []
zdatalist = []
histlist = []
for i in range(len(names)):
    xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su_S_North', ['ndarray'] )

    zdatalist += [  zdata ]
    xdatalist += [ xdata * factors[i] * 2 *1e-3 ]
    histlist += [ [datfiles[i].sensor_val_S_North_set.ndarray, datfiles[i].hist_S_North.ndarray] ]

def expdecay(x, A, x0, alpha, y0):
    return A*np.exp(-(x/x0)**alpha) + y0
##%%

zdata_newlist = []
for line_id in range(len(names)):
    zdata_new = 1- hist_auto_1d(histlist[line_id][0], histlist[line_id][1], xdatalist[line_id],  plot_thresholding=False )[0]
    zdata_newlist += [ 1-zdata_new]
#%%  17, 18, 19, 20, 26, 27, 29, 30, 31, 40, ...
# plt.figure()
# plt.plot(zdata_new)
# plt.grid()
#%%

p1s = []
T2s = []
alphas = []
plt.figure(figsize=(10,8))
for i in range(len(names)):

    zoffset = i*0.2
    xdata, zdata = xdatalist[i], zdata_newlist[i]    
    if factors[i] == 512:
        remove_pts = [17, 18, 19, 20, 26, 27, 29, 30, 31, ] + list(range(40, 51))
        mask = list( set(range(51)) - set(remove_pts) )
        xdata = xdata[mask]
        zdata = zdata[mask]
    plt.plot(xdata, zdata+zoffset,  label = f'CPMG {factors[i]}')
    plt.scatter(xdata, zdata+zoffset, s=6, label = f'CPMG {factors[i]}')

        
    if factors[i] == 512:
        from functools import partial
        bounds = ((0.2,  1, 0.5, 0.35), ( 0.5, np.inf, 4, 0.45)   )
        p0 = [0.35, 50*np.sqrt(factors[i]), 1.5, 0.36]
        _, p1 = fit_data(xdata, zdata, func=expdecay, p0=p0, plot=False, return_cov=True, bounds=bounds)
        xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
        plt.plot(xx, expdecay(xx, *p1)+ zoffset )           
    else:
        p0 = [0.25, 50*np.sqrt(factors[i]), 1.5, 0.38]
        _, p1 = fit_data(xdata, zdata, func=expdecay,p0=p0, plot=False, return_cov=True)
        xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
        plt.plot(xx, expdecay(xx, *p1)+ zoffset )    
        
    p1s += [p1]
    T2s += [p1[1]]
    alphas += [p1[2]]    
# plt.gca().axes.yaxis.set_ticklabels([])
plt.yticks([0.5, 1])
# plt.xscale('log')    
plt.xlim(1,)    
print(  ''.join(['{:.3g}, ']*len(T2s)).format(*T2s))
print(  ''.join(['{:.3g}, ']*len(alphas)).format(*alphas))



#%%
plt.figure(figsize=(4,3))
plt.scatter(factors, T2s)
plt.xscale('log')
plt.yscale('log')
plt.plot(factors, 10**np.poly1d(np.polyfit( np.log10(factors), np.log10(T2s), deg=1))(np.log10(factors)))
beta = np.polyfit( np.log10(factors), np.log10(T2s), deg=1)[0]
plt.title(r'$T_2^{CPMG}$ ~ ' + r'$n_{\pi}^{\beta}$, ' + r'$\beta = $' + f'{beta:.3g}')

# plt.plot( (1,8), ( 30* 1**0.5, 30* 8**0.5), 'k--')

# plt.plot( (50,500), ( 30* 50**(2/3), 30* 500**(2/3)), 'k--')

plt.plot( (3,500), ( 20* 3**0.5, 20* 500**0.5), 'k--')
plt.plot( (8,128), ( 25* 8**(3/4), 25* 128**(3/4)  ), 'k--')

#%%
plt.figure(figsize=(4,3))
plt.scatter(factors, alphas)
plt.xscale('log')
plt.ylim(1,)
plt.grid()









