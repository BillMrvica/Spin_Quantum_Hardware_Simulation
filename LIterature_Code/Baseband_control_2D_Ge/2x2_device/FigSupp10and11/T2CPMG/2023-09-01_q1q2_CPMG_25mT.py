# this file is modified from Supp/T2/CPMG/2023-09-01_q1q2_CPMG_25mT.py
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
plt.rcParams.update({'font.size':16})

datadir = os.getcwd()


start_time = [    
'\\2023-09-01\\' +    '20-08-06_sweep1D_waittime_q1_CPMG_1_q2_down'   ,
'\\2023-09-01\\' +    '20-09-24_sweep1D_waittime_q1_CPMG_2_q2_down'   ,
'\\2023-09-01\\' +    '20-11-08_sweep1D_waittime_q1_CPMG_4_q2_down'   ,
'\\2023-09-01\\' +    '20-13-34_sweep1D_waittime_q1_CPMG_8_q2_down'   ,
'\\2023-09-01\\' +    '20-16-06_sweep1D_waittime_q1_CPMG_16_q2_down'   ,
'\\2023-09-01\\' +    '20-19-11_sweep1D_waittime_q1_CPMG_32_q2_down'   ,
'\\2023-09-01\\' +    '20-24-09_sweep1D_waittime_q1_CPMG_64_q2_down'   ,
'\\2023-09-01\\' +    '20-33-36_sweep1D_waittime_q1_CPMG_128_q2_down'   ,
'\\2023-09-01\\' +    '20-51-49_sweep1D_waittime_q1_CPMG_256_q2_down'   ,
# '\\2023-09-01\\' +    '21-18-50_sweep1D_waittime_q1_CPMG_512_q2_down'   ,
'\\2023-09-01\\' +    '23-41-16_sweep1D_waittime_q1_CPMG_512_q2_down'   ,
              ]

factors = np.array([ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
# factors = [ 512  ]


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
    
def expdecay(x, A, x0, alpha, y0):
    return A*np.exp(-(x/x0)**alpha) + y0

#%%

zdata_newlist = []
for line_id in range(len(names)):
    zdata_new = 1- hist_auto_1d(histlist[line_id][0], histlist[line_id][1], xdatalist[line_id], plot_thresholding=False )[0]
    zdata_newlist += [zdata_new]



#%%

p1s = []
T2s = []
alphas = []
plt.figure(figsize=(10,8))
for i in range(len(factors)):
    zoffset = i*0.0
    xdata, zdata = xdatalist[i], zdata_newlist[i]
    plt.plot(xdata, zdata+zoffset,  label = f'CPMG {factors[i]}')
    plt.scatter(xdata, zdata+zoffset, s=6, label = f'CPMG {factors[i]}')
    p0 = [0.45, 50*np.sqrt(factors[i]), 1.5, 0.52]
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
plt.ylim(1,)
plt.grid()


#%%
import helper_functions
plt.rcParams.update({'font.size':18})

datadir = os.getcwd() 

start_time = [    
'\\2023-09-01\\' +    '17-33-11_sweep1D_waittime_q2_CPMG_1_q1_down'   ,
'\\2023-09-01\\' +    '18-09-59_sweep1D_waittime_q2_CPMG_2_q1_down'   ,
'\\2023-09-01\\' +    '17-34-16_sweep1D_waittime_q2_CPMG_4_q1_down'   ,
'\\2023-09-01\\' +    '18-11-33_sweep1D_waittime_q2_CPMG_8_q1_down'   ,
'\\2023-09-01\\' +    '17-36-18_sweep1D_waittime_q2_CPMG_16_q1_down'   ,
'\\2023-09-01\\' +    '18-14-21_sweep1D_waittime_q2_CPMG_32_q1_down'   ,
'\\2023-09-01\\' +    '17-39-37_sweep1D_waittime_q2_CPMG_64_q1_down'   ,
'\\2023-09-01\\' +    '18-19-37_sweep1D_waittime_q2_CPMG_128_q1_down'   ,
'\\2023-09-01\\' +    '17-49-39_sweep1D_waittime_q2_CPMG_256_q1_down'   ,
'\\2023-09-01\\' +    '18-36-22_sweep1D_waittime_q2_CPMG_512_q1_down'   ,
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
    xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su_trig1', ['ndarray'] )
    zdatalist += [  zdata ]
    xdatalist += [ xdata * factors[i] * 2 *1e-3 ]
    histlist += [ [datfiles[i].sensorval_S_North_trig1_set.ndarray, datfiles[i].hist_trig1.ndarray] ]
def expdecay(x, A, x0, alpha, y0):
    return A*np.exp(-(x/x0)**alpha) + y0


zdata_newlist = []
for line_id in range(len(names)):
    zdata_new = 1- hist_auto_1d(histlist[line_id][0], histlist[line_id][1], xdatalist[line_id],  plot_thresholding=False )[0]
    zdata_newlist += [ 1-zdata_new]
#%%

p1s = []
T2s = []
alphas = []
plt.figure(figsize=(10,8))
for i in range(len(names)):
    zoffset = i*0.2
    xdata, zdata = xdatalist[i], zdata_newlist[i]
    plt.plot(xdata, zdata+zoffset,  label = f'CPMG {factors[i]}')
    plt.scatter(xdata, zdata+zoffset, s=6, label = f'CPMG {factors[i]}')
    p0 = [0.45, 50*np.sqrt(factors[i]), 1.5, 0.52]
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


