# this file is modified from Supp/T2/CPMG/2023-08-02_q1q2_CPMG_10mT.py
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
plt.rcParams.update({'font.size':14})

datadir = os.getcwd() 

start_time = [    
# '\\2023-08-02\\' +    '18-15-43_sweep1D_waittime_q1_CPMG_1_q2_down'   ,
# '\\2023-08-02\\' +    '18-16-27_sweep1D_waittime_q1_CPMG_2_q2_down'   ,
# '\\2023-08-02\\' +    '18-17-16_sweep1D_waittime_q1_CPMG_4_q2_down'   ,
# '\\2023-08-02\\' +    '18-18-15_sweep1D_waittime_q1_CPMG_8_q2_down'   ,
'\\2023-08-02\\' +    '20-06-34_sweep1D_waittime_q1_CPMG_1_q2_down'   ,
'\\2023-08-02\\' +    '20-08-12_sweep1D_waittime_q1_CPMG_2_q2_down'   ,
'\\2023-08-02\\' +    '20-10-12_sweep1D_waittime_q1_CPMG_4_q2_down'   ,
# '\\2023-08-02\\' +    '20-12-40_sweep1D_waittime_q1_CPMG_8_q2_down'   ,
# '\\2023-08-02\\' +    '20-21-40_sweep1D_waittime_q1_CPMG_8_q2_down'   ,
'\\2023-08-02\\' +    '21-00-27_sweep1D_waittime_q1_CPMG_8_q2_down'   ,
'\\2023-08-02\\' +    '18-19-31_sweep1D_waittime_q1_CPMG_16_q2_down'   ,
'\\2023-08-02\\' +    '18-21-19_sweep1D_waittime_q1_CPMG_32_q2_down'   ,
'\\2023-08-02\\' +    '18-24-00_sweep1D_waittime_q1_CPMG_64_q2_down'   ,
'\\2023-08-02\\' +    '21-38-43_sweep1D_waittime_q1_CPMG_128_q2_down'   ,
'\\2023-08-02\\' +    '21-55-48_sweep1D_waittime_q1_CPMG_256_q2_down'   ,
'\\2023-08-02\\' +    '23-35-11_sweep1D_waittime_q1_CPMG_512_q2_down'   ,
              ]



factors = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512])

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
# for i in range(8, 10):
    zoffset = i*0.3
    xdata, zdata = xdatalist[i], zdata_newlist[i]
    plt.plot(xdata, zdata+zoffset, label = f'CPMG {factors[i]}')
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
plt.rcParams.update({'font.size':18})

datadir = os.getcwd() 

start_time = [    
'\\2023-08-02\\' +    '20-29-18_sweep1D_waittime_q2_CPMG_1_q1_down'   ,
# '\\2023-08-02\\' +    '20-30-32_sweep1D_waittime_q2_CPMG_2_q1_down'   ,
'\\2023-08-02\\' +    '20-44-34_sweep1D_waittime_q2_CPMG_2_q1_down'   ,
# '\\2023-08-02\\' +    '20-32-07_sweep1D_waittime_q2_CPMG_4_q1_down'   ,
'\\2023-08-02\\' +    '20-49-11_sweep1D_waittime_q2_CPMG_4_q1_down'   ,
# '\\2023-08-02\\' +    '10-42-38_sweep1D_waittime_q2_CPMG_8_q1_down'   ,
# '\\2023-08-02\\' +    '18-29-47_sweep1D_waittime_q2_CPMG_8_q1_down'   ,
'\\2023-08-02\\' +    '21-15-34_sweep1D_waittime_q2_CPMG_8_q1_down'   ,
'\\2023-08-02\\' +    '18-31-09_sweep1D_waittime_q2_CPMG_16_q1_down'   ,
'\\2023-08-02\\' +    '18-33-09_sweep1D_waittime_q2_CPMG_32_q1_down'   ,
'\\2023-08-02\\' +    '18-36-13_sweep1D_waittime_q2_CPMG_64_q1_down'   ,
'\\2023-08-02\\' +    '22-21-11_sweep1D_waittime_q2_CPMG_128_q1_down'   ,
'\\2023-08-02\\' +    '22-40-31_sweep1D_waittime_q2_CPMG_256_q1_down'   ,
'\\2023-08-03\\' +    '00-03-13_sweep1D_waittime_q2_CPMG_512_q1_down'   ,
              ]


factors = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512 ])
# factors = [ 64, 128, 256, 512]


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
#%%
zdata_newlist = []
for line_id in range(len(names)):
   
    if '\\2023-08-03\\' +    '00-03-13_sweep1D_waittime_q2_CPMG_512_q1_down'   in names[line_id]:

        zdata_new = 1- hist_auto_1d(histlist[line_id][0], histlist[line_id][1], xdatalist[line_id],  plot_thresholding=False, expect_pos_global=False, kernal=np.array([1,3,3,1])/8 )[0]
        zdata_newlist += [zdata_new]
    else:
        zdata_new = 1- hist_auto_1d(histlist[line_id][0], histlist[line_id][1], xdatalist[line_id],  plot_thresholding=False, expect_pos_global=False )[0]
        zdata_newlist += [zdata_new]
    

    
    
    
#%%

p1s = []
T2s = []
alphas = []
plt.figure(figsize=(10,8))
for i in range(len(names)):
# for i in range(7,10):
    zoffset = i*0.3
    xdata, zdata = xdatalist[i], zdata_newlist[i]
    plt.plot(xdata, zdata+zoffset, label = f'CPMG {factors[i]}')
    plt.scatter(xdata, zdata+zoffset, s=6, label = f'CPMG {factors[i]}')
    # plt.plot(xdata, zdata+zoffset, 'o-', label = f'CPMG {factors[i]}')
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


