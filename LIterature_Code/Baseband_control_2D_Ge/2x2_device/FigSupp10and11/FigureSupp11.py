# this file is modified from Supp/T2/CPMG/2023-09-02_q1q2_CPMG_5mT.py
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigSupp10and11','')
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
dataQ1 = {}

plt.rcParams.update({'font.size':15})

CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
MRSIZE = 2
SCSIZE = MRSIZE**2
LW = 0.5
from mpl_toolkits.axes_grid1 import Divider, Size
fig = plt.figure( figsize=(14,8))
h = [Size.Fixed(1), Size.Fixed(5), Size.Fixed(1), Size.Fixed(5), Size.Fixed(1),]
v = [Size.Fixed(0.5), Size.Fixed(5), Size.Fixed(1), Size.Fixed(5), Size.Fixed(1),]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
ax31 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=3, ny=1))



datadir = os.getcwd() + '\\T2CPMG'




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



zdata_newlist = []
for line_id in range(len(names)):
    zdata_new = 1- hist_auto_1d(histlist[line_id][0], histlist[line_id][1], xdatalist[line_id], plot_thresholding=False)[0]
    zdata_newlist += [zdata_new]





p1s = []
T2s = []
alphas = []

ax = ax11
for i in range(len(factors)):
    zoffset = i*0.8
    xdata, zdata = xdatalist[i], zdata_newlist[i]
    


    if factors[i] == 512:
        xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
        remove_pts =  list(range(40, 51))
        mask = list( set(range(51)) - set(remove_pts) )
        ax.plot(xdata[mask], zdata[mask]+zoffset,  'o', c=CQ1, markersize=MRSIZE, linestyle='') 
        ax.plot(xdata[remove_pts], zdata[remove_pts]+zoffset,  'o', c='g', markersize=MRSIZE, linestyle='') 
        dataQ1[str(i)] = dict(xdata_r=xdata[remove_pts], zdata_r=zdata[remove_pts])
        xdata = xdata[mask]
        zdata = zdata[mask]        
        from functools import partial
        bounds = ((0.2,  1, 0.5, 0.45), ( 0.5, np.inf, 4, 0.55)   )
        p0 = [0.35, 50*np.sqrt(factors[i]), 1.5, 0.47]
        _, p1 = fit_data(xdata, zdata, func=expdecay, p0=p0, plot=False, return_cov=True, bounds=bounds)
        ax.plot(xx, expdecay(xx, *p1)+ zoffset, 'k-', linewidth=LW)         
        dataQ1[str(i)].update( dict(xdata=xdata, zdata=zdata, xx=xx, zz=expdecay(xx, *p1), p1=p1, factor=factors[i])  )
    
    else:
        ax.plot(xdata, zdata+zoffset, 'o', c=CQ1, markersize=MRSIZE, linestyle='') 
        
        p0 = [0.25, 50*np.sqrt(factors[i]), 1.5, 0.38]
        _, p1 = fit_data(xdata, zdata, func=expdecay,p0=p0, plot=False, return_cov=True)
        xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
        ax.plot(xx, expdecay(xx, *p1)+ zoffset , 'k-', linewidth=LW)        
        dataQ1[str(i)] = dict(xdata=xdata, zdata=zdata, xx=xx, zz=expdecay(xx, *p1), p1=p1, factor=factors[i])
    
    p1s += [p1]
    T2s += [p1[1]]
    alphas += [p1[2]]

ax.set_xlim(1,)    






dataQ2 = {}

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


zdata_newlist = []
for line_id in range(len(names)):
    zdata_new = 1- hist_auto_1d(histlist[line_id][0], histlist[line_id][1], xdatalist[line_id],  plot_thresholding=False )[0]
    zdata_newlist += [ 1-zdata_new]


p1s = []
T2s = []
alphas = []

ax = ax31
for i in range(len(names)):

    zoffset = i*0.8
    xdata, zdata = xdatalist[i], zdata_newlist[i]    


    if factors[i] == 512:
        xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
        remove_pts = [17, 18, 19, 20, 26, 27, 29, 30, 31, ] + list(range(40, 51))
        mask = list( set(range(51)) - set(remove_pts) )
        
        ax.plot(xdata[mask], zdata[mask]+zoffset,  'o', c=CQ2, markersize=MRSIZE, linestyle='') 
        ax.plot(xdata[remove_pts], zdata[remove_pts]+zoffset,  'o', c='g', markersize=MRSIZE, linestyle='') 
        dataQ2[str(i)] = dict(xdata_r=xdata[remove_pts], zdata_r=zdata[remove_pts])
        xdata = xdata[mask]
        zdata = zdata[mask]        
        from functools import partial
        bounds = ((0.2,  1, 0.5, 0.35), ( 0.5, np.inf, 4, 0.45)   )
        p0 = [0.35, 50*np.sqrt(factors[i]), 1.5, 0.36]
        _, p1 = fit_data(xdata, zdata, func=expdecay, p0=p0, plot=False, return_cov=True, bounds=bounds)
        
        ax.plot(xx, expdecay(xx, *p1)+ zoffset, 'k-', linewidth=LW)    
        dataQ2[str(i)].update( dict(xdata=xdata, zdata=zdata, xx=xx, zz=expdecay(xx, *p1), p1=p1, factor=factors[i])  )
            
    else:
        ax.plot(xdata, zdata+zoffset,  'o', c=CQ2, markersize=MRSIZE, linestyle='') 
        
        p0 = [0.25, 50*np.sqrt(factors[i]), 1.5, 0.38]
        _, p1 = fit_data(xdata, zdata, func=expdecay,p0=p0, plot=False, return_cov=True)
        xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
        ax.plot(xx, expdecay(xx, *p1)+ zoffset, 'k-', linewidth=LW)    
        dataQ2[str(i)] = dict(xdata=xdata, zdata=zdata, xx=xx, zz=expdecay(xx, *p1), p1=p1, factor=factors[i])
            




#%%
plt.rcParams.update({'font.size':4})

CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
# MRSIZE = 2.5
SCSIZE = 5 #MRSIZE**2
from mpl_toolkits.axes_grid1 import Divider, Size
HSPACE = 0.3
VSPACE = 0.3
fig = plt.figure( figsize=(10,8))
# h = [Size.Fixed(HSP), Size.Fixed(1.5), Size.Fixed(HSP), Size.Fixed(1.5), Size.Fixed(HSP), Size.Fixed(1.5), Size.Fixed(HSP), Size.Fixed(1.5)]
# v = [Size.Fixed(1), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5), ]
h = []
for _ in range(5):
    h += [Size.Fixed(HSPACE), Size.Fixed(1.3)]
v = []
for i in range(4):
    if i in [0,3]:
        v += [Size.Fixed(VSPACE), Size.Fixed(1.3)]
    else:
        v += [Size.Fixed(VSPACE), Size.Fixed(1.3)]


divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
axs = np.zeros((4,5),dtype=object)
for i in range(axs.shape[0]):
    for j in range(axs.shape[1]):
        
        axs[axs.shape[0]-i-1,j] = fig.add_axes(divider.get_position(),
                           axes_locator=divider.new_locator(nx=(1+2*j), ny=1+2*i))

xtks = [[0,100,200,300],
        [0, 200, 400, 600],
        [0, 400, 800],
        [0, 100, 200],
        [0,100,200,300],
        [0, 400, 800, 1200],
        [0, 800, 1600, 2400],
        [0, 1000, 2000, 3000],
        [0, 2000, 4000], 
        [0, 2500, 5000]  ]



zall = [dataQ1[str(i)]['zdata'] for i in range(len(dataQ1.keys()))] + [dataQ1['9']['zdata_r']]
zmin, zmax = min(np.hstack(zall))  , max(np.hstack(zall))  
for i in range(len(dataQ1.keys())):
    xdata, zdata = dataQ1[str(i)]['xdata'], dataQ1[str(i)]['zdata']
    xx, zz = dataQ1[str(i)]['xx'], dataQ1[str(i)]['zz']
    ax = axs[i//5,i%5]
    ax.plot(xdata, zdata , 'o', c=CQ1, markersize=MRSIZE, linestyle='') 
    if dataQ1[str(i)]['factor'] == 512:
        ax.plot(dataQ1[str(i)]['xdata_r'], dataQ1[str(i)]['zdata_r'] , 'o', c='k', markersize=MRSIZE, linestyle='')     

    ax.plot(xx, zz ,  'k-', linewidth=LW) 
    if dataQ1[str(i)]['factor'] == 512:
        ax.set_yticks([0.5, 0.8])
    else:
        ax.set_yticks([0.5, 0.8])    
    
    ax.set_ylim(zmin-0.03, zmax+0.03)
    ax.set_xticks(xtks[i])


xtks = [[0,100,200],
        [0, 200, 400, ],
        [0, 200, 400, 600],
        [0, 100, 200, 300],
        [0,250,500],
        [0, 400, 800, 1200],
        [0, 1000, 2000],
        [0, 1000, 2000, 3000],
        [0, 2000, 4000], 
        [0, 2500, 5000]  ]


zall = [dataQ2[str(i)]['zdata'] for i in range(9)] 
zmin, zmax = min(np.hstack(zall))  , max(np.hstack(zall))  
for i in range(len(dataQ2.keys())):
    xdata, zdata = dataQ2[str(i)]['xdata'], dataQ2[str(i)]['zdata']
    xx, zz = dataQ2[str(i)]['xx'], dataQ2[str(i)]['zz']
    ax = axs[2+i//5,i%5]
    ax.plot(xdata, zdata , 'o', c=CQ2, markersize=MRSIZE, linestyle='') 
    if dataQ2[str(i)]['factor'] == 512:
        ax.plot(dataQ2[str(i)]['xdata_r'], dataQ2[str(i)]['zdata_r'] , 'o', c='k', markersize=MRSIZE, linestyle='')     
        
    ax.plot(xx, zz ,  'k-', linewidth=LW) 
    if dataQ2[str(i)]['factor'] == 512:
        ax.set_yticks([0, 0.2, 0.4, 0.6])
    else:
        ax.set_yticks([0.4, 0.6])
        ax.set_ylim(zmin-0.01, zmax+0.01)
    ax.set_xticks(xtks[i])

# filename = 'FigureSupp11' # 'cpmg_5mT'
# plt.savefig(filename+'.pdf')


