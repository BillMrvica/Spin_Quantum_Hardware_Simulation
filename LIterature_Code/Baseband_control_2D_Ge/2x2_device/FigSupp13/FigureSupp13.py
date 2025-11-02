# this file is modified from Fig2/Fig_2c_exchange_25mT.py

import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigSupp13','')
import sys
sys.path.insert(1, path_notebook)
import helper_functions
from helper_functions import data1d_attr, data2d_attr, data_setvaribles, hist_auto_su, hist_auto_1d


import numpy as np
import matplotlib.pyplot as plt
from projects.notebook_tools.notebook_tools import get_data_from, fit_data
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
from matplotlib import pyplot as plt
import qtt
import sys
import glob
from qcodes.data.data_set import load_data
import datetime


MRSIZE=2.5

##  modify from 2023-08-29_q1q2_ramsey_measure_exchanges.py     and     2023-08-29_hahn_residual_exchange.py


#%%  copy from above 
# q1 hahn, flip q2 in the second part of wait time. 'waittime_flip_other' means the time duration when q2 is in the flpped status
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_cz-_1100_s'] = 50 #lp.linspace(1, 101, 21, axis = 1, name = 'ramptime_q1q4', unit = 'ns')
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_s-_1100_cz'] = 50 
# device.settings['q1,q2/v_setpoints/_1100_cz/vpB12'] =  60   <--- measure exchange here with Hahn type sequence


CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
MRSIZE = 2.5
SCSIZE = MRSIZE**2
from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(10,6))
h = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5),  Size.Fixed(1),Size.Fixed(1.5),  Size.Fixed(1),Size.Fixed(4.5),  Size.Fixed(1),Size.Fixed(1.5)]
v = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
ax31 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=3, ny=1))


ax13 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=3))
ax33 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=3, ny=3))




datadir = os.getcwd()




start_time = [    
              '\\2023-08-29\\' +    '18-34-02_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # reference 
              '\\2023-08-29\\' +    '18-35-48_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # two block = 150* 23.496316 ns
              '\\2023-08-29\\' +    '18-36-43_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # three block = 300* 23.496316 ns
              
              '\\2023-08-29\\' +    '18-34-56_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # three block = 450* 23.496316 ns
              ]



names = [glob.glob(datadir + x)[0]  for x in start_time]
datfiles = []
for i in range(len(names)):
    datfiles += [load_data(names[i])]

xdatalist = []
zdatalist = []
for i in range(len(names)):
    xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su_S_North', ['ndarray'] )
    zdatalist += [  zdata ]
    xdatalist += [ xdata ]
    
    
 
from scipy.fft import fft, fftfreq



def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0

p1s = []
p1stds = []
for i in range(len(names)):
    # plt.figure() 

    zdata = zdatalist[i]
    xdata = xdatalist[i]
    zf = fft(zdata)
    xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    


    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 ]
    p1std, p1 = fit_data(xdata, zdata,func=osc_without_decay,p0=p0, plot=False, return_cov=True)
    p1s += [p1]
    p1stds += [p1std]


waittime2phase = 2*np.pi*datfiles[0].metadata['settings']['q1']['x90']['frequency']*1e-9
titlestr = ''
# plt.figure()  
ax = ax13  
for i in range(len(names)):
    name = names[i]
    for j in range(len(name)):
        try:
            tempstr = datetime.datetime.strptime(name[j:j+19], '%Y-%m-%d\\%H-%M-%S').strftime('%Y-%m-%d_%H-%M-%S') 
            break
        except ValueError as Err: 
            tempstr = None
            pass
    if i!=0:
        tempstr = tempstr[-8:]        
    xdata = xdatalist[i]
    zdata = zdatalist[i]
    
    zoffset = i*0.1
    ax.plot(xdata*waittime2phase, zdata +zoffset, marker='o' , markersize=MRSIZE, linestyle='', linewidth=0.5)
    xxx = np.linspace(xdata[0], xdata[-1], (len(xdata)-1)*5+1  )
    ax.plot(xxx*waittime2phase, osc_without_decay(xxx, *p1s[i]) +zoffset, marker='', c='k',  linestyle='-', linewidth=0.5)
    titlestr += (tempstr + ', ')




p1arr = np.array(p1s)
p1arr[:,1] = np.abs(p1arr[:,1])


waittime_flip_other = np.array([0, 150, 300, 450]) * (1e9/datfiles[0].metadata['settings']['q1']['x90']['frequency']   )
correction_idle = (20+37.75+20) * 14.7e-6 * 2*np.pi
correction_ramp = (50+50) * (14.7e-6 + 9e-6)/2 * 2*np.pi
correction_ramp_error =  (50+50) * (14.7e-6 - 9e-6)/2 * 2*np.pi

correction_phases = np.array([0, 1, 2, 3]) * (correction_idle + correction_ramp)


def func1d(x, a, b ):
    return a*x + b

y_error_0 = np.array([p1stds[i][3] for i in range(len(p1stds))])
y_error = np.sqrt( correction_ramp_error**2 + y_error_0**2)
p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error, absolute_sigma=True)

ax = ax33

for i in range(4):
    ax.errorbar(waittime_flip_other[i]/1e3, (p1arr[:,3]-correction_phases)[i],
                 yerr = y_error[i], c=f'C{i}',
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=0.5)
xxx = np.array([-2,14])*1e3
ax.plot( xxx/1e3, func1d(xxx, *p1_func1d) , 'k', linestyle='-', linewidth=0.5)



ax.set_xlim(-1, 12)
ax.set_ylim(-0.55, 0.3)
ax.set_yticks([-0.4, -0.2, 0, 0.2])
ax.set_xticks([0, 5, 10])

print( p1_func1d[0]/(2*np.pi) , ' +- ', p1_func1d_std[0]/(2*np.pi) )

# filename = 'residual_exchange_60mV'
# plt.savefig(filename+'.pdf') 


# hahn_exchange += [   [ 60,     p1_func1d[0]/(2*np.pi), p1_func1d_std[0]/(2*np.pi)  ]     ]
# 8.794349031790484e-06  +-  5.776751816971709e-07 with fit_data( ... ,func=func1d, plot=False, return_cov=True, sigma=y_error)
# 8.794349031746253e-06  +-  1.296686565820007e-06 with fit_data( ... ,func=func1d, plot=False, return_cov=True, sigma=y_error, absolute_sigma=True)



# pulse detail:
# datfiles[0].metadata['circuit']['statements']
 # 'Gxpi2:0 1,2',
 # 'Gxpi2:0 1,2',
 # 'Wait , duration=37.7500',
 # 'Wait , duration=20',
 # 'SetV 1,2, setpoint=_1100_cz, wait_after=0',
 # 'Wait , duration=3524.4474',
 # 'SetV 1,2, setpoint=_1100_s, wait_after=0',
 # 'Wait , duration=20',
 # 'Wait , duration=37.7500',
 # 'Wait , duration=20',
 # 'SetV 1,2, setpoint=_1100_cz, wait_after=0',
 # 'Wait , duration=3524.4474',
 # 'SetV 1,2, setpoint=_1100_s, wait_after=0',
 # 'Wait , duration=20',
 # 'Wait , duration=37.7500',
 # 'Wait , duration=20',
 
 # datfiles[0].metadata['settings']['q1,q2']['v_setpoints']['_ramp_times']['_1100_s-_1100_cz']
 #  50 ns
 






