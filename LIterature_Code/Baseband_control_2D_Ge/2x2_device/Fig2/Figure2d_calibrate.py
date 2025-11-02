# this file is modified from Fig2/Fig_2c_calib.py
# modified from 2023-08-20_v0_1d_cphase.py
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

import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\Fig2','')
import sys
sys.path.insert(1, path_notebook)
from helper_functions import data1d_attr, data2d_attr, data_setvaribles


    
#%%
from scipy.fft import fft, fftfreq



def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0





#%%
datadir = os.getcwd()

start_time = [    # t_Ramp =  23ns vpb12 = -76
'\\2023-08-20\\' +    '13-46-07_sweep1D_waittimecalibration_cphase_q2_down'   ,
'\\2023-08-20\\' +    '13-46-34_sweep1D_waittimecalibration_cphase_q2_up'   
              ]



names = [glob.glob(datadir + x)[0]  for x in start_time]
datfiles = []
for i in range(len(names)):
    datfiles += [load_data(names[i])]

xdatalist = []
zdatalist = []
for i in range(len(names)):
    # xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su_S_North', ['ndarray'] )
    # xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su00', ['ndarray'] )
    if 'q2'  in  names[i]:
        xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su10', ['ndarray'] )
    elif 'q1'  in  names[i]:
        xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su01', ['ndarray'] )
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
for i in range(len(names)):
    # plt.figure() 

    zdata = zdatalist[i]
    xdata = xdatalist[i]
    zf = fft(zdata)
    xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    


    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 ]
    _, p1 = fit_data(xdata, zdata,func=osc_without_decay,p0=p0, plot=False, return_cov=True)
    p1s += [p1]



titlestr = '' # 'B=-65, tramp = 46.5'
plt.figure()    
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
    plt.scatter(xdata, zdata, label=tempstr  )
    # plt.scatter(xdata, xdatalist[i], label=names[i][0:19)
    xxx = np.linspace(xdata[0], xdata[-1], (len(xdata)-1)*5+1  )
    plt.plot(xxx, osc_without_decay(xxx, *p1s[i] ) )

    titlestr += (tempstr + ', ')
plt.legend()
    
plt.title( titlestr   )




p1arr = np.array(p1s)
p1arr[:,1] = np.abs(p1arr[:,1])

for i in range(2):
    if p1arr[i,0] > 0:
        pass
    else:
        p1arr[i,3] += np.pi
        p1arr[i,0] *= -1
        
phi0, phi1 = p1arr[:,3]
f0, f1 = p1arr[:,1]
print( phi0/np.pi, phi1/np.pi )
print( phi0/np.pi - phi1/np.pi )
print( (phi0 - phi1)/(2*np.pi*(f0+f1)/2) )
print( phi0/np.pi, phi1/np.pi )
print( (2*np.pi - phi0)/(2*np.pi*(f0+f1)/2)  )
print( (1*np.pi - phi1)/(2*np.pi*(f0+f1)/2)  )

print( (1*np.pi - phi0)/(2*np.pi*(f0+f1)/2)  )
print( (2*np.pi - phi1)/(2*np.pi*(f0+f1)/2)  )


#%%
from mpl_toolkits.axes_grid1 import Divider, Size
CQ1 = '#0D77BC'
CQ2 = '#DD5E27'


fig = plt.figure( figsize=(6,4))
h = [Size.Fixed(0.5), Size.Fixed(1.0), Size.Fixed(0.7), Size.Fixed(1.0)]
v = [Size.Fixed(0.5), Size.Fixed(1.0), Size.Fixed(0.7), Size.Fixed(1.0)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)

ax11 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=1))
ax13 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=3))


waittime2phase = 2*np.pi*datfiles[0].metadata['settings']['q1']['x90']['frequency']*1e-9

p1s = [np.array([-0.41588117,  0.04214265,  0.43105793, -0.58083328]),
 np.array([ 0.4073058 ,  0.04241945,  0.43609975, -0.65035388])]

xdatalist = [np.array([ 0.,  2.,  4.,  6.,  8., 10., 12., 14., 16., 18., 20., 22., 24.,
        26., 28., 30.]),
 np.array([ 0.,  2.,  4.,  6.,  8., 10., 12., 14., 16., 18., 20., 22., 24.,
        26., 28., 30.])]

zdatalist = [np.array([0.078, 0.012, 0.062, 0.202, 0.44 , 0.598, 0.788, 0.846, 0.832,
        0.628, 0.408, 0.218, 0.08 , 0.018, 0.1  , 0.218]),
 np.array([0.76 , 0.856, 0.798, 0.646, 0.498, 0.26 , 0.108, 0.024, 0.072,
        0.214, 0.406, 0.644, 0.81 , 0.824, 0.786, 0.636])]

pointstyle = ['o', 'o']
linestyle = ['-', '--']
colors = [CQ1, CQ2]
for i in range(2):
    xdata = xdatalist[i]
    zdata = zdatalist[i]
    ax11.scatter(xdata*waittime2phase, zdata, label=tempstr  , s=10, c=colors[i] , marker= pointstyle[i] )
    xxx = np.linspace(xdata[0], xdata[-1], (len(xdata)-1)*5+1  )
    ax11.plot(xxx*waittime2phase, osc_without_decay(xxx, *p1s[i] ) , 'k'  , linewidth=0.5 )
ax11.set_ylim(0,0.9)
ax11.set_yticks([0,0.9])
# ax11.set_xticks([0, 15, 30])
ax11.set_xticks([0, 4, 8])
ax11.plot(  (14*waittime2phase,)*2, (0,1), 'r--' , linewidth=0.5 )

# filename = 'Figure2d_calibrate' # 'Fig_2c_calib'
# plt.savefig(filename+'.pdf')




