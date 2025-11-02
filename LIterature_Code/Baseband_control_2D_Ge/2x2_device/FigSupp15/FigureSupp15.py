# this file is modified from Supp\CZ\CZsearch.py
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigSupp15','')
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
datadir = os.getcwd()  

DEG = np.pi/180
from scipy.fft import fft, fftfreq

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0


def fitosc_wo_decay(xdata, zdata, plot=False):
    zf = fft(zdata)
    xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 ]
    if plot:
        plt.figure()
    _, p1 = fit_data(xdata, zdata,func=osc_without_decay,p0=p0, plot=plot, return_cov=True)
    return p1



def find_phi( zdata_x, zdata_y,  xzdata_ref=[] ):
    
    
    if xzdata_ref ==  []:

        zall = np.hstack( (zdata_x, zdata_y))
        zallsort = np.sort(zall)    
        zmin = np.average(zallsort[0:3])
        zmax = np.average(zallsort[-3:])
    else:

        p1 = fitosc_wo_decay( xzdata_ref[0], xzdata_ref[1], plot=False)

        zmax = p1[2] + abs(p1[0]) 
        zmin = p1[2] - abs(p1[0]) 
    
    zmid = (zmax + zmin)/2
    zamp = (zmax - zmin)/2

    y = (zdata_y - zmid)
    x = (-zdata_x + zmid)    
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return x, y, r, phi, zmid, zamp



#%%

start_times = dict(down='2023-08-19\\10-16-51', up='2023-08-19\\10-17-08')

qubit_id = 1



data_findCZ_ref = {}  
zdatalist_ref = {}
xdatalist_ref = {}
for up_down_str in ['down', 'up']:
   
    for xyval in [ 'ref']:     
        key = up_down_str + '_' + xyval
        
        start_time = start_times[up_down_str]
        datfiles, fnames = get_data_from(start_time, start_time, rootfolder=datadir, only_complete = False) 
        data = datfiles[0]
        xdata, zdata = data.waittime_set.ndarray, data.su_S_North.ndarray
        xdatalist_ref[key] = xdata
        zdatalist_ref[key] = zdata











#%%

datadir = os.getcwd()

start_time = '2023-08-19\\10-29-16'
end_time = '2023-08-19\\10-29-52'


datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

p_lists = [   
            dict( step= dict(name='q1,q2/cz/gate_Jon/vpB12', val=-76, tag='vpB12'),  )   
]

data_lists = {}             
p_list = p_lists[0]


steptag = p_list['step']['tag']
stepval = p_list['step']['val']
stepflabel =  '_' + steptag + f'_{stepval:.3g}'

qubit_id = 1

data_findCZ = {}  
zdatalist = {}
xdatalist = {}
for up_down_str in ['down', 'up']:
    for xyval in ['x', 'y',]:     
        key = up_down_str + '_' + xyval
        flabel_findCZ = stepflabel + '_findCZ' +  f'_q{qubit_id}_q{int(3-qubit_id)}_' + key 

        
        for k in range(len(fnames)):
            if '_'+str(int(stepval)) in fnames[k] and key in fnames[k]:
                break
        
        data_findCZ[key] = datfiles[k]
        xdata,  zdata = data_findCZ[key].t_pulse_set.ndarray, data_findCZ[key].su_S_North.ndarray
        xdatalist[key] = xdata
        zdatalist[key] = zdata
data_lists[ stepflabel ] = dict(data_findCZ=data_findCZ, xdatalist=xdatalist, zdatalist=zdatalist)

###
    
d_dn = find_phi(zdatalist['down_x'],  zdatalist['down_y'],   xzdata_ref=[ xdatalist_ref['down_ref'],  zdatalist_ref['down_ref']  ] )     
d_up = find_phi(zdatalist['up_x'],  zdatalist['up_y'],  xzdata_ref=[ xdatalist_ref['up_ref'],  zdatalist_ref['up_ref']  ] )     
    
    
xdata =     xdatalist['down_x']
    


zmid, zamp = d_dn[4:6]



#%%
CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
MRSIZE = 2.5
SCSIZE = MRSIZE**2
from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(14,8))
h = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5),  Size.Fixed(1),Size.Fixed(1.5),  Size.Fixed(1),Size.Fixed(1.5),  Size.Fixed(1),Size.Fixed(1.5)]
v = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
ax31 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=3, ny=1))
ax51 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=5, ny=1))
ax71 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=7, ny=1))
ax91 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=9, ny=1))

ax13 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=3))
ax33 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=3, ny=3))
ax53 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=5, ny=3))
ax73 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=7, ny=3))
ax93 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=9, ny=3))

zmid, zamp = d_dn[4:6]

waittime2phase = 2*np.pi*datfiles[0].metadata['settings']['q1']['x90']['frequency']*1e-9

ax = ax13
ax.scatter( (xdatalist_ref['down_ref']+100)*waittime2phase/(np.pi),  zdatalist_ref['down_ref'] , c=CQ1, s=SCSIZE)
xxx = np.linspace(min(xdatalist_ref['down_ref'])-30, max(xdatalist_ref['down_ref'])+30, len(xdata)*5)
p1 = fitosc_wo_decay( xdatalist_ref['down_ref'],  zdatalist_ref['down_ref'], plot=False)
ax.plot((xxx+100)*waittime2phase/(np.pi), osc_without_decay(xxx, *p1), 'k', linewidth=0.5 , )
ax.set_xticks([8, 9, 10, 11 ])
ax.set_xlim(8,11)
ax.set_yticks([0.1, 0.9])  

ax = ax33
ax.scatter((xdatalist_ref['up_ref']+100)*waittime2phase/(np.pi),  zdatalist_ref['up_ref']  , c=CQ2, s=SCSIZE)
xxx = np.linspace(min(xdatalist_ref['up_ref'])-30, max(xdatalist_ref['up_ref'])+30, len(xdata)*5)
p1 = fitosc_wo_decay( xdatalist_ref['up_ref'],  zdatalist_ref['up_ref'], plot=False)
ax.plot((xxx+100)*waittime2phase/(np.pi), osc_without_decay(xxx, *p1), 'k', linewidth=0.5 )
ax.set_xticks([8, 9, 10, 11 ])
ax.set_xlim(8,11)
ax.set_yticks([0.1, 0.9])  

ax = ax53
ax.plot(xdata, d_dn[0]/d_dn[5] , '^--', c=CQ1, linewidth=0.5 , markersize=MRSIZE)
ax.plot(xdata, d_dn[1]/d_dn[5] , 'v-', c=CQ1, linewidth=0.5 , markersize=MRSIZE)
ax.set_xticks([18, 23, 28])
ax.set_yticks([-1, 0, 1])  


ax = ax73
ax.plot(xdata, d_up[0]/d_up[5] , '^--', c=CQ2, linewidth=0.5 , markersize=MRSIZE)
ax.plot(xdata, d_up[1]/d_up[5] , 'v-', c=CQ2, linewidth=0.5 , markersize=MRSIZE)
ax.set_xticks([18, 23, 28])
ax.set_yticks([-1, 0, 1])  

ax = ax11
phase_dn = np.array(d_dn[3])
phase_dn[xdata>21.2 ] += 2*np.pi
ax.plot(xdata, phase_dn/np.pi , 'o-', c=CQ1, linewidth=0.5 , markersize=MRSIZE)
phase_up = np.array(d_up[3])
phase_up[xdata>25.7 ] += 2*np.pi
ax.plot(xdata, (phase_up+2*np.pi)/np.pi , 'o-', c=CQ2, linewidth=0.5 , markersize=MRSIZE)
ax.set_xticks([18, 23, 28])
ax.set_yticks([0, 2, 4])  



ax = ax31

diff_phi = phase_up - phase_dn + 2*np.pi

ax.plot(xdata, diff_phi/np.pi ,'C0o ', linewidth=0.5 , markersize=MRSIZE)    
ax.plot( (min(xdata),max(xdata)), (1,1), 'r--', linewidth=0.5)

linearfit = np.polyfit(xdata,diff_phi/np.pi,1 )
ax.plot(  xdata, np.poly1d(linearfit)(xdata) ,'k' , linewidth=0.5 )
intersect1 = (1-linearfit[1])/linearfit[0]
ax.set_xticks([18, 23, 28])
ax.set_yticks([0.7, 1, 1.3])
# ax.grid()

ax = ax51
def voltage_to_J(voltage, j2v_params=dict(j_max=0.24e6, alpha=0.059 )  ):
    j_max, alpha = j2v_params['j_max'], j2v_params['alpha']
    '''
    voltage (float): fraction of max voltage. unit mV
    J = j_max when voltage == 1.0
    '''
    return j_max*np.exp(-alpha*(voltage))  # unit Hz

# hamming window
cz_params = np.array([ [43.47, -64], 
              [38.4, -66],
              [34.5, -68],
              [30.8, -70],
              [28.23, -72],
              [25.08, -74],
              [23.00, -76],
              [20.64, -78],
              [18.88, -80],
              [16.80, -82],
              [14.95, -84]
              ])
# hann window
cz_hann = np.array([ [46.87, -64], 
              [42.23, -66],
              [38.10, -68],
              [34.31, -70],
              [31.00, -72],
              [27.81, -74],
              [25.47, -76],
              [22.74, -78],
              [20.89, -80],
              [18.74, -82],
              [16.84, -84]
              ])

ax.scatter(cz_params[:,1], cz_params[:,0],  s=SCSIZE, c='C0')
ax.plot(cz_params[:,1] , 0.5*1e9/(voltage_to_J(cz_params[:,1])*1.08), 'C0' , linewidth=0.5 )

ax.scatter(cz_hann[:,1], cz_hann[:,0],  s=SCSIZE, c='C3')
ax.plot(cz_hann[:,1] , 0.5*1e9/(voltage_to_J(cz_hann[:,1])  ), 'C3' , linewidth=0.5 )
ax.set_xticks([-84, -74, -64])
ax.set_yticks([15, 30, 45])



# filename = 'FigureSupp15' # 'CZsearch'
# plt.savefig(filename+'.pdf')




