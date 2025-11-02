# this file is modified from Fig2\CZ_T2\ramsey_auto_postprocessing_twodecay.py

import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigSupp14','')
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

datadir = os.getcwd() 
from scipy.fft import fft, fftfreq
def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0    
def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0
def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay_white(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**1)+y0

def osc_with_twodecay(t, A, f, y0, phi, T2, T2w):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2 - (t/T2w))+y0

#%%

def fit_CZ_T2(start_time, end_time, phase0=None, plotflag=False):
    # start_time = '2023-08-31\\14-55-10'
    # end_time = '2023-08-31\\14-55-57'
    
    # qubit_id, up_down_str  = 1, 'down'
    
    datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
    data_precal = datfiles[0]
    data = datfiles[1]
    
    if fnames[0][-4::]=='down':
        qubit_id, up_down_str  = int(fnames[0][-9]), 'down'
    elif fnames[0][-2::]=='up':
        qubit_id, up_down_str  = int(fnames[0][-7]), 'up'
    else:
        raise ValueError
    
    vpB12 = data.metadata['settings']['q1,q2']['v_setpoints']['_1100_cz']['vpB12']
    t_exp_start = datetime.datetime.strptime(data.metadata['loop']['ts_start'],'%Y-%m-%d %H:%M:%S')
    t_exp_end = datetime.datetime.strptime(data.metadata['loop']['ts_end'],'%Y-%m-%d %H:%M:%S')
    t_exp = (t_exp_end-t_exp_start).total_seconds()
    
    xdata_precal = data_precal.t_cz_set.ndarray
    zdata_precal = data_precal.su_S_North.ndarray
    
        
    zf = fft(zdata_precal)
    xf = fftfreq(len(xdata_precal), d=abs(xdata_precal[1]-xdata_precal[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]  
    bounds = [[0.25, -np.inf, -np.inf, -np.inf], [0.5, np.inf, np.inf, np.inf]]
    if phase0 == None:
        phase0 = 0
    p0 = [abs(np.max(zdata_precal)-np.min(zdata_precal))/2, f0, 0.5, phase0,]
    _, p1_precal = fit_data(xdata_precal, zdata_precal,func=osc_without_decay,p0=p0, plot=False, return_cov=True, bounds=bounds)
    
    f_precal = p1_precal[1]
    # print(f0)
    # print(f_precal)
    # print(p1_precal[0])
    if f_precal < 0.005 or f_precal > 0.70 or abs(p1_precal[0])<0.15:
        raise ValueError('Fitting fail ')
    
    
    flabel =   f'_q{qubit_id}_ramsey_q{int(3-qubit_id)}_' + up_down_str
    
    resultlabel =   f'q{qubit_id}_q{int(3-qubit_id)}_' 
    resultlabel += 'dn' if up_down_str == 'down' else 'up'
    # print(start_time)
    # t_exp = 0
    if '2023-08-29' in end_time and '16-16-26' in end_time:
        # print(len(data.t_cz_set.ndarray))
        xdata = data.t_cz_set.ndarray[0:120]
        zdata = data.su_S_North.ndarray[0:120]    
        # xdata = data.t_cz_set.ndarray
        # zdata = data.su_S_North.ndarray        
        t_exp = np.round((len(xdata)/len(data.t_cz_set.ndarray))  *(t_exp_end-t_exp_start).total_seconds())
        # plotflag=True
    else:
        xdata = data.t_cz_set.ndarray
        zdata = data.su_S_North.ndarray
    delta_fs = [n/abs(xdata[0]-xdata[1]) - f_precal  for n in range(1,10) ]
    delta_f = delta_fs[np.argmin(np.abs(delta_fs))]
    # print(delta_f)
    # delta_f = 0.08
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, delta_f, 0.5, 0 , 300]
    # print(delta_fs)
    _, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=False, return_cov=True)
    
    # plt.figure()
    # _, p1_white = fit_data(xdata, zdata,func=osc_with_decay_white,p0=p0, plot=False, return_cov=True)
    # print( f'T2 = {abs(p1[4]):.3g} ns, or {abs(p1_white[4]):.3g} ns' )


    
    f_finals = [n/abs(xdata[0]-xdata[1])+p1[1]  for n in range(1,10) ] + [n/abs(xdata[0]-xdata[1])-p1[1]  for n in range(1,10) ]
    f_final = f_finals[np.argmin(np.abs(f_finals-f_precal))]
    
    # plt.figure()
    _, p1_alpha = fit_data(xdata, zdata,func=osc_with_decay_alpha,p0=list(p1) + [2], plot=False, return_cov=True)
    
    # print(p1_alpha)
    f_final_alphas = [n/abs(xdata[0]-xdata[1])+p1_alpha[1]  for n in range(1,10) ] + [n/abs(xdata[0]-xdata[1])-p1_alpha[1]  for n in range(1,10) ]
    f_final_alpha = f_final_alphas[np.argmin(np.abs(f_final_alphas-f_precal))]

    _, p1_two = fit_data(xdata, zdata,func=osc_with_twodecay,p0=p0+[2e3], plot=False, return_cov=True)
    print(f'At vB12={vpB12:.3g}mV:')
    print(f'T2 = {abs(p1[4]):.3g} ns (fit to gaussian)')
    print(f'T2_gaussian = {abs(p1_two[4]):.3g} ns, T2_exp = {abs(p1_two[5]):.3g} ns (fit to gaussian * exp model) \n' )
    
    ###
    if plotflag:
        fig = plt.figure(figsize=(6,3))
        plt.subplots_adjust(top=0.8 , left=0.05, right=0.95)
        plt.suptitle(data.metadata['location'][:19] + f'  vpB12={vpB12:.3g}mV'+ '\n' +  \
                     flabel +  f', T2* = {abs(p1[4]):.3g} ns'  + '\n' + \
                     f'f_precal = {f_precal*1e3:.3g} MHz, '  +  f'f_final = {f_final*1e3:.3g} MHz' + '\n'  \
                     )
        # print(p1, f_final)
        plt.subplot(121)
        plt.scatter(xdata_precal, zdata_precal,)
        xx = np.linspace(min(xdata_precal), max(xdata_precal), (len(xdata_precal)-1)*5+1)
        plt.plot(xx, osc_without_decay(xx, *p1_precal))
            
        plt.subplot(122)
        plt.scatter(xdata, zdata)
        xx = np.linspace(min(xdata), max(xdata), (len(xdata)-1)*5+1)
        plt.plot(xx, osc_with_decay(xx, *p1))

    return  resultlabel, (vpB12, qubit_id, up_down_str, t_exp), f_precal, f_final, f_final_alpha, p1_precal, p1, p1_alpha, (xdata, zdata, ), p1_two


T2_all = dict(q1_q2_dn=[], q1_q2_up=[], 
                q2_q1_dn=[], q2_q1_up=[], )
T2alpha_all = dict(q1_q2_dn=[], q1_q2_up=[], 
                q2_q1_dn=[], q2_q1_up=[], )
Texp_all = dict(q1_q2_dn=[], q1_q2_up=[], 
                q2_q1_dn=[], q2_q1_up=[], )
Tend_all = dict(q1_q2_dn=[], q1_q2_up=[], 
                q2_q1_dn=[], q2_q1_up=[], )

results = dict()

plt.close('all')


start_end_times = [
    ['2023-08-29\\14-24-03', '2023-08-29\\14-24-44'], # -25mV
    ['2023-08-29\\14-25-30', '2023-08-29\\14-26-12'],
    ['2023-08-29\\14-26-59', '2023-08-29\\14-27-42'],
    ['2023-08-29\\14-29-05', '2023-08-29\\14-29-49'],

    ['2023-08-29\\14-16-43', '2023-08-29\\14-17-23'], # -45mV
    ['2023-08-29\\14-18-08', '2023-08-29\\14-18-49'],
    ['2023-08-29\\14-19-39', '2023-08-29\\14-20-22'],
    ['2023-08-29\\14-21-46', '2023-08-29\\14-22-31'],   

    
    ['2023-08-29\\14-01-58', '2023-08-29\\14-02-40'], # -55mV
    ['2023-08-29\\14-03-24', '2023-08-29\\14-04-07'],
    ['2023-08-29\\16-27-35', '2023-08-29\\16-28-22'],
    ['2023-08-29\\14-07-02', '2023-08-29\\14-07-47'],       
    

    ['2023-08-29\\14-58-53', '2023-08-29\\14-59-34'], # -60mV
    ['2023-08-29\\14-59-59', '2023-08-29\\15-00-41'],
    ['2023-08-29\\16-23-33', '2023-08-29\\16-24-20'],
    ['2023-08-29\\15-07-25', '2023-08-29\\15-08-13'],

    ['2023-08-29\\15-21-47', '2023-08-29\\15-22-28'], # -65mV
    ['2023-08-29\\15-23-42', '2023-08-29\\15-24-26'],
    ['2023-08-29\\16-20-03', '2023-08-29\\16-20-50'],
    ['2023-08-29\\15-29-50', '2023-08-29\\15-30-33'],    

    ['2023-08-29\\15-33-56', '2023-08-29\\15-34-36'], # -70mV
    ['2023-08-29\\15-35-26', '2023-08-29\\15-36-09'],
    ['2023-08-29\\16-15-38', '2023-08-29\\16-16-26'],
    ['2023-08-29\\15-40-09', '2023-08-29\\15-40-52'],   

    ['2023-08-29\\15-44-58', '2023-08-29\\15-45-38'], # -75mV
    ['2023-08-29\\15-46-15', '2023-08-29\\15-46-57'],
    ['2023-08-29\\15-47-58', '2023-08-29\\15-48-38'],
    ['2023-08-29\\15-52-36', '2023-08-29\\15-53-19'],   

    ['2023-08-29\\15-55-55', '2023-08-29\\16-05-35'],  # -80mV
    ['2023-08-29\\16-06-38', '2023-08-29\\16-07-22'],
    ['2023-08-29\\16-08-36', '2023-08-29\\16-09-19'],
    ['2023-08-29\\16-11-03', '2023-08-29\\16-11-47'],
    
  ]

example_plots = dict()
for i in range(len(start_end_times)):
    start_time, end_time = start_end_times[i]
    phase0 = 0
    plotflag = False
    if start_time == '2023-08-29\\14-59-59':
        phase0 = 3
    elif start_time == '2023-08-29\\14-18-08':
        phase0 = 3
        plotflag = False
        
    result  = fit_CZ_T2(start_time, end_time, phase0 = phase0, plotflag=plotflag)
    
    vpB12 = result[1][0]
    p1 = result[6]
    p1_alpha = result[7]
    f_final, f_final_alpha = result[3], result[4]
    
    results[result[0]] = dict(vpB12=vpB12, result=result, 
                              T2=p1[4], T2_alpha=p1_alpha[4], alpha=p1_alpha[5])

    T2_all[result[0]] +=  [[vpB12, p1[4], f_final]]
    T2alpha_all[result[0]] +=  [[vpB12, p1_alpha[4], f_final_alpha, p1_alpha[5]] ]
    Texp_all[result[0]] += [[vpB12, result[1][3]]]                                    
    Tend_all[result[0]] += [[vpB12, result[8][0][-1]]]  ####
    if -66 <= vpB12 <= -64:
        example_plots[result[0]] = dict(xdata=result[8][0], zdata=result[8][1], p1=p1)   
    

resultlabels = ['q1_q2_dn', 'q1_q2_up', 'q2_q1_dn', 'q2_q1_up']
for rlabel in resultlabels:
    T2_all[rlabel + '_vpB12'] = [ T2_all[rlabel][i] [0]  for i in range(len(T2_all[rlabel])) ]
    T2_all[rlabel + '_T2'] = [ np.abs(T2_all[rlabel][i] [1])  for i in range(len(T2_all[rlabel])) ]
    T2_all[rlabel + '_f'] = [ T2_all[rlabel][i] [2]  for i in range(len(T2_all[rlabel])) ]

    T2alpha_all[rlabel + '_vpB12'] = [ T2alpha_all[rlabel][i] [0]  for i in range(len(T2alpha_all[rlabel])) ]
    T2alpha_all[rlabel + '_T2'] = [np.abs( T2alpha_all[rlabel][i] [1])  for i in range(len(T2alpha_all[rlabel])) ]
    T2alpha_all[rlabel + '_f'] = [ T2alpha_all[rlabel][i] [2]  for i in range(len(T2alpha_all[rlabel])) ]
    T2alpha_all[rlabel + '_alpha'] = [ T2alpha_all[rlabel][i] [3]  for i in range(len(T2alpha_all[rlabel])) ]

    Texp_all[rlabel + '_vpB12'] = [ Texp_all[rlabel][i] [0]  for i in range(len(Texp_all[rlabel])) ]
    Texp_all[rlabel + '_Texp'] = [ Texp_all[rlabel][i] [1]  for i in range(len(Texp_all[rlabel])) ]
    Tend_all[rlabel + '_Tend'] = [ Tend_all[rlabel][i] [1]  for i in range(len(Tend_all[rlabel])) ]   ####
# q1_q2_dn=[], q1_q2_up=[], 
                # q2_q1_dn=[], q2_q1_up=[], )


#%%
CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
MRSIZE = 10 #2.5
SCSIZE = MRSIZE**2

x = T2_all['q1_q2_dn'+ '_vpB12']

plt.figure(figsize=(18,6))
plt.subplot(141)
y = T2_all['q1_q2_dn'+ '_T2']
plt.plot(x, y, marker='v', c=CQ1, markersize=MRSIZE,  linestyle='-')

y = T2_all['q1_q2_up'+ '_T2']
plt.plot(x, y, marker='^', c=CQ1, markersize=MRSIZE, linestyle='-')

y = T2_all['q2_q1_dn'+ '_T2']
plt.plot(x, y, marker='v', c=CQ2, markersize=MRSIZE, linestyle='-')

y = T2_all['q2_q1_up'+ '_T2']
plt.plot(x, y, marker='^', c=CQ2, markersize=MRSIZE, linestyle='-')



plt.subplot(142)
y = T2_all['q1_q2_dn'+ '_f']
plt.plot(x, y, marker='v', c=CQ1, markersize=MRSIZE,  linestyle='-')

y = T2_all['q1_q2_up'+ '_f']
plt.plot(x, y, marker='^', c=CQ1, markersize=MRSIZE, linestyle='-')

y = T2_all['q2_q1_dn'+ '_f']
plt.plot(x, y, marker='v', c=CQ2, markersize=MRSIZE, linestyle='-')

y = T2_all['q2_q1_up'+ '_f']
plt.plot(x, y, marker='^', c=CQ2, markersize=MRSIZE, linestyle='-')


plt.subplot(143)
Jq1 = np.array(T2_all['q1_q2_up'+ '_f']) - np.array(T2_all['q1_q2_dn'+ '_f'])
plt.plot(x, Jq1, marker='o', c=CQ1, markersize=MRSIZE,  linestyle='-')

rlabel = 'q2_q1_dn'
Jq2 = np.array(T2_all['q2_q1_up'+ '_f']) - np.array(T2_all['q2_q1_dn'+ '_f'])
plt.plot(x, Jq2, marker='o', c=CQ2, markersize=MRSIZE, linestyle='-')

J = (Jq1+Jq2)/2



plt.subplot(144)
rlabel = 'q1_q2_dn'
x = T2_all[rlabel+ '_vpB12']
y = T2_all[rlabel+ '_T2']
plt.plot(x, y*J*2, marker='v', c=CQ1, markersize=MRSIZE,  linestyle='-')

rlabel = 'q1_q2_up'
y = T2_all[rlabel+ '_T2']
plt.plot(x, y*J*2, marker='^', c=CQ1, markersize=MRSIZE, linestyle='-')

rlabel = 'q2_q1_dn'
y = T2_all[rlabel+ '_T2']
plt.plot(x, y*J*2, marker='v', c=CQ2, markersize=MRSIZE, linestyle='-')

rlabel = 'q2_q1_up'
y = T2_all[rlabel+ '_T2']
plt.plot(x, y*J*2, marker='^', c=CQ2, markersize=MRSIZE, linestyle='-')





#%%


CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
MRSIZE = 2.5
SCSIZE = MRSIZE**2
LW = 0.5
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


x = T2_all['q1_q2_dn'+ '_vpB12']
 

ax = ax13
ax.grid(linewidth=0.5) 
y = np.array(T2_all['q1_q2_dn'+ '_T2'])
ax.plot(x, y/1e3, marker='v', c=CQ1, markersize=MRSIZE,  linestyle='-', linewidth=LW)
y = np.array(T2_all['q1_q2_up'+ '_T2'])
ax.plot(x, y/1e3, marker='^', c=CQ1, markersize=MRSIZE, linestyle='--', linewidth=LW)
y = np.array(T2_all['q2_q1_dn'+ '_T2'])
ax.plot(x, y/1e3, marker='v', c=CQ2, markersize=MRSIZE, linestyle='-', linewidth=LW)
y = np.array(T2_all['q2_q1_up'+ '_T2'])
ax.plot(x, y/1e3, marker='^', c=CQ2, markersize=MRSIZE, linestyle='--', linewidth=LW)
ax.set_xticks([-80, -60, -40, -20 ])
# ax.set_xlim(-90, -10 )
ax.set_yticks([0, 1, 2, 3, 4])  
ax.set_ylim(0, 4.5 )




ax = ax11
xdata = example_plots['q1_q2_dn']['xdata']
zdata = example_plots['q1_q2_dn']['zdata']
p1 = example_plots['q1_q2_dn']['p1']
ax.plot(xdata, zdata, marker='v', c=CQ1, markersize=MRSIZE,  linestyle='', linewidth=LW)
xxx = np.linspace(min(xdata), max(xdata), len(xdata)*5)
ax.plot(xxx, osc_with_decay(xxx, *p1), 'k-', linewidth=LW)
# ax.set_xticks([-80, -60, -40, -20 ])
ax.set_yticks([ 0.2, 0.8])  
ax.set_ylim(0.1, 0.95 )


ax = ax31
xdata = example_plots['q1_q2_up']['xdata']
zdata = example_plots['q1_q2_up']['zdata']
p1 = example_plots['q1_q2_up']['p1']
ax.plot(xdata, zdata, marker='^', c=CQ1, markersize=MRSIZE,  linestyle='', linewidth=LW)
xxx = np.linspace(min(xdata), max(xdata), len(xdata)*5)
ax.plot(xxx, osc_with_decay(xxx, *p1), 'k--', linewidth=LW)
# ax.set_xticks([-80, -60, -40, -20 ])
ax.set_yticks([ 0.2, 0.8])  
ax.set_ylim(0.1, 0.95 )


ax = ax51
xdata = example_plots['q2_q1_dn']['xdata']
zdata = example_plots['q2_q1_dn']['zdata']
p1 = example_plots['q2_q1_dn']['p1']
ax.plot(xdata, zdata, marker='v', c=CQ2, markersize=MRSIZE,  linestyle='', linewidth=LW)
xxx = np.linspace(min(xdata), max(xdata), len(xdata)*5)
ax.plot(xxx, osc_with_decay(xxx, *p1), 'k-', linewidth=LW)
ax.set_xticks([0, 1500, 3000 ])
ax.set_yticks([ 0.2, 0.8])  
ax.set_ylim(0.1, 0.95 )


ax = ax71
xdata = example_plots['q2_q1_up']['xdata']
zdata = example_plots['q2_q1_up']['zdata']
p1 = example_plots['q2_q1_up']['p1']
ax.plot(xdata, zdata, marker='^', c=CQ2, markersize=MRSIZE,  linestyle='', linewidth=LW)
xxx = np.linspace(min(xdata), max(xdata), len(xdata)*5)
ax.plot(xxx, osc_with_decay(xxx, *p1), 'k--', linewidth=LW)
ax.set_xticks([0, 250, 500 ])
ax.set_yticks([ 0.2, 0.8])  
ax.set_ylim(0.1, 0.95 )





# filename = 'exchange_on_T2'
# plt.savefig(filename+'.pdf')


# ramsey experiment time (Seconds)
 # 'q1_q2_dn_Texp': np.array([42.0, 41.0, 40.0, 21.0, 39.0, 20.0, 18.0, 21.0]),
 # 'q1_q2_up_Texp': np.array([44.0, 47.0, 50.0, 28.0, 58.0, 31.0, 32.0, 33.0]),
 # 'q2_q1_dn_Texp': np.array([80.0, 80.0, 133.0, 125.0, 126.0, 108.0, 39.0, 44.0]),
 # 'q2_q1_up_Texp': np.array([84.0, 88.0, 90.0, 108.0, 50.0, 52.0, 56.0, 38.0])





ax = ax91
ax.grid(linewidth=0.5) 
y = np.array(T2_all['q1_q2_dn'+ '_T2'])
ax.plot(x, y/1e3, marker='v', c=CQ1, markersize=MRSIZE,  linestyle='', linewidth=LW)
y = np.array(T2_all['q1_q2_up'+ '_T2'])
ax.plot(x, y/1e3, marker='^', c=CQ1, markersize=MRSIZE, linestyle='', linewidth=LW)
y = np.array(T2_all['q2_q1_dn'+ '_T2'])
ax.plot(x, y/1e3, marker='v', c=CQ2, markersize=MRSIZE, linestyle='', linewidth=LW)
y = np.array(T2_all['q2_q1_up'+ '_T2'])
ax.plot(x, y/1e3, marker='^', c=CQ2, markersize=MRSIZE, linestyle='', linewidth=LW)
ax.set_xticks([-80, -60, -40, -20 ])
# ax.set_xlim(-90, -10 )
ax.set_yticks([0, 1, 2, 3, 4])  
ax.set_ylim(0, 4.5 )



def twoQlevels(xxx, tso_c, tso_theta, Jxxx0, L1_0, L1_1, L2_0, L2_1):
    sq2 = np.sqrt(2)
    U = 2.56e-3*1.602e-19 / (1e9*6.62e-34)
    tso_phi = 0 #0.9
    tso_dir = np.array([ np.sin(tso_theta)*np.cos(tso_phi), 
               np.sin(tso_theta)*np.sin(tso_phi)  , 
               np.cos(tso_theta)])
    L1 = L1_0, L1_1
    L2 = L2_0, L2_1
    evs = np.zeros((6,len(xxx)))
    evecs = np.zeros((6,6,len(xxx)),dtype=complex)          
    for i in range(len(xxx)):
        Jxxx = Jxxx0*np.exp(-0.059*(xxx[i]-10) ) *1e-3
        t_all = np.sqrt(U * Jxxx /4)
        dEz = (np.poly1d(L2)(xxx[i]) - np.poly1d(L1)(xxx[i]))/2e3
        sEz = (np.poly1d(L2)(xxx[i]) + np.poly1d(L1)(xxx[i]))/2e3
        e = 0
        tc, tx, ty, tz = t_all* np.array([ np.cos(tso_c/2)   ] + list(np.sin(tso_c/2)*tso_dir)  )
        H = np.array([[       U+e,           0, -ty+1j*tx,  sq2*tc, -1j*sq2*tz, -ty-1j*tx],
                      [         0,         U-e, -ty+1j*tx,  sq2*tc, -1j*sq2*tz, -ty-1j*tx],
                      [ -ty-1j*tx,   -ty-1j*tx,       sEz,       0,          0,         0],
                      [    sq2*tc,      sq2*tc,         0,       0,        dEz,         0],
                      [ 1j*sq2*tz,   1j*sq2*tz,         0,     dEz,          0,         0],
                      [ -ty+1j*tx,   -ty+1j*tx,         0,       0,          0,      -sEz],
        ]) 
    
        ev, evec = np.linalg.eigh(H)
        evs[:,i] = np.real(ev)
        evecs[:,:,i] = evec
    evs *= 1000    
    return evs

def twoQfreq(xxx_in, tso_c, tso_theta, Jxxx0, L1_0, L1_1, L2_0, L2_1):
    Lxxx = round(len(xxx_in)/4)
    collect = []
    for N in range(4):
        xxx = xxx_in[N*Lxxx: N*Lxxx+Lxxx]
        evs = twoQlevels(xxx, tso_c, tso_theta, Jxxx0, L1_0, L1_1, L2_0, L2_1)
        if N == 0:
            collect += [  evs[1,:]-evs[0,:]  ]
        elif N==1:
            collect += [  evs[3,:]-evs[2,:]  ]
        elif N==2:
            collect += [  evs[2,:]-evs[0,:]  ]
        else:
            collect += [  evs[3,:]-evs[1,:]  ]
    return np.hstack(collect)


p1 = np.array([ 1.26609665e+00,  5.60351262e-01,  3.04382405e-01, -1.05262210e-01,
                        4.98468800e+01, -8.97941921e-02,  9.38339262e+01])
xxx0 = np.linspace(-80.5, -24.5, 1001)



xdata = np.hstack([xxx0, xxx0, xxx0, xxx0] )  

Ndiv = len(xxx0)
zdata_div = {}
for i in range(4):
    zdata_div[str(i)] = twoQfreq(xdata, *p1) [i*Ndiv:i*Ndiv+Ndiv] *1e-3 # in GHz



def meas2ref(T2, Texp, Tevolve):
    Texp_ref = 1000
    Tevolve_ref = 1e-6
    return T2*np.sqrt(  (np.log( 0.401/(Tevolve/Texp)  ))/(np.log( 0.401/(Tevolve_ref/Texp_ref)  ))  )
    

def ref2meas(T2_ref, Texp, Tevolve):
    Texp_ref = 1000
    Tevolve_ref = 1e-6
    return T2_ref*np.sqrt(  (np.log( 0.401/(Tevolve_ref/Texp_ref)  ))/(np.log( 0.401/(Tevolve/Texp)  ))  )


from scipy.interpolate import CubicSpline
from scipy.optimize import minimize


def tominimize(x_param):
    delta_vpB12, delta_L1, delta_L2 = x_param
    

    deltaRate = np.array([])
    for i, key in enumerate(['q1_q2_dn', 'q1_q2_up', 'q2_q1_dn', 'q2_q1_up']):
        dEdvpB12_1_interp_data = (zdata_div[str(i)][2:]-zdata_div[str(i)][:-2])/(xxx0[2]-xxx0[0])
        dEdvpB12_2_interp_data = np.diff(zdata_div[str(i)],2)/((xxx0[1]-xxx0[0])**2)
        dEdvpB12_1_interp = CubicSpline(xxx0[1:-1], dEdvpB12_1_interp_data)
        dEdvpB12_2_interp = CubicSpline(xxx0[1:-1], dEdvpB12_2_interp_data)
        
        x = T2_all[key + '_vpB12']
        dEdvpB12_1 = dEdvpB12_1_interp(x)
        dEdvpB12_2 = dEdvpB12_2_interp(x)
        dEdvpB12_2 = 0
        T2 = abs( np.sqrt(2)/(2*np.pi*(dEdvpB12_1*delta_vpB12+0.5*dEdvpB12_2*delta_vpB12*delta_vpB12))   )
        
        
        
        p1_dL1 = np.array([0,0,0,0, delta_L1,0,0])
        evs = (twoQlevels(x, *(p1+p1_dL1)) - twoQlevels(x, *p1) )*1e-3
        if i == 0:
            tempvar =  evs[1,:]-evs[0,:]  
        elif i==1:
            tempvar =  evs[3,:]-evs[2,:] 
        elif i==2:
            tempvar =  evs[2,:]-evs[0,:] 
        else:
            tempvar =  evs[3,:]-evs[1,:] 
        T2_L1 = abs( np.sqrt(2)/(2*np.pi*tempvar  )   )
        
        p1_dL2 = np.array([0,0,0,0,0,0, delta_L2])
        evs = (twoQlevels(x, *(p1+p1_dL2)) - twoQlevels(x, *p1) )*1e-3
        if i == 0:
            tempvar =  evs[1,:]-evs[0,:]  
        elif i==1:
            tempvar =  evs[3,:]-evs[2,:] 
        elif i==2:
            tempvar =  evs[2,:]-evs[0,:] 
        else:
            tempvar =  evs[3,:]-evs[1,:]         
        T2_L2 = abs( np.sqrt(2)/(2*np.pi*tempvar  )   )
        T2 = 1/np.sqrt(1/T2**2 + 1/T2_L1**2 + 1/T2_L2**2) 

            
        y = np.array(T2_all[key + '_T2'])
        Texp = np.array(Texp_all[key+ '_Texp'])
        Tend = np.array(Tend_all[key+ '_Tend'])
        # deltaRate = np.concatenate((deltaRate, 1/T2 - 1/meas2ref(y, Texp*(y/Tend), y*1e-9) )) 
        deltaRate = np.concatenate((deltaRate, 1/ref2meas(T2, T2*(Texp/Tend), T2*1e-9)  - 1/y )) 
        # deltaRate = np.concatenate((deltaRate, 1/T2 - 1/y ) ) 
        
    return np.sum( deltaRate**2  )
# plt.plot(valvP4, T2,label='vP4 dephasing',color='k')
res_minimized = minimize(tominimize, x0=np.array([0.7, 50e-3, 50e-3])   ) #, method='Nelder-Mead', tol=1e-6)

delta_vpB12, delta_L1, delta_L2 = res_minimized['x']
# delta_vpB12 = 0.78517063602
# delta_L1 = 0.05023395514319382
# delta_L2 = 0.0634475149594602

print(f'T2_L1 = { abs( np.sqrt(2)/(2*np.pi*delta_L1  )   ):.4g}'  )
print(f'T2_L2 = { abs( np.sqrt(2)/(2*np.pi*delta_L2  )   ):.4g}'  )

print(f'S_1/f_L1 = {   delta_L1**2/np.log(0.401/(1e-6/1000)) :.4g} MHz^2'  )
print(f'S_1/f_L2 = {  delta_L2**2/np.log(0.401/(1e-6/1000))   :.4g}'  )
print(f'S_1/f_vpB12 = { delta_vpB12**2/np.log(0.401/(1e-6/1000))   :.4g}'  )

# plt.figure()
# ax=plt.gca()

ax = ax91



for i, key in enumerate(['q1_q2_dn', 'q1_q2_up', 'q2_q1_dn', 'q2_q1_up']):
    dEdvpB12_1 = (zdata_div[str(i)][2:]-zdata_div[str(i)][:-2])/(xxx0[2]-xxx0[0])
    dEdvpB12_2 = np.diff(zdata_div[str(i)],2)/((xxx0[1]-xxx0[0])**2)
    dEdvpB12_2 = 0
    T2 = abs( np.sqrt(2)/(2*np.pi*(dEdvpB12_1*delta_vpB12+0.5*dEdvpB12_2*delta_vpB12*delta_vpB12))   )
    
    x = xxx0[1:-1]
    p1_dL1 = np.array([0,0,0,0, delta_L1,0,0])
    evs = (twoQlevels(x, *(p1+p1_dL1)) - twoQlevels(x, *p1) )*1e-3
    if i == 0:
        tempvar =  evs[1,:]-evs[0,:]  
    elif i==1:
        tempvar =  evs[3,:]-evs[2,:] 
    elif i==2:
        tempvar =  evs[2,:]-evs[0,:] 
    else:
        tempvar =  evs[3,:]-evs[1,:] 
    T2_L1 = abs( np.sqrt(2)/(2*np.pi*tempvar  )   )
    
    p1_dL2 = np.array([0,0,0,0,0,0, delta_L2])
    evs = (twoQlevels(x, *(p1+p1_dL2)) - twoQlevels(x, *p1) )*1e-3
    if i == 0:
        tempvar =  evs[1,:]-evs[0,:]  
    elif i==1:
        tempvar =  evs[3,:]-evs[2,:] 
    elif i==2:
        tempvar =  evs[2,:]-evs[0,:] 
    else:
        tempvar =  evs[3,:]-evs[1,:]         
    T2_L2 = abs( np.sqrt(2)/(2*np.pi*tempvar  )   )
    
    T2 = 1/np.sqrt(1/T2**2 + 1/T2_L1**2 + 1/T2_L2**2) 

        

    x = T2_all[key+ '_vpB12'].copy()
    Texp = Texp_all[key+ '_Texp'].copy()
    Tend = Tend_all[key+ '_Tend'].copy()
    if min(xxx0) < min(x):
        x += [min(xxx0) ]
        Texp += [Texp[-1]]
        Tend += [Tend[-1]]
    if max(xxx0) > max(x):
        x = [max(xxx0) ] + x
        Texp = [Texp[0]] + Texp
        Tend = [Tend[0]] + Tend    
    x = np.array(x)
    Texp = np.array(Texp)
    Tend = np.array(Tend)
        
    tempvar = CubicSpline( x[::-1] , (Texp/Tend)[::-1]  )
    T2 = ref2meas(T2, T2 *tempvar(xxx0[1:-1]), T2*1e-9) 
    Lstyle = dict(q1_q2_dn='C0-', q1_q2_up='C0--', q2_q1_dn='C1-', q2_q1_up='C1--',  )
    ax.plot(xxx0[1:-1], T2*1e-3   , Lstyle[key] ,linewidth=LW  )        

ax.set_xlim(-83, -20 )


# filename = 'FigSupp14bd' #  'exchange_on_T2'
# plt.savefig(filename+'.pdf')





# from Fig_2c_pulse.py
vgate_raw = np.array([-14.        , -33.19103993, -34.07529296, -36.46701613,
       -39.79035263, -43.4992084 , -47.231775  , -50.79479297,
       -54.10208114, -57.12518634, -59.86383871, -62.33022131,
       -64.54107231, -66.51388733, -68.26519338, -69.80984705,
       -71.16082699, -72.32925436, -73.32450833, -74.15437145,
       -74.82517308, -75.34191688, -75.70838623, -75.92722594,
       -76.        , -75.92722594, -75.70838623, -75.34191688,
       -74.82517308, -74.15437145, -73.32450833, -72.32925436,
       -71.16082699, -69.80984705, -68.26519338, -66.51388733,
       -64.54107231, -62.33022131, -59.86383871, -57.12518634,
       -54.10208114, -50.79479297, -47.231775  , -43.4992084 ,
       -39.79035263, -36.46701613, -34.07529296, -14.        ]) + 11.249341

vgate = CubicSpline(np.arange(len(vgate_raw)), vgate_raw)

t = np.linspace(0, 46, 1001)


from scipy.interpolate import CubicSpline



Texp_2qRB = 26790
delta_vpB12 = delta_vpB12*np.sqrt(np.log(0.401/(0.108e-6/Texp_2qRB))  /  np.log(0.401/(1e-6/1000)))
delta_L1 = delta_L1*np.sqrt(np.log(0.401/(0.108e-6/Texp_2qRB))  /  np.log(0.401/(1e-6/1000)))
delta_L2 = delta_L2*np.sqrt(np.log(0.401/(0.108e-6/Texp_2qRB))  /  np.log(0.401/(1e-6/1000)))
dE = dict(dvpB12=[], dL1=[], dL2=[])
for i in range(4):
    dEdvpB12_1 = (zdata_div[str(i)][2:]-zdata_div[str(i)][:-2])/(xxx0[2]-xxx0[0])
    dEdvpB12_2 = np.diff(zdata_div[str(i)],2)/((xxx0[1]-xxx0[0])**2)
    dEdvpB12_2 = 0
    tempvar = dEdvpB12_1*delta_vpB12+0.5*dEdvpB12_2*delta_vpB12*delta_vpB12  # unit: MHz
    dE['dvpB12'] += [   CubicSpline(xxx0[1:-1], tempvar)     ]
    
    x = xxx0[1:-1]
    p1_dL1 = np.array([0,0,0,0, delta_L1,0,0])
    evs = (twoQlevels(x, *(p1+p1_dL1)) - twoQlevels(x, *p1) )*1e-3
    if i == 0:
        tempvar =  evs[1,:]-evs[0,:]  
    elif i==1:
        tempvar =  evs[3,:]-evs[2,:] 
    elif i==2:
        tempvar =  evs[2,:]-evs[0,:] 
    else:
        tempvar =  evs[3,:]-evs[1,:] 
    dE['dL1'] += [   CubicSpline(xxx0[1:-1], tempvar)     ]
    
    p1_dL2 = np.array([0,0,0,0,0,0, delta_L2])
    evs = (twoQlevels(x, *(p1+p1_dL2)) - twoQlevels(x, *p1) )*1e-3
    if i == 0:
        tempvar =  evs[1,:]-evs[0,:]  
    elif i==1:
        tempvar =  evs[3,:]-evs[2,:] 
    elif i==2:
        tempvar =  evs[2,:]-evs[0,:] 
    else:
        tempvar =  evs[3,:]-evs[1,:]         
    dE['dL2'] += [   CubicSpline(xxx0[1:-1], tempvar)     ]



# plt.figure()
# plt.plot( t, vgate(t)   )
sigma_phis = {}
for key in dE.keys():
    sigma_phis[key] = []
    for i in range(4):
        tempvar = 2*np.pi* np.sum(   dE[key][i]( vgate(t)  )*(t[1]-t[0])       )    
        sigma_phis[key] += [  tempvar   ]
        max_freq_diff =  np.max( np.abs(dE[key][i]( vgate(t)  ))) 
        avg_freq_diff = np.mean(   dE[key][i]( vgate(t)  )  )

    
    sigma_phis[key] = np.array(sigma_phis[key])


sigma_freq_dot1 = delta_L1*1e-3*np.sqrt(np.log(0.401/(0.108e-6/Texp_2qRB))  /  np.log(0.401/(1e-6/1000)))
sigma_q1 = 2*np.pi*sigma_freq_dot1*(16 + 15 + 15 + 16 )


sigma_freq_dot2 = delta_L2*1e-3*np.sqrt(np.log(0.401/(0.108e-6/Texp_2qRB))  /  np.log(0.401/(1e-6/1000)))
sigma_q2 = 2*np.pi*sigma_freq_dot2*(16 + 15 + 15 + 16  ) 

s0 = np.eye(2)
sx = np.array([[0, 1],
               [1, 0]])
sy = np.array([[0, -1j],
               [1j, 0]])
sz = np.array([[1, 0],
               [0,-1]])

for key in sigma_phis.keys():
    print('\n' + 'fluctuation in ', key)
    print('phase accumulation fluctuation for four energy levels = ' + ''.join(['{:.3g}, ']*len(sigma_phis[key])).format(*sigma_phis[key]) + ' rad')
    UiH =  np.diag([1,1,1,-1])
    tempvar = sigma_phis[key]
    R_exchange = np.exp(1j*np.diag([0, tempvar[2],tempvar[0], tempvar[0]+tempvar[3]]))
    R_exchange[3] *= -1
    
    if key in ['dL1','dL2']:
        R_singleQ = np.kron(  s0*np.cos(sigma_q1/2 ) + 1j*(sz )*np.sin(sigma_q1/2), 
                              s0*np.cos(sigma_q2/2 ) + 1j*(sz )*np.sin(sigma_q2/2)  )
        R_exchange = R_singleQ @ R_exchange
    
    d = 4
    print('Incoherent error contribution = ', 1-(np.abs(np.trace( R_exchange @ UiH  ))**2 + d )/(d  *(d+1))   )




