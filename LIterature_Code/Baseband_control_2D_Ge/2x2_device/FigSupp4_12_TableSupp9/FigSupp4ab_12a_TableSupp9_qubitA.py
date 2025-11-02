# this file is modified from Supp/fit_tc_angle/AC14_v0.py

#%%
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigSupp4_12_TableSupp9','')
#%%
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from matplotlib import colors

sys.path.insert(1, path_notebook)
# import notebook_tools
# from notebook_tools import get_data_from
#%%
import numpy as np
import matplotlib.pyplot as plt
from projects.notebook_tools.notebook_tools import get_data_from, fit_data
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
from matplotlib import pyplot as plt
import qtt
import sys

# from notebook_tools import fit_data
from helper_functions import data1d_attr, data2d_attr, data_setvaribles
# from thresholding_2023_03_02 import Gauss2,_estimate_double_gaussian_parameters,thresholded_data,thresholded_2d_data
# from thresholding import Gauss2,_estimate_double_gaussian_parameters,thresholded_data,thresholded_2d_data
#%%
start_time = '2023-08-29\\20-54-19'



end_time = '2023-08-29\\21-45-08'

datadir = os.getcwd()  
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

datfile = datfiles[0]
#%%
from scipy.fft import fft, fftfreq

def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_with_decay_1(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2))+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0


xdata_coarse = np.hstack([datfiles[0].wait_q4_set.ndarray, datfiles[3].wait_q4_set.ndarray])
ydata_coarse = datfiles[0].vP1_Q1Q4_AC_set.ndarray
zdata_coarse = np.hstack([datfiles[0].su_S_North.ndarray, datfiles[3].su_S_North.ndarray])





xdata_fine = np.vstack([datfiles[1].wait_q4_set.ndarray, datfiles[2].wait_q4_set.ndarray])
ydata_fine = np.hstack([datfiles[1].vP1_Q1Q4_AC_set.ndarray, datfiles[2].vP1_Q1Q4_AC_set.ndarray])
zdata_fine = np.vstack([datfiles[1].su_S_North.ndarray, datfiles[2].su_S_North.ndarray])


#%%

fig = plt.figure(figsize=(13,6))
plt.subplots_adjust(bottom=0.15, top=0.8)
# plt.plot(xf, zf)
plt.subplot(121)

img = plt.pcolormesh(xdata_coarse, ydata_coarse, zdata_coarse  )
# plt.xlabel('wait_q2' + ' (' + 'ns' + ')')
# plt.ylabel('wait_q3' + ' (' + 'ns' + ')')

plt.title(start_time)
cb = fig.colorbar(img, orientation='vertical', aspect=10, shrink=0.8)
cb.set_label('spin up probability')


plt.subplot(122)

img = plt.pcolormesh(xdata_fine, ydata_fine, zdata_fine  )


#%%
f_ramsey = dict(q1_coarse=[], q2_coarse=[], q1_fine=[] )


T2_ramsey = dict(q1_coarse=[], q2_coarse=[], q1_fine=[] )



plt.close('all')



for i in range(len(ydata_coarse)):
    xdata = xdata_coarse[i,:]
    zdata = zdata_coarse[i,:]
        
    zf = fft(zdata[0:])
    xf = fftfreq(len(xdata[0:]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure(figsize=(10,5))
    plt.title(ydata_coarse[i])
    # plt.plot(xf, np.abs(zf))
    
    if i in [7,8]:
        zdata = zdata[xdata<=90]
        xdata = xdata[xdata<=90]
    
    
    if i == 13:
        f0 = 0.086
    elif i==17:
        f0 = 0.089
    # else:
    #     f0 = p1[1]
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)
    # bounds = ((-np.inf,-np.inf,-np.inf,-np.inf,0), (np.inf,np.inf,np.inf,np.inf,np.inf))
    # p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay_1,p0=p0, plot=True, return_cov=True, bounds=bounds)
        

    f_ramsey['q1_coarse'] += [ [ydata_coarse[i], p1[1],p1std[1] ]    ]
    T2_ramsey['q1_coarse'] += [ [ydata_coarse[i], p1[4],p1std[4] ]    ]
    


#%%

plt.close('all')



for i in range(len(ydata_fine)):
    xdata = xdata_fine[i,:]
    zdata = zdata_fine[i,:]
        
    zf = fft(zdata[0:])
    xf = fftfreq(len(xdata[0:]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure(figsize=(10,5))
    plt.title(ydata_fine[i])
    # plt.plot(xf, np.abs(zf))
    
    # if i in [7,8]:
    #     zdata = zdata[xdata<=90]
    #     xdata = xdata[xdata<=90]
    
    
    # if i == 13:
    #     f0 = 0.086
    # elif i==17:
    #     f0 = 0.089
    # else:
    #     f0 = p1[1]
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)
        
    # bounds = ((-np.inf,-np.inf,-np.inf,-np.inf,-np.inf), (np.inf,np.inf,np.inf,np.inf,np.inf))
    # p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay_1,p0=p0, plot=True, return_cov=True)
    
    f_ramsey['q1_fine'] += [ [ydata_fine[i], p1[1],p1std[1] ]    ]
    T2_ramsey['q1_fine'] += [ [ydata_fine[i], p1[4],p1std[4] ]    ]


    
#%%
vP4_coarse = np.array([f_ramsey['q1_coarse'][i][0] for i in range(len(f_ramsey['q1_coarse'])) ]  )
f_coarse = np.array([f_ramsey['q1_coarse'][i][1] for i in range(len(f_ramsey['q1_coarse'])) ]  )
f_coarse_err = np.array([f_ramsey['q1_coarse'][i][2] for i in range(len(f_ramsey['q1_coarse'])) ]  )
T2_coarse = np.array([T2_ramsey['q1_coarse'][i][1] for i in range(len(T2_ramsey['q1_coarse'])) ]  )

plt.close('all')
# plt.figure()
# plt.scatter( vP4_coarse, f_coarse    )

vP4_fine = np.array([f_ramsey['q1_fine'][i][0] for i in range(len(f_ramsey['q1_fine'])) ]  )
f_fine = np.array([f_ramsey['q1_fine'][i][1] for i in range(len(f_ramsey['q1_fine'])) ]  )
f_fine_err = np.array([f_ramsey['q1_fine'][i][2] for i in range(len(f_ramsey['q1_fine'])) ]  )
T2_fine = np.array([T2_ramsey['q1_fine'][i][1] for i in range(len(T2_ramsey['q1_fine'])) ]  )
T2_fine_err = np.array([T2_ramsey['q1_fine'][i][2] for i in range(len(T2_ramsey['q1_fine'])) ]  )
  
# plt.scatter( vP4_fine, f_fine    )
    

#%%
min_vP4_RamseyQ1Q4, max_vP4_RamseyQ1Q4 = 42, 62
min_vP1_RamseyQ1Q4, max_vP1_RamseyQ1Q4 = 35, 55

Lever_arm_1=0.094
Lever_arm_4=0.078
Lever_arm_detuning_Q1Q4=0.5*Lever_arm_4+0.5*np.abs((min_vP1_RamseyQ1Q4-max_vP1_RamseyQ1Q4)/(min_vP4_RamseyQ1Q4-max_vP4_RamseyQ1Q4))*Lever_arm_1



#%%

h_planck=6.63e-34
e=1.602e-19 



def fres_evolution(vP, tc, vP_0, fA_0, fA_1, fB_0, fB_1,theta_deg):
    Lever_arm_detuning = Lever_arm_detuning_Q1Q4
    det=Lever_arm_detuning*(vP-vP_0)*1e-3*e/(h_planck*1e9)
    theta=theta_deg*np.pi/180
    fA=fA_0+fA_1*(vP-vP_0)
    fB=fB_0+fB_1*(vP-vP_0)
    Delta=np.sqrt(det**2+tc**2)
    num=(2*det**2+tc**2)*(fA**2+fB**2)+2*det*(fB**2-fA**2)*Delta+2*fA*fB*tc**2*np.cos(theta)
    return (np.sqrt(num))/(2*Delta)



#%%
vP_tot=np.concatenate([vP4_coarse,vP4_fine])
f_tot=np.concatenate([f_coarse,f_fine])
f_tot=f_tot[np.argsort(vP_tot)]
vP_tot=np.sort(vP_tot)

f_tot_err=np.concatenate([f_coarse_err,f_fine_err])
f_tot_err=f_tot_err[np.argsort(vP_tot)]



def fres_evolution_v1(vP, tc, vP_0, theta_deg):
    maskA = vP_tot<=49.5
    fA_1, fA_0 = np.polyfit( vP_tot[maskA], f_tot[maskA],  deg=1)
    maskB = vP_tot>=57.5
    fB_1, fB_0 = np.polyfit( vP_tot[maskB], f_tot[maskB],  deg=1)
    Lever_arm_detuning = Lever_arm_detuning_Q1Q4
    det=Lever_arm_detuning*(vP-vP_0)*1e-3*e/(h_planck*1e9)
    theta=theta_deg*np.pi/180
    fA=fA_0+fA_1*(vP-0)
    fB=fB_0+fB_1*(vP-0)
    Delta=np.sqrt(det**2+tc**2)
    num=(2*det**2+tc**2)*(fA**2+fB**2)+2*det*(fB**2-fA**2)*Delta+2*fA*fB*tc**2*np.cos(theta)
    return (np.sqrt(num))/(2*Delta)

def H4x4_v1(vP, tc, vP_0, theta_deg):
    maskA = vP_tot<=49.5
    fA_1, fA_0 = np.polyfit( vP_tot[maskA], f_tot[maskA],  deg=1)
    maskB = vP_tot>=57.5
    fB_1, fB_0 = np.polyfit( vP_tot[maskB], f_tot[maskB],  deg=1)
    Lever_arm_detuning = Lever_arm_detuning_Q1Q4
    det=Lever_arm_detuning*(vP-vP_0)*1e-3*e/(h_planck*1e9)
    theta=theta_deg*np.pi/180
    fA=fA_0+fA_1*(vP-0)
    fB=fB_0+fB_1*(vP-0)
    s_z = np.array([[1,0],[0,-1]])
    s_x = np.array([[0,1],[1, 0]])
    s_y = np.array([[0,-1j],[1j, 0]])
    s_0 = np.eye(2)    
    Hsxs0 = np.kron(s_x, s_0)
    Hszs0 = np.kron(s_z, s_0)
    H_Ez = np.array([[ fA,   0,              0,                 0],
                        [  0,  -fA,             0,                 0],
                        [  0,   0, fB*np.cos(theta),  fB*np.sin(theta)],
                        [  0,   0, fB*np.sin(theta), -fB*np.cos(theta)]]) /2
    H4x4 = Hsxs0*tc + Hszs0*det + H_Ez
    ev, evec = np.linalg.eigh( H4x4  ) 
    return ev, evec


def fres_evolution_v2(vP, tc, vP_0):
    theta_deg = 41.5
    maskA = vP_tot<=49.5
    fA_1, fA_0 = np.polyfit( vP_tot[maskA], f_tot[maskA],  deg=1)
    maskB = vP_tot>=57
    fB_1, fB_0 = np.polyfit( vP_tot[maskB], f_tot[maskB],  deg=1)
    Lever_arm_detuning = Lever_arm_detuning_Q1Q4
    det=Lever_arm_detuning*(vP-vP_0)*1e-3*e/(h_planck*1e9)
    theta=theta_deg*np.pi/180
    fA=fA_0+fA_1*(vP-0)
    fB=fB_0+fB_1*(vP-0)
    Delta=np.sqrt(det**2+tc**2)
    num=(2*det**2+tc**2)*(fA**2+fB**2)+2*det*(fB**2-fA**2)*Delta+2*fA*fB*tc**2*np.cos(theta)
    return (np.sqrt(num))/(2*Delta)




#%%
from functools import partial



CQ1 = '#0D77BC'
CQ2 = '#DD5E27'

fig = plt.figure(figsize=(8,8))
plt.subplots_adjust(bottom=0.15, top=0.8)
plt.subplot(221)

img = plt.pcolormesh(xdata_coarse, ydata_coarse, zdata_coarse  )

cb = fig.colorbar(img, orientation='vertical', aspect=10, shrink=0.8)
cb.set_label('spin up probability')

plt.subplot(222)

img = plt.pcolormesh(xdata_fine, ydata_fine, zdata_fine  )


plt.subplot(223)

plt.scatter( vP4_coarse, f_coarse    )
plt.scatter( vP4_fine, f_fine    )
plt.xlim(50, 60)


vP_fit=np.linspace(vP_tot[0], vP_tot[-1],2001)

plt.subplot(224)
plt.plot(vP_tot,f_tot,color=CQ1,marker='v',markersize=4, linestyle='none',label='Data')

init_guess=[9, 54.5,   35]

pstd, popt = fit_data(vP_tot,f_tot, p0=init_guess,func=fres_evolution_v1,plot=False,return_cov=True)
plt.plot(vP_fit, fres_evolution_v1(vP_fit,*popt),label='Fit',color='k')

print('popt = ' + ''.join(['{:.3g}, ']*len(popt)).format(*popt)   )
# popt = 26.9, 54.8, 64.6, 
print('pstd = ' + ''.join(['{:.3g}, ']*len(pstd)).format(*pstd)   )
# pstd = 0.592, 0.0309, 1.79, 
fit_param = popt




#%%
valvP4_ = np.linspace(51,59,1001)
evs = np.zeros((len(valvP4_), 4))
transition_mat = np.zeros((len(valvP4_),4,4))
f_arr = np.zeros((len(valvP4_),4,4))
for i in range(len(valvP4_)):
    vP = valvP4_[i]
    ev, U = H4x4_v1(vP,*popt)
    dvP = 0.001
    dUdvP = (H4x4_v1(vP+dvP,*popt)[1] - H4x4_v1(vP-dvP,*popt)[1])/(2*dvP)    
    for j in range(4):
        for k in range(4):
            Udag_dUdvP_jk = (U.conj().T @ dUdvP)[j,k] 
            f_arr[i,j,k] = ev[k] - ev[j]
            transition_mat[i,j,k]  =  np.pi * f_arr[i,j,k]*1e9 *1.054e-34 *  Udag_dUdvP_jk 
transition_mat = np.array(transition_mat)




#%%
k_B = 1.38e-23
T_env = 700 # unit: Kalvin
S_vPvP = 2 * 50 * k_B *T_env  # unit: V^2/Hz
plt.figure()

j = 0
for k in [1,]:
    T1 = 1/(S_vPvP*1e6 * (transition_mat[:,j,k]/1.054e-34)**2)
    plt.plot(valvP4_, 0.18+0.35 - np.exp(-300e-6*(2/T1))*0.35 )
    plt.plot(valvP4_, 0.89-0.35 + np.exp(-300e-6*(2/T1))*0.35 )


datfiles = get_data_from('2023-08-29\\23-37-56', '2023-08-29\\23-39-29', rootfolder=datadir, only_complete = False) [0]
xdata, ydata = datfiles[0].vP1_Q1Q4_AC_set.ndarray,  datfiles[0].su_S_North.ndarray
plt.scatter(xdata-0.5, ydata )
xdata, ydata = datfiles[1].vP1_Q1Q4_AC_set.ndarray,  datfiles[1].su_S_North.ndarray 
plt.scatter(xdata-0.5, ydata )
plt.xlabel('vP4 (mV)')
plt.ylabel('Spin up probability after idling for 300 us at vP4')

decay_per_ramp = 2e-9*(abs(valvP4_[0]-valvP4_[-1])/25.4)*np.mean(1/T1)
print('\n' + f'At photon temperature {T_env}K the 50Ohm coax produces the noise',   S_vPvP/50 , ' Watt/Hz' \
      ,f'or {10*np.log10(1e3*S_vPvP/50) :.5g} dBm/Hz' )
print('thermalization contribute to avg infidelity ', decay_per_ramp*4* 2/3  )
print('minimum decay time = ', min(T1),' seconds')
total_infidelity = decay_per_ramp*4* 2/3


#%%
CMAP = 'viridis'
CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5),  Size.Fixed(0.7),Size.Fixed(1.5)]
v = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
ax31 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=3, ny=1))
ax51 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=5, ny=1))


ax = ax11
x = xdata_coarse
y = ydata_coarse
z = zdata_coarse
img = ax.pcolormesh( x, y, z , shading='auto', cmap=CMAP,rasterized=True)
cb=fig.colorbar(img, ax=ax, location='top', ticks=[0.2, 0.8],shrink=0.05,aspect=4,anchor=(0.3, -0.25))
ax.set_xlabel('t (ns)')
ax.set_ylabel(r'$\rm vP_4$ (mV)')
ax.set_xticks( np.linspace(np.min(x), np.max(x), 3) )
ax.set_yticks( np.linspace(np.min(y), np.max(y), 3)  )


ax = ax31
x = xdata_fine
y = ydata_fine
z = zdata_fine
img = ax.pcolormesh( x, y, z , shading='auto', cmap=CMAP,rasterized=True)
cb=fig.colorbar(img, ax=ax, location='top', ticks=[0.2, 0.8] ,shrink=0.05,aspect=4,anchor=(0.9, -0.25))
ax.set_xlabel('t (ns)')
ax.set_ylabel(r'$\rm vP_4$ (mV)')
ax.set_xticks( np.linspace(np.min(x), np.max(x), 3) )
ax.set_yticks( np.linspace(np.min(y), np.max(y), 3)  )


ax = ax51
vP_fit=np.linspace(vP_tot[0], vP_tot[-1],201)
ax.plot(vP_tot,f_tot*1e3,color=CQ1,marker='v',markersize=2, linestyle='none',label='Data')
ax.plot(vP_fit, fres_evolution_v1(vP_fit,*popt)*1e3,label='Fit',color='k', linewidth=0.5)
ylimit = np.array(ax.get_ylim())
ax.set_ylim(ylimit)
ax.set_xlabel(r'$\rm vP_4$ (mV)')
ax.set_ylabel('frequency (MHz)')
ax.set_xticks( np.linspace(np.min(vP_tot), np.max(vP_tot), 3) )

# filename = 'Supp_fit_tc_q1'
# plt.savefig(filename+'.pdf')

#%%
# data points that show obvious decay: vP3_coarse from 53 to 55.5
mT2fit = np.logical_and(vP4_fine>52.99, vP4_fine<55.51)

CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
from mpl_toolkits.axes_grid1 import Divider, Size
fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5),  Size.Fixed(0.7),Size.Fixed(1.5)]
v = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))

x, y, yerr = vP4_fine[mT2fit], np.abs(T2_fine)[mT2fit], np.abs(T2_fine_err)[mT2fit]
for i in range(len(x)):
    plt.errorbar(x[i], y[i],
                 yerr = yerr[i], c=CQ1,
                 fmt ='o', markersize=3, linestyle='', linewidth=0.5)
    
valvP4_ = np.linspace(42,62,10001)
E = fres_evolution_v1(valvP4_, *fit_param)
dEdvP4 = np.diff(E)/np.diff(valvP4_)
dEdvP4_1 = np.diff(E)/np.diff(valvP4_)[0]
dEdvP4_1 = (dEdvP4_1[:-1] + dEdvP4_1[1:])/2
dEdvP4_2 = np.diff(E,2)/(np.diff(valvP4_)[0]**2)
valvP4 = valvP4_[1:-1] 

from scipy.interpolate import CubicSpline
dEdvP4_1_interp = CubicSpline(valvP4[np.logical_and(valvP4>52.99, valvP4<55.51)], dEdvP4_1[np.logical_and(valvP4>52.99, valvP4<55.51)]  )
dEdvP4_2_interp = CubicSpline(valvP4[np.logical_and(valvP4>52.99, valvP4<55.51)], dEdvP4_2[np.logical_and(valvP4>52.99, valvP4<55.51)]  )


from scipy.optimize import minimize
def tominimize(delta_vP4):
    x, y = vP4_fine[mT2fit], np.abs(T2_fine)[mT2fit]
    dEdvP4_1 = dEdvP4_1_interp(x)
    dEdvP4_2 = dEdvP4_2_interp(x)
    return np.sum( (1/abs( np.sqrt(2)/(2*np.pi*(dEdvP4_1*delta_vP4+0.5*dEdvP4_2*delta_vP4*delta_vP4))   ) - 1/y)**2  )

res_minimized = minimize(tominimize, x0=np.array([0.22])   ) 

delta_vP4 = res_minimized['x'][0]
T2 = abs( np.sqrt(2)/(2*np.pi*(dEdvP4_1*delta_vP4+0.5*dEdvP4_2*delta_vP4*delta_vP4))   )
ax11.plot(valvP4, T2, label='vP4 dephasing',color='k')
print(f'delta_vP4 = {delta_vP4}')
# delta_vP4 = 0.19381304639146746


plt.xticks([53, 54, 55, 56])
plt.xlim(52.5, 56.5)
plt.ylim(0, 200)
plt.yticks([0, 100, 200])

# estimate the required barrier gate noise delta_vB14 to produce the same order of dephasing time 
kappa = 0.059
tempvar_0 = fres_evolution_v1(valvP4, 26.9*np.exp(-0.5*kappa*0.01), *fit_param[1:])
tempvar_1 = fres_evolution_v1(valvP4, 26.9*np.exp(-0.5*kappa*0.0), *fit_param[1:])
tempvar_2 = fres_evolution_v1(valvP4, 26.9*np.exp( 0.5*kappa*0.01), *fit_param[1:])
dEdvB14_1 = (tempvar_0 - tempvar_2 ) /(2*0.01)
dEdvB14_2 = (tempvar_0 -2*tempvar_1 + tempvar_2 ) /(0.01**2)
delta_vB14 = 8  
T2 = abs( np.sqrt(2)/(2*np.pi*(dEdvB14_1*delta_vB14+0.5*dEdvB14_2*delta_vB14*delta_vB14))   )
ax31 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=3, ny=1))
for i in range(len(x)):
    plt.errorbar(x[i], y[i],
                 yerr = yerr[i], c=CQ1,
                 fmt ='o', markersize=3, linestyle='', linewidth=0.5)
ax31.plot(valvP4, T2,label='vpB14 dephasing',color='b')

plt.xticks([53, 54, 55, 56])
plt.xlim(52.5, 56.5)
plt.ylim(0, 200)
plt.yticks([0, 100, 200])


# filename = 'Supp_T2_ac14'
# plt.savefig(filename+'.pdf')





#%%  Nielsen2002. equation 17  ( 	Phys. Lett. A 303.  arXiv:quant-ph/0205035  )
import scipy
import random
s0 = np.eye(2)
sx = np.array([[0, 1],
               [1, 0]])
sy = np.array([[0, -1j],
               [1j, 0]])
sz = np.array([[1, 0],
               [0,-1]])
P1q = [s0, sx, sy, sz]

dtq2 = 1.9489
dtq3 = 1.9411
DEG = np.pi/180
fq1 = 3.38324710e-02*5/4 
fq4     = 7.15275643e-02*5/4   # qubit frequency
theta = 41.517 * DEG       # q3 angle from sigma_z axis (toward sigma_x axis)

wait_q4_1= 3.747 
wait_q1 =  19.33 
wait_q4_2 =  10.17 
add_wait_q1 =  9.4 


def Q1Xgate(fq1, fq4, dt_es, t_hres=[0,0,0,0], theta=41.517 * DEG):
    H1_0 = (1/2) *  sz
    H4_0 = (1/2) * (np.cos(theta)*sz + np.sin(theta)*sx   )
    Tq1 = 23.645922876871747
    output = scipy.linalg.expm(  -2*np.pi*1j* H1_0 *fq1* ( Tq1 + dt_es + t_hres[0])  )
    output = scipy.linalg.expm(  -2*np.pi*1j* H4_0 *fq4* (wait_q4_1+dtq3  - 2*dt_es - t_hres[0] +t_hres[1])     )  @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H1_0 *fq1* (wait_q1+dtq2 +2*dt_es -t_hres[1] + t_hres[2])  ) @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H4_0 *fq4* (wait_q4_2+dtq3 -2*dt_es -t_hres[2] + t_hres[3])  ) @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H1_0 *fq1* (add_wait_q1+dtq2 + Tq1 +dt_es -t_hres[3] )  ) @ output
    U = scipy.linalg.expm(  -1j* 0.5*sz* (-109.87825962836682)*np.pi/180    )
    return U.T.conj() @  output @ U

temp_var = 5

mat = Q1Xgate(fq1, fq4, 0)
u0 = np.real( np.trace( mat @ s0 )/2  )
ux = np.imag( np.trace( mat @ sx )/2  )
uy = np.imag( np.trace( mat @ sy )/2  )
uz = np.imag( np.trace( mat @ sz )/2  )

U = s0* u0 +  (ux*sx + uy*sy + uz*sz) 
uvec = np.array([ux, uy, uz])
if np.sqrt(np.sum(uvec**2)) <= 1e-14:
    uvec = np.array([0,0,1])
else:
    uvec = uvec/np.sqrt(np.sum(uvec**2))
rot_angle = 2*np.arccos( u0)
gate_polar_angle = np.arctan2( np.sqrt(uvec[0]**2 + uvec[1]**2), uvec[2]   ) 
gate_azu_angle = np.arctan2( uvec[1], uvec[0] )
print('\n'+'X90 gate rotation: ', rot_angle/DEG,     gate_polar_angle /DEG , gate_azu_angle/DEG  )


#%%  infidelity in single qubit space

noise_correlated = False

Nrand = 10_000
T2_dot1 = 7000 # unit ns
T2_dot4 = 4500 # unit ns
delta_vP4 = res_minimized['x'][0]  
# if detuning noise 1meV = 241.69 GHz
# es[0] = -337.43 GHz
# es[600] = 226.43 GHz
# -> dt_es = tramp* 241.69/(337.43+226.43) = 0.857269 ns 
print('\n'+'based on Pedersen2007-equation(5)')

for i_noise_source, name in enumerate(['Dot1 Larmor frequency fluctuation', 
                                 'Dot4 Larmor frequency fluctuation', 
                                 'Detuning noise']):
    idxvec = np.zeros((3,))
    idxvec[i_noise_source] = 1
    sigma_noise = np.array([  -idxvec[0]*np.sqrt(2)/(T2_dot1*2*np.pi) ,
                              -idxvec[1]*np.sqrt(2)/(T2_dot4*2*np.pi) ,
                               idxvec[2]*delta_vP4 *Lever_arm_detuning_Q1Q4* 0.857269  ]) 

    
    
    noise = np.empty((Nrand,3))
    for i_rand in range(Nrand):
        if noise_correlated:
            tempvar = random.gauss(0, 1)  
            for i in range(3):
                noise[i_rand, i] = tempvar * sigma_noise[i]
        else:
            for i in range(3):
                noise[i_rand, i] = random.gauss(0, 1) * sigma_noise[i]
    
    
    # Ui = scipy.linalg.expm(  -2*np.pi*1j* sx/8    ) 
    Ui = Q1Xgate(fq1, fq4, 0 )
    UiH = Ui.conj().T
    
    
    
    tr2sum = 0
    for i_rand in range(Nrand):
        dfq1 = noise[i_rand,0]
        dfq4 = noise[i_rand,1]
        dt_es = noise[i_rand,2]
        dtheta = 0
        R = Q1Xgate(fq1+dfq1, fq4+dfq4, dt_es, theta=theta + dtheta )
        RH = R.conj().T
        tr2sum += np.abs(np.trace( RH @ Ui  ))**2/Nrand
    d = 2
    F =  (tr2sum + d )/(d  *(d+1))
    print(name + ' incoherent error = ', 1-F)
    total_infidelity += 1-F
    # not consider 1/f integral:
    # 2.553131377835971e-05
    # 5.016118646672396e-06
    # 7.24794470664536e-05    
    


#%%

import xarray as xr
def get_awg_output(t, samples, analogue_shift, digital_filter_mode = 1):

    fname = os.path.dirname(__file__) + f'/keysight_pulse_response_{digital_filter_mode}.hdf5'
    pulse_response = xr.open_dataset(fname, engine='h5netcdf')
    t_response = pulse_response.coords['t'].data
    response = pulse_response['y'].data / 0.77
    sr = round(1/(t_response[1]-t_response[0]))

    t = np.linspace(t[0], t[-1]+1, len(t)*sr, endpoint=False)
    d = np.zeros(len(samples)*sr)
    d[::sr] = samples

    d = np.convolve(d, response)
    n_before = round(-t_response[0]*sr)
    n_after = round(t_response[-1]*sr)
    return t+analogue_shift, d[n_before: -n_after]


def render_ramp(t_start, t_stop, v_start, v_stop, samples):

    dv_dt = (v_stop - v_start) / (t_stop - t_start)
    i_start, dt_start = divmod(t_start, 1.0)
    i_stop, dt_stop = divmod(t_stop, 1.0)
    i_start = int(i_start)
    i_stop = int(i_stop)

    samples[i_start] += v_start*(1-dt_start)
    samples[i_stop] += v_stop*dt_stop - dv_dt*dt_stop
    samples[i_start+1:i_stop] += np.linspace(v_start + dv_dt*(1-dt_start),
                                             v_stop - dv_dt*dt_stop,
                                             i_stop-i_start-1, endpoint=False)

def quantize_amplitude(wave):
    min_vstep = 0.00038 # 0.38 mV
    wave = np.round(wave/min_vstep) * min_vstep
    return wave


def plot_samples(t, samples, color = 'k', linestyle = 'solid'):
    t2 = np.repeat(t, 2)[1:-1]
    samples2 = np.repeat(samples[:-1], 2)
    return plt.plot(t2, samples2, color)


t_max = 110
n_pulses = 1
t_ramp = 2.0
amplitude = 0.2
analogue_shift = 0
digital_filter_mode = 1
t = np.arange(t_max)


sarr = np.linspace(0, 1, 101)[:-1]
deviation = np.empty((len(sarr),4))
for n, s in  enumerate(sarr):
    wave = np.zeros(t_max)
    
    points = [(0.0, 0.0)]
    Tq1 = 23.645922876871747
    t0 = Tq1 + s
    # pt.title(f"first ramp starts at {t0:.2f} ns")

    for i in range(n_pulses):
        
        # add 1 ramped pulse
        points += [(t0, 0.0)]
        t0 += t_ramp
        points += [(t0, amplitude)]
        t0 += wait_q4_1
        points += [(t0, amplitude)]
        t0 += t_ramp
        points += [(t0, 0.0)]
        t0 += wait_q1
        points += [(t0, 0.0)]
        t0 += t_ramp
        points += [(t0, amplitude)]
        t0 += wait_q4_2  
        points += [(t0, amplitude)]
        t0 += t_ramp
        points += [(t0, 0.0)]
        t0 += add_wait_q1 + Tq1       
    points += [(t_max-1, 0.0)]

    xy = np.array(points).T

    for i in range(len(points)-1):
        t_start, v_start = points[i]
        t_stop, v_stop = points[i+1]
        if t_stop > t_start:
            render_ramp(t_start, t_stop, v_start, v_stop, wave)

    wave = quantize_amplitude(wave)
    ta, out = get_awg_output(t, wave, analogue_shift, digital_filter_mode)
    
    for i in range(2):
        mx = np.logical_and(ta>xy[0][4*i+1], ta<xy[0][4*i+2]   )
        ta_mid = np.interp(amplitude/2, out[mx], ta[mx])
        deviation[n,2*i] =  ta_mid - (xy[0][4*i+1]+xy[0][4*i+2])/2  
    
        mx = np.logical_and(ta>xy[0][4*i+3], ta<xy[0][4*i+4]   )
        ta_mid = np.interp(amplitude/2, out[mx][::-1], ta[mx][::-1])
        deviation[n,2*i+1] =  ta_mid - (xy[0][4*i+3]+xy[0][4*i+4])/2  

collect = []
for idx in range(deviation.shape[0]):

    Ui = Q1Xgate(fq1, fq4, 0, t_hres= deviation[idx,:])
    Nrand = deviation.shape[0]
    tr2sum = 0
    tr2arr = []
    for i_rand in range(Nrand):
    
        R = Q1Xgate(fq1, fq4, 0, t_hres= deviation[i_rand,:])
        RH = R.conj().T
        tr2sum += np.abs(np.trace( RH @ Ui  ))**2/Nrand
        tr2arr += [np.abs(np.trace( RH @ Ui  ))**2]
    d = 2
    F =  (tr2sum + d )/(d  *(d+1))
    collect += [1-F]
    # print('based on Pedersen2007-equation(5), total infidelity = ', 1-F)
    tr2arr = np.array(tr2arr)

print('\n'+'waveform uncertainty contribute to incoherent error (avg. gate infidelity): ')
print('max ', np.max(collect))
print('min ', np.min(collect))
print('mean ', np.mean(collect))
# max  0.00014627532007083222
# min  3.990454197244642e-05
# mean  7.859483809478651e-05
plt.hist(collect)

print('\n'+'qubit 1 X90 gate summary')
print(f'total_infidelity minimum = {total_infidelity + np.min(collect) :.5g}',  )
print(f'total_infidelity maximum= {total_infidelity + np.max(collect) :.5g}',  )

