# this file is modified from SUpp/fit_tc_angle/AC23_v0.py


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
start_time = '2023-08-29\\22-26-44'



end_time = '2023-08-29\\22-50-21'

datadir = os.getcwd()  
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

datfile = datfiles[0]
#%%
from scipy.fft import fft, fftfreq

def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0

def osc_with_twodecay(t, A, f, y0, phi, T2, T2w):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2 - (t/T2w))+y0

xdata_coarse = np.hstack([datfiles[0].waittime_set.ndarray, datfiles[1].waittime_set.ndarray])
ydata_coarse = datfiles[0].vP3_Q2Q3_AC_set.ndarray
zdata_coarse = np.hstack([datfiles[0].su_S_North.ndarray, datfiles[1].su_S_North.ndarray])




#%%

fig = plt.figure(figsize=(13,6))
plt.subplots_adjust(bottom=0.15, top=0.8)
# plt.plot(xf, zf)
plt.subplot(121)

img = plt.pcolormesh(xdata_coarse, ydata_coarse, zdata_coarse  )


plt.title(start_time)
cb = fig.colorbar(img, orientation='vertical', aspect=10, shrink=0.8)
cb.set_label('spin up probability')

#%%
f_ramsey = dict(q1_coarse=[], q2_coarse=[], q1_fine=[] )

T2_ramsey = dict(q1_coarse=[], q2_coarse=[], q1_fine=[] )
T2w_ramsey = dict(q1_coarse=[], q2_coarse=[], q1_fine=[] )



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

    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)

    f_ramsey['q2_coarse'] += [ [ydata_coarse[i], p1[1],p1std[1] ]    ]
    T2_ramsey['q2_coarse'] += [ [ydata_coarse[i], p1[4],p1std[4] ]    ]
    
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200, 800]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_twodecay,p0=p0, plot=True, return_cov=True)
    T2w_ramsey['q2_coarse'] += [ [ydata_coarse[i], p1[4],p1std[4] , p1[5],p1std[5] ]    ]    



    
#%%
vP3_coarse = np.array([f_ramsey['q2_coarse'][i][0] for i in range(len(f_ramsey['q2_coarse'])) ]  )
f_coarse = np.array([f_ramsey['q2_coarse'][i][1] for i in range(len(f_ramsey['q2_coarse'])) ]  )
f_coarse_err = np.array([f_ramsey['q2_coarse'][i][2] for i in range(len(f_ramsey['q2_coarse'])) ]  )
T2_coarse = np.array([T2_ramsey['q2_coarse'][i][1] for i in range(len(T2_ramsey['q2_coarse'])) ]  )
T2_coarse_err = np.array([T2_ramsey['q2_coarse'][i][2] for i in range(len(T2_ramsey['q2_coarse'])) ]  )

plt.close('all')
# plt.figure()
# plt.scatter( vP3_coarse, f_coarse    )



#%%

min_vP3_RamseyQ2Q3, max_vP3_RamseyQ2Q3 = 32, 42
min_vP2_RamseyQ2Q3, max_vP2_RamseyQ2Q3 = 49, 59

Lever_arm_2=0.084
Lever_arm_3=0.08
Lever_arm_detuning_Q2Q3=0.5*Lever_arm_3+0.5*np.abs((min_vP2_RamseyQ2Q3-max_vP2_RamseyQ2Q3)/(min_vP3_RamseyQ2Q3-max_vP3_RamseyQ2Q3))*Lever_arm_2


#%%

h_planck=6.63e-34
e=1.602e-19 



def fres_evolution(vP, tc, vP_0, fA_0,  fB_0, fA_1, fB_1, theta_deg,  ):
    Lever_arm_detuning = Lever_arm_detuning_Q2Q3
    det=Lever_arm_detuning*(vP-vP_0)*1e-3*e/(h_planck*1e9)
    theta=theta_deg*np.pi/180
    fA=fA_0+fA_1*(vP-vP_0)
    fB=fB_0+fB_1*(vP-vP_0)
    Delta=np.sqrt(det**2+tc**2)
    num=(2*det**2+tc**2)*(fA**2+fB**2)+2*det*(fB**2-fA**2)*Delta+2*fA*fB*tc**2*np.cos(theta)
    return (np.sqrt(num))/(2*Delta)



#%%
vP_tot=np.concatenate([vP3_coarse])
f_tot=np.concatenate([f_coarse])
f_tot=f_tot[np.argsort(vP_tot)]
vP_tot=np.sort(vP_tot)

f_tot_err=np.concatenate([f_coarse_err])
f_tot_err=f_tot_err[np.argsort(vP_tot)]


def fres_evolution_v1(vP, tc, vP_0, theta_deg):
    
    maskA = vP_tot<=33
    fA_1, fA_0 = np.polyfit( vP_tot[maskA], f_tot[maskA],  deg=1)
    maskB = vP_tot>=40
    fB_1, fB_0 = np.polyfit( vP_tot[maskB], f_tot[maskB],  deg=1)
    # fA_1, fA_0 = 0, f_tot[np.argmin(vP_tot)]
    # fB_1, fB_0 = 0, f_tot[np.argmax(vP_tot)]
    
    # print(fA_1, fA_0)
    Lever_arm_detuning = Lever_arm_detuning_Q2Q3
    det=Lever_arm_detuning*(vP-vP_0)*1e-3*e/(h_planck*1e9)
    theta=theta_deg*np.pi/180
    fA=fA_0+fA_1*(vP-0)
    fB=fB_0+fB_1*(vP-0)
    Delta=np.sqrt(det**2+tc**2)
    num=(2*det**2+tc**2)*(fA**2+fB**2)+2*det*(fB**2-fA**2)*Delta+2*fA*fB*tc**2*np.cos(theta)
    return (np.sqrt(num))/(2*Delta)

def H4x4_v1(vP, tc, vP_0, theta_deg):
    maskA = vP_tot<=33
    fA_1, fA_0 = np.polyfit( vP_tot[maskA], f_tot[maskA],  deg=1)
    maskB = vP_tot>=40
    fB_1, fB_0 = np.polyfit( vP_tot[maskB], f_tot[maskB],  deg=1)
 
    Lever_arm_detuning = Lever_arm_detuning_Q2Q3
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
    theta_deg = 44.7
    maskA = vP_tot<=33
    fA_1, fA_0 = np.polyfit( vP_tot[maskA], f_tot[maskA],  deg=1)
    maskB = vP_tot>=40
    fB_1, fB_0 = np.polyfit( vP_tot[maskB], f_tot[maskB],  deg=1)

    
    # print(fA_1, fA_0)
    Lever_arm_detuning = Lever_arm_detuning_Q2Q3
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

fig = plt.figure(figsize=(12,12))
plt.subplots_adjust(bottom=0.15, top=0.8)
# plt.plot(xf, zf)
plt.subplot(221)

img = plt.pcolormesh(xdata_coarse, ydata_coarse, zdata_coarse  )

cb = fig.colorbar(img, orientation='vertical', aspect=10, shrink=0.8)
cb.set_label('spin up probability')

plt.subplot(222)



plt.subplot(223)

plt.scatter( vP3_coarse, f_coarse    )



vP_fit=np.linspace(vP_tot[0], vP_tot[-1],2001)

plt.subplot(224)
plt.plot(vP_tot,f_tot,color=CQ2,marker='v',markersize=4, linestyle='none',label='Data')

init_guess=[40, 37.4,   50]
pstd, popt = fit_data(vP_tot,f_tot, p0=init_guess,func=fres_evolution_v1,plot=False,return_cov=True)
plt.plot(vP_fit, fres_evolution_v1(vP_fit,*popt),label='Fit',color='k')
# plt.plot(vP_fit, fres_evolution_v2(vP_fit,*popt),label='Fit',color='k')
# plt.plot(vP_fit, fres_evolution_v1(vP_fit,popt[0], 36.5, 44.7),label='Fit',color='k')
print('popt = ' + ''.join(['{:.3g}, ']*len(popt)).format(*popt)   )
# popt = 19.9, 36.6, 51, 
print('pstd = ' + ''.join(['{:.3g}, ']*len(pstd)).format(*pstd)   )
# pstd = 1.18, 0.0545, 1.55
fit_param = popt.copy()


#%%
valvP4_ = np.linspace(35,40,1001)

transition_mat = np.zeros((len(valvP4_),4,4))
f_arr = np.zeros((len(valvP4_),4,4))
for i in range(len(valvP4_)):
    vP = valvP4_[i]

    ev, U = H4x4_v1(vP,*popt)
    dvP = 0.001
    dUdvP = (H4x4_v1(vP+dvP,*popt)[1] - H4x4_v1(vP-dvP,*popt)[1])/(2*dvP)    
    for j in range(4):
        for k in range(4):
            state0, state1 = U[:,j], U[:,k]
        
            Udag_dUdvP_jk = (U.conj().T @ dUdvP)[j,k] 
            f_arr[i,j,k] = ev[k] - ev[j]
            transition_mat[i,j,k]  =  np.pi * f_arr[i,j,k]*1e9 *1.054e-34 *  Udag_dUdvP_jk 
transition_mat = np.array(transition_mat)
#%%
k_B = 1.38e-23
T_env = 700 # unit Kalvin
S_vPvP = 2 * 50 * k_B *T_env  # unit: V^2/Hz

plt.figure()

j = 0
for k in [1,]:
    T1 = 1/(S_vPvP*1e6 * (transition_mat[:,j,k]/1.054e-34)**2)

    plt.plot(valvP4_-0.5, 0.18+0.35 - np.exp(-300e-6*(2/T1))*0.35 )
    plt.plot(valvP4_-0.5, 0.89-0.35 + np.exp(-300e-6*(2/T1))*0.35 )


datfiles = get_data_from('2023-08-29\\23-21-57', '2023-08-29\\23-26-18', rootfolder=datadir, only_complete = False) [0]
xdata, ydata = datfiles[0].vP3_Q2Q3_AC_set.ndarray,  datfiles[0].su_S_North.ndarray
plt.scatter(xdata-0.5, ydata )
xdata, ydata = datfiles[1].vP3_Q2Q3_AC_set.ndarray,  datfiles[1].su_S_North.ndarray 
plt.scatter(xdata-0.5, ydata )
plt.xlabel('vP3 (mV)')
plt.ylabel('Spin up probability after idling for 300 us at vP3')

decay_per_ramp = 2e-9*(abs(valvP4_[0]-valvP4_[-1])/19)*np.mean(1/T1)

print('\n' + f'At photon temperature {T_env}K the 50Ohm coax produces the noise',   S_vPvP/50 , ' Watt/Hz' \
      ,f'or {10*np.log10(1e3*S_vPvP/50) :.5g} dBm/Hz' )
print('thermalization contribute to avg infidelity ', decay_per_ramp*2* 2/3  )
print('minimum decay time = ', min(T1) ,' seconds')
total_infidelity = decay_per_ramp*2* 2/3

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
ax.set_ylabel(r'$\rm vP_3$ (mV)')
ax.set_xticks( np.linspace(np.min(x), np.max(x), 3) )
ax.set_yticks( np.linspace(np.min(y), np.max(y), 3)  )


ax = ax31
vP_fit=np.linspace(vP_tot[0], vP_tot[-1],201)
ax.plot(vP_tot,f_tot*1e3,color=CQ2,marker='v',markersize=2, linestyle='none',label='Data')
ax.plot(vP_fit, fres_evolution_v1(vP_fit,*popt)*1e3,label='Fit',color='k', linewidth=0.5)
ylimit = np.array(ax.get_ylim())
ax.set_ylim(ylimit)
ax.set_xlabel(r'$\rm vP_3$ (mV)')
ax.set_ylabel('frequency (MHz)')

ax.set_xticks( np.linspace(np.min(vP_tot), np.max(vP_tot), 3) )

# filename = 'Supp_fit_tc_q2'
# plt.savefig(filename+'.pdf')


#%%
# data points that show obvious decay: vP3_coarse = [38.5, 38. , 37.5, 37. , 36.5]
mT2fit = np.logical_and(vP3_coarse>36.4, vP3_coarse<38.6)


CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
from mpl_toolkits.axes_grid1 import Divider, Size
fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5),  Size.Fixed(0.7),Size.Fixed(1.5)]
v = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
x, y, yerr = vP3_coarse[mT2fit], np.abs(T2_coarse)[mT2fit], np.abs(T2_coarse_err)[mT2fit]
for i in range(len(x)):
    plt.errorbar(x[i], y[i],
                 yerr = yerr[i], c=CQ2,
                 fmt ='o', markersize=3, linestyle='', linewidth=0.5)

valvP3_ = np.linspace(27, 46,60001)
E = fres_evolution_v1(valvP3_, *fit_param)
dEdvP3 = np.diff(E)/np.diff(valvP3_)
dEdvP3_1 = np.diff(E)/np.diff(valvP3_)[0]
dEdvP3_1 = (dEdvP3_1[:-1] + dEdvP3_1[1:])/2
dEdvP3_2 = np.diff(E,2)/(np.diff(valvP3_)[0]**2)
valvP3 = valvP3_[1:-1] 

from scipy.interpolate import CubicSpline
dEdvP3_1_interp = CubicSpline(valvP3[np.logical_and(valvP3>36.49, valvP3<38.51)], dEdvP3_1[np.logical_and(valvP3>36.49, valvP3<38.51)]  )
dEdvP3_2_interp = CubicSpline(valvP3[np.logical_and(valvP3>36.49, valvP3<38.51)], dEdvP3_2[np.logical_and(valvP3>36.49, valvP3<38.51)]  )
# delta_vP4 = 0.22
# T2 = abs( np.sqrt(2)/(2*np.pi*dEdvP4_1*delta_vP4))
# plt.plot(valvP4, T2,label='vP4 dephasing',color='k')

from scipy.optimize import minimize
def tominimize(delta_vP3):
    x, y = vP3_coarse[mT2fit], np.abs(T2_coarse)[mT2fit]
    dEdvP3_1 = dEdvP3_1_interp(x)
    dEdvP3_2 = dEdvP3_2_interp(x)
    # print(dEdvP4_1)
    # print(dEdvP4_2)
    # print(abs( np.sqrt(2)/(2*np.pi*(dEdvP4_1*delta_vP4+0.5*dEdvP4_2*delta_vP4*delta_vP4))   ) )
    return np.sum( (1/abs( np.sqrt(2)/(2*np.pi*(dEdvP3_1*delta_vP3+0.5*dEdvP3_2*delta_vP3*delta_vP3))   ) - 1/y)**2  )
res_minimized = minimize(tominimize, x0=np.array([0.15])   ) #, method='Nelder-Mead', tol=1e-6)

delta_vP3 = res_minimized['x'][0]
T2 = abs( np.sqrt(2)/(2*np.pi*(dEdvP3_1*delta_vP3+0.5*dEdvP3_2*delta_vP3*delta_vP3))   )
plt.plot(valvP3, T2, label='vP3 dephasing',color='k')
print(f'delta_vP3 = {delta_vP3}')

plt.xlim(36, 39)
plt.ylim(0, 450)
plt.yticks([0, 200, 400])


# estimate the required barrier gate noise delta_vB14 to produce the same order of dephasing time 
kappa = 0.059
tempvar_0 = fres_evolution_v1(valvP3, 26.9*np.exp(-0.5*kappa*0.01), *fit_param[1:])
tempvar_1 = fres_evolution_v1(valvP3, 26.9*np.exp(-0.5*kappa*0.0), *fit_param[1:])
tempvar_2 = fres_evolution_v1(valvP3, 26.9*np.exp( 0.5*kappa*0.01), *fit_param[1:])
dEdvB14_1 = (tempvar_0 - tempvar_2 ) /(2*0.01)
dEdvB14_2 = (tempvar_0 -2*tempvar_1 + tempvar_2 ) /(0.01**2)
delta_vB14 = 6
T2 = abs( np.sqrt(2)/(2*np.pi*(dEdvB14_1*delta_vB14+0.5*dEdvB14_2*delta_vB14*delta_vB14))   )
ax31 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=3, ny=1))
for i in range(len(x)):
    plt.errorbar(x[i], y[i],
                 yerr = yerr[i], c=CQ1,
                 fmt ='o', markersize=3, linestyle='', linewidth=0.5)
ax31.plot(valvP3, T2,label='vpB14 dephasing',color='b')



plt.xlim(36, 39)
plt.ylim(0, 450)
plt.yticks([0, 200, 400])

# filename = 'Supp_T2_ac23'
# plt.savefig(filename+'.pdf')





#%%  Nielsen2002. equation 17  ( 	Phys. Lett. A 303.  arXiv:quant-ph/0205035  )
import scipy, scipy.linalg
import random
s0 = np.eye(2)
sx = np.array([[0, 1],
               [1, 0]])
sy = np.array([[0, -1j],
               [1j, 0]])
sz = np.array([[1, 0],
               [0,-1]])
P1q = [s0, sx, sy, sz]

dtq2 = 2.28847028
dtq3 = 1.54215906
DEG = np.pi/180
fq2 = 0.070926 *5/4 
fq3     = 0.06203767*5/4   # qubit frequency
theta = 44.718 * DEG       # q3 angle from sigma_z axis (toward sigma_x axis)

wait_q3 = 4.86  
wait_q2 = 3.42  



def Q2Xgate(fq2, fq3, dt_es, t_hres=[0,0], theta=44.718 * DEG):
    H2_0 = (1/2) *  sz
    H3_0 = (1/2) * (np.cos(theta)*sz + np.sin(theta)*sx   )
    Tq2 = 11.279361588134112
    output = scipy.linalg.expm(  -2*np.pi*1j* H2_0 *fq2* ( Tq2 + dt_es + t_hres[0])  )
    output = scipy.linalg.expm(  -2*np.pi*1j* H3_0 *fq3* (wait_q3+dtq3  - 2*dt_es - t_hres[0] + t_hres[1])     )  @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H2_0 *fq2* (wait_q2+dtq2 + Tq2 +dt_es  - t_hres[1])  ) @ output
    U = scipy.linalg.expm(  -1j* 0.5*sz* ( -88.902233307162 )*np.pi/180    )
    return U.T.conj() @  output @ U
temp_var = 5

mat = Q2Xgate(fq2, fq3, 0)
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
print('\n'+'X90 gate rotation: ',  rot_angle/DEG,     gate_polar_angle /DEG , gate_azu_angle/DEG  )


#%%  infidelity in single qubit space

noise_correlated = False

Nrand = 10_000
T2_dot2 = 4500  # unit ns
T2_dot3 = 4500  # unit ns
delta_vP3 = res_minimized['x'][0]  # 0.14384030863490807
# if detuning noise 1meV = 241.69 GHz
# es[0] = -190.49 GHz
# es[600] = 186.52 GHz
# des = 0.24169 GHz
# -> dt = tramp* 241.69/(190.49+186.52) = 1.28214 ns
print('\n'+'based on Pedersen2007-equation(5)')

for i_noise_source, name in enumerate(['Dot2 Larmor frequency fluctuation', 
                                 'Dot3 Larmor frequency fluctuation', 
                                 'Detuning noise']):
    idxvec = np.zeros((3,))
    idxvec[i_noise_source] = 1
    
    sigma_noise = np.array([  -idxvec[0]*np.sqrt(2)/(T2_dot2*2*np.pi)   ,
                              -idxvec[1]*np.sqrt(2)/(T2_dot3*2*np.pi)   ,
                               idxvec[2]*delta_vP3 *Lever_arm_detuning_Q2Q3*1.28214  ] ) 
    

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
    Ui = Q2Xgate(fq2, fq3, 0 )
    UiH = Ui.conj().T
    
    
    tr2sum = 0
    for i_rand in range(Nrand):
        dfq2 = noise[i_rand,0]
        dfq3 = noise[i_rand,1]
        dt_es = noise[i_rand,2]
        dtheta = 0
        R = Q2Xgate(fq2+dfq2, fq3+dfq3, dt_es, theta=theta + dtheta  )
        RH = R.conj().T
        tr2sum += np.abs(np.trace( RH @ Ui  ))**2/Nrand
    d = 2
    F =  (tr2sum + d )/(d  *(d+1))
    print(name + ' incoherent error = ', 1-F)
    total_infidelity += 1-F
    # not consider 1/f integral:
    # 6.900268236798013e-06
    # 6.74077409890117e-07
    # 1.2729936340605263e-06

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


t_max = 50
n_pulses = 1
t_ramp = 2.0
amplitude = 0.2
analogue_shift =  0
digital_filter_mode = 1
t = np.arange(t_max)

# sarr = np.linspace(0, 1, 1001)[:-1]
sarr = np.linspace(0, 1, 101)[:-1]
deviation = np.empty((len(sarr),2))
for n, s in  enumerate(sarr):
    wave = np.zeros(t_max)
    
    points = [(0.0, 0.0)]
    Tq2 = 11.279361588134112
    t0 = Tq2 + s
    # pt.title(f"first ramp starts at {t0:.2f} ns")

    for i in range(n_pulses):
        
        # add 1 ramped pulse
        points += [(t0, 0.0)]
        t0 += t_ramp
        points += [(t0, amplitude)]
        t0 += wait_q3
        points += [(t0, amplitude)]
        t0 += t_ramp
        points += [(t0, 0.0)]
        t0 += wait_q2 + Tq2       
    points += [(t_max-1, 0.0)]

    xy = np.array(points).T

    for i in range(len(points)-1):
        t_start, v_start = points[i]
        t_stop, v_stop = points[i+1]
        if t_stop > t_start:
            render_ramp(t_start, t_stop, v_start, v_stop, wave)

    wave = quantize_amplitude(wave)
    ta, out = get_awg_output(t, wave, analogue_shift, digital_filter_mode)
    
    for i in range(1):
        mx = np.logical_and(ta>xy[0][4*i+1], ta<xy[0][4*i+2]   )
        ta_mid = np.interp(amplitude/2, out[mx], ta[mx])
        deviation[n,2*i] =  ta_mid - (xy[0][4*i+1]+xy[0][4*i+2])/2  
    
        mx = np.logical_and(ta>xy[0][4*i+3], ta<xy[0][4*i+4]   )
        ta_mid = np.interp(amplitude/2, out[mx][::-1], ta[mx][::-1])
        deviation[n,2*i+1] =  ta_mid - (xy[0][4*i+3]+xy[0][4*i+4])/2  

collect = []
for idx in range(deviation.shape[0]):

    Ui = Q2Xgate(fq2, fq3, 0, t_hres= deviation[idx,:])
    Nrand = deviation.shape[0]
    tr2sum = 0
    tr2arr = []
    for i_rand in range(Nrand):
    
        R = Q2Xgate(fq2, fq3, 0, t_hres= deviation[i_rand,:])
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
# max  0.0001719630042238407
# min  5.1050957200593494e-05
# mean  0.00010120311863615617
plt.figure()
plt.hist(collect)


print('\n'+'qubit 2 X90 gate summary')
print(f'total_infidelity minimum = {total_infidelity + np.min(collect) :.5g}',  )
print(f'total_infidelity maximum= {total_infidelity + np.max(collect) :.5g}',  )




