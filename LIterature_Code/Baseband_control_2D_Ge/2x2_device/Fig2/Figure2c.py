# this file is modified from Fig2/Fig_2c_exchange_25mT.py

import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\Fig2','')
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
#%%
f_ramsey = dict(fq1_q2dn=[], fq1_q2up=[], fq2_q1dn=[], fq2_q1up=[],  )
example_plots = dict()

#%%
start_time = '2023-08-29\\19-13-56'
end_time = start_time
datadir = os.getcwd()
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


from scipy.fft import fft, fftfreq

xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  


from scipy.fft import fft, fftfreq

def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0

plt.close('all')

guess = [0.0378, ]
fq1_q2dn = []

for i in range(len(ydatas)):
    
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()
    
    if i ==0:
        f0 = 0.0378
    else:
        f0 = p1[1]
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)
    
        
    fq1_q2dn += [p1[1]]
    f_ramsey['fq1_q2dn'] += [ [ydatas[i], p1[1],p1std[1] ]    ]

    
fq1_q2dn = np.array(fq1_q2dn)    

# plt.close('all')
#%%
start_time = '2023-08-29\\19-24-28'
end_time = start_time
datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  

plt.close('all')


fq1_q2up = []

for i in range(len(ydatas)):
    
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()

    
    if i == 0:
        f0 = 0.07
    elif i==1:
        f0 = 0.07
    else:
        f0 = p1[1]
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)
        
    fq1_q2up += [p1[1]]
    f_ramsey['fq1_q2up'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    
fq1_q2up = np.array(fq1_q2up)    

    
#%%
plt.figure()
plt.subplot(121)
plt.plot(ydatas, fq1_q2up*1e3)
plt.plot(ydatas, fq1_q2dn*1e3)
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel('fq1 (MHz)')

plt.subplot(122)
plt.plot(ydatas, (fq1_q2up - fq1_q2dn)*1e3   )
plt.yscale('log')
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel('fq1 (MHz)')
plt.grid()        
# plt.ylim(0.1, 30)



#%%
start_time = '2023-08-29\\19-35-13'
end_time = start_time
datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


from scipy.fft import fft, fftfreq


ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  


from scipy.fft import fft, fftfreq

def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0


plt.close('all')

fq2_q1dn = []

for i in range(len(ydatas)):
    xdata = datfile.t_cz_set.ndarray[0,:] 
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()

    if i >= 15:
        f0 = 0.09
    
    if i == 6:
        zdata = zdata[xdata<=80]
        xdata = xdata[xdata<=80]
                
    elif i == 7:
        zdata = zdata[xdata<=60]
        xdata = xdata[xdata<=60]  
                 
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)
    plt.plot(xdata, zdata)
    # 
        
    fq2_q1dn += [p1[1]]
    f_ramsey['fq2_q1dn'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    
fq2_q1dn = np.array(fq2_q1dn)    


#%%
start_time = '2023-08-29\\19-45-52'
end_time = start_time
datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  



plt.close('all')

fq2_q1up = []

for i in    range(len(ydatas)):
    
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()

    bounds = ((-np.inf, )*5, (np.inf, )*5   )
    phi = 0.4
    T2 = 200
    if i < 2:
        bounds = ((0.25, -np.inf, -np.inf, -np.inf, 100), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
        # f0 = 0.095
        phi = -0.2
        T2 = 800
    elif i==2:
        bounds = ((0.3, -np.inf, -np.inf, -np.inf, 100), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
        f0 = 0.118
        phi = -2
        T2 = 400        
    elif i==6:
        bounds = ((0.1, -np.inf, -np.inf, -np.inf, 400), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
        f0 = 0.14
        phi = 1.5
        T2 = 600       
    elif i==7:
        bounds = ((0.1, -np.inf, -np.inf, -np.inf, 400), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
        # f0 = 0.14
        phi = 1.5
        T2 = 600         
        
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), phi , T2]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True, bounds = bounds )

        
    fq2_q1up += [p1[1]]
    f_ramsey['fq2_q1up'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    
fq2_q1up = np.array(fq2_q1up)    


    
#%%
plt.figure()
plt.subplot(121)
plt.plot(ydatas, fq2_q1up*1e3)
plt.plot(ydatas, fq2_q1dn*1e3)
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel('fq2 (MHz)')

plt.subplot(122)
plt.plot(ydatas, (fq2_q1up - fq2_q1dn)*1e3   )
plt.yscale('log')
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel('fq2 (MHz)')
plt.grid()        
# plt.ylim(0.1, 30)




























#%%
start_time = '2023-08-29\\12-01-48'
end_time = start_time
datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


from scipy.fft import fft, fftfreq

xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  


from scipy.fft import fft, fftfreq

def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0

plt.close('all')

guess = [0.0378, ]
fq1_q2dn = []

for i in range(len(ydatas)):
    
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()

    
    if i ==0:
        f0 = 0.0378
    else:
        f0 = p1[1]
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)
    
        
    fq1_q2dn += [p1[1]]
    f_ramsey['fq1_q2dn'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    
    if -66 <= ydatas[i] <= -64:
        example_plots['fq1_q2dn'] = dict(xdata=xdata, zdata=zdata, p1=p1, func=osc_with_decay)
    
fq1_q2dn = np.array(fq1_q2dn)    

# plt.close('all')
#%%
start_time = '2023-08-29\\12-10-12'
end_time = start_time
datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  

plt.close('all')


fq1_q2up = []

for i in range(len(ydatas)):
    
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()

    if i < 2:
        f0 = 0.055
    else:
        f0 = p1[1]
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)

        
    fq1_q2up += [p1[1]]
    f_ramsey['fq1_q2up'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    if -66 <= ydatas[i] <= -64:
        example_plots['fq1_q2up'] = dict(xdata=xdata, zdata=zdata, p1=p1, func=osc_with_decay)
        
fq1_q2up = np.array(fq1_q2up)    

    
#%%
plt.figure()
plt.subplot(121)
plt.plot(ydatas, fq1_q2up*1e3)
plt.plot(ydatas, fq1_q2dn*1e3)
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel('fq1 (MHz)')

plt.subplot(122)
plt.plot(ydatas, (fq1_q2up - fq1_q2dn)*1e3   )
plt.yscale('log')
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel('fq1 (MHz)')
plt.grid()        
# plt.ylim(0.1, 30)




#%%
start_time = '2023-08-29\\12-19-00'
end_time = start_time
datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


from scipy.fft import fft, fftfreq

xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  


from scipy.fft import fft, fftfreq

def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0


plt.close('all')

fq2_q1dn = []

for i in range(len(ydatas)):
    
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()

    
    if i >= 15:
        f0 = 0.09
    # else:
    #     f0 = p1[1]
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)
    
        
    fq2_q1dn += [p1[1]]
    f_ramsey['fq2_q1dn'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    
fq2_q1dn = np.array(fq2_q1dn)    


#%%
start_time = '2023-08-29\\12-27-29'
end_time = start_time
datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  



plt.close('all')

fq2_q1up = []

for i in    range(len(ydatas)):
    
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()
    bounds = ((-np.inf, )*5, (np.inf, )*5   )
    phi = 0.4
    T2 = 200
    if i < 2:
        bounds = ((0.25, -np.inf, -np.inf, -np.inf, 100), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
        # f0 = 0.095
        phi = -0.2
        T2 = 800
    elif i==2:
        bounds = ((0.3, -np.inf, -np.inf, -np.inf, 100), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
        f0 = 0.118
        phi = -2
        T2 = 400        
    elif i==6:
        bounds = ((0.25, -np.inf, -np.inf, -np.inf, 400), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
        f0 = 0.11
        phi = 1.5
        T2 = 600       
        
        
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), phi , T2]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True, bounds = bounds )
        

    
        
    fq2_q1up += [p1[1]]
    f_ramsey['fq2_q1up'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    
fq2_q1up = np.array(fq2_q1up)    


    
#%%
plt.figure()
plt.subplot(121)
plt.plot(ydatas, fq2_q1up*1e3)
plt.plot(ydatas, fq2_q1dn*1e3)
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel('fq2 (MHz)')

plt.subplot(122)
plt.plot(ydatas, (fq2_q1up - fq2_q1dn)*1e3   )
plt.yscale('log')
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel('fq2 (MHz)')
plt.grid()        
# plt.ylim(0.1, 30)












#%%
start_time = '2023-08-29\\09-20-06'
end_time = start_time
datadir = os.getcwd()
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


from scipy.fft import fft, fftfreq

xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  


from scipy.fft import fft, fftfreq

def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0

plt.close('all')

guess = [0.0378, ]
fq1_q2dn = []

for i in range(len(ydatas)):

    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()

    
    if i ==0:
        f0 = 0.0378
    else:
        f0 = p1[1]
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)
    
    # plt.figure()
    # plt.subplots_adjust(top=0.7)
    # plt.scatter(xdata, zdata)
    # xx = np.linspace(min(xdata), max(xdata), len(xdata)*10-9)
    # plt.plot(xx, osc_without_decay(xx, *p1))
    
    # plt.title(start_time + '\n' + \
    #         'A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0   \n' + \
    #           ''.join(['{:.3g}, ']*len(p1)).format(*p1)  )
        
    fq1_q2dn += [ [ydatas[i], p1[1],p1std[1] ]    ]
    f_ramsey['fq1_q2dn'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    
fq1_q2dn = np.array(fq1_q2dn)    

# plt.close('all')
#%%
start_time = '2023-08-29\\09-50-25'
end_time = start_time
datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  

plt.close('all')


fq1_q2up = []

for i in range(len(ydatas)):
    
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()
    # plt.plot(xf, np.abs(zf))
    
    if i < 2:
        f0 = 0.055
    else:
        f0 = p1[1]
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)
    

        
    fq1_q2up += [ [ydatas[i], p1[1],p1std[1] ]    ]
    f_ramsey['fq1_q2up'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    
fq1_q2up = np.array(fq1_q2up)    

    
#%%
plt.figure()
plt.subplot(121)
plt.plot(ydatas, fq1_q2up*1e3)
plt.plot(ydatas, fq1_q2dn*1e3)
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel('fq1 (MHz)')

plt.subplot(122)
plt.plot(ydatas, (fq1_q2up - fq1_q2dn)*1e3   )
plt.yscale('log')
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel('fq1 (MHz)')
plt.grid()        
# plt.ylim(0.1, 30)



    
#%%
start_time = '2023-08-29\\10-22-03'
end_time = start_time
datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


from scipy.fft import fft, fftfreq

xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  


from scipy.fft import fft, fftfreq

def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0


plt.close('all')

fq2_q1dn = []

for i in range(len(ydatas)):
    
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()

    
    if i >= 15:
        f0 = 0.09
    # else:
    #     f0 = p1[1]
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 200]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)
    
        
    fq2_q1dn += [ [ydatas[i], p1[1],p1std[1] ]    ]
    f_ramsey['fq2_q1dn'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    
fq2_q1dn = np.array(fq2_q1dn)    


#%%
start_time = '2023-08-29\\10-51-54'
end_time = start_time
datadir = os.getcwd() 
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


xdata = datfile.t_cz_set.ndarray[0,:] 
ydatas = datfile._1100_cz_vpB12_set.ndarray 
zdatas = datfile.su_S_North.ndarray  



plt.close('all')

fq2_q1up = []

for i in     range(len(ydatas)):
    
    zdata = zdatas[i,:]
        
    zf = fft(zdata[0:21])
    xf = fftfreq(len(xdata[0:21]), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    plt.figure()
    bounds = ((-np.inf, )*5, (np.inf, )*5   )
    phi = 0.4
    T2 = 400
    if i < 2:
        bounds = ((0.35, -np.inf, -np.inf, -np.inf, 600), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
        # f0 = 0.095
        phi = -0.2
        T2 = 800
    elif i==2:
        bounds = ((0.35, -np.inf, -np.inf, -np.inf, 300), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
        f0 = 0.097
        phi = -2
        T2 = 400        
        
    elif 5 <= i <= 7:
        f0 = 0.1
    elif i >= 12 and i<15:
        f0 = 0.105
    # elif i in [8,]:
    #     f0 = 0.1
    #     phi = 1.8
    elif i >= 15:
        f0 = 0.092
        phi = 0.2
        T2 = 600
    elif 8 <= i <= 11:    
        bounds = ((0.35, -np.inf, -np.inf, -np.inf, 600), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
        phi = 0.5
        T2 = 600
        
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), phi , T2]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True, bounds = bounds )
        
        
    fq2_q1up += [ [ydatas[i], p1[1],p1std[1] ]    ]
    f_ramsey['fq2_q1up'] += [ [ydatas[i], p1[1],p1std[1] ]    ]
    
fq2_q1up = np.array(fq2_q1up)    


    




#%%

f_ramsey_long = dict(fq1_q2dn=[], fq1_q2up=[], fq2_q1dn=[], fq2_q1up=[],  )



#%%

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0


start_time_list = [  '2023-08-29\\16-36-40', 
                   '2023-08-29\\16-44-50', 
                   '2023-08-29\\16-54-36', 
                   '2023-08-29\\17-02-45', ]
ydatas = [ -10, 0, 10, 20 ]

plt.close('all')


for i in  range(len(ydatas)):
    start_time = start_time_list[i]
    end_time = start_time
    datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
    
    data = datfiles[0]
    
    
    xdata = data.t_cz_set.ndarray
    zdata = data.su_S_North.ndarray#[:60]

        
    zf = fft(zdata)
    xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
    if i in [0,2]:
        bounds = ((-np.inf, -np.inf, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf, np.inf)   )
    elif i==3:
        bounds = ((0.3, 0.0003, -np.inf, -np.inf, 2e3), (0.5, 0.002, np.inf, np.inf, 10e3)   )
    elif i==1:
        bounds = ((0.3, 0.00045, 0.4, -np.inf, 0), (0.5, 0.0007, 0.6, np.inf, 10e3)   )
        
    plt.figure()
    plt.subplots_adjust(top=0.8)
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 5000,]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True,  bounds = bounds)
    # print(p1)
    # plt.title(data.metadata['location'][:19] + '\n' +\
    #                   'A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0   \n' + \
    #           ''.join(['{:.3g}, ']*len(p1)).format(*p1)  )
    # plt.xlim(0, np.max(xdata))
    # plt.ylim(0,1)
    # plt.yticks([0, 0.25, 0.5, 0.75, 1])
    # plt.plot(xdata, osc_with_decay(xdata, p1[0], p1[1], p1[2], p1[3] ,p1[4]  ) )
    # plt.plot(xdata, osc_with_decay(xdata, p1[0], device.settings[f'q{qubit_id}/x90/frequency']/1e9, p1[2], p1[3] ,p1[4]  ) )
    
    # print( 2*1e3/abs(xdata[1]-xdata[0])- p1[1]*1e3  )
    # print(p1[4])
    # print( f'{p1[4]:.4g} +- {p1std[4]:.4g}')
    
    
    
    f_ramsey_long['fq1_q2dn'] += [ [ ydatas[i] , 2/abs(xdata[1]-xdata[0])- p1[1], p1std[1] ]    ]




#%%

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0


start_time_list = [   '2023-08-29\\16-39-13', 
                   '2023-08-29\\16-50-07', 
                   '2023-08-29\\16-56-27', 
                   '2023-08-29\\17-00-31', 
 ]
ydatas = [ -10, 0, 10, 20 ]

plt.close('all')


for i in  range(len(ydatas)):
    start_time = start_time_list[i]
    end_time = start_time
    datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
    
    data = datfiles[0]
    
    
    xdata = data.t_cz_set.ndarray
    zdata = data.su_S_North.ndarray#[:60]

        
    zf = fft(zdata)
    xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    
    
        
    plt.figure()
    plt.subplots_adjust(top=0.8)
    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 5000,]
    p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay,p0=p0, plot=True, return_cov=True)#,  bounds = bounds)
    # print(p1)
    # plt.title(data.metadata['location'][:19] + '\n' +\
    #                   'A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0   \n' + \
    #           ''.join(['{:.3g}, ']*len(p1)).format(*p1)  )
    # plt.xlim(0, np.max(xdata))
    # plt.ylim(0,1)
    # plt.yticks([0, 0.25, 0.5, 0.75, 1])
    # plt.plot(xdata, osc_with_decay(xdata, p1[0], p1[1], p1[2], p1[3] ,p1[4]  ) )
    # plt.plot(xdata, osc_with_decay(xdata, p1[0], device.settings[f'q{qubit_id}/x90/frequency']/1e9, p1[2], p1[3] ,p1[4]  ) )
    
    # print( 2*1e3/abs(xdata[1]-xdata[0])- p1[1]*1e3  )
    # print(p1[4])
    # print( f'{p1[4]:.4g} +- {p1std[4]:.4g}')
    
    
    
    f_ramsey_long['fq1_q2up'] += [ [ ydatas[i] , 2/abs(xdata[1]-xdata[0])- p1[1], p1std[1] ]    ]


















#%%

hahn_exchange = []


#%%  q1 hahn, flip q2 in the second part of wait time. 'waittime_flip_other' means the time duration when q2 is in the flpped status


datadir = os.getcwd()


start_time = [    
              '\\2023-08-29\\' +    '17-08-27_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # waittime_flip_other = 10* 23.496316 ns
              '\\2023-08-29\\' +    '17-15-24_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # waittime_flip_other = 110* 23.496316 ns
              '\\2023-08-29\\' +    '17-15-01_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # waittime_flip_other = 210* 23.496316 ns
              '\\2023-08-29\\' +    '17-15-48_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # waittime_flip_other = 310* 23.496316 ns
'\\2023-08-29\\' +    '17-08-57_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # waittime_flip_other = 410* 23.496316 ns
              ]



names = [glob.glob(datadir + x)[0]  for x in start_time]
datfiles = []
for i in range(len(names)):
    datfiles += [load_data(names[i])]

xdatalist = []
zdatalist = []
for i in range(len(names)):
    xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su_S_North', ['ndarray'] )
    # xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su00', ['ndarray'] )
    # xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su10', ['ndarray'] )
    # xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su01', ['ndarray'] )
    zdatalist += [  zdata ]
    xdatalist += [ xdata ]
    
    
 
from scipy.fft import fft, fftfreq



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


titlestr = ''
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





p1arr = np.array(p1s)
p1arr[:,1] = np.abs(p1arr[:,1])


        
p1arr[0,3] -= np.pi



waittime_flip_other = np.array([10,110, 210, 310, 410]) * 23.496316

def func1d(x, a, b ):
    return a*x + b
p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3],func=func1d, plot=False, return_cov=True, sigma=[p1stds[i][3] for i in range(len(p1stds))])

plt.figure()
plt.scatter( waittime_flip_other, p1arr[:,3]  )
plt.plot( waittime_flip_other, func1d(waittime_flip_other, *p1_func1d)  )
plt.title( titlestr )
plt.xlabel(' waittime_flip_other')
plt.ylabel('fitted phase')

print( p1_func1d[0]/(2*np.pi)  )

hahn_exchange += [   [ 80,     p1_func1d[0]/(2*np.pi), p1_func1d_std[0]/(2*np.pi)  ]     ]

###
correction_phases = 0
y_error = np.array([p1stds[i][3] for i in range(len(p1stds))])
p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error)
# p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error, absolute_sigma=True)
v_setpoint = datfiles[0].metadata['settings']['q1,q2']['v_setpoints']['_1100_s']
vpB12_setp, vP1_setp, vP2_setp = [ v_setpoint[k] for k in ['vpB12', 'vP1', 'vP2'] ]
print('\n'+ f'vB12 at {vpB12_setp:.3g} mV, (vP1, vP2) at ({vP1_setp:.3g}, {vP2_setp:.3g}) mV:')
print(f'J = {1e6*p1_func1d[0]/(2*np.pi):.4g} +- {1e6*p1_func1d_std[0]/(2*np.pi):.3g} kHz' )
# print( p1_func1d[0]/(2*np.pi) , ' +- ', p1_func1d_std[0]/(2*np.pi) )
# 1.479663968405374e-05  +-  1.2371176179841088e-06
# 1.479663968405374e-05  +-  1.0155352774452808e-06 if add absolute_sigma=True
plt.figure()
for i in range(5):
    plt.errorbar(waittime_flip_other[i]/1e3, (p1arr[:,3]-correction_phases)[i],
                 yerr = y_error[i], c=f'C{i}',
                 fmt ='o', markersize=8, linestyle='', linewidth=2)

xxx = np.array([-1,11])*1e3
plt.plot( xxx/1e3, func1d(xxx, *p1_func1d) , 'k', linestyle='-', linewidth=1)



#%%  q1 hahn, flip q2 in the second part of wait time. 'waittime_flip_other' means the time duration when q2 is in the flpped status
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_cz-_1100_s'] = 50 
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_s-_1100_cz'] = 50 
# device.settings['q1,q2/v_setpoints/_1100_cz/vpB12'] =  20   <--- measure exchange here with Hahn type sequence




datadir = os.getcwd()




start_time = [    
              '\\2023-08-29\\' +    '18-03-16_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # reference 
              '\\2023-08-29\\' +    '18-04-07_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # one block = 100* 23.496316 ns
              '\\2023-08-29\\' +    '18-05-03_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # two block = 200* 23.496316 ns
              '\\2023-08-29\\' +    '18-05-59_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # three block = 300* 23.496316 ns
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
    _, p1 = fit_data(xdata, zdata,func=osc_without_decay,p0=p0, plot=False, return_cov=True)
    p1s += [p1]
    p1stds += [p1std]


titlestr = ''
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

p1arr[2,3] += np.pi
p1arr[3,3] += np.pi

waittime_flip_other = np.array([0,100, 200, 300]) * (1e9/datfiles[0].metadata['settings']['q1']['x90']['frequency']   )
correction_idle = (20+37.75+20) * 14.7e-6 * 2*np.pi
correction_ramp = (50+50) * (14.7e-6 + 107e-6)/2 * 2*np.pi
correction_ramp_error =  (50+50) * (14.7e-6 - 107e-6)/2 * 2*np.pi

correction_phases = np.array([0, 1, 2, 3]) * (correction_idle + correction_ramp)


def func1d(x, a, b ):
    return a*x + b

y_error_0 = np.array([p1stds[i][3] for i in range(len(p1stds))])
y_error = np.sqrt( correction_ramp_error**2 + y_error_0**2)
# p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error)
p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error, absolute_sigma=True)


_, ax = plt.subplots()
for i in range(4):
    ax.errorbar(waittime_flip_other[i]/1e3, (p1arr[:,3]-correction_phases)[i],
                 yerr = y_error[i], c=f'C{i}',
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=0.5)

# plt.plot( waittime_flip_other, np.poly1d(np.polyfit(waittime_flip_other,p1arr[:,3],1))(waittime_flip_other)  )
# ax.plot( waittime_flip_other, func1d(waittime_flip_other, *p1_func1d) ,  )
# ax.plot( waittime_flip_other/1e3, func1d(waittime_flip_other, *p1_func1d) , 'k', linestyle='-', linewidth=0.5)
xxx = np.array([-2,14])*1e3
ax.plot( xxx/1e3, func1d(xxx, *p1_func1d) , 'k', linestyle='-', linewidth=0.5)



v_setpoint = datfiles[0].metadata['settings']['q1,q2']['v_setpoints']['_1100_cz']
vpB12_setp, vP1_setp, vP2_setp = [ v_setpoint[k] for k in ['vpB12', 'vP1', 'vP2'] ]
print('\n'+ f'vB12 at {vpB12_setp:.3g} mV, (vP1, vP2) at ({vP1_setp:.3g}, {vP2_setp:.3g}) mV:')
print(f'J = {1e6*p1_func1d[0]/(2*np.pi):.4g} +- {1e6*p1_func1d_std[0]/(2*np.pi):.3g} kHz' )
# print( p1_func1d[0]/(2*np.pi) , ' +- ', p1_func1d_std[0]/(2*np.pi) )
# 0.00010737814389439737  +-  9.910729296842075e-07
# 0.00010737814389439737  +-  2.0583290074391397e-06 if set absolute_sigma=True
hahn_exchange += [   [ 20,     p1_func1d[0]/(2*np.pi), p1_func1d_std[0]/(2*np.pi)  ]     ]












#%%  q1 hahn, flip q2 in the second part of wait time. 'waittime_flip_other' means the time duration when q2 is in the flpped status
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_cz-_1100_s'] = 50
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_s-_1100_cz'] = 50 
# device.settings['q1,q2/v_setpoints/_1100_cz/vpB12'] =  40   <--- measure exchange here with Hahn type sequence


datadir = os.getcwd()




start_time = [    
              '\\2023-08-29\\' +    '18-28-46_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # reference 
              '\\2023-08-29\\' +    '18-29-55_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # one block = 100* 23.496316 ns
              '\\2023-08-29\\' +    '18-30-55_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # two block = 200* 23.496316 ns
              '\\2023-08-29\\' +    '18-31-41_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # three block = 300* 23.496316 ns
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


    zdata = zdatalist[i]
    xdata = xdatalist[i]
    zf = fft(zdata)
    xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    


    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 ]
    _, p1 = fit_data(xdata, zdata,func=osc_without_decay,p0=p0, plot=False, return_cov=True)
    p1s += [p1]
    p1stds += [p1std]


titlestr = ''
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
    xxx = np.linspace(xdata[0], xdata[-1], (len(xdata)-1)*5+1  )
    plt.plot(xxx, osc_without_decay(xxx, *p1s[i] ) )

    titlestr += (tempstr + ', ')
plt.legend()
    
plt.title( titlestr   )




p1arr = np.array(p1s)
p1arr[:,1] = np.abs(p1arr[:,1])



waittime_flip_other = np.array([0, 100, 200, 300]) * (1e9/datfiles[0].metadata['settings']['q1']['x90']['frequency']   )
correction_idle = (20+37.75+20) * 14.7e-6 * 2*np.pi
correction_ramp = (50+50) * (14.7e-6 + 31e-6)/2 * 2*np.pi
correction_ramp_error =  (50+50) * (14.7e-6 - 31e-6)/2 * 2*np.pi

correction_phases = np.array([0, 1, 2, 3]) * (correction_idle + correction_ramp)


def func1d(x, a, b ):
    return a*x + b

y_error_0 = np.array([p1stds[i][3] for i in range(len(p1stds))])
y_error = np.sqrt( correction_ramp_error**2 + y_error_0**2)
# p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error)
p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error, absolute_sigma=True)

# plt.figure()
_, ax = plt.subplots()
# plt.scatter( waittime_flip_other, p1arr[:,3]  )
# ax.plot(waittime_flip_other/1e3,  p1arr[:,3]-correction_phases, marker='o' , markersize=MRSIZE, linestyle='', linewidth=0.5)


# ax.errorbar(waittime_flip_other/1e3, p1arr[:,3]-correction_phases,
#              yerr = y_error, 
#              fmt ='o', markersize=MRSIZE, linestyle='', linewidth=0.5)
for i in range(4):
    ax.errorbar(waittime_flip_other[i]/1e3, (p1arr[:,3]-correction_phases)[i],
                 yerr = y_error[i], c=f'C{i}',
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=0.5)

# plt.plot( waittime_flip_other, np.poly1d(np.polyfit(waittime_flip_other,p1arr[:,3],1))(waittime_flip_other)  )
# ax.plot( waittime_flip_other, func1d(waittime_flip_other, *p1_func1d) ,  )
# ax.plot( waittime_flip_other/1e3, func1d(waittime_flip_other, *p1_func1d) , 'k', linestyle='-', linewidth=0.5)
xxx = np.array([-2,14])*1e3
ax.plot( xxx/1e3, func1d(xxx, *p1_func1d) , 'k', linestyle='-', linewidth=0.5)




v_setpoint = datfiles[0].metadata['settings']['q1,q2']['v_setpoints']['_1100_cz']
vpB12_setp, vP1_setp, vP2_setp = [ v_setpoint[k] for k in ['vpB12', 'vP1', 'vP2'] ]
print('\n'+ f'vB12 at {vpB12_setp:.3g} mV, (vP1, vP2) at ({vP1_setp:.3g}, {vP2_setp:.3g}) mV:')
print(f'J = {1e6*p1_func1d[0]/(2*np.pi):.4g} +- {1e6*p1_func1d_std[0]/(2*np.pi):.3g} kHz' )
# print( p1_func1d[0]/(2*np.pi) , ' +- ', p1_func1d_std[0]/(2*np.pi) )
# 3.131426299613079e-05  +-  1.1113817803159844e-06
# 3.131426299613079e-05  +-  1.8679437027282455e-06 if set absolute_sigma=True
hahn_exchange += [   [ 40,     p1_func1d[0]/(2*np.pi), p1_func1d_std[0]/(2*np.pi)  ]     ]





#%%  q1 hahn, flip q2 in the second part of wait time. 'waittime_flip_other' means the time duration when q2 is in the flpped status
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_cz-_1100_s'] = 50 
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_s-_1100_cz'] = 50 
# device.settings['q1,q2/v_setpoints/_1100_cz/vpB12'] =  60   <--- measure exchange here with Hahn type sequence

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
    _, p1 = fit_data(xdata, zdata,func=osc_without_decay,p0=p0, plot=False, return_cov=True)
    p1s += [p1]
    p1stds += [p1std]



titlestr = ''
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
    xxx = np.linspace(xdata[0], xdata[-1], (len(xdata)-1)*5+1  )
    plt.plot(xxx, osc_without_decay(xxx, *p1s[i] ) )

    titlestr += (tempstr + ', ')
plt.legend()
    
plt.title( titlestr   )



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
# p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error)
p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error, absolute_sigma=True)

# plt.figure()
_, ax = plt.subplots()
# plt.scatter( waittime_flip_other, p1arr[:,3]  )
# ax.plot(waittime_flip_other/1e3,  p1arr[:,3]-correction_phases, marker='o' , markersize=MRSIZE, linestyle='', linewidth=0.5)


# ax.errorbar(waittime_flip_other/1e3, p1arr[:,3]-correction_phases,
#              yerr = y_error, 
#              fmt ='o', markersize=MRSIZE, linestyle='', linewidth=0.5)
for i in range(4):
    ax.errorbar(waittime_flip_other[i]/1e3, (p1arr[:,3]-correction_phases)[i],
                 yerr = y_error[i], c=f'C{i}',
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=0.5)

# plt.plot( waittime_flip_other, np.poly1d(np.polyfit(waittime_flip_other,p1arr[:,3],1))(waittime_flip_other)  )
# ax.plot( waittime_flip_other, func1d(waittime_flip_other, *p1_func1d) ,  )
# ax.plot( waittime_flip_other/1e3, func1d(waittime_flip_other, *p1_func1d) , 'k', linestyle='-', linewidth=0.5)
xxx = np.array([-2,14])*1e3
ax.plot( xxx/1e3, func1d(xxx, *p1_func1d) , 'k', linestyle='-', linewidth=0.5)


ax.set_xlim(-1, 12)
ax.set_ylim(-0.55, 0.3)
ax.set_yticks([-0.4, -0.2, 0, 0.2])
ax.set_xticks([0, 5, 10])

v_setpoint = datfiles[0].metadata['settings']['q1,q2']['v_setpoints']['_1100_cz']
vpB12_setp, vP1_setp, vP2_setp = [ v_setpoint[k] for k in ['vpB12', 'vP1', 'vP2'] ]
print('\n'+ f'vB12 at {vpB12_setp:.3g} mV, (vP1, vP2) at ({vP1_setp:.3g}, {vP2_setp:.3g}) mV:')
print(f'J = {1e6*p1_func1d[0]/(2*np.pi):.4g} +- {1e6*p1_func1d_std[0]/(2*np.pi):.3g} kHz' )
# print( p1_func1d[0]/(2*np.pi) , ' +- ', p1_func1d_std[0]/(2*np.pi) )
# 8.794337756526609e-06  +-  5.83986508353492e-07
# 8.794337756526609e-06  +-  1.2415211781628862e-06 if set absolute_sigma=True
hahn_exchange += [   [ 60,     p1_func1d[0]/(2*np.pi), p1_func1d_std[0]/(2*np.pi)  ]     ]



#%%  q1 hahn, flip q2 in the second part of wait time. 'waittime_flip_other' means the time duration when q2 is in the flpped status
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_cz-_1100_s'] = 50 
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_s-_1100_cz'] = 50 
# device.settings['q1,q2/v_setpoints/_1100_cz/vpB12'] =  50   <--- measure exchange here with Hahn type sequence



datadir = os.getcwd()




start_time = [    
              '\\2023-08-29\\' +    '18-55-24_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # reference 
              '\\2023-08-29\\' +    '18-58-37_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # two block = 150* 23.496316 ns
              '\\2023-08-29\\' +    '18-59-29_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # three block = 300* 23.496316 ns
              
              '\\2023-08-29\\' +    '18-56-30_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # three block = 450* 23.496316 ns
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

    zdata = zdatalist[i]
    xdata = xdatalist[i]
    zf = fft(zdata)
    xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    


    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 ]
    _, p1 = fit_data(xdata, zdata,func=osc_without_decay,p0=p0, plot=False, return_cov=True)
    p1s += [p1]
    p1stds += [p1std]


titlestr = ''
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


waittime_flip_other = np.array([0, 150, 300, 450]) * (1e9/datfiles[0].metadata['settings']['q1']['x90']['frequency']   )
correction_idle = (20+37.75+20) * 14.7e-6 * 2*np.pi
correction_ramp = (50+50) * (14.7e-6 + 18e-6)/2 * 2*np.pi
correction_ramp_error =  (50+50) * (14.7e-6 - 18e-6)/2 * 2*np.pi

correction_phases = np.array([0, 1, 2, 3]) * (correction_idle + correction_ramp)


def func1d(x, a, b ):
    return a*x + b

y_error_0 = np.array([p1stds[i][3] for i in range(len(p1stds))])
y_error = np.sqrt( correction_ramp_error**2 + y_error_0**2)
p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error)
# p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error, absolute_sigma=True)
# plt.figure()
_, ax = plt.subplots()

for i in range(4):
    ax.errorbar(waittime_flip_other[i]/1e3, (p1arr[:,3]-correction_phases)[i],
                 yerr = y_error[i], c=f'C{i}',
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=0.5)


xxx = np.array([-2,14])*1e3
ax.plot( xxx/1e3, func1d(xxx, *p1_func1d) , 'k', linestyle='-', linewidth=0.5)



v_setpoint = datfiles[0].metadata['settings']['q1,q2']['v_setpoints']['_1100_cz']
vpB12_setp, vP1_setp, vP2_setp = [ v_setpoint[k] for k in ['vpB12', 'vP1', 'vP2'] ]
print('\n'+ f'vB12 at {vpB12_setp:.3g} mV, (vP1, vP2) at ({vP1_setp:.3g}, {vP2_setp:.3g}) mV:')
print(f'J = {1e6*p1_func1d[0]/(2*np.pi):.4g} +- {1e6*p1_func1d_std[0]/(2*np.pi):.3g} kHz' )
# print( p1_func1d[0]/(2*np.pi) , ' +- ', p1_func1d_std[0]/(2*np.pi) )
# 1.775589023033733e-05  +-  1.8946512814938308e-06
# 1.775589023033733e-05  +-  1.2411709983072937e-06 if set absolute_sigma=True
# hahn_exchange += [   [ 60,     p1_func1d[0]/(2*np.pi), p1_func1d_std[0]/(2*np.pi)  ]     ]


hahn_exchange += [   [ 50,     p1_func1d[0]/(2*np.pi), p1_func1d_std[0]/(2*np.pi)  ]     ]





#%%  q1 hahn, flip q2 in the second part of wait time. 'waittime_flip_other' means the time duration when q2 is in the flpped status
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_cz-_1100_s'] = 50 
# device.settings['q1,q2/v_setpoints/_ramp_times/_1100_s-_1100_cz'] = 50 
# device.settings['q1,q2/v_setpoints/_1100_cz/vpB12'] =  30   <--- measure exchange here with Hahn type sequence



datadir = os.getcwd()




start_time = [    
              '\\2023-08-29\\' +    '19-02-09_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # reference 
              '\\2023-08-29\\' +    '19-03-01_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # one block = 100* 23.496316 ns
              '\\2023-08-29\\' +    '19-04-05_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # two block = 200* 23.496316 ns
              '\\2023-08-29\\' +    '19-05-10_sweep1D_waittime_end_Hahn_residual_exchange'   ,    # three block = 300* 23.496316 ns
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


    zdata = zdatalist[i]
    xdata = xdatalist[i]
    zf = fft(zdata)
    xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
    f0 = xf[np.argmax(abs(zf[1:])) + 1]    


    p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 ]
    _, p1 = fit_data(xdata, zdata,func=osc_without_decay,p0=p0, plot=False, return_cov=True)
    p1s += [p1]
    p1stds += [p1std]


titlestr = ''
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
    xxx = np.linspace(xdata[0], xdata[-1], (len(xdata)-1)*5+1  )
    plt.plot(xxx, osc_without_decay(xxx, *p1s[i] ) )

    titlestr += (tempstr + ', ')
plt.legend()
    
plt.title( titlestr   )



p1arr = np.array(p1s)
p1arr[:,1] = np.abs(p1arr[:,1])



waittime_flip_other = np.array([0, 100, 200, 300]) * (1e9/datfiles[0].metadata['settings']['q1']['x90']['frequency']   )
correction_idle = (20+37.75+20) * 14.7e-6 * 2*np.pi
correction_ramp = (50+50) * (14.7e-6 + 58e-6)/2 * 2*np.pi
correction_ramp_error =  (50+50) * (14.7e-6 - 58e-6)/2 * 2*np.pi

correction_phases = np.array([0, 1, 2, 3]) * (correction_idle + correction_ramp)


def func1d(x, a, b ):
    return a*x + b

y_error_0 = np.array([p1stds[i][3] for i in range(len(p1stds))])
y_error = np.sqrt( correction_ramp_error**2 + y_error_0**2)
p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error)
# p1_func1d_std, p1_func1d = fit_data(waittime_flip_other, p1arr[:,3]-correction_phases,func=func1d, plot=False, return_cov=True, sigma=y_error, absolute_sigma=True)

# plt.figure()
_, ax = plt.subplots()

for i in range(4):
    ax.errorbar(waittime_flip_other[i]/1e3, (p1arr[:,3]-correction_phases)[i],
                 yerr = y_error[i], c=f'C{i}',
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=0.5)

xxx = np.array([-2,14])*1e3
ax.plot( xxx/1e3, func1d(xxx, *p1_func1d) , 'k', linestyle='-', linewidth=0.5)



v_setpoint = datfiles[0].metadata['settings']['q1,q2']['v_setpoints']['_1100_cz']
vpB12_setp, vP1_setp, vP2_setp = [ v_setpoint[k] for k in ['vpB12', 'vP1', 'vP2'] ]
print('\n'+ f'vB12 at {vpB12_setp:.3g} mV, (vP1, vP2) at ({vP1_setp:.3g}, {vP2_setp:.3g}) mV:')
print(f'J = {1e6*p1_func1d[0]/(2*np.pi):.4g} +- {1e6*p1_func1d_std[0]/(2*np.pi):.3g} kHz' )
# 5.765974046464719e-05  +-  2.189138831991015e-06
# 5.765974046464719e-05  +-  1.9065555924832138e-06 if set absolute_sigma = True
hahn_exchange += [   [ 30,     p1_func1d[0]/(2*np.pi), p1_func1d_std[0]/(2*np.pi)  ]     ]

















#%%
plt.close('all')
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size':15})



updn_upup = np.array([  f_ramsey['fq2_q1up'][i][1] for i in range(len(f_ramsey['fq2_q1up']))    ])*1e3
dndn_dnup = np.array([  f_ramsey['fq2_q1dn'][i][1] for i in range(len(f_ramsey['fq2_q1dn']))    ])*1e3
dnup_upup = np.array([  f_ramsey['fq1_q2up'][i][1] for i in range(len(f_ramsey['fq1_q2up']))    ])*1e3
dnup_upup_err = np.array([  f_ramsey['fq1_q2up'][i][2] for i in range(len(f_ramsey['fq1_q2up']))    ])*1e3
dndn_updn = np.array([  f_ramsey['fq1_q2dn'][i][1] for i in range(len(f_ramsey['fq1_q2dn']))    ])*1e3
dndn_updn_err = np.array([  f_ramsey['fq1_q2dn'][i][2] for i in range(len(f_ramsey['fq1_q2dn']))    ])*1e3
vpB12 = np.array([  f_ramsey['fq1_q2dn'][i][0] for i in range(len(f_ramsey['fq1_q2dn']))    ])


dnup_upup_long = np.array([  f_ramsey_long['fq1_q2up'][i][1] for i in range(len(f_ramsey_long['fq1_q2up']))    ])*1e3
dnup_upup_long_err = np.array([  f_ramsey_long['fq1_q2up'][i][2] for i in range(len(f_ramsey_long['fq1_q2up']))    ])*1e3
dndn_updn_long = np.array([  f_ramsey_long['fq1_q2dn'][i][1] for i in range(len(f_ramsey_long['fq1_q2dn']))    ])*1e3
dndn_updn_long_err = np.array([  f_ramsey_long['fq1_q2dn'][i][2] for i in range(len(f_ramsey_long['fq1_q2dn']))    ])*1e3
vpB12_long = np.array([  f_ramsey_long['fq1_q2dn'][i][0] for i in range(len(f_ramsey_long['fq1_q2dn']))    ])

df_hahn = np.array( [hahn_exchange[i][1] for i in range(len(hahn_exchange))] )*1e3
df_hahn_err = np.array( [hahn_exchange[i][2] for i in range(len(hahn_exchange))] )*1e3
vpB12_hahn = np.array( [hahn_exchange[i][0] for i in range(len(hahn_exchange))] )



MRSIZE = 2
from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(5,3))
h = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.0)]
v = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.0)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)

ax11 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=1))

plt.grid(linewidth=0.5)  
mask = vpB12<-10
# plt.scatter(vpB12[mask], (dnup_upup - dndn_updn)[mask]   , c='b', s=4,  label=r'$\Delta fq1$')
# plt.scatter(vpB12[mask], (updn_upup - dndn_dnup)[mask]   , c='b', label=r'$\Delta fq2$')
# plt.scatter(vpB12_long, (dnup_upup_long - dndn_updn_long)    , c='b', s=4,  label=r'$\Delta fq1$')


# plt.scatter(  vpB12_hahn[1:],  df_hahn[1:], c='r',  s=4, )
# plt.scatter(  vpB12_hahn[0],  df_hahn[0], c='r', marker='^',  s=4, )

ax = plt.gca()
#### ax.scatter(vpB12[mask], (dnup_upup - dndn_updn)[mask]   , c='b', s=SCSIZE,  label=r'$\Delta fq1$')
for i in range(len(vpB12[mask])):
    ax.errorbar(vpB12[mask][i], (dnup_upup - dndn_updn)[mask][i],
                 yerr = np.sqrt(dnup_upup_err**2+dndn_updn_err**2)[i], c='b',
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=1,  label=r'$\Delta fq1$')
# plt.scatter(vpB12[mask], (updn_upup - dndn_dnup)[mask]   , c='b', label=r'$\Delta fq2$')



### ax.scatter(vpB12_long, (dnup_upup_long - dndn_updn_long)    , c='b', s=SCSIZE,  label=r'$\Delta fq1$')
for i in range(len(vpB12_long)):
    ax.errorbar(vpB12_long[i], (dnup_upup_long - dndn_updn_long)[i],
                 yerr = np.sqrt(dnup_upup_long_err**2+dndn_updn_long_err**2)[i], c='b',
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=1,  label=r'$\Delta fq1$')
    
#### ax.scatter(  vpB12_hahn[1:],  df_hahn[1:], c='r',  s=5, )
for i in range(1,len( vpB12_hahn)):
    ax.errorbar( vpB12_hahn[i], df_hahn[i],
                 yerr = df_hahn_err[i], c='r',
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=1 )
#### ax.scatter(  vpB12_hahn[0],  df_hahn[0], c='r', marker='^',  s=5, )
for i in [0]:
    ax.errorbar( vpB12_hahn[i], df_hahn[i],
                 yerr = df_hahn_err[i], c='r',
                 fmt ='^', markersize=MRSIZE, linestyle='', linewidth=1 )


plt.yscale('log')
plt.xlabel('_1100_cz_vpB12 (mV)')
plt.ylabel(r'$\Delta fq$ (MHz)')
      
plt.ylim(0.005, 80)
plt.yticks([0.01, 0.1, 1, 10])
# plt.xlim(-90, 90)
# plt.xticks([-80, 0, 80])
# plt.legend()
# plt.plot(vpB12, 0.24*np.exp(-0.059*(vpB12-10) )  , 'C2', linewidth=0.5  )
xxx = np.linspace(-90,90,2)
plt.plot(xxx, 0.24*np.exp(-0.059*(xxx-10) )  , 'C2', linewidth=0.5  )
plt.xlim(-99, 90)
plt.xticks([-80, 0, 80])

# filename = 'Figure2c' # 'Fig_2c_exchange_25mT'
# plt.savefig(filename+'.pdf')




