# this file is modifid from Fig1/Fig_1d_25mT.py
# also from the following:
#  2023-09-01_q1q2_CPMG_25mT.py
# 2023-08-30-18-25-06_q1q2_Hahn_25mT.py
# 2023-08-30_and_09-03_T2star_25mT.py
# import os
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\Fig1','')
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
# from notebook_tools import fit_data
import glob
from datetime import datetime
import glob

datadir = os.getcwd()




#%%
start_time = '2023-08-30\\15-31-53'
end_time = start_time
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

data = datfiles[0]


xdata = data.waittime_set.ndarray[0,:]*1e-3
zdatas = data.su_S_North.ndarray#[:60]


from scipy.fft import fft, fftfreq

def osc_with_decay_alpha(t, A, f, y0, phi, T2, alpha):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0

def osc_with_decay(t, A, f, y0, phi, T2):
    return A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**2)+y0

def osc_without_decay(t, A, f, y0, phi):
    return A*np.cos(2*np.pi*t*f+phi)+y0





    
zdata = np.average(zdatas, axis=0)
    
zf = fft(zdata)
xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
f0 = xf[np.argmax(abs(zf[1:])) + 1]    

plt.figure()
plt.subplots_adjust(top=0.8)
p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 5, 1.5]
p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay_alpha,p0=p0, plot=True, return_cov=True)
# plt.title(data.metadata['location'][:19] + '\n' +\
#                   'A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0   \n' + \
#           ''.join(['{:.3g}, ']*len(p1)).format(*p1)  )
# plt.xlim(0, np.max(xdata))
plt.ylim(0,1)
plt.yticks([0, 0.25, 0.5, 0.75, 1])


print( f'T2* = {p1[4]:.4g} +- {p1std[4]:.4g} us')

xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
q1_t2star = dict(xdata=xdata, zdata=zdata, xx=xx, yy=osc_with_decay_alpha(xx,*p1))


start_time = '2023-08-30\\15-58-03'
end_time = start_time
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

data = datfiles[0]


xdata = data.waittime_set.ndarray[0,:] * 1e-3
zdatas = data.su_S_North.ndarray#[:60]




    
zdata = np.average(zdatas, axis=0)
    
zf = fft(zdata)
xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
f0 = xf[np.argmax(abs(zf[1:])) + 1]    

plt.figure()
plt.subplots_adjust(top=0.8)
p0 = [abs(np.max(zdata)-np.min(zdata))/2, f0, np.mean(zdata), 0 , 5, 1.5]
p1std, p1 = fit_data(xdata, zdata,func=osc_with_decay_alpha,p0=p0, plot=True, return_cov=True)
# plt.title(data.metadata['location'][:19] + '\n' +\
#                   'A*np.cos(2*np.pi*t*f+phi)*np.exp(-(t/T2)**alpha)+y0   \n' + \
#           ''.join(['{:.3g}, ']*len(p1)).format(*p1)  )
plt.ylim(0,1)
plt.yticks([0, 0.25, 0.5, 0.75, 1])    
    

print( f'T2* = {p1[4]:.4g} +- {p1std[4]:.4g} us')
xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
q2_t2star = dict(xdata=xdata, zdata=zdata, xx=xx, yy=osc_with_decay_alpha(xx,*p1))














#%%
start_time = '2023-08-30\\18-25-06'
end_time = start_time


datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]
    

def expdecay(x, A, x0, alpha, y0):
    return A*np.exp(-(x/x0)**alpha) + y0
x, y = 2*datfile.waittime_set.ndarray/1e3, datfile.su_S_North.ndarray



# plt.figure()
p0 = [0.45, 50, 1.5, 0.5]
p1std, p1 = fit_data(x, y, func=expdecay,p0=p0, plot=False, return_cov=True)

plt.figure(figsize=(5,5))
plt.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top = 0.9)
# plt.subplot(121)
plt.scatter(x, y, s=2)
xx = np.linspace(min(x), max(x), 10*len(x))
plt.plot(xx, expdecay(xx, *p1))
plt.xlabel('time '+ r'$2\tau$'+'(us)')
plt.ylabel('Hahn echo amplitude ')

# plt.subplot(122)

plt.title(start_time +'\n' + f'T2H = {p1[1]:.3g} us, alpha={p1[2]:.3g}',fontsize=12)

q1_hahn = dict(x=x, y=y, xx=xx, yy=expdecay(xx, *p1))

print( f'T2Hahn = {p1[1]:.4g} +- {p1std[1]:.4g} us')

#%%
start_time = '2023-08-30\\18-28-23'
end_time = start_time


datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

datfile = datfiles[0]


def expdecay(x, A, x0, alpha, y0):
    return A*np.exp(-(x/x0)**alpha) + y0
x, y = 2*datfile.waittime_set.ndarray/1e3, datfile.su_S_North.ndarray


# plt.figure()
p0 = [0.45, 30, 1.5, 0.5]
p1std, p1 = fit_data(x, y, func=expdecay,p0=p0, plot=False, return_cov=True)


plt.figure(figsize=(5,5))
plt.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top = 0.9)
# plt.subplot(121)
plt.scatter(x, y, s=2)
xx = np.linspace(min(x), max(x), 10*len(x))
plt.plot(xx, expdecay(xx, *p1))
plt.xlabel('time '+ r'$2\tau$'+'(us)')
plt.ylabel('Hahn echo amplitude ')

# plt.subplot(122)

plt.title(start_time +'\n' + f'T2H = {p1[1]:.3g} us, alpha={p1[2]:.3g}',fontsize=12)
    
q2_hahn = dict(x=x, y=y, xx=xx, yy=expdecay(xx, *p1))

print( f'T2Hahn = {p1[1]:.4g} +- {p1std[1]:.4g} us')
#%%
# plt.rcParams.update({'font.size':18})

datadir = os.getcwd()



start_time = [    
'\\2023-09-01\\' +    '21-18-50_sweep1D_waittime_q1_CPMG_512_q2_down' , 
              ]


factors = np.array([ 512])
# factors = [ 64, 128, 256, 512]


names = [glob.glob(datadir + x)[0]  for x in start_time]
datfiles = []
for i in range(len(names)):
    datfiles += [load_data(names[i])]

xdatalist = []
zdatalist = []
histlist = []
for i in range(len(names)):
    # xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su_S_North', ['ndarray'] )
    xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su_trig0', ['ndarray'] )
    zdatalist += [  zdata ]
    xdatalist += [ xdata * factors[i] * 2 *1e-3 ]
    # histlist += [ [datfiles[i].sensor_val_S_North_set.ndarray, datfiles[i].hist_S_North.ndarray] ]
    histlist += [ [datfiles[i].sensorval_S_North_trig0_set.ndarray, datfiles[i].hist_trig0.ndarray] ]
    # histlist += [ [datfiles[i].sensorval_S_North_trig1_set.ndarray, datfiles[i].hist_trig1.ndarray] ]
    
    
def expdecay(x, A, x0, alpha, y0):
    return A*np.exp(-(x/x0)**alpha) + y0


zdata_newlist = []
for line_id in range(len(names)):
    zdata_new = 1- hist_auto_1d(histlist[line_id][0], histlist[line_id][1], xdatalist[line_id], expect_pos_global=False, plot_thresholding=False )[0]
    zdata_newlist += [zdata_new]



q1_cpmg = {}
p1s = []
T2s = []
alphas = []
plt.figure(figsize=(5,4))
for i in range(len(factors)):
# for i in range(7,10):
    zoffset = i*0.0
    xdata, zdata = xdatalist[i], zdata_newlist[i]
    q1_cpmg = dict(xdata=xdata, zdata=zdata,  )
    plt.scatter(xdata, zdata+zoffset, s=2, label = f'CPMG {factors[i]}')
    if factors[i] == 512:
        remove_pts =  list(range(38, 51))
        mask = list( set(range(51)) - set(remove_pts) )
        xdata = xdata[mask]
        zdata = zdata[mask]
        
        
    # plt.scatter(xdata, zdata+zoffset, s=2, label = f'CPMG {factors[i]}')
    p0 = [0.45, 50*np.sqrt(factors[i]), 1.5, 0.52]
    p1std, p1 = fit_data(xdata, zdata, func=expdecay,p0=p0, plot=False, return_cov=True)
    xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
    plt.plot(xx, expdecay(xx, *p1)+ zoffset )
    p1s += [p1]
    T2s += [p1[1]]
    alphas += [p1[2]]
    
    q1_cpmg.update( dict( xx=xx, yy=expdecay(xx, *p1), p1=p1 ) )
    print( f'{p1[1]:.4g} +- {p1std[1]:.4g}')
# plt.gca().axes.yaxis.set_ticklabels([])
plt.yticks([0.5, 1])
# plt.xscale('log')    
plt.xlim(10,)    
print(  ''.join(['{:.3g}, ']*len(T2s)).format(*T2s))


#%%
datadir = os.getcwd()

start_time = [    
'\\2023-09-01\\' +    '18-36-22_sweep1D_waittime_q2_CPMG_512_q1_down'    ,
              ]


factors = np.array([ 512])


names = [glob.glob(datadir + x)[0]  for x in start_time]
datfiles = []
for i in range(len(names)):
    datfiles += [load_data(names[i])]

xdatalist = []
zdatalist = []
histlist = []
for i in range(len(names)):
    # xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su_S_North', ['ndarray'] )
    xdata,  zdata = data1d_attr(datfiles[i],  ['ndarray'], 'su_trig1', ['ndarray'] )
    zdatalist += [  zdata ]
    xdatalist += [ xdata * factors[i] * 2 *1e-3 ]
    # histlist += [ [datfiles[i].sensor_val_S_North_set.ndarray, datfiles[i].hist_S_North.ndarray] ]
    histlist += [ [datfiles[i].sensorval_S_North_trig1_set.ndarray, datfiles[i].hist_trig1.ndarray] ]
def expdecay(x, A, x0, alpha, y0):
    return A*np.exp(-(x/x0)**alpha) + y0


zdata_newlist = []
for line_id in range(len(names)):
    zdata_new = 1- hist_auto_1d(histlist[line_id][0], histlist[line_id][1], xdatalist[line_id], expect_pos_global=False, plot_thresholding=False )[0]
    zdata_newlist += [1-zdata_new]

q2_cpmg = {}
p1s = []
T2s = []
alphas = []
plt.figure(figsize=(5,4))
for i in range(len(factors)):
# for i in range(7,10):
    zoffset = i*0.0
    xdata, zdata = xdatalist[i], zdata_newlist[i]
    plt.scatter(xdata, zdata+zoffset, s=6, label = f'CPMG {factors[i]}')
    p0 = [0.45, 50*np.sqrt(factors[i]), 1.5, 0.52]
    p1std, p1 = fit_data(xdata, zdata, func=expdecay,p0=p0, plot=False, return_cov=True)
    xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
    plt.plot(xx, expdecay(xx, *p1)+ zoffset )
    p1s += [p1]
    T2s += [p1[1]]
    alphas += [p1[2]]    
    
    q2_cpmg = dict(xdata=xdata, zdata=zdata, xx=xx, yy=expdecay(xx, *p1), p1=p1 )   

plt.yticks([0.5, 1])
# plt.xscale('log')    
plt.xlim(1,)    
print(  ''.join(['{:.3g}, ']*len(T2s)).format(*T2s))





#%%
CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5)]
v = [Size.Fixed(0.5), Size.Fixed(1.0), Size.Fixed(0.7), Size.Fixed(1.0)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
ax13 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=3))
ax31 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=3, ny=1))
ax33 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=3, ny=3))
ax51 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=5, ny=1))
ax53 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=5, ny=3))

ax = ax13
ax.scatter(q1_t2star['xdata'], q1_t2star['zdata'], s=3, c=CQ1)
ax.plot(q1_t2star['xx'], q1_t2star['yy'], 'k', linewidth=0.5)
ax.set_xlim(0, max(q1_t2star['xdata']) )
# ax.set_ylim(0,1)
# ax.set_yticks([0, 1])    
ax.set_yticks([0.1, 0.9])    

ax = ax11
ax.scatter(q2_t2star['xdata'], q2_t2star['zdata'], s=3, c=CQ2)
ax.plot(q2_t2star['xx'], q2_t2star['yy'], 'k', linewidth=0.5)
ax.set_xlim(0, max(q2_t2star['xdata']) )
# ax.set_ylim(0,1)
# ax.set_yticks([0, 1])    
ax.set_yticks([0.1, 0.9])    

ax = ax33
ax.scatter(q1_hahn['x'], q1_hahn['y'], s=3, c=CQ1)
ax.plot(q1_hahn['xx'], q1_hahn['yy'], 'k', linewidth=0.5)
ax.set_xlim(0, max(q1_hahn['x']) )
ax.set_ylim(0.45, 0.95)
ax.set_yticks([0.5,  0.9])    


ax = ax31
ax.scatter(q2_hahn['x'], q2_hahn['y'], s=3, c=CQ2)
ax.plot(q2_hahn['xx'], q2_hahn['yy'], 'k', linewidth=0.5)
ax.set_xlim(0, max(q2_hahn['x']) )
ax.set_ylim(0.45, 0.95)
ax.set_yticks([0.5,  0.9])   



ax = ax53
ax.scatter(q1_cpmg['xdata']/1000, q1_cpmg['zdata'], s=3, c=CQ1)
ax.plot(q1_cpmg['xx']/1000, q1_cpmg['yy'], 'k', linewidth=0.5)
ax.set_xlim(0, max(q1_cpmg['xdata']/1000) )
ax.set_ylim(0.3, 0.9)
ax.set_yticks([0.4, 0.8])  
# ax.set_xlim(0.1, max(q1_cpmg['xdata']/1000) ) 
# ax.set_xscale('log')


ax = ax51
ax.scatter(q2_cpmg['xdata']/1000, q2_cpmg['zdata'], s=3, c=CQ2)
ax.plot(q2_cpmg['xx']/1000, q2_cpmg['yy'], 'k', linewidth=0.5)
ax.set_xlim(0, max(q2_cpmg['xdata']/1000) )
ax.set_ylim(0.4, 0.9)
ax.set_yticks([0.4, 0.8])   
# ax.set_xlim(0.1, max(q2_cpmg['xdata']/1000) )
# ax.set_xscale('log')


# filename = 'Figure1e'    # 'Fig_1d_25mT'
# plt.savefig(filename+'.pdf')








