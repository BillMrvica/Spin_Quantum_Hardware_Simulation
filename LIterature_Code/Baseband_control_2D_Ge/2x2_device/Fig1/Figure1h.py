# this file is modified from Fig1/Fig_1e.py
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
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
from matplotlib import pyplot as plt
import qtt
import sys

# from notebook_tools import fit_data

#%%

start_time = '2023-08-14\\12-28-17'

end_time = '2023-08-14\\12-53-59'

import os
datadir = os.getcwd()
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

datfile = datfiles[0]
datfile_1 = datfiles[1]
# datfile_2 = datfiles[2]

from scipy.fft import fft, fftfreq


xdata = datfile.su_S_North.ndarray 
xlabel = 'n_shuttles'
xunit = datfile.su_S_North.unit
zdata = datfile.n_shuttle_gate_set.ndarray


xdata_1 = datfile_1.su_S_North.ndarray 
xunit = datfile_1.su_S_North.unit
zdata_1 = datfile_1. n_shuttle_gate_set.ndarray


# xdata_2 = datfile_2 .su_S_North.ndarray 

# zdata_2 = datfile_2. n_shuttle_gate_set.ndarray


# zf = fft(zdata)
# xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
# zf = zf[xf>0]
# xf = xf[xf>0]

plt.figure(figsize=(9,4))
plt.subplots_adjust(bottom=0.15, top=0.8)
# plt.plot(xf, zf)
# plt.subplot(121)
# plt.plot(xf*1e3, np.log(np.abs(zf)))
# plt.ylabel('log FFT')
# plt.xlabel('freq (MHz)')

# plt.subplot(122)
plt.plot(zdata, xdata, 'bo--' ,label = 'q1x90 when q2 is down')
plt.plot(zdata_1, xdata_1, 'ro--', label = 'q1x90 when q2 is up')
# plt.plot(zdata_2, xdata_2, 'go--', label = 'no pulse on q1')
plt.legend()
# plt.plot(zdata_2, xdata_2, 'g')
plt.xlabel(xlabel + ' (' + xunit + ')')
plt.title(start_time+', shuttle gate')

q1_x90 = dict(y=xdata, x=zdata)




#%%

start_time = '2023-08-14\\14-35-23'

end_time = '2023-08-14\\14-38-20'

datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

datfile = datfiles[0]
datfile_1_ = datfiles[1]


from scipy.fft import fft, fftfreq


xdata = datfile.su_S_North.ndarray 
xlabel = 'n_shuttles'
xunit = datfile.su_S_North.unit
zdata = datfile.n_shuttle_gate_set.ndarray


xdata_1 = datfile_1_.su_S_North.ndarray
xunit = datfile_1.su_S_North.unit
zdata_1 = datfile_1_.n_shuttle_gate_set.ndarray


# xdata_2 = datfile_2 .su_S_North.ndarray 

# zdata_2 = datfile_2. n_shuttle_gate_set.ndarray


# zf = fft(zdata)
# xf = fftfreq(len(xdata), d=abs(xdata[1]-xdata[0])  )
# zf = zf[xf>0]
# xf = xf[xf>0]

plt.figure(figsize=(9,4))
plt.subplots_adjust(bottom=0.15, top=0.8)
# plt.plot(xf, zf)
# plt.subplot(121)
# plt.plot(xf*1e3, np.log(np.abs(zf)))
# plt.ylabel('log FFT')
# plt.xlabel('freq (MHz)')

# plt.subplot(122)
plt.plot(zdata, xdata, 'bo--' ,label = 'q2x90 when q1 is down')
plt.plot(zdata_1, xdata_1, 'ro--', label = 'q2x90 when q1 is up')
# plt.plot(zdata_2, xdata_2, 'go--', label = 'no pulse on q1')
plt.legend()
# plt.plot(zdata_2, xdata_2, 'g')
plt.xlabel(xlabel + ' (' + xunit + ')')
plt.title(start_time+', shuttle gate')

q2_x90 = dict(y=xdata, x=zdata)


#%%
from mpl_toolkits.axes_grid1 import Divider, Size
CQ1 = '#0D77BC'
CQ2 = '#DD5E27'

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(3.1), ]
v = [Size.Fixed(0.5), Size.Fixed(1.0), Size.Fixed(0.7), Size.Fixed(1.0)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))

ax13 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=3))


ax = ax13
# ax.scatter(q1_t2star['xdata'], q1_t2star['zdata'], s=3)
ax.plot(q1_x90['x'], q1_x90['y'], CQ1 , marker = 'o', markersize=3, linewidth=0.5)
ax.set_xlim(0, 80 )
ax.set_xticks([0, 40, 80]) 
ax.set_ylim(0,1)
ax.set_yticks([0, 1])    
# ax_twinx = ax.twinx()  
# ax_twinx.set_xticks([])
secax_y = ax.secondary_xaxis('top', functions=( lambda x: x*4, lambda x: x/4))   
secax_y.set_xlabel('number of shuttles')
secax_y.set_xticks([0, 160, 320]) 

ax = ax11
# ax.scatter(q2_x90['xdata'], q2_x90['zdata'], s=3)
ax.plot(q2_x90['x'], q2_x90['y'], CQ2 , marker='o', markersize=3, linewidth=0.5)
ax.set_xlim(0, 80 )
ax.set_xticks([0, 40, 80]) 
ax.set_ylim(0,1)
ax.set_yticks([0, 1])    
secax_y = ax.secondary_xaxis('top', functions=( lambda x: x*2, lambda x: x/2))   
secax_y.set_xlabel('number of shuttles')
secax_y.set_xticks([0, 80, 160]) 



# filename = 'Figure1h'
# plt.savefig(filename+'.pdf')







