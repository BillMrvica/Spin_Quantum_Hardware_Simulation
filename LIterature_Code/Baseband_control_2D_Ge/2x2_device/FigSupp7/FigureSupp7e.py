# this file is modified from Supp/unequal_time/equal_time_repeat_x90.py
import numpy as np
import matplotlib.pyplot as plt
from projects.notebook_tools.notebook_tools import get_data_from, fit_data
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
from matplotlib import pyplot as plt
import qtt
import sys


#%%
import os
start_time = '2023-07-17\\09-13-53'

end_time = start_time

datadir = os.getcwd()
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

datfile = datfiles[0]

from scipy.fft import fft, fftfreq


xdata = datfile.su_S_North.ndarray 
xlabel = 'n_shuttles'
xunit = datfile.su_S_North.unit
zdata = datfile.n_shuttle_gate_set.ndarray



q1_x90 = dict(y=xdata, x=zdata)




#%%
start_time = '2023-07-17\\09-44-53'

end_time = start_time

datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

datfile = datfiles[0]


from scipy.fft import fft, fftfreq


xdata = datfile.su_S_North.ndarray 
xlabel = 'n_shuttles'
xunit = datfile.su_S_North.unit
zdata = datfile.n_shuttle_gate_set.ndarray



q2_x90 = dict(y=xdata, x=zdata)


#%%
from mpl_toolkits.axes_grid1 import Divider, Size
CQ1 = '#0D77BC'
CQ2 = '#DD5E27'

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(4.5), ]
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
ax.set_title('equal time')


ax = ax11
# ax.scatter(q2_x90['xdata'], q2_x90['zdata'], s=3)
ax.plot(q2_x90['x'], q2_x90['y'], CQ2 , marker='o', markersize=3, linewidth=0.5)
ax.set_xlim(0, 80 )
ax.set_xticks([0, 40, 80]) 
ax.set_ylim(0,1)
ax.set_yticks([0, 1])    
secax_y = ax.secondary_xaxis('top', functions=( lambda x: x*2, lambda x: x/2))   
secax_y.set_xlabel('number of shuttles')
secax_y.set_xticks([0, 160, 320]) 



# filename = 'FigureSupp7e' # 'equal_time_repeat_x90'
# plt.savefig(filename+'.pdf')







