# this file is modified from Fig2/readout_fromGST/template_2d_histogram.py
import numpy as np
import  matplotlib.pyplot as plt
##%%
import pickle

hist_dict = {'min_vals': [120, 100], 
              'max_vals': [210, 194], 
              'thresholds': [164, 145  ], 
              'Nbins': [41, 41], 
              'invert_result': [False, True]  }
        
        
x_bins = np.linspace(hist_dict['min_vals'][0], hist_dict['max_vals'][0], hist_dict['Nbins'][0])
y_bins = np.linspace(hist_dict['min_vals'][1], hist_dict['max_vals'][1], hist_dict['Nbins'][1])


idxlist = [1, 4, 7, 16]
datapkl = pickle.load(  open( 'dataset_XYICphase_max_max_lengths_82023-08-20_15-35-39.pickle',  'rb')  )

from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.0), Size.Fixed(0.7), Size.Fixed(1.0)]
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

axs = [ax11, ax13, ax31, ax33  ]

for i in range(4):
    idx = idxlist[i]
    data = datapkl[idx]['shots'][0]
    read0 = data[0] 
    read1 = data[1]  
    
    
    print(data[2:])
    
    
    ax = axs[i]
    h, xedges, yedges, img = ax.hist2d(read0, read1, bins =[x_bins, y_bins], cmap='Blues', vmin=0, vmax=21)
    print(np.max(h))

    ax.set_xticks([120, 160, 200])
    ax.set_yticks([100, 140, 180])

    if i == 0:
        cb=plt.colorbar(img, ax=ax, location='top',ticks=[0, 15], shrink=0.05,aspect=4,anchor=(0.61,-2.55))
        cb.ax.tick_params(labelsize=12) 
# plt.clim(0,40)
# filename = 'Figure2b'  #'Fig_2b_histogram'
# plt.savefig(filename+'.pdf')




#%%

read0s = []
read1s = []
for i in range(4):
    idx = idxlist[i]
    data = datapkl[idx]['shots'][0]
    read0s += [data[0]]
    read1s += [data[1]]  
read0s = np.hstack(read0s)
read1s = np.hstack(read1s)

#%%
vth0 = 164
vth1 = 145

import random
colorlist = ['C0', ]*500 + ['C1', ]*500 + ['C2', ]*500 + ['C3', ]*500 
idxlist = list(range(len(colorlist)))
random.shuffle( idxlist )

read0s_ = [read0s[i] for i in idxlist]
read1s_ = [read1s[i] for i in idxlist]
colorlist_ = [colorlist[i] for i in idxlist]



plotfirst = []
for i in range(len(read0s_)):
    condition00 = colorlist_[i] == 'C0' and read0s_[i] < vth0  and read1s_[i] > vth1
    condition01 = colorlist_[i] == 'C1' and read0s_[i] < vth0  and read1s_[i] < vth1
    condition10 = colorlist_[i] == 'C2' and read0s_[i] > vth0  and read1s_[i] > vth1
    condition11 = colorlist_[i] == 'C3' and read0s_[i] > vth0  and read1s_[i] < vth1
    if  condition00 or condition01 or condition10 or condition11:
        plotfirst += [i]


plotsecond = list(set(range(len(read0s_))) - set(plotfirst))


fig = plt.figure( figsize=(6,6))    
plt.scatter(np.array(read0s_)[plotfirst], np.array(read1s_)[plotfirst], c=np.array(colorlist_)[plotfirst], s=20)

plt.scatter(np.array(read0s_)[plotsecond], np.array(read1s_)[plotsecond], c=np.array(colorlist_)[plotsecond], s=20)








