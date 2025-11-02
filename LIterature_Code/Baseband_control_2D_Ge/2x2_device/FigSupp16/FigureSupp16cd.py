# this file is modified from Supp\CZ\CZ_leakage\Hamming_tramp_vpB12_v1.py
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigSupp16','')
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


#%%
datadir = os.getcwd()



datfiles = []
start_times = [     '2023-08-18\\22-47-23',  # init state qA up qB down, CZ^8
              ]
for start_time in start_times:
    end_time = start_time
    datfiles += [ *get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) [0]  ]

datfiles_hamming = list(datfiles)


#%%

cz_hamming = np.array([ [43.47, -64], 
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



#%%

zdata_min = []
zdata_max = []
for datfiles in [datfiles_hamming, datfiles_hamming] :
    for i_dat, datfile in enumerate(datfiles):
        for i_z, zname in enumerate(['su00', 'su01', 'su10', 'su11']):
            _, _, zdata = data2d_attr(datfile,  ['ndarray'], zname, ['ndarray'] )
            zdata_min += [np.min(zdata)]
            zdata_max += [np.max(zdata)]
zmin = min(zdata_min)
zmax = max(zdata_max) 

CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
# MRSIZE = 2.5
SCSIZE = 5 #MRSIZE**2
LW = 0.5
from mpl_toolkits.axes_grid1 import Divider, Size
HSPACE = 0.3
VSPACE = 0.3
fig = plt.figure( figsize=(40,8))
# h = [Size.Fixed(HSP), Size.Fixed(1.5), Size.Fixed(HSP), Size.Fixed(1.5), Size.Fixed(HSP), Size.Fixed(1.5), Size.Fixed(HSP), Size.Fixed(1.5)]
# v = [Size.Fixed(1), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5), ]
h = []
for _ in range(4):
    h += [Size.Fixed(HSPACE), Size.Fixed(1.3)]
v = []
for i in range(6):
    if i in [0,3]:
        v += [Size.Fixed(VSPACE), Size.Fixed(0.9)]
    else:
        v += [Size.Fixed(VSPACE), Size.Fixed(1.3)]


divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
axs = np.zeros((4,4),dtype=object)
for i in range(axs.shape[0]):
    for j in range(axs.shape[1]):
        
        axs[axs.shape[0]-i-1,j] = fig.add_axes(divider.get_position(),
                           axes_locator=divider.new_locator(nx=(1+2*j), ny=1+2*i))


plt.rcParams.update({'font.size':8})

# fig, axs = plt.subplots(6, 4, figsize=(12,18)) #, sharex=True, sharey=True)
# plt.subplots_adjust( hspace=0.3, wspace=0.3)


        
datfiles = datfiles_hamming
cz_params = cz_hamming
for i_dat, datfile in enumerate(datfiles):
    for i_z, zname in enumerate(['su00', 'su01', 'su10', 'su11']):
        xdata, ydata, zdata = data2d_attr(datfile,  ['ndarray'], zname, ['ndarray'] )
        xlabel, ylabel, _ = data2d_attr(datfile,  ['name'],  )

        # ax = axs[i_dat*3-3+1, i_z]
        ax = axs[-2, i_z]
        img = ax.pcolormesh(xdata, ydata, zdata, vmin=zmin, vmax=zmax, rasterized=True, )
        # if i_dat != 3:
        #     ax.set_xticks([])
        ax.scatter(cz_params[:,0], cz_params[:,1],  s=SCSIZE, c='C1')
        # ax.plot( cz_params[:,1], 0.5*1e9/( 0.24e6*np.exp(-0.059*cz_params[:,1])*1.08),   'r' , linewidth=0.5 )
        ax.set_yticks([-85, -75, -65])
        ax.set_xticks([12, 22, 32, 42, 52])

cb=fig.colorbar(img, location='top',ticks=[0.1, 0.8], shrink=0.03, aspect=5,anchor=(0.5,0) )
        

# axs[0,4].set_title(f'{ydata[idx]:.3g}mV')



colorA = 'C2'
idx = 6
# axs[0,5].set_title(f'{ydata[idx]:.3g}mV')
yt = [[0.1, 0.2], [0, 0.2], [0.65, 0.85], [0, 0.1]]
yl = [[0.07,0.23], [0,0.23], [0.6, 0.9], [0.0,0.12]]
for datfiles in [datfiles_hamming, ] :
    for i_dat, datfile in enumerate(datfiles):
        for i_z, zname in enumerate(['su00', 'su01', 'su10', 'su11']):
            xdata, _, zdata = data2d_attr(datfile,  ['ndarray'], zname, ['ndarray'] )
            ax = axs[3, i_z]
            ax.plot( xdata,  zdata[idx,:], colorA, linewidth=LW)
            ax.set_xticks([12, 22, 32, 42, 52])
            ax.set_yticks(yt[i_z])
            ax.set_ylim(*yl[i_z])
            ax.set_xlim([12, 52])
            
            ax = axs[-2, i_z]
            ax.plot( (47,52),(ydata[idx],)*2, colorA  )

colorB = 'C0'            
idx = 5
for datfiles in [ datfiles_hamming] :
    for i_dat, datfile in enumerate(datfiles):
        for i_z, zname in enumerate(['su00', 'su01', 'su10', 'su11']):
            xdata, _, zdata = data2d_attr(datfile,  ['ndarray'], zname, ['ndarray'] )
            ax = axs[3, i_z]
            ax.plot( xdata,  zdata[idx,:], colorB, linewidth=LW)
            ax.set_xticks([12, 22, 32, 42, 52])
            ax.set_xlim([12, 52])
            
            ax = axs[-2, i_z]
            ax.plot( (47,52),(ydata[idx],)*2, colorB  )            

colorC = 'r'
idx = 4
# axs[0,5].set_title(f'{ydata[idx]:.3g}mV')
# yt = [[0, 0.4], [0, 0.3], [0.8, 0.4], [0, 0.2]]
# yl = [[0,0.45], [0,0.4], [0.28, 0.9], [0,0.25]]
for datfiles in [datfiles_hamming, ] :
    for i_dat, datfile in enumerate(datfiles):
        for i_z, zname in enumerate(['su00', 'su01', 'su10', 'su11']):
            xdata, _, zdata = data2d_attr(datfile,  ['ndarray'], zname, ['ndarray'] )
            ax = axs[3, i_z]
            ax.plot( xdata,  zdata[idx,:], colorC, linewidth=LW)
            ax.set_xticks([12, 22, 32, 42, 52])
            ax.set_yticks(yt[i_z])
            ax.set_ylim(*yl[i_z])
            ax.set_xlim([12, 52])
            
            ax = axs[-2, i_z]
            ax.plot( (47,52),(ydata[idx],)*2, colorC  )



# filename = 'FigureSupp16cd' # 'hamming_-73mV_-77mV_init01'
# plt.savefig(filename+'.pdf')


