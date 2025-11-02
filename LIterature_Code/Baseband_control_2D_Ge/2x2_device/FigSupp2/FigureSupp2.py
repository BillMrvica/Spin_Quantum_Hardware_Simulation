# this file is mofidied from Supp/fit_tc_angle/Supp_CSD_shuttle_gates.py
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\FigureSupp2','')
#%%

import sys
import numpy as np
import matplotlib.pyplot as plt
import qcodes
# from qcodes.data.data_set import load_data
from qcodes_loop.data.data_set import load_data
from projects.notebook_tools.notebook_tools import get_data_from, fit_data
# from qcodes.data import data_set
from qcodes_loop.data import data_set
sys.path.insert(1,path_notebook)
# import notebook_tools
# from notebook_tools import get_data_from
from matplotlib import colors
# from notebook_tools import fit_data
from mpl_toolkits.axes_grid1 import Divider, Size

#%%

W = 2.68
FSx = 4
FSy = 4 


MKR = 5

_1100_s = {
  'vP1': 30,
  'vP2': 45,
  'vP3': 46,
  'vP4': 65, }
_1010_s = {
  'vP1': 30,
  'vP2': 64, 
  'vP3': 27, 
  'vP4': 65, }
_0101_s = {
  'vP1': 58.5,
  'vP2': 45,
  'vP3': 46,
  'vP4': 39.6 }

collect_zdata = []

#%%

start_time = '2023-07-10\\' + '17-05-12' 

end_time = start_time 
datadir =  path 
datfiles, fnames = get_data_from(start_time, end_time, num = 1, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]


x = datfile.vP1_set.ndarray
y = datfile.vP4_set.ndarray
z = datfile.ch3.ndarray
collect_zdata += [ z.reshape(-1,)]

fig=plt.figure( figsize=(FSx,FSy))
plt.clf()
h = [Size.Fixed(1), Size.Fixed(W)]
v = [Size.Fixed(0.5), Size.Fixed(W)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
img = plt.pcolormesh( x, y, z, shading='auto', cmap='Blues',rasterized=True)


plt.plot( _1100_s['vP1'], _1100_s['vP4'], 'ko', markersize=MKR)
plt.plot( _0101_s['vP1'], _0101_s['vP4'], 'ko', markersize=MKR)
# plt.plot( _0110_s['vP1'], _0110_s['vP4'], 'ko', markersize=MKR)
plt.plot(  (_1100_s['vP1'], _0101_s['vP1']) , (_1100_s['vP4'],_0101_s['vP4'])   )

cb=plt.colorbar(location='top',ticks=[-200,-120],shrink=0.2,aspect=8,anchor=(0.9, -0.25))


cb.set_label(label=r'$V_{\rm sensor}$ (mV)',position=(-1,1),labelpad=-25, fontsize=12, rotation=0)

plt.xlabel(r'$\rm vP_1$ (mV)')
plt.ylabel(r'$\rm vP_4$ (mV)')


#%%

start_time = '2023-07-10\\' + '17-08-32'
end_time = start_time 
datadir =  path
datfiles, fnames = get_data_from(start_time, end_time, num = 1, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]

x = datfile.vP2_set.ndarray
y = datfile.vP3_set.ndarray
z = datfile.ch3.ndarray
collect_zdata += [ z.reshape(-1,)]

fig=plt.figure( figsize=(FSx,FSy))
plt.clf()

h = [Size.Fixed(1), Size.Fixed(W)]
v = [Size.Fixed(0.5), Size.Fixed(W)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=1))
img = plt.pcolormesh( x, y[0,:], z, shading='auto', cmap='Blues',rasterized=True)


plt.plot( _1100_s['vP2'], _1100_s['vP3'], 'ko', markersize=MKR)
plt.plot( _1010_s['vP2'], _1010_s['vP3'], 'ko', markersize=MKR)
# plt.show()
plt.plot(  (_1100_s['vP2'], _1010_s['vP2']) , (_1100_s['vP3'],_1010_s['vP3'])   )

# cb=plt.colorbar(location='top',ticks=[-260,-60],shrink=0.2,aspect=8,anchor=(0.9, -0.25))
cb=plt.colorbar(location='top' ,shrink=0.2,aspect=8,anchor=(0.9, -0.25))
# cb.set_label(label=r'$V_{\rm sensor}$ (mV)',position=(-1,1),labelpad=-25, fontsize=12, rotation=0)

plt.xlabel(r'$\rm vP_2$ (mV)')
plt.ylabel(r'$\rm vP_3$ (mV)')

#%%


start_time = '2023-07-10\\' + '17-05-37'
end_time = start_time 
datadir =  path
datfiles, fnames = get_data_from(start_time, end_time, num = 1, rootfolder=datadir, only_complete = False) 
datfile = datfiles[0]

x = datfile.vP1_set.ndarray
y = datfile.vP2_set.ndarray
z = datfile.ch3.ndarray
collect_zdata += [ z.reshape(-1,)]

fig=plt.figure( figsize=(FSx,FSy))
plt.clf()

h = [Size.Fixed(1), Size.Fixed(W)]
v = [Size.Fixed(0.5), Size.Fixed(W)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=1))
img = plt.pcolormesh( x, y, z, shading='auto', cmap='Blues',rasterized=True)

# ax.plot(x, z[:,0])

plt.plot( _1100_s['vP1'], _1100_s['vP2'], 'ko', markersize=MKR)
plt.plot( _0101_s['vP1'], _0101_s['vP2'], 'ko', markersize=MKR)
plt.plot( _1010_s['vP1'], _1010_s['vP2'], 'ko', markersize=MKR)

# plt.show()


# cb=plt.colorbar(location='top',ticks=[-260,-60],shrink=0.2,aspect=8,anchor=(0.9, -0.25))
cb=plt.colorbar(location='top' ,shrink=0.2,aspect=8,anchor=(0.9, -0.25))
# cb.set_label(label=r'$V_{\rm sensor}$ (mV)',position=(-1,1),labelpad=-25, fontsize=12, rotation=0)

plt.xlabel(r'$\rm vP_1$ (mV)')
plt.ylabel(r'$\rm vP_2$ (mV)')




#%%
collect_zdata_ = np.hstack(collect_zdata)

zmin = np.min(collect_zdata_)
zmax = np.max(collect_zdata_)

#%%
CMAP = 'viridis'
CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5)]
v = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(0.7), Size.Fixed(1.5)]
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




ax = ax11

start_time = '2023-07-10\\' + '17-05-37'
datfile = get_data_from(start_time, start_time, num = 1, rootfolder=path, only_complete = False) [0][0]
x = datfile.vP1_set.ndarray
y = datfile.vP2_set.ndarray
z = datfile.ch3.ndarray

img = ax.pcolormesh( x, y, z, shading='auto', cmap=CMAP,rasterized=True, vmin=zmin, vmax=zmax)
ax.plot( _1100_s['vP1'], _1100_s['vP2'], 'ko', markersize=MKR)
ax.plot( _0101_s['vP1'], _0101_s['vP2'], 'k^', markersize=MKR)
ax.plot( _1010_s['vP1'], _1010_s['vP2'], 'kv', markersize=MKR)
cb=plt.colorbar(img, ax=ax, location='top',ticks=[180, 240],shrink=0.05,aspect=4,anchor=(0.9, -0.25))
# cb.set_label(label=r'$V_{\rm sensor}$ (mV)',position=(-1,1),labelpad=-25, fontsize=12, rotation=0)
ax.set_xlabel(r'$\rm vP_1$ (mV)')
ax.set_ylabel(r'$\rm vP_2$ (mV)')
# ax.set_xlim(np.min(x), np.max(x))
# ax.set_ylim(np.min(y), np.max(y))
ax.set_xticks( np.linspace(np.min(x), np.max(x), 3) )
ax.set_yticks( np.linspace(np.min(y), np.max(y), 3)  )


ax = ax31

start_time = '2023-07-10\\' + '17-05-12' 
datfile = get_data_from(start_time, start_time, num = 1, rootfolder=path, only_complete = False) [0][0]
x = datfile.vP1_set.ndarray
y = datfile.vP4_set.ndarray
z = datfile.ch3.ndarray

img = ax.pcolormesh( x, y, z, shading='auto', cmap=CMAP,rasterized=True, vmin=zmin, vmax=zmax )
ax.plot( _1100_s['vP1'], _1100_s['vP4'], 'ko', markersize=MKR)
ax.plot( _0101_s['vP1'], _0101_s['vP4'], 'k^',  markersize=MKR)
# cb=plt.colorbar(img, ax=ax, location='top',ticks=[-200,-120],shrink=0.2,aspect=8,anchor=(0.9, -0.25))
# cb.set_label(label=r'$V_{\rm sensor}$ (mV)',position=(-1,1),labelpad=-25, fontsize=12, rotation=0)
ax.set_xlabel(r'$\rm vP_1$ (mV)')
ax.set_ylabel(r'$\rm vP_4$ (mV)')
ax.set_xticks( np.linspace(np.min(x), np.max(x), 3) )
ax.set_yticks( np.linspace(np.min(y), np.max(y), 3)  )


ax = ax51

start_time = '2023-07-10\\' + '17-08-32'
datfile = get_data_from(start_time, start_time, num = 1, rootfolder=path, only_complete = False) [0][0]
x = datfile.vP2_set.ndarray
y = datfile.vP3_set.ndarray
z = datfile.ch3.ndarray

img = ax.pcolormesh( x, y[0,:],  z, shading='auto', cmap=CMAP,rasterized=True, vmin=zmin, vmax=zmax)
ax.plot( _1100_s['vP2'], _1100_s['vP3'], 'ko', markersize=MKR)
ax.plot( _1010_s['vP2'], _1010_s['vP3'], 'kv',  markersize=MKR)
# cb=plt.colorbar(img, ax=ax, location='top',ticks=[-200,-120],shrink=0.2,aspect=8,anchor=(0.9, -0.25))
# cb.set_label(label=r'$V_{\rm sensor}$ (mV)',position=(-1,1),labelpad=-25, fontsize=12, rotation=0)
ax.set_xlabel(r'$\rm vP_2$ (mV)')
ax.set_ylabel(r'$\rm vP_3$ (mV)')
ax.set_xticks( np.linspace(np.min(x), np.max(x), 3) )
ax.set_yticks( np.linspace(np.min(y), np.max(y), 3)  )



 

# filename = 'FigureSupp2' # 'Supp_CSD'
# plt.savefig(filename+'.pdf')













