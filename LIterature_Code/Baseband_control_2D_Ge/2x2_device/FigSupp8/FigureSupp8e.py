# this file is modified from Supp/unequal_time/plotAllXY.py
import os
import matplotlib.pyplot as plt
import numpy as np


from projects.notebook_tools.notebook_tools import get_data_from, fit_data

datadir = os.getcwd()
start_time = '2023-08-16\\22-41-30'

end_time = start_time
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

data = datfiles[0]



MKRSIZE = 5

CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(10,8))
h = [Size.Fixed(0.5), Size.Fixed(1.5*2), Size.Fixed(0.7), Size.Fixed(1.5*2), ]
v = [Size.Fixed(0.5), Size.Fixed(1.0*2), Size.Fixed(0.7), Size.Fixed(1.0*2)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
ax13 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=3))
ax31 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=3, ny=1))
ax33 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=3, ny=3))

ax = ax13
ax.plot(   data.gate_sequence_index_set , data.arrays['1_counts'] , markersize=MKRSIZE, c=CQ1, marker='o' , linewidth=0.5)
ymax, ymin = np.average(data.arrays['1_counts'][0:5]), np.average(data.arrays['1_counts'][17:21])
yideal = np.array([1,]*5+[0.5,]*12+[0,]*4) * (ymax-ymin) + ymin
ax.plot( data.gate_sequence_index_set, yideal, 'k', linewidth=0.5)

ax.set_title('unequal time')
ax.set_ylim(0,1)
ax.set_yticks([0,0.5,1])

seqs = [ 'II', r'$\rm X^2X^2$', r'$\rm Y^2Y^2$', r'$\rm X^2Y^2$', r'$\rm Y^2X^2$', 'XI', 'YI', 'XY', 'YX'
        , r'$\rm XY^2$', r'$\rm YX^2$', r'$\rm X^2Y$', r'$\rm Y^2X$', r'$\rm XX^2$', r'$\rm X^2X$'
        , r'$\rm YY^2$', r'$\rm Y^2Y$', r'$\rm X^2I$', r'$\rm Y^2I$', 'XX', 'YY']
ax.set_xticks(range(21), seqs, rotation='vertical', fontsize=6)   






start_time = '2023-08-16\\22-42-48'
end_time = start_time
datfiles, fnames = get_data_from(start_time, end_time, rootfolder=datadir, only_complete = False) 

data = datfiles[0]




ax = ax11
ax.plot(   data.gate_sequence_index_set , data.arrays['1_counts'] , markersize=MKRSIZE, c=CQ2, marker='o' , linewidth=0.5)
ymax, ymin = np.average(data.arrays['1_counts'][0:5]), np.average(data.arrays['1_counts'][17:21])
yideal = np.array([1,]*5+[0.5,]*12+[0,]*4) * (ymax-ymin) + ymin
ax.plot( data.gate_sequence_index_set, yideal, 'k', linewidth=0.5)

ax.set_title('unequal time')
ax.set_ylim(0,1)
ax.set_yticks([0,0.5,1])

seqs = [ 'II', r'$\rm X^2X^2$', r'$\rm Y^2Y^2$', r'$\rm X^2Y^2$', r'$\rm Y^2X^2$', 'XI', 'YI', 'XY', 'YX'
        , r'$\rm XY^2$', r'$\rm YX^2$', r'$\rm X^2Y$', r'$\rm Y^2X$', r'$\rm XX^2$', r'$\rm X^2X$'
        , r'$\rm YY^2$', r'$\rm Y^2Y$', r'$\rm X^2I$', r'$\rm Y^2I$', 'XX', 'YY']
ax.set_xticks(range(21), seqs, rotation='vertical', fontsize=6)   







 

# filename = 'FigureSupp8e' # 'Supp_allxy_unequal'
# plt.savefig(filename+'.pdf')













