# this file is modified from Supp/shuttle_model/q1q4_comparison_H4x4_vs_H2x2_vs_simple_loop.py
import matplotlib.pyplot as plt
import numpy as np
import pickle

#%%
import pickle
# results = pickle.load( open('FigSupp5_qutip_simulation_results_QuTip5_23.pickle', 'rb'))
results = pickle.load( open('FigSupp5_qutip_simulation.pickle', 'rb'))

pair = results['pair']
#%%


finalstate_q2_f = results['print_residual_False']['finalstate_q2_f']
finalstate_q2_reduction_f = results['print_residual_False']['finalstate_q2_reduction_f']
finalstate_q2_simple_f = results['print_residual_False']['finalstate_q2_simple_f']
waittimes = results['print_residual_False']['waittimes']

from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5),  Size.Fixed(1),Size.Fixed(1.5)]
v = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
ax31 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=3, ny=1))
ax51 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=5, ny=1))






ax = ax11
ax.plot(waittimes, finalstate_q2_f[0,:,0] ,label='Fit',color='C0', linewidth=0.5)
ax.plot(waittimes, finalstate_q2_reduction_f[0,:,0] ,color='C1', linewidth=0.5)
ax.plot(waittimes, finalstate_q2_simple_f[0,:,0] ,color='C2', linewidth=0.5)
ax.set_xlabel('waittime_' +  pair + ' (ns)')
ax.set_ylabel('P_down')
ax.set_xticks( np.linspace(np.min(waittimes), np.max(waittimes), 3) )

ax = ax31
# ax.plot(waittimes, finalstate_q2_simple_f[0,:,0] - finalstate_q2_f[0,:,0] ,label='Fit',color='C0', linewidth=0.5)
ax.plot(waittimes, finalstate_q2_reduction_f[0,:,0] - finalstate_q2_f[0,:,0] ,color='C0', linewidth=0.5)
ax.set_xlabel('waittime_' +  pair + ' (ns)')
ax.set_ylabel('delta P_down')
ax.set_xticks( np.linspace(np.min(waittimes), np.max(waittimes), 3) )

ax = ax51
ax.plot(waittimes, finalstate_q2_simple_f[0,:,0] - finalstate_q2_f[0,:,0] ,label='Fit',color='C0', linewidth=0.5)
# ax.plot(waittimes, finalstate_q2_reduction_f[0,:,0] - finalstate_q2_f[0,:,0] ,color='C1', linewidth=0.5)
ax.set_xlabel('waittime_' +  pair + ' (ns)')
ax.set_ylabel('delta P_down')
ax.set_xticks( np.linspace(np.min(waittimes), np.max(waittimes), 3) )




# filename = 'Supp_shuttle_model_waittimes_q1'
# plt.savefig(filename+'.pdf')
    
#%%
finalstate_q2 = results['print_residual_True']['finalstate_q2']
finalstate_q2_reduction = results['print_residual_True']['finalstate_q2_reduction']
finalstate_q2_simple = results['print_residual_True']['finalstate_q2_simple']
ts = results['print_residual_True']['ts']
from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5),  Size.Fixed(1),Size.Fixed(1.5)]
v = [Size.Fixed(0.5), Size.Fixed(1.5), Size.Fixed(1), Size.Fixed(1.5)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
ax31 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=3, ny=1))
ax51 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=5, ny=1))






ax = ax11
ax.plot(ts, finalstate_q2[0,0,:,0]*(-2)+1 ,label='Fit',color='C0', linewidth=0.5)
ax.plot(ts, finalstate_q2_reduction[0,0,:,0]*(-2)+1 ,color='C1', linewidth=0.5)
ax.plot(ts, finalstate_q2_simple[0,0,:,0]*(-2)+1 ,color='C2', linewidth=0.5)
ax.set_xlabel('waittime_' +  pair + ' (ns)')
ax.set_ylabel('instantaneous sigma_z')
ax.set_xticks( np.linspace(np.min(ts), np.max(ts), 3) )

ax = ax31
ax.plot(ts, finalstate_q2_reduction[0,0,:,0]*(-2) - finalstate_q2[0,0,:,0]*(-2) ,color='C0', linewidth=0.5)
ax.set_xlabel('waittime_' +  pair + ' (ns)')
ax.set_ylabel('delta sigma_z')
ax.set_xticks( np.linspace(np.min(ts), np.max(ts), 3) )

ax = ax51
ax.plot(ts, finalstate_q2_simple[0,0,:,0]*(-2) - finalstate_q2[0,0,:,0]*(-2) ,label='Fit',color='C0', linewidth=0.5)
ax.set_xlabel('waittime_' +  pair + ' (ns)')
ax.set_ylabel('delta sigma_z')
ax.set_xticks( np.linspace(np.min(ts), np.max(ts), 3) )

plt.show()

# filename = 'Supp_shuttle_model_ts_q1'
# plt.savefig(filename+'.pdf')










