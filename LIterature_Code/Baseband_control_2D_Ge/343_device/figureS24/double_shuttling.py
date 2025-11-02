# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 17:02:42 2023

@author: vjohn
"""

# %% imports
import numpy as np
import matplotlib.pylab as plt

from array_343.figures.hopping_spins_paper.figureS24.QuantumDotArray import QuantumDotArray


from array_343.figures.hopping_spins_paper.figureS24.helpers import (shuttle_spin, time_vs_time_scan, time_scan_bloch,
                     plot_bloch_shuttle_seq, plot_time_vs_time)

# %% define quantum dot array
array = QuantumDotArray(2)

array.Q1.g_factor = 0.05
array.Q1.phi = 0
array.Q1.theta = 0

array.Q2.g_factor = 0.05
array.Q2.phi = -0.5*np.pi
array.Q2.theta = 0.3*np.pi

# %% Global input parameters

B = 20e-3
initial_position = (0, 0, 1)

shuttle_sequence = [1, 2, 1, 2, 1]

bloch_time_sweeps = [0, 2]

time_sweep_dot1 = ['Q1', 0, 300, 151]
time_sweep_dot2 = ['Q2', 0, 300, 151]

n_time_sweep_1 = 1
n_time_sweep_2 = 1
shuttle_back = False

array.Q1.shuttle_wait = array.Q2.lamor_period(B)*1e9/4
array.Q2.shuttle_wait = array.Q2.lamor_period(B)*1e9/4

# %%

z_values_2D = time_vs_time_scan(array, shuttle_sequence,
                                B, initial_position,
                                time_sweep_dot1, time_sweep_dot2,
                                n_trig_dot1=n_time_sweep_1,
                                n_trig_dot2=n_time_sweep_2,
                                back=shuttle_back)

z_values_2D_norm = (z_values_2D + 1) / 2

# %%

# Plotting the results as a colormap with normalized scale
fig = plt.figure(figsize=(10, 5))

ax0 = fig.add_subplot(121)
ax1 = fig.add_subplot(122, projection='3d')

plot_time_vs_time(z_values_2D_norm, time_sweep_dot1,
                  time_sweep_dot2, ax0, fig=fig,
                  show_colorbar=True)

plot_bloch_shuttle_seq(array, shuttle_sequence, initial_position,
                       bloch_time_sweeps, B, back=shuttle_back,
                       n_time_sweep_1=n_time_sweep_1,
                       n_time_sweep_2=n_time_sweep_2,
                       fig=fig, axes=ax1)

plt.show()

# %%

shuttle_sequence = [1, 2, 1, 2]
shuttle_back = False

bloch_time_sweeps = [0, 0]

time_sweep_dot1 = ['Q1', 0, 300, 151]
time_sweep_dot2 = ['Q2', 0, 300, 151]

bloch_point_size = 1
font_size = 12

fig = plt.figure(figsize=(2*3.9, 2*2.7))
# fig.suptitle(
#     f'array.Q3.phi = {array.Q2.phi/np.pi:.2f}$\pi$, shuttle_back = {shuttle_back}')

ax0 = fig.add_subplot(251, projection='3d')
ax1 = fig.add_subplot(252, projection='3d')
ax2 = fig.add_subplot(253, projection='3d')
ax3 = fig.add_subplot(254, projection='3d')

ax4 = fig.add_subplot(151, projection='3d')
ax5 = fig.add_subplot(152, projection='3d')
ax6 = fig.add_subplot(153, projection='3d')
ax7 = fig.add_subplot(154, projection='3d')

axs = [ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7]

bloch_spheres = []

for n, shuttle_wait in enumerate(np.linspace(0, array.Q1.lamor_period(B)*1e9, 8)):
    array.Q1.shuttle_wait = shuttle_wait

    z_values_2D = time_vs_time_scan(array, shuttle_sequence,
                                    B, initial_position,
                                    time_sweep_dot1, time_sweep_dot2,
                                    back=shuttle_back)

    bloch = plot_bloch_shuttle_seq(array, shuttle_sequence, initial_position,
                                   bloch_time_sweeps, B, fig=fig, axes=axs[n],
                                   back=shuttle_back, bloch_point_size=bloch_point_size,
                                   font_size=font_size,
                                   vector_width=1, view=[45, 10])
    bloch_spheres.append(bloch)
    bloch.view = [0, 0]
    bloch.render()

    axs[n].legend().set_visible(False)

plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.1,
                    hspace=0.8)

plt.savefig('double_shuttling_bloch_sphere.pdf', dpi=300)
plt.show()


# %%

plt.figure(figsize=(2.8, 2.3))

time_sweep_1D = ['Q1', 0, 3*array.Q1.lamor_period(B)*1e9, 211]

time = np.linspace(time_sweep_1D[1],
                   time_sweep_1D[2],
                   time_sweep_1D[3])

array.print_readable_snapshot()
result = time_scan_bloch(array, shuttle_sequence, B, initial_position,
                         time_sweep_1D)
result_z = result[:, 2]/2+1/2

plt.plot(time, result_z)
plt.xlabel(r'$t_{\mathrm{sweep}}$ (ns)')
plt.ylabel(r'$P_{\downarrow}$')
plt.ylim(0, 1)

for n in np.linspace(0, 70, 8):
    n = int(n)
    plt.scatter([time[n]], [result_z[n]], color='black', marker='x')
# plt.xlim(0, 200)

plt.tight_layout()
plt.savefig('double_shuttling_1D.pdf', dpi=300)
plt.show()
