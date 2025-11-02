# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 09:15:18 2023

@author: vjohn
"""

from mpl_toolkits.axes_grid1 import make_axes_locatable
#from QuantumDotArray import QuantumDotArray
from array_343.figures.hopping_spins_paper.figureS24.QuantumDotArray import QuantumDotArray
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import math

from qutip import Bloch, Qobj
import matplotlib.lines as mlines
import matplotlib.colors as mcolors

# %%
# Constants
mu_B = 9.274009994e-24  # Bohr magneton in J/T
hbar = 1.054571817e-34  # Reduced Planck's constant in J*s


def normalize(vector):
    """Normalize a vector to unit length."""
    norm = np.linalg.norm(vector)
    if norm == 0:
        return vector
    return vector / norm


def vector_norm(theta, phi, r=1):
    """
    Convert spherical coordinates to cartesian coordinates.

    :param r: radius
    :param theta: azimuthal angle (angle in xy-plane from the x-axis), in radians
    :param phi: polar (or zenith) angle (angle from the z-axis), in radians
    :return: tuple containing cartesian coordinates (x, y, z)
    """
    x = r * math.sin(phi) * math.cos(theta)
    y = r * math.sin(phi) * math.sin(theta)
    z = r * math.cos(phi)

    return x, y, z


def larmor_precession(g, B, axis, dt, unit='ns'):
    """Compute the Larmor precession rotation matrix for a given g-factor, magnetic field, quantization axis, and time."""
    if unit == 'ns':
        dt = dt*1e-9
    elif unit == 's':
        dt = dt
    else:
        raise ValueError("unit can only be 'ns' or 's'.")
    theta = g * mu_B * B / hbar * dt
    n_x, n_y, n_z = axis
    rotation_matrix = np.array([
        [np.cos(theta) + n_x**2 * (1 - np.cos(theta)), n_x*n_y*(1 - np.cos(theta)) -
         n_z*np.sin(theta), n_x*n_z*(1 - np.cos(theta)) + n_y*np.sin(theta)],
        [n_y*n_x*(1 - np.cos(theta)) + n_z*np.sin(theta), np.cos(theta) + n_y **
         2 * (1 - np.cos(theta)), n_y*n_z*(1 - np.cos(theta)) - n_x*np.sin(theta)],
        [n_z*n_x*(1 - np.cos(theta)) - n_y*np.sin(theta), n_z*n_y*(1 - np.cos(theta)
                                                                   ) + n_x*np.sin(theta), np.cos(theta) + n_z**2 * (1 - np.cos(theta))]
    ])
    return [rotation_matrix] if np.isscalar(dt) else [rotation_matrix] * len(dt)


def shuttle_spin(B, g_factors, quantization_axes, times, initial_position):
    """Modified function to simulate the shuttling of a spin through an array of quantum dots."""
    position = np.array(initial_position)
    final_positions = []

    # Iterate through each quantum dot
    for idx, (g, axis, dt) in enumerate(zip(g_factors, quantization_axes, times)):
        if isinstance(dt, np.ndarray):  # For dots with varying time
            for t in dt:
                rotation_matrix = larmor_precession(g, B, axis, t)[0]
                new_position = np.dot(rotation_matrix, position)
                final_positions.append(new_position)
        else:
            rotation_matrix = larmor_precession(g, B, axis, dt)[0]
            position = np.dot(rotation_matrix, position)
            final_positions.append(position)

    return final_positions

# %% QuantumDotArray functions


def time_scan_bloch(array: QuantumDotArray,
                    shuttle_sequence, B, initial_position,
                    time_sweep_dot):
    waiting_times = []
    sweep_integer = None
    for int_sweep, dot_number in enumerate(shuttle_sequence):
        dot_name = f'Q{dot_number}'
        waiting_times.append(getattr(getattr(array, dot_name), 'shuttle_wait'))
        if dot_name == time_sweep_dot[0]:
            sweep_integer = int_sweep

    # Simulation for varying times in dots 2 and 3
    time_range_dot = np.linspace(time_sweep_dot[1],
                                 time_sweep_dot[2],
                                 time_sweep_dot[3])

    vectors = []

    for i, t2 in enumerate(time_range_dot):
        times_current = waiting_times.copy()
        times_current
        times_current[sweep_integer] = t2
        final_position_current = shuttle_spin(B,
                                              array.g_factors(
                                                  shuttle_sequence),
                                              array.axes(shuttle_sequence),
                                              times_current,
                                              initial_position)
        vectors.append(final_position_current)
    return np.array(vectors)[:, -1]


def time_vs_time_scan(array: QuantumDotArray,
                      shuttle_sequence, B, initial_position,
                      time_sweep_dot1, time_sweep_dot2, back=False,
                      n_trig_dot1=1, n_trig_dot2=1):
    waiting_times = []
    sweep_integers = [None, None]
    n_sweep_dot1 = 1
    n_sweep_dot2 = 1

    if back:
        full_shuttle_seq = shuttle_sequence[1:] + shuttle_sequence[::-1][1:]
    else:
        full_shuttle_seq = shuttle_sequence[1:]

    for int_sweep, dot_number in enumerate(full_shuttle_seq):
        dot_name = f'Q{dot_number}'
        waiting_times.append(getattr(getattr(array, dot_name), 'shuttle_wait'))
        if dot_name == time_sweep_dot1[0]:
            if n_sweep_dot1 == n_trig_dot1:
                sweep_integers[0] = int_sweep
            n_sweep_dot1 += 1
        elif dot_name == time_sweep_dot2[0]:
            if n_sweep_dot2 == n_trig_dot2:
                sweep_integers[1] = int_sweep
            n_sweep_dot2 += 1

    # Simulation for varying times in dots 2 and 3
    time_range_dot1 = np.linspace(time_sweep_dot1[1],
                                  time_sweep_dot1[2],
                                  time_sweep_dot1[3])
    time_range_dot2 = np.linspace(time_sweep_dot2[1],
                                  time_sweep_dot2[2],
                                  time_sweep_dot2[3])

    z_values_2D = np.zeros((len(time_range_dot2), len(time_range_dot1)))

    for i, t2 in enumerate(time_range_dot1):
        for j, t3 in enumerate(time_range_dot2):
            times_current = waiting_times.copy()
            times_current[sweep_integers[0]] = t2
            times_current[sweep_integers[1]] = t3
            final_position_current = shuttle_spin(B,
                                                  array.g_factors(
                                                      full_shuttle_seq),
                                                  array.axes(full_shuttle_seq),
                                                  times_current,
                                                  initial_position)
            z_values_2D[j, i] = final_position_current[-1][2]
    return z_values_2D


def process_dot(dot, counter, trigger_count, wait_time, B, array,
                time_step=None):
    if counter == trigger_count:
        wait_time += getattr(getattr(array,
                             f'Q{dot}'), 'lamor_period')(B) * 1e9
        if time_step is None:
            time_step = int(wait_time/1)
        return [f'Q{dot}', 0, wait_time, int(time_step)]
    else:
        if time_step is None:
            time_step = int(wait_time/1)
        return [f'Q{dot}', 0, wait_time, int(time_step)]


def map_integers_to_colors(integers):
    """Map a list of integers to a list of colors."""
    # Find unique integers and sort them for consistent color mapping
    unique_integers = sorted(set(integers))

    # Get a color palette with enough colors for each unique integer
    colors = list(mcolors.TABLEAU_COLORS.values())

    # Create a dictionary that maps integers to colors
    color_map = dict(zip(unique_integers, colors))

    # Map the integers in the original list to colors
    colors_for_integers = [color_map[i] for i in integers]

    return colors_for_integers


def plot_bloch_shuttle_seq(array, shuttle_sequence, initial_position,
                           time_sweeps, B, back: bool, fig=None, axes=None,
                           n_time_sweep_1=1, n_time_sweep_2=1,
                           time_sweep_back=False,
                           bloch_point_size=10, font_size=15,
                           vector_width=1, view=[-60, 30]):
    if back:
        full_shuttle_seq = shuttle_sequence + shuttle_sequence[::-1][1:]
    else:
        full_shuttle_seq = shuttle_sequence

    b = Bloch()
    b.axes = axes
    b.fig = fig
    b.font_size = font_size
    b.point_marker = 'o'
    b.sphere_alpha = 0.05
    b.frame_alpha = 0
    b.frame_width = 0.1
    b.vector_width = vector_width
    b.point_size = list(bloch_point_size*np.ones(len(full_shuttle_seq)))
    b.zlabel = [r'$|\!\downarrow\!\rangle$',
                r'$|\!\uparrow\!\rangle$']
    b.xlabel = ['', '']
    b.ylabel = ['', '']

    colors = map_integers_to_colors(full_shuttle_seq)
    b.point_color = colors[1:]
    b.vector_color = colors

    handles = []

    start_dot = shuttle_sequence[0]
    quant_axis = getattr(getattr(array, f'Q{start_dot}'), 'axis')
    b.add_vectors(quant_axis)

    color = b.vector_color[0]
    line = mlines.Line2D([], [], color=color,
                         markersize=15, label=f'Dot {start_dot}')
    handles.append(line)

    counters = {time_sweeps[0]: 0, time_sweeps[1]: 0}
    trigger_counts = {
        time_sweeps[0]: n_time_sweep_1, time_sweeps[1]: n_time_sweep_2}

    for n, dot in enumerate(full_shuttle_seq[1:]):
        wait_time = getattr(getattr(array, f'Q{dot}'), 'shuttle_wait')

        if dot in time_sweeps:
            counters[dot] += 1
            time_sweep = process_dot(dot, counters[dot], trigger_counts[dot],
                                     wait_time, B, array)
        else:
            time_step = int(wait_time/1)+1
            time_sweep = [f'Q{dot}', 0, wait_time, time_step]

        bloch_time_evolution = time_scan_bloch(array, full_shuttle_seq[1:][:n+1],
                                               B, initial_position,
                                               time_sweep)
        points = [bloch_time_evolution[:, 0],
                  bloch_time_evolution[:, 1],
                  bloch_time_evolution[:, 2]]
        b.add_points(points, meth='s')

    for n, dot in enumerate(shuttle_sequence[1:]):
        quant_axis = getattr(getattr(array, f'Q{dot}'), 'axis')
        b.add_vectors(quant_axis)

        # Create a dummy line with the same color for legend
        color = b.vector_color[n+1]
        line = mlines.Line2D([], [], color=color,
                             markersize=15, label=f'Dot {dot}')
        handles.append(line)

    # b.view = view
    b.render()

    b.axes.view_init(elev=view[1], azim=view[0])

    # Add legend to the Bloch sphere using the underlying Axes
    b.axes.legend(handles=handles, loc='upper right')

    return b


def plot_time_vs_time(z_values_2D, time_sweep_dot1, time_sweep_dot2,
                      ax, fig=None, show_colorbar=True):
    pcm = ax.imshow(z_values_2D, aspect='auto', origin='lower',
                    extent=[0, time_sweep_dot1[-2],
                            0, time_sweep_dot2[-2]],
                    cmap='viridis', vmin=0, vmax=1)

    if show_colorbar and fig is not None:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(pcm, label='spin down probability', boundaries=np.linspace(0, 1, 11),
                     cax=cax, orientation='vertical')
    ax.set_xlabel(f'Time in {time_sweep_dot1[0]} (ns)')
    ax.set_ylabel(f'Time in {time_sweep_dot2[0]} (ns)')
    ax.set_title(
        f'Varying Times in {time_sweep_dot1[0]} and {time_sweep_dot2[0]}')
