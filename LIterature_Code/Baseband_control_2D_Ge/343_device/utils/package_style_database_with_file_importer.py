# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:44:55 2023

@author: Francesco
"""

# %% importing core-tools
from os.path import isfile, join
from os import listdir
from numbers import Number
import logging
import inspect
from scipy.stats import linregress
import imageio
import scipy.optimize as optimize
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from matplotlib import gridspec as gridspec
from matplotlib.gridspec import GridSpec
from functools import partial
from mpl_interactions import zoom_factory, panhandler
from matplotlib.patches import Polygon
from mpl_point_clicker import clicker
from copy import copy
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
                                                   mark_inset)
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.colors as colors
from scipy import ndimage, integrate
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt, butter
import json
from numpy import savetxt, loadtxt
from numpy import save, savez, load
from numpy import asarray
from numpy.linalg import inv
import numpy as np
import pylab
import importlib
import time
import scipy.constants as constant
from core_tools.drivers.hardware.utility import load_virtual_gate_matrix_from_ds
from core_tools.utility.powerpoint import addPPTslide
from core_tools.drivers.hardware.hardware import hardware
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm, Normalize

from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib_inline
import matplotlib as mpl
import core_tools as ct


#%% loading from file
from core_tools.data.ds.reader import load_by_uuid, load_by_id, set_data_location #loading from files

# %% checking user
import sys
import os

#%%
path_datasets = r'C:\Code\spin-projects\stations\vaughan\array_343\figures\hopping_spins_paper\datasets'
set_data_location(None) #loading from the database
set_data_location(path_datasets) #loading from the database


# %% testing connection
print('testing data import...')
dat = load_by_uuid(1693861160775283691)

# %% defining style
plt.rcParams.update({'font.size': 7})
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams["figure.figsize"] = (2.28, 1.73)
plt.rcParams['figure.dpi'] = 150
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['mathtext.default'] = 'it'
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 'small'
plt.rcParams['legend.scatterpoints'] = 1
plt.rcParams['axes.labelpad'] = 4
fontsize_text = 6

qubits = dict()
qubits_colors = [
    "#1f77b4",  # muted blue
    "#ff7f0e",  # safety orange
    "#2ca02c",  # cooked asparagus green
    "#d62728",  # brick red
    "#9467bd",  # muted purple
    "#8c564b",  # chestnut brown
    "#e377c2",  # raspberry yogurt pink
    "#7f7f7f",  # middle gray
    "#bcbd22",  # curry yellow-green
    "#17becf"   # blue-teal
]

for i in range(1, 11):
    qubits[f'{i}'] = qubits_colors[i-1]


save_fig = False
titles = False
save_image = False

field_factor = 0.69

ueVtoGHz = 0.2418
Plank_constant = 4.136  # ueV/GHz
kT = 11.9  # ueV

letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
           'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']

Blues = cm.get_cmap('Blues', 256)
New_Blues = ListedColormap(Blues(np.linspace(0.0, 0.90, int(256*0.90))))
Reds = cm.get_cmap('Reds', 256)
New_Reds = ListedColormap(Reds(np.linspace(0.0, 0.90, int(256*0.90))))


def cl():
    plt.close('all')

