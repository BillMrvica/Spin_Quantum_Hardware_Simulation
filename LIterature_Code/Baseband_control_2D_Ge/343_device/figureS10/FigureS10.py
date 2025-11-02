# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 13:45:34 2024

@author: Francesco
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import array_343.figures.hopping_spins_paper.utils.analysis_tools as tools

#%% g-factor
dir_name = os.path.dirname(__file__)    
# Construct the file name and save the figure
file_name = f"variability_gfactor.csv"
file_path = os.path.join(dir_name, file_name)


# Assuming there are no headers in the file, we'll set them manually
column_headers = ['phi', 'theta', 'mean', 'standard_deviation']
data = pd.read_csv(file_path, header=None, names=column_headers)

# Display the first few rows to verify
print(data.head())
data.columns = column_headers

pivot_table = data.pivot("phi", "theta", "mean")
phi = np.array(pivot_table.columns) * 180/np.pi
theta = np.array(pivot_table.index)* 180/np.pi
mean = np.array(pivot_table)

pivot_table = data.pivot("phi", "theta", "standard_deviation")
phi = np.array(pivot_table.columns) * 180/np.pi
theta = np.array(pivot_table.index)* 180/np.pi
standard_deviation = np.array(pivot_table)

# Redoing the plot with a different approach to ensure the colorbar displays correctly
#plt.figure(figsize=cm2inch(10, 8))
fig, axs = plt.subplots(1,2, figsize = tools.cm2inch(15, 8))

cp = axs[0].contourf(phi, theta, mean, cmap='viridis')
#axs[0].contour(phi, theta, mean, colors='k')
fig.colorbar(cp, ax=axs[0], label='Average $g$-factor')
axs[0].set_xlabel('$\\theta \, (\deg)$ ')
axs[0].set_ylabel('$\phi \, (\deg)$ ')

cp = axs[1].contourf(phi, theta, standard_deviation, cmap='viridis')
#axs[1].contour(phi, theta, standard_deviation, colors='k')
fig.colorbar(cp, ax=axs[1], label='Standard deviation of the $g$-factor')
axs[1].set_xlabel('$\\theta \, (\deg)$ ')
axs[1].set_ylabel('$\phi \, (\deg)$ ')



plt.tight_layout()
plt.show()
if save_fig:
    # Get the directory name of the current script
    dir_name = os.path.dirname(__file__)    
    # Construct the file name and save the figure
    file_name = f"figureS10_gfactor.pdf"
    file_path = os.path.join(dir_name, file_name)
    plt.savefig(file_path
                , format = 'pdf'
                , dpi = 300
                , transparent = True)  


#%% quantisation axis
dir_name = os.path.dirname(__file__)    
# Construct the file name and save the figure
file_name = f"variability_quant.csv"
file_path = os.path.join(dir_name, file_name)


# Assuming there are no headers in the file, we'll set them manually
column_headers = ['phi', 'theta', 'mean', 'standard_deviation']
data = pd.read_csv(file_path, header=None, names=column_headers)

# Display the first few rows to verify
print(data.head())
data.columns = column_headers

pivot_table = data.pivot("phi", "theta", "mean")

phi = np.array(pivot_table.columns)
theta = np.array(pivot_table.index)
mean = np.array(pivot_table)

pivot_table = data.pivot("phi", "theta", "standard_deviation")
phi = np.array(pivot_table.columns) 
theta = np.array(pivot_table.index)
standard_deviation = np.array(pivot_table)

# Redoing the plot with a different approach to ensure the colorbar displays correctly
#plt.figure(figsize=cm2inch(10, 8))
fig, axs = plt.subplots(1,2, figsize = tools.cm2inch(15, 8))

cp = axs[0].contourf(phi, theta, mean, cmap='viridis')
#axs[0].contour(phi, theta, mean, colors='k')
fig.colorbar(cp, ax=axs[0], label='quantization axis difference, $\Delta \Phi \, (\deg)$')
axs[0].set_xlabel('$\\theta \, (\deg)$ ')
axs[0].set_ylabel('$\phi \, (\deg)$ ')

cp = axs[1].contourf(phi, theta, standard_deviation, cmap='viridis')
#axs[1].contour(phi, theta, standard_deviation, colors='k')
fig.colorbar(cp, ax=axs[1], label='Standard deviation of the quantization axis difference $\, (\deg)$')
axs[1].set_xlabel('$\\theta \, (\deg)$ ')
axs[1].set_ylabel('$\phi \, (\deg)$ ')

plt.tight_layout()
plt.show()
if save_fig:
    # Get the directory name of the current script
    dir_name = os.path.dirname(__file__)    
    # Construct the file name and save the figure
    file_name = f"figureS10_quantizationaxis.pdf"
    file_path = os.path.join(dir_name, file_name)
    plt.savefig(file_path
                , format = 'pdf'
                , dpi = 300
                , transparent = True)  