# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 11:45:50 2023

@author: Francesco
"""

import array_343.figures.hopping_spins_paper.utils.analysis_tools as tools

uuid_maps = dict()
uuid_maps['73'] = 1690356693883283691
uuid_maps['710'] = 1690356720952283691
uuid_maps['106'] = 1690356775835283691
uuid_maps['26']= 1690356856203283691
uuid_maps['25']= 1690356883143283691
uuid_maps['14']= 1690356936770283691
uuid_maps['48']= 1690356964739283691
uuid_maps['15']= 1690356991668283691
uuid_maps['85']= 1690357018477283691
uuid_maps['96']= 1690357099743283691
uuid_maps['95']= 1690357126532283691



channel_maps = dict()
channel_maps['73'] = [2]
channel_maps['710'] = [4, 2]
channel_maps['106'] = [1, 4]
channel_maps['26'] = [1, 4]
channel_maps['25'] = [4, 1]
channel_maps['14'] = [1, 3] 
channel_maps['48'] = [4, 3]
channel_maps['15'] = [1, 4]
channel_maps['85'] = 4
channel_maps['96'] = [4, 1] 
channel_maps['95'] = 4

#%%
N = 2

plt.close('all')

for uuid_map in uuid_maps:
    
    print('----------')
    print(f'Map {uuid_map}')
    uuid = uuid_maps[uuid_map]
    
    dat = load_by_uuid(uuid)
    
    channels = channel_maps[uuid_map]
    if not isinstance(channels, list):  # If a single channel is provided, convert it to a list
        channels = [channels]

    x, y = dat.m1_2.y(), dat.m1_2.x()
    x_label, y_label, z_label = dat.m1_2.y.label, dat.m1_2.x.label, dat.m1_2.label
    x_unit, y_unit, z_unit = dat.m1_2.y.unit, dat.m1_2.x.unit, dat.m1_2.unit
    
    i, j = 3, len(x)  # Slicing only in the x direction
    
    # Initialize sum of derivatives
    sum_dz = np.zeros((len(y), j-i))

    for channel in channels:
        z = getattr(dat, f'm1_{channel}').z()  
        
        # Cutting and processing data in the x axis
        z = z[:, i:j]
        dx, dy, dz = tools.derivative(x[i:j], y, z, N, method='gaussian')

        # Accumulate the sum of derivatives
        sum_dz += dz
    
    # Plotting
    fig, axs = plt.subplots(1, 1, figsize=tools.cm2inch(4.5, 4.5), constrained_layout=True)
    sum_dz = (sum_dz - np.min(sum_dz))/(np.max(sum_dz) - np.min(sum_dz) )
    c = axs.pcolor(dx, dy, sum_dz, shading='auto')  # Changed to shading='auto' for better color mapping
    axs.set_xlabel('$\Delta$' + x_label + ' (' + x_unit + ')')
    axs.set_ylabel('$\Delta$' + y_label + ' (' + y_unit + ')')
    #fig.colorbar(c, ax=axs, orientation='horizontal')  # Added colorbar to the plot
    
    axs.scatter(0, 0, marker = 's', color = 'white', s = 10, edgecolor = 'black', linewidths = 0.5)
    #axs.grid(True)
    
    plt.show()

    if save_fig:
        # Get the directory name of the current script
        dir_name = os.path.dirname(__file__)
        
        # Construct the file name and save the figure
        file_name = f"{uuid_map}_{'_'.join(map(str, channels))}.pdf"
        file_path = os.path.join(dir_name, file_name)
        plt.savefig(file_path, format='pdf', dpi=300, transparent=True)

#%% plotting the colorbar


name = 'viridis'

 

fig, axes = plt.subplots(figsize = tools.cm2inch(4, 1.8), nrows=1)

gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))
axes.imshow(gradient, aspect='auto', cmap=name)

axes.set_yticks([])
axes.set_xticks([0, 256])
axes.set_xticklabels(['0','1'])

axes.set_xlabel('Norm. sensor signal (a.u.)')


plt.tight_layout()
plt.show()

if save_fig:
    # Get the directory name of the current script
    dir_name = os.path.dirname(__file__)    
    # Construct the file name and save the figure
    file_name = f"colormap.pdf"
    file_path = os.path.join(dir_name, file_name)
    plt.savefig(file_path, format='pdf', dpi=300, transparent=True)






    
#%% 
# Coordinates for 3-4-3 configuration
coordinates = [
    (0, 0), (2, 0), (4, 0),
    (-1., -1), (1, -1), (3, -1), (5, -1),
    (0, -2), (2, -2), (4, -2)
]

# Generating random quantities for demonstration
quantities = np.zeros(10)

# Normalize quantities for coloring
norm = mcolors.Normalize(vmin=quantities.min(), vmax=quantities.max(), clip=True)
mapper = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.Blues) 

fig, ax = plt.subplots(figsize = tools.cm2inch(2.5, 2.5))
for coord, quantity in zip(coordinates, quantities):
    circle = plt.Circle(coord, 0.5, facecolor=mapper.to_rgba(quantity)
                        , edgecolor = 'black', linewidth = 0.5 )
    ax.add_artist(circle)

# Setting aspect ratio and limits for visualization
ax.set_aspect('equal', 'box')
ax.set_xlim(-2, 6)
ax.set_ylim(-3, 1)
ax.axis('off')


plt.tight_layout()
plt.show()
if save_fig:
    dir_name = os.path.dirname(__file__)
    file_name = f"343_architecture.pdf"
    file_path = os.path.join(dir_name, file_name)
    plt.savefig(file_path, format='pdf', dpi=300, transparent=True)
