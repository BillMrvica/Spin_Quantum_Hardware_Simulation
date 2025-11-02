# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 13:36:01 2023

@author: Francesco
"""

#%%extracting virtual gate matrix  from data id

import array_343.figures.hopping_spins_paper.utils.analysis_tools as tools

dat = load_by_uuid(1690357153432283691)

gates = dat.snapshot['station']['instruments']['hardware']['virtual_gates']['exchanges']['real_gate_names']

virtual_gates = dat.snapshot['station']['instruments']['hardware']['virtual_gates']['exchanges']['virtual_gate_names']

matrix =  dat.snapshot['station']['instruments']['hardware']['virtual_gates']['exchanges']['virtual_gate_matrix']

matrix = json.loads(matrix)

#from list to array
matrix = np.array(matrix)
#39x39

#inverting the matrix
matrix_inv = inv(matrix)
#%% matrix 

fig = plt.figure(figsize = tools.cm2inch(14, 14)) #12.5 , 12.5
ax = fig.add_subplot(111)
im = ax.matshow(matrix, cmap = 'viridis_r') #viridis

ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
ax.xaxis.set_ticks(range(len(gates)))
ax.set_xticklabels(gates,  rotation = 70, fontsize = 7)


ax.yaxis.set_ticks(range(len(virtual_gates)))
ax.set_yticklabels(virtual_gates, fontsize = 7)

#fig.colorbar(cax, shrink = 0.5, location = 'bottom', orientation = 'horizontal')

cax3_inset = inset_axes(ax, width='3%', height='30%', loc='upper right')
cbar = fig.colorbar(im, cax=cax3_inset, orientation='vertical')
tick_locator = ticker.MaxNLocator(nbins=3)
cbar.locator = tick_locator
cbar.update_ticks()
cax3_inset.yaxis.set_ticks_position('left')
cax3_inset.yaxis.set_label_position('left')
cax3_inset.yaxis.set_tick_params(direction='in', labelcolor='black')

plt.tight_layout()
plt.show()


#%% matrix inv

fig = plt.figure(figsize = tools.cm2inch(14, 14)) #12.5 , 12.5
ax = fig.add_subplot(111)
im = ax.matshow(matrix_inv, cmap = 'viridis') #viridis

ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
ax.xaxis.set_ticks(range(len(virtual_gates)))
ax.set_xticklabels(virtual_gates,  rotation = 70, fontsize = 7)


ax.yaxis.set_ticks(range(len(gates)))
ax.set_yticklabels(gates, fontsize = 7)

#fig.colorbar(cax, shrink = 0.5, location = 'bottom', orientation = 'horizontal')

cax3_inset = inset_axes(ax, width='3%', height='25%', loc='upper right')
cbar = fig.colorbar(im, cax=cax3_inset, orientation='vertical')
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
cax3_inset.yaxis.set_ticks_position('left')
cax3_inset.yaxis.set_label_position('left')
cax3_inset.yaxis.set_tick_params(direction='in', labelcolor='black')

plt.tight_layout()
plt.show()

    
#%%
#%% cropping columns
matrix_inv_del = np.delete(matrix_inv, [25, 24, 23, 22, 21, 20 , 19, 18, 17, 16, 15, 14, 13, 12],  1 )

#%% cropping rows
matrix_inv_del_del = np.delete(matrix_inv_del,  [25, 24, 23, 22], 0 )

#%% cropping non-fast gates and virtual gates.
virtual_gates_del = np.delete(virtual_gates, [25, 24, 23, 22, 21, 20 , 19, 18, 17, 16, 15, 14, 13, 12], 0)
gates_del = np.delete(gates, [25, 24, 23, 22], 0)


#%% final figure
fig = plt.figure(figsize = tools.cm2inch(12, 12)) #12.5 , 12.5
ax = fig.add_subplot(111)
im = ax.matshow(matrix_inv_del_del, cmap = 'viridis') #viridis

ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
ax.xaxis.set_ticks(range(len(virtual_gates_del)))
ax.set_xticklabels(virtual_gates_del,  rotation = 70, fontsize = 7)


ax.yaxis.set_ticks(range(len(gates_del)))
ax.set_yticklabels(gates_del, fontsize = 7)

fig.colorbar(im, shrink = 0.3,  orientation = 'vertical', fraction = 0.1)

for i in range(len(gates_del)):
    for j in range(len(virtual_gates_del)):
        matrix_element = np.round(matrix_inv_del_del[i, j], 2)
        
        if np.round(matrix_element,2) != 0:
            if np.round(matrix_element,2) >0:
                text = ax.text(j, i, f'{matrix_element:.2f}',
                               ha="center", va="center", color="black", fontsize = 5.5)
            if np.round(matrix_element,2) <0: 
                text = ax.text(j, i, f'{matrix_element:.2f}',
                               ha="center", va="center", color="white", fontsize = 5.5)
            
            
ax.spines[:].set_visible(False)
ax.set_xticks(np.arange(matrix_inv_del_del.shape[1]+1)-.5, minor=True)
ax.set_yticks(np.arange(matrix_inv_del_del.shape[0]+1)-.5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth= 1)
ax.tick_params(which="minor", top=False, left=False)

plt.tight_layout()
plt.show()

if save_fig:
    # Get the directory name of the current script
    dir_name = os.path.dirname(__file__)    
    # Construct the file name and save the figure
    file_name = f"virtual_matrix_exchanges.pdf"
    file_path = os.path.join(dir_name, file_name)
    plt.savefig(file_path, format='pdf', dpi=300, transparent=True)



plt.tight_layout()
plt.show()





