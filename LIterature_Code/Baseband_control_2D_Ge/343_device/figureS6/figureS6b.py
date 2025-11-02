# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 14:25:29 2023

@author: Francesco
"""

 
from array_343.figures.utils.package_style_database import thresholded_data

#%%


uuid_maps = dict()

uuid_maps['10'] = 1692884839064283691
uuid_maps['AC'] = 1692884848085283691
uuid_maps['01'] = 1692884857289283691



#%%

for map in uuid_maps:


    print('----------')
    print(f'Shuttling {map}')
    uuid = uuid_maps[map]
    
    dat = load_by_uuid( uuid)
    
        
    # Sample data
    x, x_label = dat.m1_3.y(), dat.m1_3.y.label
    y, y_label = dat.m1_3.x(), dat.m1_3.x.label
    z = dat.m1_3()  
    
    
    fig, axs = plt.subplots(1, 1, figsize = cm2inch(5,5))

    c = axs.pcolor(x,y,z, shading = 'auto')
    axs.set_ylabel('$\mathrm{U_{4,8}}$ (mV)')
    axs.set_xlabel(r'$\mathrm{\epsilon_{4,8}}$ (mV)')
    
    axs.scatter(0, 0, marker = 's', color = 'white', s = 10, edgecolor = 'black', linewidths = 0.5)   
        
    #fig.colorbar(c, ax=axs[0], label = '$P_\mathrm{up}$', shrink = 0.5, orientation = 'horizontal', location = 'top' )
    #fig.colorbar(c, ax=axs[1], label = '$A_\mathrm{FFT}$', shrink = 0.5, orientation = 'horizontal', location = 'top')
    
    fig.tight_layout()
    plt.show()
    if save_fig:        
        dir_name = os.path.dirname(__file__)           
        file_name = f"CSD_{map}.pdf"
        file_path = os.path.join(dir_name, file_name)
        plt.savefig(file_path, format='pdf', dpi=300, transparent=True)

    
   
  