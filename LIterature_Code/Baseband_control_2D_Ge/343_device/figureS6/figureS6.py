# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 14:25:29 2023

@author: Francesco
"""


import array_343.figures.hopping_spins_paper.utils.analysis_tools as tools

#%% 

uuid_detuning_time_maps = dict()

#uuid_detuning_time_maps['4,8'] = 1691762567742283691
uuid_detuning_time_maps['8,5'] = 1691771577064283691
#uuid_detuning_time_maps['3,7'] = 1690556783281283691
#uuid_detuning_time_maps['6,10'] = 1689772467030283691


#%%

for map in uuid_detuning_time_maps:


    print('----------')
    print(f'Shuttling {map}')
    uuid = uuid_detuning_time_maps[map]
    
    dat = load_by_uuid(uuid)
    
        
    # Sample data
    x, x_label = dat.m1_3.j(), dat.m1_3.j.label
    y, y_label = dat.m1_3.i(), dat.m1_3.i.label
    z = dat.m1_3()
    
    # adjust the split depending on the dataset 
    # '4,8': 0.5
    # '8,5': 0.5
    # '3,7': 0.6
    # '6,10': 0.6
    
    split = 0.5
    max_diff = 16
    pop = []     
    pop = [list(tools.thresholded_data(z_2d, x, split = split, max_diff = max_diff, plot = False, sensor = None)[0]) for z_2d in z]
    pop = np.array(pop)
    
    
    
    fig, axs = plt.subplots(1, 2, figsize = tools.cm2inch(9,6))

    c = axs[0].pcolor(x,y,pop, shading = 'auto')
    axs[0].set_xlabel('time (ns)')
    axs[0].set_ylabel(r'$\mathrm{\epsilon_{' + str(map) + r'}}$ (mV)')
    
   
    z_fft = np.fft.fft(pop, axis=1)
    
    n = x.size
    timestep = x[1] - x[0]
    freq = np.fft.fftfreq(n, d=timestep)
    
    half_n = n // 2
    fft_norm = np.abs(z_fft)[:, 1:half_n]
    
    # Normalize by the global maximum
    global_max = fft_norm.max()
    fft_norm /= global_max
    
    f = axs[1].pcolor(freq[1:half_n], y, fft_norm, shading = 'auto')
    axs[1].set_xlabel('$f$ (GHz)')
    axs[1].set_ylabel(r'$\mathrm{\epsilon_{' + str(map) + r'}}$ (mV)')   
    
            
    fig.colorbar(c, ax=axs[0], label = '$P_\mathrm{\\uparrow}$', shrink = 0.5, orientation = 'horizontal', location = 'top' )
    fig.colorbar(c, ax=axs[1], label = '$A_\mathrm{FFT}$', shrink = 0.5, orientation = 'horizontal', location = 'top')
    
    fig.tight_layout()
    plt.show()
    if save_fig:        
        dir_name = os.path.dirname(__file__)           
        file_name = f"oscillations_{map}.pdf"
        file_path = os.path.join(dir_name, file_name)
        plt.savefig(file_path, format='pdf', dpi=300, transparent=True)

    
   
  