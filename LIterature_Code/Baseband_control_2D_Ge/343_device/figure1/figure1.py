# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 14:43:29 2023

@author: Francesco
"""


import array_343.figures.hopping_spins_paper.utils.analysis_tools as tools

#%% plotting magnetic field dependence Q6

uuid_field = dict()
uuid_field['Q6'] = 1693861160775283691

for qubit in uuid_field:

    print('----------')
    print(f'Qubit {qubit}')
    uuid = uuid_field[qubit]
    
    dat = load_by_uuid(uuid)
    
    x, y, z = dat.m1_2.y(), dat.m1_2.x(), dat.m1_3() #2
    x_label, y_label, z_label = dat.m1_2.y.label, dat.m1_2.x.label, dat.m1_2.label
    x_unit, y_unit, z_unit = dat.m1_2.y.unit, dat.m1_2.x.unit, dat.m1_2.unit
    
    
    i, j = 0, y.size 
    
    y = y * field_factor
    
    #cutting data in the y axis
    x = x
    y = y[i:j]
    z = z[i:j, :]
    
    split = 0.48 
    max_diff = 16
    pop = []     
    pop = [list(tools.thresholded_data(z_2d, x, split = split, max_diff = max_diff, plot = False, sensor = None)[0]) for z_2d in z]
    pop = np.array(pop)
    
    z = pop
    
    fig, axs = plt.subplots(2, 2, figsize=tools.cm2inch(9, 5),
                            gridspec_kw={'width_ratios': [1.0, 1.0], 'height_ratios': [0.2, 1]})

    c = axs[1, 0].pcolor(x, y, z, shading='auto', cmap='Blues')
    axs[1, 0].set_xlabel(f'$t_\\mathrm{{D{qubit[1:]}}}$ (ns)')
    axs[1, 0].set_ylabel('$B$ (T)')
     
    # Span the colorbar across the height of the heatmap subplot
    colorbar1 = fig.colorbar(c 
                             , ax=[axs[0, 0]]
                             , shrink = 0.6
                             , fraction = 0.7
                             , orientation = 'horizontal'
                             , location = 'top'
                             , pad = 0.00                         
                              )
     
    axs[0, 0].remove()
    axs[0, 1].remove()

   
    
    freq_values, amp_values, tau_values, C_values, phi_values = [], [], [], [], []

    for q, data in enumerate(z):
        print('-----')
        print(f'B = {y[q]:.4f} T')
        freq, fft_norm, peak_freq, peak_amplitude = tools.find_major_peak_frequency(x, data, height_threshold=0.05)
        try:
            freq_guess = peak_freq[np.argmax(peak_amplitude)]
        except:
            freq_guess = 0.01
        
        fit_params = tools.fit_Ramsey_data(x, data, freq_guess)
        if fit_params is not None:
            A_fit, f_fit, phi_fit, tau_fit, C_fit = fit_params[0]
            freq_guess = f_fit
            
            amp_values.append(A_fit)
            freq_values.append(f_fit)
            phi_values.append(phi_fit)
            tau_values.append(tau_fit)        
            C_values.append(C_fit)
        else:
            # Handle the case when fitting fails
            amp_values.append(np.nan)
            freq_values.append(np.nan)
            tau_values.append(np.nan)
            phi_values.append(np.nan)
            C_values.append(np.nan)



    z_fft = []
    for zx in z:
        zfft = np.fft.fft(zx)
        z_fft.append(zfft)
        
    z_fft = np.array(z_fft)

    q = x.size
    timestep = x[1]-x[0]
    freq = np.fft.fftfreq(q, d=timestep)
    fft_norm = np.abs(z_fft)[:,1:round(len(freq)/2)]/ np.abs(zfft).max()
    
    freq, fft_norm, peak_freq, peak_amplitude = tools.find_major_peak_frequency(x, data, height_threshold=0.02)
    
    try: 
        freq_guess = peak_freq[np.argmax(peak_amplitude)]
    except:
        freq_guess = 0.001
        
    fit_params = tools.fit_Ramsey_data(x, data, freq_guess)

    # Fit a linear function to frequency vs. magnetic field
    valid_indices = ~np.isnan(freq_values)  # Find non-NaN indices
    valid_frequencies = np.array(freq_values)[valid_indices]
    valid_fields =  np.array(y)[valid_indices]

    popt, _ = curve_fit(tools.linear_func, valid_fields, valid_frequencies,
                        bounds=([0, -0.001], [2, 0.001]))

    pks = valid_frequencies
    indices = np.where(np.abs(np.diff(pks)) > 0.005)[0] + 1
    sequences_pks = np.split(pks, indices)
    y_pks = np.split(valid_fields, indices)
    gs = []

    for sequence_pk, y_pk in zip(sequences_pks, y_pks):        
        if len(sequence_pk) > 15:
            popt, _ = curve_fit(tools.linear_func, y_pk, sequence_pk,
                                bounds=([0, -0.005], [5, 0.005]))
            alf = popt[0]
            axs[1, 1].scatter(sequence_pk*1e3, y_pk, s=6, c='black')
            x_fit = np.linspace(y_pk[0], y_pk[-1], 10)
            y_fit = popt[0]*x_fit + popt[1]
            g = 1e9 * tools.g_factor(1/alf) #check that g_factor is actually a function
            gs.append(g)   
            axs[1, 1].plot(y_fit*1e3,x_fit, label=f'$g$ = {np.round(g, 3)}')   
    
    
    axs[1, 1].set_xlabel('f (MHz)')
    axs[1, 1].set_xlim(xmin=0)    
    axs[1, 1].set_ylim((0,0.06*field_factor))
    
    g_av = np.average(gs)  
    print(f'Average g-factor: {g_av:.3f}')  
    
    axs[1, 1].legend(loc = 'lower right')   


    
    plt.tight_layout()    
    plt.show()
    if save_fig:
        dir_name = os.path.dirname(__file__)    
        file_name = f"Q6_data.pdf"
        file_path = os.path.join(dir_name, file_name)
        plt.savefig(file_path
                    , format = 'pdf'
                    , dpi = 300
                    , transparent = True)  


#%% plotting T2* at 41.4 mT

uuid_long_time_traces = dict()
uuid_long_time_traces['Q6'] =  1692171914624283691 

for qubit in uuid_long_time_traces:


    print('----------')
    print(f'Qubit {qubit}')
    uuid = uuid_long_time_traces[qubit]
    

    
    dat = load_by_uuid( uuid)
    
    # Sample data
    x, x_label = dat.m1_3.y(), dat.m1_3.y.label
    y, y_label  = dat.m1_3.x(), dat.m1_3.x.label
    z = dat.m1_3.z()
    
    # Check Nan data
    x_valid = ~np.isnan(x)
    y_valid = ~np.isnan(y)
    z = z[np.ix_(y_valid, x_valid)]

    x = x[x_valid]
    y = y[y_valid]
    
    qubit_number = int(y_label[-1])
    color = qubits[f'{qubit[1:]}']
    
    #add threshold
    z, _ = tools.thresholded_data(z, y, split = 0.5, max_diff = 10, plot = False, sensor = None)
    
    fig, axs = plt.subplots(1,1, figsize = tools.cm2inch(9.5, 4.0)) #9, 4.5
    axs.set_xlabel(f'$t_\\mathrm{{D{qubit[1:]}}}$ (ns)')
    #axs.set_ylabel('$P_\mathrm{up}$')   

    
    axs.scatter(y, z, s = 3)    
    
    #finding main frequency
    freq, fft_norm, peak_freq, peak_amplitude = tools.find_major_peak_frequency(y, z, height_threshold=0.035) 
    
    try:
        freq_guess = peak_freq[np.argmax(peak_amplitude)]
    except:
        freq_guess = 0.01
    #print(f'freq guess {freq_guess:.3f}')
    
    
    
    #fitting the time trace considering the main frequency
    fit_params = tools.fit_Ramsey_data(y, z, freq_guess)
    
    if fit_params is not None:
        
        fit_params = fit_params[0]
        A_fit, f_fit, phi_fit, tau_fit, C_fit = fit_params
        
        g_factor_number = tools.g_factor_from_frequency(f_fit * 1e9, B = field_factor * 0.06)
        print(f'G_factor = {g_factor_number:.4f}')
        
        y_fit = np.linspace(np.min(y), np.max(y), 1802)
    
        axs.plot(y_fit, tools.Ramsey(y_fit, *fit_params)
                       , color = 'black'
                       , linestyle = '--'
                       , linewidth = 0.5
                       , label = f'$\\mathrm{{D{qubit[1:]}}}: f = $ {fit_params[1]*1e3:.1f} MHz, $T^*_2$ = {fit_params[3]*1e-3:.3} $\\mu$s'
                       )
    
    axs.set_yticks([0.2, 0.4, 0.6])
    axs.set_xlim(-10, 1600)
    
    plt.legend(loc = 'upper right', fontsize = 5.5)
    
    plt.tight_layout()
    plt.show()

    if save_fig:
        dir_name = os.path.dirname(__file__)    
        file_name = f"Q6_T2star.pdf"
        file_path = os.path.join(dir_name, file_name)
        plt.savefig(file_path
                    , format = 'pdf'
                    , dpi = 300
                    , transparent = True) 






#%%
# fitted from B-field dependence in file: figureS2/figureS2_v3.py
g_factors = dict()
g_factors['Q1'] = 0.018
g_factors['Q2'] = 0.042
g_factors['Q3'] = 0.008
g_factors['Q4'] = 0.012
g_factors['Q5'] = 0.104
g_factors['Q6'] = 0.062
g_factors['Q7'] = 0.013
g_factors['Q8'] = 0.041
g_factors['Q9'] = 0.066
g_factors['Q10'] = 0.048

for key, val in zip(g_factors.keys(), g_factors.values()):
    print(f'{key}: {val:.3f}')
    
g_factor_average = np.average(list(g_factors.values()))
g_factor_std = np.std(list(g_factors.values()))
    

# fitted in file: figureS1/figureS1
T2 = dict()
T2['Q1'] = 1.69
T2['Q2'] = 1.51
T2['Q3'] = 0.99
T2['Q4'] = 1.65
T2['Q5'] = 1.59
T2['Q6'] = 1.12
T2['Q7'] = 0.88
T2['Q8'] = 1.70
T2['Q9'] =  0.29
T2['Q10'] = 1.49

for key, val in zip(T2.keys(), T2.values()):
    print(f'{key}: {val:.3f}')
    
T2_average = np.average(list(T2.values()))
T2_std = np.std(list(T2.values()))



#%% g factors
# Coordinates for 3-4-3 configuration
coordinates = [
    (0, 0), (2, 0), (4, 0),
    (-1., -1), (1, -1), (3, -1), (5, -1),
    (0, -2), (2, -2), (4, -2)
]

# Generating random quantities for demonstration
quantities = np.array(list(g_factors.values()))


norm = mcolors.Normalize(vmin=0, vmax=0.105, clip=True)
mapper = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.Blues) 

fig, ax = plt.subplots(figsize = tools.cm2inch(4.5, 4.5)) #5, 4.5
for coord, quantity in zip(coordinates, quantities):
    circle = plt.Circle(coord, 0.5, facecolor=mapper.to_rgba(quantity)
                        , edgecolor = 'black', linewidth = 0.5 )
    ax.add_artist(circle)

# Setting aspect ratio and limits for visualization
ax.set_aspect('equal', 'box')
ax.set_xlim(-2, 6)
ax.set_ylim(-3, 1)
ax.axis('off')

# Displaying the colorbar
cbar = plt.colorbar(mapper, ax=ax, orientation = 'vertical', location = 'right'
                , shrink = 0.3, aspect = 10
                , ticks = [0, 0.05, 0.10])
cbar.set_label('$g$-factor')

plt.tight_layout()
plt.show()
if save_fig:
    dir_name = os.path.dirname(__file__)    
    file_name = f"g-factors.pdf"
    file_path = os.path.join(dir_name, file_name)
    plt.savefig(file_path
                , format = 'pdf'
                , dpi = 300
                , transparent = True)
#%% decay_times
# Coordinates for 3-4-3 configuration
coordinates = [
    (0, 0), (2, 0), (4, 0),
    (-1., -1), (1, -1), (3, -1), (5, -1),
    (0, -2), (2, -2), (4, -2)
]


quantities = np.array([1688.3040097 , 1513.75626497,  992.11270191, 1650.59444553,
       1592.42909085, 1118.08458152,  878.30891318, 1698.73153599,
        294.978, 1493.55736558])

for quantity in quantities:
    print(np.round(quantity*1e-3, 2))


# Normalize quantities for coloring
norm = mcolors.Normalize(vmin=0, vmax=1700, clip=True)
mapper = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.Blues) 

fig, ax = plt.subplots(figsize = tools.cm2inch(4.5, 4.5)) #5, 4.5
for coord, quantity in zip(coordinates, quantities):
    circle = plt.Circle(coord, 0.5, facecolor=mapper.to_rgba(quantity)
                        , edgecolor = 'black', linewidth = 0.5 )
    ax.add_artist(circle)

# Setting aspect ratio and limits for visualization
ax.set_aspect('equal', 'box')
ax.set_xlim(-2, 6)
ax.set_ylim(-3, 1)
ax.axis('off')

# Displaying the colorbar
cbar = plt.colorbar(mapper, ax=ax, orientation = 'vertical', location = 'right'
                , shrink = 0.3, aspect = 10
                , ticks = [0, 850, 1700])
cbar.set_label('$T^*_2$ (ns)')

plt.tight_layout()
plt.show()
if save_fig:
    dir_name = os.path.dirname(__file__)    
    file_name = f"decay-times.pdf"
    file_path = os.path.join(dir_name, file_name)
    plt.savefig(file_path
                , format = 'pdf'
                , dpi = 300
                , transparent = True)
    


