# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:34:16 2024

@author: Francesco
"""

def normalize_array(x):
    y = (x - np.min(x)) / (np.max(x) - np.min(x))
    return y

from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt, butter
import array_343.figures.hopping_spins_paper.utils.analysis_tools as tools

#%% determining the splitter value per qubit dataset
uuid_field = dict()

uuid_field['Q1'] = 1693558701888283691
uuid_field['Q2'] = 1693855737612283691
uuid_field['Q3'] = 1689846576231283691
uuid_field['Q4'] = 1693398973179283691
uuid_field['Q5'] = 1693847413281283691
uuid_field['Q6'] = 1693861160775283691
uuid_field['Q7'] = 1689780917716283691 
uuid_field['Q8'] = 1693841659702283691
uuid_field['Q9'] = 1693836731202283691
uuid_field['Q10'] = 1693868872715283691

frequency_guess = dict()
frequency_guess['Q1'] = 0.01
frequency_guess['Q2'] = 0.01
frequency_guess['Q3'] = 0.001 
frequency_guess['Q4'] = 0.01
frequency_guess['Q5'] = 0.01
frequency_guess['Q6'] = 0.01
frequency_guess['Q7'] = 0.0057
frequency_guess['Q8'] = 0.01
frequency_guess['Q9'] = 0.0057
frequency_guess['Q10'] = 0.01


#%% updated guesses 

frequency_guess = dict()
frequency_guess['Q1'] = 0.01
frequency_guess['Q2'] = 0.01
frequency_guess['Q3'] = 0.0021
frequency_guess['Q4'] = 0.01
frequency_guess['Q5'] = 0.01
frequency_guess['Q6'] = 0.01
frequency_guess['Q7'] = 0.0023
frequency_guess['Q8'] = 0.01
frequency_guess['Q9'] = 0.01
frequency_guess['Q10'] = 0.01



#%% defining the guess for the threshold separating the two histrogram peaks

splitter = dict()

for qubit in uuid_field:
    uuid = uuid_field[qubit]
    dat = load_by_uuid( uuid)


    field, time, sensor_val, z_full = dat.m1_3.i(), dat.m1_3.j(), dat.m1_3.k(), dat.m1_3()
    x_label, y_label, z_label = dat.m1_2.y.label, dat.m1_2.x.label, dat.m1_2.label
    x_unit, y_unit, z_unit = dat.m1_2.y.unit, dat.m1_2.x.unit, dat.m1_2.unit
    
    z = z_full[0]
    sensor_val = normalize_array(sensor_val)
    
    fig, axs = plt.subplots(2, 1, figsize = tools.cm2inch(5, 10))
    axs[0].set_ylabel('time (ns)')
    axs[0].pcolor(sensor_val, time, z, shading = 'auto')
    averaged_histograms = np.average(z, axis = 0)
    axs[1].plot(sensor_val, averaged_histograms)
    axs[1].set_xlabel('Normalized sensor value')
    
    # Find the largest peaks in the smoothed gradient
    minimal_distance = 5
    number_largest_peaks = 2
    peaks, _ = find_peaks(averaged_histograms, distance=minimal_distance)  
    largest_peaks = peaks[np.argsort(averaged_histograms[peaks])[-number_largest_peaks:]]
    
    axs[1].scatter(sensor_val[largest_peaks], averaged_histograms [largest_peaks]
                    , color='red', zorder=5
                    , label=f'{number_largest_peaks} Largest Peaks')
    
    center = np.average(sensor_val[largest_peaks])
    splitter[qubit] = center
    
    plt.tight_layout()
    plt.show()
    
#%% re-threshold all the  (takes some minutes)
dataset = dict()

for qubit in uuid_field:
    print(f'QUBIT {qubit}')
    uuid = uuid_field[qubit]
    dat = load_by_uuid( uuid)
    
    y, x, sensor_val, z = dat.m1_3.i(), dat.m1_3.j(), dat.m1_3.k(), dat.m1_3()
    x_label, y_label, z_label = dat.m1_2.y.label, dat.m1_2.x.label, dat.m1_2.label
    x_unit, y_unit, z_unit = dat.m1_2.y.unit, dat.m1_2.x.unit, dat.m1_2.unit
    
    i, j = 0, y.size   
    
    #cutting data in the y axis
    x = x
    y = y[i:j]
    z = z[i:j, :]
    
    split = splitter[qubit]
    max_diff = 16
    pop = []     
    pop = [list(tools.thresholded_data(z_2d, x, split = split, max_diff = max_diff, plot = False, sensor = None)[0]) for z_2d in z]
    pop = np.array(pop)
    z = pop
    dataset[qubit] = z
    
#%% reassessing the plots

g_fac = dict()

for qubit in uuid_field:
    print('----------')
    print(f'Qubit {qubit}')
    uuid = uuid_field[qubit]
    
    dat = load_by_uuid( uuid)
    
    x, y, z = dat.m1_2.y(), dat.m1_2.x(), dat.m1_3() #2
    x_label, y_label, z_label = dat.m1_2.y.label, dat.m1_2.x.label, dat.m1_2.label
    x_unit, y_unit, z_unit = dat.m1_2.y.unit, dat.m1_2.x.unit, dat.m1_2.unit
    
    y = y * field_factor
    
    i, j = 0, y.size    
    B = 0.01
    l = 15 
    
    #cutting data in the y axis
    x = x
    y = y[i:j]
    z = dataset[qubit] #thresholded z
    
    fig, axs = plt.subplots(1, 2, figsize=tools.cm2inch(7.5, 5), sharey = True)

    # Plot the heatmap in the first subplot
    c = axs[0].pcolor(x, y, z, shading='auto', cmap =  'Blues')
    axs[0].set_xlabel(f'$t_\\mathrm{{D{qubit[1:]}}}$ (ns)')
    axs[0].set_ylabel('$B$ (T)')
    
    # Create colorbar axis
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("top", size="5%", pad=0.1)
    
    cb = fig.colorbar(c, cax=cax, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    
    # Create colorbar axis
    divider = make_axes_locatable(axs[1])
    cax1 = divider.append_axes("top", size="5%", pad=0.1)
    cax1.remove()
    
    freq_values, amp_values, tau_values, C_values, phi_values = [], [], [], [], []

    for q, data in enumerate(z):
        print('-----')
        print(f'B = {y[q]:.4f} T')
        
        if qubit == 'Q3':
            freq, fft_norm, peak_freq, peak_amplitude = tools.find_major_peak_frequency(x, data, height_threshold=0.01)
        else:
            freq, fft_norm, peak_freq, peak_amplitude = tools.find_major_peak_frequency(x, data, height_threshold=0.05)
        
        try:
            freq_guess = peak_freq[np.argmax(peak_amplitude)]
        except:
            freq_guess = frequency_guess[qubit] 
        
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


    n = tools.find_nearest_index(y, B)
    data =  z[n]
    
    z_fft = []
    for zx in z:
        zfft = np.fft.fft(zx)
        z_fft.append(zfft)
        
    z_fft = np.array(z_fft)

    # q = x.size
    # timestep = x[1]-x[0]
    # freq = np.fft.fftfreq(q, d=timestep)
    # fft_norm = np.abs(z_fft)[:,1:round(len(freq)/2)]/ np.abs(zfft).max()
    
    # freq, fft_norm, peak_freq, peak_amplitude = tools.find_major_peak_frequency(x, data, height_threshold=0.02)
    
    # try: 
    #     freq_guess = peak_freq[np.argmax(peak_amplitude)]
    # except:
    #     freq_guess = frequency_guess[qubit]
        
    # fit_params = tools.fit_Ramsey_data(x, data, freq_guess)
    
    # Fit a linear function to frequency vs. magnetic field
    valid_indices = ~np.isnan(freq_values)  # Find non-NaN indices
    valid_frequencies = np.array(freq_values)[valid_indices]
    valid_fields =  np.array(y)[valid_indices]
    
    # Filter valid_fields and adjust valid_frequencies correspondingly
    if qubit == 'Q7' or  qubit == 'Q1' or qubit == 'Q4' or qubit == 'Q3':
        threshold = 0.02
        filtered_indices = valid_fields > threshold
        valid_fields = valid_fields[filtered_indices]
        valid_frequencies = valid_frequencies[filtered_indices]

    popt, _ = curve_fit(tools.linear_func, valid_fields, valid_frequencies,
                        bounds=([0, -0.001], [5, 0.001]))

    pks = valid_frequencies
    indices = np.where(np.abs(np.diff(pks)) >  0.005)[0] + 1
    sequences_pks = np.split(pks, indices)
    y_pks = np.split(valid_fields, indices)
    gs = []

    for sequence_pk, y_pk in zip(sequences_pks, y_pks):        
        if len(sequence_pk) > l:
            popt, _ = curve_fit(tools.linear_func, y_pk, sequence_pk,
                                #bounds=([0, -0.0001], [5, 0.0001]))
                                bounds=([0, -0.005], [5, 0.005]))
            alf = popt[0]
            axs[1].scatter(sequence_pk, y_pk, s=6, c='black')
            x_fit = np.linspace(y_pk[0], y_pk[-1], 10)
            y_fit = popt[0]*x_fit + popt[1]
            g = 1e9 * tools.g_factor(1/alf)
            gs.append(g)   
            axs[1].plot(y_fit,x_fit, label=f'g = {np.round(g, 3)}', linestyle = '--')   
    
    
    axs[1].set_xlabel('$f$ (GHz)')
    
    axs[0].set_ylim((0,np.max(y)))
    axs[1].set_xlim(left = 0)
    g_av = np.average(gs)  
    print(f'Average g-factor: {g_av:.3f}')
    g_fac[f'{qubit}'] = g_av
    
    for ax in np.ravel(axs):
        ax.legend(loc = 'lower right')
    
    plt.tight_layout()
    plt.show()


    if save_fig:
        dir_name = os.path.dirname(__file__)    
        file_name = f"panel_q{qubit[1:]}.pdf"  #pdf
        file_path = os.path.join(dir_name, file_name)
        plt.savefig(file_path
                    , format = 'pdf' #pdf
                    , dpi = 300
                    , transparent = True)  
        
#%% printing the values
for (key, val) in zip(g_fac.keys(), g_fac.values()):
    print(key, np.round(val,3))
    

