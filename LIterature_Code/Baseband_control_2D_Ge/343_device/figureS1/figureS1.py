# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 14:43:29 2023

@author: Francesco
"""



import array_343.figures.hopping_spins_paper.utils.analysis_tools as tools

#%% using a modified version of fit Ramsey data with different boundaries for the estimations.

def fit_Ramsey_data(x, data, freq_guess, short_measurement=False, tau_guess=10e3):
    """
    Fit Ramsey data and return best-fit parameters and their errors.

    Parameters:
    x (array): Independent variable.
    data (array): Dependent variable (data to be fitted).
    freq_guess (float): Initial guess for the frequency.
    short_measurement (bool): Flag for short measurement.
    tau_guess (float): Initial guess for decay time.

    Returns:
    tuple: Best-fit parameters and their standard deviations, or None if fitting fails.
    """

    # Guess amplitude based on measurement type
    amplitude_guess = 2 * np.std(data) * 1.5
    if short_measurement:
        amplitude_guess = 0.5 * np.abs(np.max(data) - np.min(data))

    offset_guess = np.average(data)
    initial_guess = [amplitude_guess, freq_guess, 0.0, tau_guess, offset_guess]
    
    
    #note that here the frequency is in GHz, therefore the range is given also in GHz, and not in Hz.
    try:
        best_fit_params, cov_matrix = curve_fit(
            tools.Ramsey, x, data, p0=initial_guess,
            # bounds=([amplitude_guess - 0.20, freq_guess - 10e6, -np.pi, 100, -1.0],
            #         [amplitude_guess + 0.20, freq_guess + 10e6, np.pi, 10e3, 1.0])
            bounds=([amplitude_guess - 0.20, freq_guess - 0.01, -2*np.pi, 100, 0],
                    [amplitude_guess + 0.20, freq_guess + 0.01, 2*np.pi, 10e3, 1.0])

        )

        # Calculate standard deviations (square root of the diagonal of the covariance matrix)
        perr = np.sqrt(np.diag(cov_matrix))

        # Print the best-fit parameters and their errors
        print("Best-fit parameters and their errors:")
        param_names = ["Amplitude (A)", "Frequency", "Phase", "Decay Time (tau)", "Constant Offset (C)"]
        for i, (param, error) in enumerate(zip(best_fit_params, perr)):
            print(f"{param_names[i]}: {param:.3f} \u00B1 {error:.3f}")

        return best_fit_params, perr

    except Exception as e:
        print(f"Error fitting data: {e}")
        return None
    
#%%


uuid_long_time_traces = dict()

# uuid_long_time_traces['Q1'] =  1693557950033283691
# uuid_long_time_traces['Q2'] =  1692020384925283691 
# uuid_long_time_traces['Q3'] =  1692282513469283691  
# uuid_long_time_traces['Q4'] =  1693415065370283691
# uuid_long_time_traces['Q5'] =  1692012726066283691
# uuid_long_time_traces['Q6'] =  1692171914624283691 
# uuid_long_time_traces['Q7'] =  1692729232124283691
# uuid_long_time_traces['Q8'] = 1692013520094283691 
uuid_long_time_traces['Q9'] =  1693834872772283691 
# uuid_long_time_traces['Q10'] = 1692628795405283691
 
#when fitting this datasets use as guess for the amplitude the opportune method

g_factors = dict()
decay_times = dict()
visibility = dict()
time_shift = dict()

g_factors_off = dict()
decay_times_off = dict()
visibility_off = dict()
time_shift_off = dict()

#number of shuttling steps from the original dots D1 or D4 to Dn
n_shuttle = dict()
n_shuttle['Q1'] = 4
n_shuttle['Q2'] = 6
n_shuttle['Q3'] = 10
n_shuttle['Q4'] = 4
n_shuttle['Q5'] = 4
n_shuttle['Q6'] = 8
n_shuttle['Q7'] = 12
n_shuttle['Q8'] = 2
n_shuttle['Q9'] = 10
n_shuttle['Q10'] = 10

n_shuttle_array = np.array(list(n_shuttle.values()))


#%% plotting the individual T2* decay


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
    
    fig, axs = plt.subplots(1,1, figsize = tools.cm2inch(6, 4))
    axs.set_xlabel(f'$t_\\mathrm{{D{qubit[1:]}}}$ (ns)')
    axs.set_ylabel('$P_{\\uparrow}$')   


    
    axs.scatter(y, z, s = 3, c = color)    
    
    #for Q9: improve fitting by choosing to exclude the first 70 ns
    if qubit == 'Q9':
        y = y[20:]
        z = z[20:]
    
    #finding main frequency
    freq, fft_norm, peak_freq, peak_amplitude = tools.find_major_peak_frequency(y, z, height_threshold=0.035) 
    
    try:
        freq_guess = peak_freq[np.argmax(peak_amplitude)]
    except:
        freq_guess = 0.01
   
    if qubit == 'Q7' or qubit == 'Q3':
        fit_params = fit_Ramsey_data(y, z, freq_guess, short_measurement=True)
    else:
        fit_params = fit_Ramsey_data(y, z, freq_guess, short_measurement=False)
        
    
    if fit_params is not None:
        A_fit, f_fit, phi_fit, tau_fit, C_fit = fit_params[0]
       
        g_factor_val = tools.g_factor_from_frequency(f_fit * 1e9, B = field_factor* 0.06)
        print(f'G_factor = {g_factor_val:.4f}')
        
        

        
        g_factors[qubit] = g_factor_val
        decay_times[qubit] = tau_fit
        visibility[qubit] = 2*A_fit
        
        y_fit = y #np.linspace(0, np.max(y), 1000)
        axs.plot(y_fit, tools.Ramsey(y_fit, *fit_params[0])
                       , color = 'black'
                       , linestyle = '--'
                       , linewidth = 0.5
                       , label = f'$\\mathrm{{D{qubit[1:]}}}: f = $ {fit_params[0][1]*1e3:.2f} MHz, $T^*_2$ = {fit_params[0][3]*1e-3:.3} $\\mu$s'
                       )
    
    axs.set_yticks([0.2, 0.5, 0.8])
    axs.set_yticklabels(['0.2', '0.5', '0.8'])
    axs.set_ylim(0.0, 1.0)
    
    plt.legend(loc = 'upper right')
    
    plt.tight_layout()
    plt.show()
    if save_fig:
        dir_name = os.path.dirname(__file__)    
        file_name = f'long_panel_q{qubit[1:]}.pdf'
        file_path = os.path.join(dir_name, file_name)
        plt.savefig(file_path
                    , format = 'pdf'
                    , dpi = 300
                    , transparent = True)  





g_array = np.array(list(g_factors.values()))
decay_times_array = np.array(list(decay_times.values()))
visibility_array = np.array(list(visibility.values()))




#%%
plt.close('all')
#%%
fig, axs = plt.subplots(3,1, figsize = tools.cm2inch(10, 12))
axs[0].scatter(n_shuttle_array, visibility_array, c = 'blue',  label = 'without time offset')
axs[0].set_ylabel('vis.')
axs[0].set_ylim(-0.05, 1.05)
axs[0].set_xticks([2, 4, 6, 8, 10, 12])
axs[0].set_xlabel('$n_{shuttle}$')

axs[1].scatter(n_shuttle_array, decay_times_array, c = 'blue' )
axs[1].set_xticks([2, 4, 6, 8, 10, 12])
axs[1].set_ylim(0, 2000)
axs[1].set_xlabel('$n_{shuttle}$')
axs[1].set_ylabel('$T^{\\ast}_{2}$ (ns)')


axs[2].set_xticks([2, 4, 6, 8, 10, 12])
#axs[1].set_ylim(0, 2000)
axs[2].set_xlabel('$n_{shuttle}$')
axs[2].set_ylabel('$t_{0}$ (ns)')


plt.tight_layout()
plt.show()


