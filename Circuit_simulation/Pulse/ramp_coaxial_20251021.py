import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lsim
from scipy.optimize import curve_fit

def cable_loss_model(f, k1, k2):
    """
    Standard cable loss model: sqrt(f) for skin effect, f for dielectric loss.
    """
    f_safe = np.maximum(f, 1e-9)
    return k1 * np.sqrt(f_safe) + k2 * f_safe

def plot_final_corrected_response(ax, V_in_mV, ramp_time_ns, loss_params):
    """
    Plots the final, corrected transient response including the DAC signal.
    """
    # --- Parameters ---
    V_in = V_in_mV / 1000.0
    T_ramp = ramp_time_ns * 1e-9
    f_lpf = 400e6
    dac_period = 1e-9

    # Cable & Load Parameters
    Z0 = 50.0; L = 3.0
    C_load_values = [10e-12, 30e-12, 50e-12]
    c = 299792458.0
    v_prop = c * 0.70
    t_delay = L / v_prop

    # --- Simulation Time ---
    sim_time_end = max(t_delay + T_ramp + 15 * (Z0 * 50e-12), 60e-9)
    dt = 0.02e-9
    t_sim = np.arange(0, sim_time_end, dt)

    # --- 1. Generate Ideal & DAC Step Signal ---
    V_ideal = np.zeros_like(t_sim)
    if T_ramp > 1e-12:
        ramp_indices = (t_sim > 0) & (t_sim <= T_ramp)
        V_ideal[ramp_indices] = (V_in / T_ramp) * t_sim[ramp_indices]
        V_ideal[t_sim > T_ramp] = V_in
    else: # Ideal step
        V_ideal[t_sim > 0] = V_in
    
    V_dac_samples = np.interp(np.arange(0, sim_time_end, dac_period), t_sim, V_ideal)
    sample_index_for_t = (t_sim / dac_period).astype(int)
    V_dac = V_dac_samples[np.minimum(sample_index_for_t, len(V_dac_samples)-1)]
    
    ax.plot(t_sim * 1e9, V_ideal * 1000, 'k--', label='Ideal Input', alpha=0.5)
    ax.plot(t_sim * 1e9, V_dac * 1000, color='gray', label='DAC Output', alpha=0.7)

    # --- 2. Filter DAC signal with 400 MHz LPF ---
    tau_lpf = 1 / (2 * np.pi * f_lpf)
    lpf_sys = ([1], [tau_lpf, 1])
    _, V_cable_input, _ = lsim(lpf_sys, U=V_dac, T=t_sim)
    
    ax.plot(t_sim * 1e9, V_cable_input * 1000, 'm:', 
            label='Signal into Cable (Post-LPF)', lw=2)

    # --- 3. Apply Cable Attenuation (Time-Domain Method) ---
    atten_db_at_400MHz = cable_loss_model(400e6, *loss_params)
    total_db_loss = atten_db_at_400MHz * L
    loss_scaler = 10**(-total_db_loss / 20.0)
    
    V_at_load_input = V_cable_input * loss_scaler
    
    ax.plot((t_sim + t_delay) * 1e9, V_at_load_input * 1000, '-.', color='darkorange',
            label='Signal After Lossy Coax', lw=2)

    # --- 4. Calculate Final Response of Load Capacitors ---
    for C_load in C_load_values:
        tau_rc = Z0 * C_load
        load_sys = ([1], [tau_rc, 1])
        _, V_load, _ = lsim(load_sys, U=V_at_load_input, T=t_sim)
        
        ax.plot((t_sim + t_delay) * 1e9, V_load * 1000, 
                label=f'C={C_load*1e12:.0f}pF')

    ax.set_title(f'Amp = {V_in_mV} mV, Ramp = {ramp_time_ns} ns')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Voltage (mV)')
    ax.grid(True)
    ax.set_ylim(-0.05 * V_in_mV, 1.15 * V_in_mV)
    ax.set_xlim(-5, 55)
    ax.legend(fontsize='small', loc='lower right')
    
# --- Main execution ---
freq_data = np.array([100, 400, 900, 1000, 1800, 2400, 3000, 5000, 6000]) * 1e6
atten_data = np.array([0.29, 0.58, 0.88, 0.93, 1.27, 1.46, 1.65, 2.14, 2.34])
popt, _ = curve_fit(cable_loss_model, freq_data, atten_data)

amplitudes = [5, 15]
ramp_times = [0, 5, 10]

fig, axes = plt.subplots(len(amplitudes), len(ramp_times), figsize=(20, 11))
fig.suptitle('Final System Response Including DAC Output', fontsize=18, y=1.0)
fig.tight_layout(pad=4.5)

for i, amp in enumerate(amplitudes):
    for j, ramp_t in enumerate(ramp_times):
        plot_final_corrected_response(axes[i, j], amp, ramp_t, popt)

# Manually add delay line for clarity
delay_val = (3.0 / (299792458.0 * 0.70)) * 1e9
for ax_row in axes:
    for ax in ax_row:
        ax.axvline(x=delay_val, color='r', linestyle=':', lw=1)
axes[0,0].plot([], [], 'r:', label=f'Delay ({delay_val:.2f} ns)')[0]
axes[0,0].legend(fontsize='small', loc='lower right')

plt.show()