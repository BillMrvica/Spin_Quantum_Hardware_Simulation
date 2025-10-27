import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

# --- CORE MODELING FUNCTIONS ---

def create_opamp_model(dc_gain_db=100, unity_gain_freq_hz=10e6):
    """
    Creates a behavioral model of a two-stage op-amp based on target specs.
    Returns: A dictionary of intrinsic parameters.
    """
    A0 = 10**(dc_gain_db / 20)
    f_unity = unity_gain_freq_hz
    w_unity = 2 * np.pi * f_unity
    p1 = w_unity / A0
    p2_at_2pf = 2.2 * w_unity
    z1 = 10 * w_unity
    C_L_baseline = 2e-12
    R_out = 1 / (p2_at_2pf * C_L_baseline)

    print("--- Op-Amp Model Derived ---")
    print(f"DC Gain (A0): {A0:.1e} V/V ({dc_gain_db} dB)")
    print(f"Unity Gain Freq: {f_unity/1e6:.1f} MHz")
    print(f"Calculated Output Resistance (R_out): {R_out:.1f} Ω")
    print("-" * 30 + "\n")
    
    return {'A0': A0, 'p1': p1, 'z1': z1, 'R_out': R_out}

def get_tline_abcd_complex(w, line_length=1.0, Z0=50.0, R_per_meter=0.1):
    """
    Calculates the complex ABCD parameters of a transmission line at given frequencies 'w'.
    This uses the exact hyperbolic solution, which is more stable than cascading.
    """
    prop_delay_ns_per_m = 5.0
    alpha = R_per_meter / (2 * Z0) # Attenuation constant (simplified)
    beta = w / (3e8 / np.sqrt(2.25)) # Phase constant
    
    gamma = alpha + 1j * beta
    
    A = np.cosh(gamma * line_length)
    B = Z0 * np.sinh(gamma * line_length)
    C = (1/Z0) * np.sinh(gamma * line_length)
    D = np.cosh(gamma * line_length)
    
    return A, B, C, D

def get_bode_pm(freq, mag_db, phase_deg):
    """Calculates Phase Margin from Bode plot data arrays."""
    try:
        idx_unity = np.where(mag_db < 0)[0][0]
        f_unity = freq[idx_unity]
        phase_margin = 180 + phase_deg[idx_unity]
        return phase_margin, f_unity
    except IndexError:
        return np.nan, np.nan

def run_transient_sim_fft(H_freq, w_freq, t_array, u_t):
    """
    Calculates the time-domain response using FFT/IFFT from frequency-domain data.
    """
    dt = t_array[1] - t_array[0]
    
    # FFT the input signal
    U_w = np.fft.fft(u_t)
    # Get the frequency bins corresponding to the FFT
    fft_freq = np.fft.fftfreq(len(t_array), d=dt)
    w_fft = 2 * np.pi * fft_freq
    
    # Interpolate the system's frequency response onto the FFT frequency bins
    # We must interpolate the real and imaginary parts separately
    H_real = np.interp(w_fft, w_freq, np.real(H_freq))
    H_imag = np.interp(w_fft, w_freq, np.imag(H_freq))
    H_interp = H_real + 1j * H_imag
    
    # Perform the multiplication in the frequency domain
    Y_w = H_interp * U_w
    
    # IFFT back to the time domain
    y_t = np.fft.ifft(Y_w)
    
    return np.real(y_t)

# -- MAIN SCRIPT --
if __name__ == "__main__":
    opamp = create_opamp_model()
    R_out, A0, p1, z1 = opamp['R_out'], opamp['A0'], opamp['p1'], opamp['z1']
    
    w_range = np.logspace(np.log10(2*np.pi*1e4), np.log10(2*np.pi*1e9), 2000)
    freq_range = w_range / (2 * np.pi)
    s = 1j * w_range

    # Intrinsic Op-Amp Gain, A_v(s), calculated numerically
    Av = A0 * (1 - s/z1) / (1 + s/p1)
    
    # --- SCENARIO A: Baseline (2pF Load, Stable) ---
    print("--- Analyzing Scenario A: Direct 2pF Load ---")
    C_load_A = 2e-12
    Z_load_A = 1 / (s * C_load_A)
    Loop_Gain_A_freq = Av * (Z_load_A / (R_out + Z_load_A))
    mag_lg_A = 20 * np.log10(np.abs(Loop_Gain_A_freq))
    phase_lg_A = np.angle(Loop_Gain_A_freq, deg=True)
    H_A_freq = Loop_Gain_A_freq / (1 + Loop_Gain_A_freq)

    # --- SCENARIO B: Coax + Scope (Unstable) ---
    print("--- Analyzing Scenario B: 2m Coax + Scope Load ---")
    A_t, B_t, C_t, D_t = get_tline_abcd_complex(w_range, line_length=1.0, Z0=50.0)
    R_scope, C_scope = 1e6, 15e-12
    Z_L_scope = 1 / (1/R_scope + s*C_scope)
    Z_load_B_freq = (A_t * Z_L_scope + B_t) / (C_t * Z_L_scope + D_t)
    Loop_Gain_B_freq = Av * (Z_load_B_freq / (R_out + Z_load_B_freq))
    mag_lg_B = 20 * np.log10(np.abs(Loop_Gain_B_freq))
    phase_lg_B = np.angle(Loop_Gain_B_freq, deg=True)
    H_B_freq = Loop_Gain_B_freq / (1 + Loop_Gain_B_freq)
    
    # --- Plotting Loop Gains (The Core Stability Analysis) ---
    fig_lg, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig_lg.suptitle("Loop Gain Analysis (Stability)", fontsize=16)
    
    pm_A, fu_A = get_bode_pm(freq_range, mag_lg_A, phase_lg_A)
    pm_B, fu_B = get_bode_pm(freq_range, mag_lg_B, phase_lg_B)
    
    ax_mag.semilogx(freq_range, mag_lg_A, label=f'Baseline (2pF Load)\nPM = {pm_A:.1f}°, f_unity = {fu_A/1e6:.2f} MHz')
    ax_mag.semilogx(freq_range, mag_lg_B, label=f'With Coax + Scope\nPM = {pm_B:.1f}° (Unstable!)', color='C2')
    ax_mag.set_ylabel("Magnitude (dB)"); ax_mag.grid(True, which='both'); ax_mag.legend(); ax_mag.set_title("Loop Gain Magnitude")
    
    ax_phase.semilogx(freq_range, phase_lg_A)
    ax_phase.semilogx(freq_range, phase_lg_B, color='C2')
    ax_phase.set_ylabel("Phase (degrees)"); ax_phase.set_xlabel("Frequency (Hz)"); ax_phase.set_title("Loop Gain Phase"); ax_phase.grid(True, which='both')
    
    # --- Plotting Transient Responses ---
    fig_t, (ax_tA, ax_tB) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig_t.suptitle("Closed-Loop Transient Response (Unity Gain)", fontsize=16)
    t = np.linspace(0, 2e-6, 8192) # Use a power of 2 for FFT efficiency
    u = np.ones_like(t) * 0.5; u[t < 0.1e-6] = 0 # 0.5V step
    
    # Run transient simulations using the FFT-based method
    yA = run_transient_sim_fft(H_A_freq, w_range, t, u)
    ax_tA.plot(t * 1e6, u, 'r--', label='Input Step')
    ax_tA.plot(t * 1e6, yA, 'b', label='Output (2pF Load)')
    ax_tA.set_title("Scenario A: Stable Response"); ax_tA.legend(); ax_tA.grid(True); ax_tA.set_ylabel("Voltage (V)")

    yB = run_transient_sim_fft(H_B_freq, w_range, t, u)
    ax_tB.plot(t * 1e6, u, 'r--', label='Input Step')
    ax_tB.plot(t * 1e6, yB, 'g', label='Output (Coax + Scope)')
    ax_tB.set_title(f"Scenario B: Oscillatory Response (Direct Simulation via IFFT)")
    ax_tB.legend(); ax_tB.grid(True); ax_tB.set_xlabel("Time (μs)"); ax_tB.set_ylabel("Voltage (V)")
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()