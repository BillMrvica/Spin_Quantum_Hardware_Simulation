import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

# --- Helper Functions for Polynomial Transfer Function Algebra ---
# (These are for the stable transient simulation only)

def mul_tf(tf1, tf2):
    return np.polymul(tf1[0], tf2[0]), np.polymul(tf1[1], tf2[1])

def add_tf(tf1, tf2):
    num = np.polyadd(np.polymul(tf1[0], tf2[1]), np.polymul(tf2[0], tf1[1]))
    den = np.polymul(tf1[1], tf2[1])
    return num, den

def div_tf(tf1, tf2):
    return np.polymul(tf1[0], tf2[1]), np.polymul(tf1[1], tf2[0])


# --- CORE MODELING FUNCTIONS ---

def create_opamp_model_from_specs(dc_gain_db=75.55, gbw_hz=7.45e6, pm_deg=65.0, C_L_spec=2e-12):
    """
    Creates a two-pole op-amp model from standard datasheet specifications.
    The second pole is determined by the op-amp's output resistance and the specified load.
    """
    # 1. Calculate DC Gain (A0) and Dominant Pole (p1)
    A0 = 10**(dc_gain_db / 20.0)
    w_gbw = 2 * np.pi * gbw_hz
    # From GBW = A0 * p1
    p1 = w_gbw / A0
    
    # 2. Determine the location of the second pole (p2) from the phase margin spec
    # PM = 90 - atan(w_unity / p2), where w_unity is the unity-gain frequency.
    # We approximate w_unity ≈ w_gbw, which is a very good approximation for high-gain amps.
    phase_from_p2_deg = 90.0 - pm_deg
    # tan(phase_p2) = w_unity / p2  =>  p2 = w_unity / tan(phase_p2)
    p2 = w_gbw / np.tan(np.deg2rad(phase_from_p2_deg))
    
    # 3. Calculate the op-amp's intrinsic output resistance from p2 and C_L_spec
    # p2 = 1 / (R_out * C_L_spec)
    R_out = 1 / (p2 * C_L_spec)

    print("--- Op-Amp Model Derived from Specifications ---")
    print(f"Specified DC Gain: {dc_gain_db:.2f} dB  => A0 = {A0:.2e} V/V")
    print(f"Specified GBW: {gbw_hz/1e6:.2f} MHz")
    print(f"Specified PM: {pm_deg:.1f}° at {C_L_spec*1e12:.1f} pF")
    print("-" * 20)
    print(f"Calculated Dominant Pole (p1): {p1/(2*np.pi):.2f} Hz")
    print(f"Calculated Second Pole (p2): {p2/(2*np.pi)/1e6:.2f} MHz (with {C_L_spec*1e12:.1f}pF load)")
    print(f"==> Calculated Output Resistance (R_out): {R_out:.2f} Ω <==")
    print("-" * 45 + "\n")
    
    return {'A0': A0, 'p1': p1, 'R_out': R_out}

def get_tline_abcd_complex(w, line_length=1.0, Z0=50.0, R_per_meter=0.1):
    """Calculates the complex ABCD parameters using the exact hyperbolic solution."""
    # Simplified loss model, sufficient for this analysis
    prop_velocity = 2e8 # m/s, for typical PE/PTFE dielectric
    alpha = R_per_meter / (2 * Z0) 
    beta = w / prop_velocity
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
    """Calculates the time-domain response using FFT/IFFT."""
    dt = t_array[1] - t_array[0]
    U_w = np.fft.fft(u_t)
    fft_freq = np.fft.fftfreq(len(t_array), d=dt)
    w_fft = 2 * np.pi * fft_freq
    H_real = np.interp(w_fft, w_freq, np.real(H_freq))
    H_imag = np.interp(w_fft, w_freq, np.imag(H_freq))
    H_interp = H_real + 1j * H_imag
    Y_w = H_interp * U_w
    y_t = np.fft.ifft(Y_w)
    return np.real(y_t)

# -- MAIN SCRIPT --
if __name__ == "__main__":
    opamp = create_opamp_model_from_specs()
    R_out, A0, p1 = opamp['R_out'], opamp['A0'], opamp['p1']
    
    w_range = np.logspace(np.log10(2*np.pi*1e3), np.log10(2*np.pi*5e8), 2000)
    freq_range = w_range / (2 * np.pi)
    s = 1j * w_range

    # Intrinsic Op-Amp Gain, A_v(s) = A0 / (1 + s/p1)
    Av_intrinsic = A0 / (1 + s/p1)
    
    # --- SCENARIO A: Baseline (2pF Load, Stable) ---
    print("--- Verifying Model with Scenario A: Direct 2pF Load ---")
    C_load_A = 2e-12
    Z_load_A = 1 / (s * C_load_A)
    Loop_Gain_A_freq = Av_intrinsic * (Z_load_A / (R_out + Z_load_A))
    mag_lg_A = 20 * np.log10(np.abs(Loop_Gain_A_freq))
    phase_lg_A = np.angle(Loop_Gain_A_freq, deg=True)
    H_A_freq = Loop_Gain_A_freq / (1 + Loop_Gain_A_freq)

    # --- SCENARIO B: Coax + Scope (Unstable) ---
    print("--- Analyzing Scenario B: 1m Coax + Scope Load ---")
    A_t, B_t, C_t, D_t = get_tline_abcd_complex(w_range, line_length=2.0, Z0=50.0)
    R_scope, C_scope = 1e6, 15e-12
    Z_L_scope = 1 / (1/R_scope + s*C_scope)
    Z_load_B_freq = (A_t * Z_L_scope + B_t) / (C_t * Z_L_scope + D_t)
    Loop_Gain_B_freq = Av_intrinsic * (Z_load_B_freq / (R_out + Z_load_B_freq))
    mag_lg_B = 20 * np.log10(np.abs(Loop_Gain_B_freq))
    phase_lg_B = np.angle(Loop_Gain_B_freq, deg=True)
    H_B_freq = Loop_Gain_B_freq / (1 + Loop_Gain_B_freq)
    
    # --- Plotting Loop Gains (The Core Stability Analysis) ---
    fig_lg, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig_lg.suptitle("Loop Gain Analysis (Stability)", fontsize=16)
    
    pm_A, fu_A = get_bode_pm(freq_range, mag_lg_A, phase_lg_A)
    pm_B, fu_B = get_bode_pm(freq_range, mag_lg_B, phase_lg_B)
    
    ax_mag.semilogx(freq_range, mag_lg_A, label=f'Baseline (2pF Load)\nCalculated PM = {pm_A:.1f}° (Matches Spec!)')
    ax_mag.semilogx(freq_range, mag_lg_B, label=f'With Coax + Scope\nCalculated PM = {pm_B:.1f}° (Unstable!)', color='C2')
    ax_mag.set_ylabel("Magnitude (dB)"); ax_mag.grid(True, which='both'); ax_mag.legend(); ax_mag.set_title("Loop Gain Magnitude")
    
    ax_phase.semilogx(freq_range, phase_lg_A)
    ax_phase.semilogx(freq_range, phase_lg_B, color='C2')
    ax_phase.set_ylabel("Phase (degrees)"); ax_phase.set_xlabel("Frequency (Hz)"); ax_phase.set_title("Loop Gain Phase"); ax_phase.grid(True, which='both')
    
    # --- Plotting Transient Responses ---
    fig_t, (ax_tA, ax_tB) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig_t.suptitle("Closed-Loop Transient Response (Unity Gain)", fontsize=16)
    t = np.linspace(0, 4e-6, 8192)
    u = np.ones_like(t) * 0.5; u[t < 0.1e-6] = 0
    
    # Run transient simulation for Scenario A
    yA = run_transient_sim_fft(H_A_freq, w_range, t, u)
    ax_tA.plot(t * 1e6, u, 'r--', label='Input Step')
    ax_tA.plot(t * 1e6, yA, 'b', label='Output (2pF Load)')
    ax_tA.set_title("Scenario A: Stable Response"); ax_tA.legend(); ax_tA.grid(True); ax_tA.set_ylabel("Voltage (V)")

    # Run transient simulation for Scenario B
    yB = run_transient_sim_fft(H_B_freq, w_range, t, u)
    ax_tB.plot(t * 1e6, u, 'r--', label='Input Step')
    ax_tB.plot(t * 1e6, yB, 'g', label='Output (Coax + Scope)')
    ax_tB.set_title(f"Scenario B: Oscillatory Response (Direct Simulation via IFFT)")
    ax_tB.legend(); ax_tB.grid(True); ax_tB.set_xlabel("Time (μs)"); ax_tB.set_ylabel("Voltage (V)")
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()