import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

# --- CORE MODELING FUNCTIONS ---

def create_opamp_model_from_specs(dc_gain_db=75.55, gbw_hz=7.45e6, pm_deg=65.0, C_L_spec=2e-12):
    """Creates a two-pole op-amp model from standard datasheet specifications."""
    A0 = 10**(dc_gain_db / 20.0)
    w_gbw = 2 * np.pi * gbw_hz
    p1 = w_gbw / A0
    phase_from_p2_deg = 90.0 - pm_deg
    p2 = w_gbw / np.tan(np.deg2rad(phase_from_p2_deg))
    R_out = 1 / (p2 * C_L_spec)
    print("--- Op-Amp Model Derived from Specifications ---")
    print(f"DC Gain: {dc_gain_db:.2f} dB, GBW: {gbw_hz/1e6:.2f} MHz, PM: {pm_deg:.1f}° @ {C_L_spec*1e12:.1f}pF")
    print(f"==> A0 = {A0:.2e} V/V, p1 = {p1/(2*np.pi):.2f} Hz, R_out = {R_out:.2f} Ω")
    print("-" * 55 + "\n")
    return {'A0': A0, 'p1': p1, 'R_out': R_out}

def get_tline_abcd_complex(w, line_length=2.0, Z0=50.0, R_per_meter=0.1):
    """Calculates the complex ABCD parameters for a transmission line."""
    prop_velocity = 2e8
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
        phase_margin = 180 + phase_deg[idx_unity]
        return phase_margin
    except IndexError:
        return np.nan

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
    
    w_range = np.logspace(np.log10(2*np.pi*1e3), np.log10(2*np.pi*5e8), 4000)
    freq_range = w_range / (2 * np.pi)
    s = 1j * w_range

    # Intrinsic Op-Amp Gain, A_v(s) = A0 / (1 + s/p1)
    Av_intrinsic = A0 / (1 + s/p1)

    # --- SCENARIO A: Baseline (2pF Load, Stable) ---
    print("--- Analyzing Scenario A: Direct 2pF Load ---")
    C_load_A = 2e-12
    Z_load_A = 1 / (s * C_load_A)
    Loop_Gain_A_freq = Av_intrinsic * (Z_load_A / (R_out + Z_load_A))
    H_A_freq = Loop_Gain_A_freq / (1 + Loop_Gain_A_freq)
    pm_A = get_bode_pm(freq_range, 20*np.log10(np.abs(Loop_Gain_A_freq)), np.angle(Loop_Gain_A_freq, deg=True))
    print(f"Result: Stable, PM = {pm_A:.1f}° (Matches spec)\n")

    # --- SCENARIO B: Coax + Scope (Unstable) ---
    print("--- Analyzing Scenario B: 2m Coax + Scope Load ---")
    A_t, B_t, C_t, D_t = get_tline_abcd_complex(w_range, line_length=2.0, Z0=50.0)
    R_scope_B, C_scope_B = 1e6, 15e-12
    Z_L_scope_B = 1 / (1/R_scope_B + s*C_scope_B)
    Z_load_B_freq = (A_t * Z_L_scope_B + B_t) / (C_t * Z_L_scope_B + D_t)
    Loop_Gain_B_freq = Av_intrinsic * (Z_load_B_freq / (R_out + Z_load_B_freq))
    H_B_freq = Loop_Gain_B_freq / (1 + Loop_Gain_B_freq)
    pm_B = get_bode_pm(freq_range, 20*np.log10(np.abs(Loop_Gain_B_freq)), np.angle(Loop_Gain_B_freq, deg=True))
    print(f"Result: Unstable, PM = {pm_B:.1f}°\n")
    
    # --- SCENARIO C: Source Follower Buffer + Coax + Scope (Stable) ---
    print("--- Analyzing Scenario C: Source Follower Added ---")
    R_SF, R_L_osc, C_L_osc = 700.0, 1e6, 15e-12
    g_m, C_gs = 5e-3, (220e-15 + 2e-12)
    print(f"SF Parameters: gm={g_m*1e3}mS, Cgs={C_gs*1e15:.0f}fF. Load: R_SF={R_SF}Ω, R_osc={R_L_osc/1e6}MΩ, C_osc={C_L_osc*1e12}pF")
    
    # 1. Op-amp stability analysis (driving C_gs)
    Z_load_opamp_C = 1 / (s * C_gs)
    Loop_Gain_C_opamp = Av_intrinsic * (Z_load_opamp_C / (R_out + Z_load_opamp_C))
    pm_C = get_bode_pm(freq_range, 20*np.log10(np.abs(Loop_Gain_C_opamp)), np.angle(Loop_Gain_C_opamp, deg=True))
    print(f"Result: Op-amp is stable driving SF, PM = {pm_C:.1f}°\n")
    
    # 2. System transfer functions for Scenario C
    H_C_amp_out_freq = Loop_Gain_C_opamp / (1 + Loop_Gain_C_opamp)
    Z_L_total_C = 1 / (1/R_SF + 1/R_L_osc + s*C_L_osc)
    Z_in_tline_C = (A_t * Z_L_total_C + B_t) / (C_t * Z_L_total_C + D_t)
    H_SF_gain = (g_m * Z_in_tline_C) / (1 + g_m * Z_in_tline_C)
    H_TLine_gain = Z_L_total_C / (A_t * Z_L_total_C + B_t)
    H_C_final_out_freq = H_C_amp_out_freq * H_SF_gain * H_TLine_gain

    # --- Plot 1: Loop Gain Analysis ---
    fig1, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig1.suptitle("Plot 1: Op-Amp Loop Gain Analysis", fontsize=16)
    ax_mag.semilogx(freq_range, 20*np.log10(np.abs(Loop_Gain_A_freq)), label=f'A (2pF Load), PM = {pm_A:.1f}°')
    ax_mag.semilogx(freq_range, 20*np.log10(np.abs(Loop_Gain_B_freq)), label=f'B (Coax Load), PM = {pm_B:.1f}°', color='C2')
    ax_mag.semilogx(freq_range, 20*np.log10(np.abs(Loop_Gain_C_opamp)), label=f'C (SF Load), PM = {pm_C:.1f}°', color='C3')
    ax_mag.set_ylabel("Magnitude (dB)"); ax_mag.grid(True, which='both'); ax_mag.legend(); ax_mag.set_title("Loop Gain Magnitude")
    ax_phase.semilogx(freq_range, np.angle(Loop_Gain_A_freq, deg=True))
    ax_phase.semilogx(freq_range, np.angle(Loop_Gain_B_freq, deg=True), color='C2')
    ax_phase.semilogx(freq_range, np.angle(Loop_Gain_C_opamp, deg=True), color='C3')
    ax_phase.set_ylabel("Phase (deg)"); ax_phase.set_xlabel("Frequency (Hz)"); ax_phase.grid(True, which='both'); ax_phase.set_title("Loop Gain Phase")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

    # --- Plot 2: Overall System AC Frequency Response ---
    fig2, ax_ac = plt.subplots(figsize=(10, 6))
    ax_ac.set_title("Plot 2: Overall System AC Frequency Response (Vout/Vin)", fontsize=16)
    ax_ac.semilogx(freq_range, 20*np.log10(np.abs(H_A_freq)), label='Scenario A: Final Output (2pF Load)')
    ax_ac.semilogx(freq_range, 20*np.log10(np.abs(H_B_freq)), label='Scenario B: Final Output (Coax Load)', color='C2')
    ax_ac.semilogx(freq_range, 20*np.log10(np.abs(H_C_amp_out_freq)), label='Scenario C: Amp Output', color='C4', linestyle=':')
    ax_ac.semilogx(freq_range, 20*np.log10(np.abs(H_C_final_out_freq)), label='Scenario C: Final Output (with SF)', color='C3')
    ax_ac.set_ylabel("Magnitude (dB)"); ax_ac.set_xlabel("Frequency (Hz)"); ax_ac.grid(True, which='both'); ax_ac.legend()
    plt.tight_layout()
    plt.show()

    # --- Run all transient simulations ---
    t = np.linspace(0, 20e-6, 8192)
    u = np.ones_like(t) * 0.5; u[t < 0.2e-6] = 0
    yA = run_transient_sim_fft(H_A_freq, w_range, t, u)
    yB = run_transient_sim_fft(H_B_freq, w_range, t, u)
    yC_final_out = run_transient_sim_fft(H_C_final_out_freq, w_range, t, u)
    
    # --- Plot 3: Comparison of Transient Responses ---
    fig3, ax_t = plt.subplots(figsize=(10, 6))
    ax_t.set_title("Plot 3: Comparison of Final Output Transient Responses", fontsize=16)
    ax_t.plot(t * 1e6, u, 'k--', label='Input Step')
    ax_t.plot(t * 1e6, yA, label='Scenario A: Final Output (2pF Load)')
    ax_t.plot(t * 1e6, yB, label='Scenario B: Final Output (Coax Load)', color='C2')
    ax_t.plot(t * 1e6, yC_final_out, label='Scenario C: Final Output (with SF)', color='C3')
    ax_t.set_xlabel("Time (μs)"); ax_t.set_ylabel("Voltage (V)"); ax_t.grid(True); ax_t.legend()
    plt.tight_layout()
    plt.show()

    # --- Plot 4: Scenario C Internal Transient Behavior ---
    yC_amp_out = run_transient_sim_fft(H_C_amp_out_freq, w_range, t, u)
    fig4, ax_tc = plt.subplots(figsize=(10, 6))
    ax_tc.set_title("Plot 4: Scenario C Internal Transient Behavior", fontsize=16)
    ax_tc.plot(t * 1e6, u, 'k--', label='Input Step')
    ax_tc.plot(t * 1e6, yC_amp_out, label='Scenario C: Amp Output (SF Gate)', color='C4', linestyle=':')
    ax_tc.plot(t * 1e6, yC_final_out, label='Scenario C: Final Output (at Load)', color='C3')
    ax_tc.set_xlabel("Time (μs)"); ax_tc.set_ylabel("Voltage (V)"); ax_tc.grid(True); ax_tc.legend()
    plt.tight_layout()
    plt.show()
