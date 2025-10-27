import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

# -- CORE MODELING FUNCTIONS (Provided by you, with one addition) --

def get_tline_abcd_matrix(
    Z0=50.0,
    line_length=1.0,
    num_segments=100,
    R_per_meter=0.1
    ):
    """
    Calculates the cascaded ABCD matrix for the transmission line itself.
    The matrix elements A, B, C, D are returned as numpy polynomials.
    """
    # Based on standard coax properties: delay per meter is sqrt(epsilon_r)/c
    # For common dielectrics (PE, PTFE), this is ~5.0 ns/m
    prop_delay_ns_per_m = 5.0
    
    segment_length = line_length / num_segments
    segment_delay = (prop_delay_ns_per_m * 1e-9) * segment_length
    delta_R = R_per_meter * segment_length
    
    # L/C parameters from Z0 and delay
    delta_L = Z0 * segment_delay
    delta_C = segment_delay / Z0

    # Represent 's' as a polynomial [1, 0]
    s = np.array([1.0, 0.0])

    # Symmetrical Pi-model for each segment
    Z_series = np.polyadd([delta_L], [delta_R]) # Corrected to use full R and L
    Y_shunt_half = np.polyadd(np.polymul([delta_C/2], s), [0]) # No G term for simplicity

    M_z = [[np.array([1.0]), Z_series], [np.array([0.0]), np.array([1.0])]]
    M_y_half = [[np.array([1.0]), np.array([0.0])], [Y_shunt_half, np.array([1.0])]]

    def multiply_abcd_poly(m1, m2):
        A = np.polyadd(np.polymul(m1[0][0], m2[0][0]), np.polymul(m1[0][1], m2[1][0]))
        B = np.polyadd(np.polymul(m1[0][0], m2[0][1]), np.polymul(m1[0][1], m2[1][1]))
        C = np.polyadd(np.polymul(m1[1][0], m2[0][0]), np.polymul(m1[1][1], m2[1][0]))
        D = np.polyadd(np.polymul(m1[1][0], m2[0][1]), np.polymul(m1[1][1], m2[1][1]))
        return [[A, B], [C, D]]

    # A single segment matrix using a symmetrical Pi model
    M_T = multiply_abcd_poly(M_y_half, multiply_abcd_poly(M_z, M_y_half))

    # Cascade all segments
    M_line = [[np.array([1.0]), np.array([0.0])], [np.array([0.0]), np.array([1.0])]]
    for _ in range(num_segments):
        M_line = multiply_abcd_poly(M_line, M_T)
        
    return M_line[0][0], M_line[0][1], M_line[1][0], M_line[1][1]

def create_transimpedance_lti(C, D, RL, CL):
    """
    Creates a transimpedance LTI system H(s) = V_load / I_in from the block's
    C, D matrix elements and the load impedance.
    """
    s = np.array([1.0, 0.0])

    # Z_L = 1 / Y_L. Y_L = 1/RL + s*CL
    if RL > 1e-9:
        num_yl = np.polyadd(np.polymul([CL], s), [1/RL])
        den_yl = np.array([1.0])
        Z_L_num, Z_L_den = den_yl, num_yl
    else: # Purely capacitive load
        Z_L_num = np.array([1.0])
        Z_L_den = np.polymul([CL], s)

    # Denominator of H(s) is C*Z_L + D
    # Den_H = C_num/C_den * Z_L_num/Z_L_den + D_num/D_den
    Den_H_num_part1 = np.polymul(C, Z_L_num)
    Den_H_den_part1 = Z_L_den # C_den and D_den are [1.0]
    Den_H_num_part2 = np.polymul(D, Z_L_den)
    
    Den_H_num = np.polyadd(Den_H_num_part1, Den_H_num_part2)
    Den_H_den = Den_H_den_part1
    
    # H(s) = Z_L / Den_H = (Z_L_num/Z_L_den) / (Den_H_num/Den_H_den)
    Num_H = np.polymul(Z_L_num, Den_H_den)
    Den_H = np.polymul(Z_L_den, Den_H_num)
    
    # Clean up near-zero terms that can result from polynomial math
    Num_H[np.abs(Num_H) < 1e-25] = 0.0
    Den_H[np.abs(Den_H) < 1e-25] = 0.0
    
    return signal.lti(Num_H, Den_H)

# -- MAIN PROGRAM --
if __name__ == "__main__":
    
    # --- System Parameters ---
    Z0_LINE = 50.0
    LINE_LENGTH = 2.0
    N_SEGMENTS = 100
    R_PER_METER = 0.1
    
    # --- Source and Signal Parameters ---
    CURRENT_STEP_AMPS = 1e-12 # 1 pA
    prop_delay_ns_per_m = 1/np.sqrt(1/(3e8)**2 * (1.602e-19)**2) * 1e9 # From Z0 and C_per_meter in your previous example
    prop_delay_ns_per_m = 5.0 # ns/m, typical for PE dielectric
    total_delay_s = (prop_delay_ns_per_m * 1e-9) * LINE_LENGTH
    SIMULATION_DURATION = total_delay_s * 10 # Simulate for 10x the delay

    # --- Load Scenarios ---
    loads = {
        'Matched (50Ω)': {'R': 50.0, 'C': 0},
        'High-Z Scope (1MΩ || 15pF)': {'R': 1e6, 'C': 15e-12},
        'Open Circuit (~1TΩ)': {'R': 1e12, 'C': 0}
    }

    # --- Setup Plots ---
    fig, (ax_t, ax_f) = plt.subplots(2, 1, figsize=(12, 10))
    fig.suptitle(f'Response of {LINE_LENGTH}m, {Z0_LINE}Ω Coax to a {CURRENT_STEP_AMPS*1e12:.1f} pA Current Step', fontsize=16)

    # --- Script Execution ---
    print("1. Building transmission line ABCD matrix...")
    A_line, B_line, C_line, D_line = get_tline_abcd_matrix(
        Z0=Z0_LINE, line_length=LINE_LENGTH, num_segments=N_SEGMENTS, R_per_meter=R_PER_METER
    )
    print("...Line model created.\n")

    # 2. Loop through each load scenario
    for name, params in loads.items():
        print(f"2. Simulating for Load: {name}...")
        RL, CL = params['R'], params['C']

        # 3. Create the transimpedance LTI system for this load
        system = create_transimpedance_lti(C_line, D_line, RL, CL)
        
        # 4. Perform Transient Analysis
        t = np.linspace(0, SIMULATION_DURATION, 4000)
        u_current = np.ones_like(t) * CURRENT_STEP_AMPS
        u_current[t < total_delay_s/100] = 0.0 # Create clean step
        
        tout, V_load_t, x = signal.lsim(system, U=u_current, T=t)
        
        # Plot transient response, scaling voltage to picoVolts
        ax_t.plot(tout * 1e9, V_load_t * 1e12, label=f'Load: {name}', lw=2)

        # 5. Perform Frequency Analysis
        w = np.logspace(np.log10(2*np.pi*1e3), np.log10(2*np.pi*1e9), 500)
        w, mag, phase = system.bode(w)
        freq_hz = w / (2 * np.pi)
        
        # Plot frequency response (magnitude is already in dB, but represents dB-Ohms)
        ax_f.semilogx(freq_hz, mag, label=f'Load: {name}', lw=2)
        print("...Simulation complete.")

    # --- Final Plot Formatting ---
    ax_t.set_title('Transient Response at Load')
    ax_t.set_xlabel('Time (ns)')
    ax_t.set_ylabel('Load Voltage (pV)')
    ax_t.axvline(total_delay_s * 1e9, color='red', ls=':', lw=1.5, label=f'Propagation Delay ({total_delay_s*1e9:.2f} ns)')
    ax_t.grid(True, which='both', linestyle='--')
    ax_t.legend()

    ax_f.set_title('Frequency Response (Transimpedance)')
    ax_f.set_xlabel('Frequency (Hz)')
    ax_f.set_ylabel('Gain (dBΩ)')
    ax_f.grid(True, which='both', linestyle='--')
    ax_f.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()