import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

# -- CORE MODELING FUNCTIONS (No changes here) --

def get_tline_abcd_matrix(
    Z0=50.0,
    line_length=1.0,
    num_segments=100,
    delta_R=0.1,
    prop_delay_ns_per_m=5.0
    ):
    """
    Calculates the cascaded ABCD matrix for the transmission line itself.
    The matrix elements A, B, C, D are returned as numpy polynomials.
    """
    segment_delay = (prop_delay_ns_per_m * 1e-9) * (line_length / num_segments)
    delta_L = Z0 * segment_delay
    delta_C = segment_delay / Z0
    s = np.array([1.0, 0.0])
    Z_series_half = np.polyadd(np.polymul([delta_L/2], s), [delta_R/2])
    Y_shunt = np.polymul([delta_C], s)
    M_z_half = [[np.array([1.0]), Z_series_half], [np.array([0.0]), np.array([1.0])]]
    M_y = [[np.array([1.0]), np.array([0.0])], [Y_shunt, np.array([1.0])]]

    def multiply_abcd_poly(m1, m2):
        A = np.polyadd(np.polymul(m1[0][0], m2[0][0]), np.polymul(m1[0][1], m2[1][0]))
        B = np.polyadd(np.polymul(m1[0][0], m2[0][1]), np.polymul(m1[0][1], m2[1][1]))
        C = np.polyadd(np.polymul(m1[1][0], m2[0][0]), np.polymul(m1[1][1], m2[1][0]))
        D = np.polyadd(np.polymul(m1[1][0], m2[0][1]), np.polymul(m1[1][1], m2[1][1]))
        return [[A, B], [C, D]]

    M_T = multiply_abcd_poly(M_z_half, multiply_abcd_poly(M_y, M_z_half))
    M_line = [[np.array([1.0]), np.array([0.0])], [np.array([0.0]), np.array([1.0])]]
    for _ in range(num_segments):
        M_line = multiply_abcd_poly(M_line, M_T)
    return M_line[0][0], M_line[0][1], M_line[1][0], M_line[1][1], multiply_abcd_poly

def create_lti_from_abcd(A, B, C, D, Rs, RL, Cs, CL):
    """
    Creates a full LTI system H(s) = V_load / V_source from the block's ABCD matrix
    and source/load impedances.
    """
    s = np.array([1.0, 0.0])
    Yl_num = np.polyadd(np.polymul([CL], s), [1/RL if RL > 1e-9 else 0])
    Yl_den = np.array([1.0])
    Term1_num = np.polyadd(np.polymul(A, Yl_den), np.polymul(B, Yl_num))
    Term2_num = np.polyadd(np.polymul([Rs*Cs],s), [1.0])
    Term3_num = np.polymul([Rs], np.polyadd(np.polymul(C, Yl_den), np.polymul(D, Yl_num)))
    Den_H_num = np.polyadd(np.polymul(Term1_num, Term2_num), Term3_num)
    Num_H = np.array([1.0])
    Den_H_num[np.abs(Den_H_num) < 1e-18] = 0.0
    return signal.lti(Num_H, Den_H_num)

def convert_abcd_to_s_params(A, B, C, D, Z0):
    """Converts ABCD polynomial parameters to S-parameter polynomials."""
    term_B_div_Z0 = np.polymul(B, [1.0/Z0])
    term_C_mul_Z0 = np.polymul(C, [Z0])
    den_s = np.polyadd(np.polyadd(A, D), np.polyadd(term_B_div_Z0, term_C_mul_Z0))
    num_s11 = np.polyadd(np.polyadd(A, term_B_div_Z0), np.polysub(np.array([0.0]), np.polyadd(term_C_mul_Z0, D)))
    num_s21 = np.array([2.0])
    num_s12 = np.polymul([2.0], np.polysub(np.polymul(A, D), np.polymul(B, C)))
    num_s22 = np.polyadd(np.polysub(term_B_div_Z0, A), np.polysub(D, term_C_mul_Z0))
    return {'S11': (num_s11, den_s), 'S21': (num_s21, den_s), 'S12': (num_s12, den_s), 'S22': (num_s22, den_s)}

# -- PLOTTING FUNCTIONS (No changes here) --

def plot_frequency_response(system, title):
    w_min, w_max = 2 * np.pi * 1e3, 2 * np.pi * 2e9; w = np.logspace(np.log10(w_min), np.log10(w_max), 500)
    w, mag, phase = system.bode(w)
    plt.figure(figsize=(12, 6)); plt.suptitle(title, fontsize=16)
    plt.subplot(2, 1, 1); plt.semilogx(w / (2 * np.pi), mag); plt.grid(True, which='both', linestyle='--'); plt.ylabel("Magnitude (dB)"); plt.title("Frequency Response")
    plt.subplot(2, 1, 2); plt.semilogx(w / (2 * np.pi), phase); plt.grid(True, which='both', linestyle='--'); plt.xlabel("Frequency (Hz)"); plt.ylabel("Phase (degrees)")
    plt.tight_layout(rect=[0, 0, 1, 0.96])

def plot_transient_response(system, f_square, sim_duration, title):
    t = np.linspace(0, sim_duration, 4000, endpoint=False); u = signal.square(2 * np.pi * f_square * t, duty=0.5) * 0.5 + 0.5
    tout, y, x = signal.lsim(system, U=u, T=t)
    plt.figure(figsize=(12, 6)); plt.plot(tout * 1e9, u, 'r--', label='Input Signal (at source)', linewidth=1.5); plt.plot(tout * 1e9, y, 'b-', label='Output Signal (at load)', linewidth=2)
    plt.grid(True, which='both', linestyle='--'); plt.title(title, fontsize=16); plt.xlabel("Time (ns)"); plt.ylabel("Voltage (V)"); plt.legend(); plt.ylim(-0.2, 1.2); plt.tight_layout()

def plot_s_parameters(s_params, Z0, title):
    freq = np.logspace(3, 9.3, 500); w = 2 * np.pi * freq; s_db = {}
    for name, (num, den) in s_params.items():
        s_complex = np.polyval(num, 1j * w) / np.polyval(den, 1j * w); s_db[name] = 20 * np.log10(np.abs(s_complex))
    plt.figure(figsize=(12, 6)); plt.suptitle(title, fontsize=16); plt.semilogx(freq, s_db['S11'], label='S11 (Return Loss)'); plt.semilogx(freq, s_db['S21'], label='S21 (Insertion Loss)')
    plt.grid(True, which='both', linestyle='--'); plt.xlabel("Frequency (Hz)"); plt.ylabel("Magnitude (dB)"); plt.legend(); plt.ylim(-60, 5); plt.title(f"S-Parameters (Z0 = {Z0}Ω)"); plt.tight_layout(rect=[0, 0, 1, 0.96])

# -- MAIN PROGRAM --
if __name__ == "__main__":
    
    # --- You can adjust all parameters here ---
    Z0 = 50.0; LINE_LENGTH = 1.0; N_SEGMENTS = 100
    DELTA_R_PER_SEGMENT = 0.01
    SOURCE_RESISTANCE = 50.0; LOAD_RESISTANCE = 1e6
    SOURCE_CAPACITANCE = 5e-12; LOAD_CAPACITANCE = 50e-12
    
    # *** NEW PARAMETER: Series resistor at the output ***
    # To disable, set to 0.0
    SERIES_OUTPUT_RESISTOR = 0 # 1 kOhm

    SQUARE_WAVE_FREQ = 1e6; SIMULATION_DURATION = 5 / SQUARE_WAVE_FREQ 

    # --- Script Execution ---
    
    # 1. Get the ABCD matrix for the transmission line itself
    print("Building transmission line ABCD matrix...")
    A_line, B_line, C_line, D_line, multiply_abcd_poly = get_tline_abcd_matrix(
        Z0=Z0, line_length=LINE_LENGTH, num_segments=N_SEGMENTS, delta_R=DELTA_R_PER_SEGMENT
    )
    
    # *** MODIFIED: Cascade the series output resistor if it exists ***
    A_sys, B_sys, C_sys, D_sys = A_line, B_line, C_line, D_line
    if SERIES_OUTPUT_RESISTOR > 1e-9:
        print(f"Cascading {SERIES_OUTPUT_RESISTOR} Ohm series resistor at the output...")
        M_line = [[A_line, B_line], [C_line, D_line]]
        # ABCD matrix for a series resistor R is [[1, R], [0, 1]]
        M_resistor = [[np.array([1.0]), np.array([SERIES_OUTPUT_RESISTOR])], 
                      [np.array([0.0]), np.array([1.0])]]
        
        M_combined = multiply_abcd_poly(M_line, M_resistor)
        A_sys, B_sys, C_sys, D_sys = M_combined[0][0], M_combined[0][1], M_combined[1][0], M_combined[1][1]

    # 2. Create the full LTI system using the (potentially new) system ABCD matrix
    print("Creating LTI system model...")
    lti_system = create_lti_from_abcd(
        A_sys, B_sys, C_sys, D_sys, # Use the combined system matrix
        SOURCE_RESISTANCE, LOAD_RESISTANCE, SOURCE_CAPACITANCE, LOAD_CAPACITANCE
    )
    
    # 3. Convert the *original line's* ABCD matrix to S-parameters (this is unchanged)
    print("Converting line-only ABCD to S-parameters...")
    s_parameters = convert_abcd_to_s_params(A_line, B_line, C_line, D_line, Z0)

    # 4. Generate and display all plots
    # *** MODIFIED: Update title string to include the new resistor ***
    resistor_str = f", R_out_series={SERIES_OUTPUT_RESISTOR}Ω" if SERIES_OUTPUT_RESISTOR > 1e-9 else ""
    sim_title_str = f"Rs={SOURCE_RESISTANCE}Ω, RL={LOAD_RESISTANCE}Ω, CL={LOAD_CAPACITANCE*1e12:.1f}pF{resistor_str}"
    line_title_str = f"Line: {LINE_LENGTH}m, {N_SEGMENTS} segments, ΔR_seg={DELTA_R_PER_SEGMENT}Ω"

    plot_frequency_response(lti_system, f"System Frequency Response\n({sim_title_str})")
    plot_transient_response(lti_system, SQUARE_WAVE_FREQ, SIMULATION_DURATION, f"Transient Response\n({sim_title_str})")
    plot_s_parameters(s_parameters, Z0, f"S-Parameter Plot of T-Line ONLY\n({line_title_str})")
    
    plt.show()