import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# --- 1. Define the Physical Parameters ---
# These parameters are based on the paper's supplementary information.
# Note: These are example values and may need to be adjusted to precisely
# match the experimental results shown in the paper.
w0 = 7.0 * 2 * np.pi  # Resonator bare frequency (GHz)
Lr = 1.0             # Resonator inductance (arbitrary units)
Cr = 1 / (w0**2 * Lr) # Resonator capacitance
Zr = np.sqrt(Lr / Cr)  # Resonator impedance
E0 = 2.7 * 1e3       # SQUID Josephson energy (GHz)
delta = 0.0526       # SQUID junction asymmetry
Phi0 = 1.0           # Flux quantum (normalized to 1)
h_bar = 1.0          # Reduced Planck constant (normalized to 1)

# --- 2. Define the Time-Dependent Functions for Magnetic Flux ---
def flux_step(t, t_rise, phi_initial, phi_final):
    """Defines a magnetic flux step."""
    if t < 0:
        return phi_initial
    elif t < t_rise:
        return phi_initial + (phi_final - phi_initial) * t / t_rise
    else:
        return phi_final

def flux_impulse_train(t, t_impulse, t_interval, num_impulses, phi_base, phi_peak):
    """Defines a train of magnetic flux impulses."""
    for i in range(num_impulses):
        t_start = i * t_interval
        t_end = t_start + t_impulse
        if t_start <= t < t_end:
            return phi_peak
    return phi_base

# --- 3. Define the ODE System for the Simulation Parameters ---
def odes(t, y, flux_func, t_args):
    """
    System of ODEs for the parameters A, B, C, D, E, F.
    y = [A, B, C, D, E, F] (complex values)
    """
    A, B, C, D, E, F = y

    # Calculate time-dependent parameters based on the magnetic flux
    phi_ext = flux_func(t, **t_args)
    tan_arg = np.pi * phi_ext / Phi0
    
    # Avoid singularity in tan function
    if np.isclose(np.cos(tan_arg), 0):
        tan_val = 1e15 # A large number to approximate infinity
    else:
        tan_val = np.tan(tan_arg)

    theta = np.arctan(delta * tan_val)
    Ej = E0 * np.cos(tan_arg) * np.sqrt(1 + delta**2 * tan_val**2)
    
    # Avoid division by zero for Lj
    if np.isclose(Ej, 0):
        Lj = 1e15 # A large inductance
    else:
        Lj = (Phi0 / (2 * np.pi))**2 / Ej

    L_tot = Lr + Lj
    Z_t = np.sqrt(L_tot / Cr)

    v = w0 * (Lr / L_tot + 1) / 2
    f = w0 * (Lr / L_tot - 1) / 4
    g = (Phi0 * theta / (2 * np.pi * L_tot)) * np.sqrt(Lr * w0 / (2 * h_bar))

    # ODEs from Supplementary Note 2
    dAdt = -1j * (np.conj(f) * (2 * C + B**2) + np.conj(g) * B)
    dBdt = -1j * ((4 * C * np.conj(f) + v) * B + 2 * C * np.conj(g) + g)
    dCdt = -1j * (4 * np.conj(f) * C**2 + 2 * v * C + f)
    dDdt = -1j * (4 * C * np.conj(f) + v) * (D + 1)
    dEdt = -1j * (D + 1) * (2 * np.conj(f) * B + np.conj(g))
    dFdt = -1j * np.conj(f) * (D + 1)**2

    return [dAdt, dBdt, dCdt, dDdt, dEdt, dFdt]

# --- 4. Function to Run Simulation and Plot Results ---
def run_and_plot(flux_func, t_args, t_span, t_eval, y0, title):
    """Runs the simulation and plots input, time, and frequency domain signals."""
    sol = solve_ivp(odes, t_span, y0, args=(flux_func, t_args), t_eval=t_eval, method='RK45')

    A, B, C, D, E, F = sol.y

    # Calculate expectation value of the annihilation operator (time-domain signal)
    # Based on equations from Supplementary Note 2
    r_prime = np.arctanh(2 * np.abs(C))
    psi_prime = np.angle(C)
    alpha = B / np.cosh(r_prime)
    a_t = np.cosh(r_prime) * alpha + np.sinh(r_prime) * np.exp(1j * psi_prime) * np.conj(alpha)

    # --- Plotting ---
    fig = plt.figure(figsize=(12, 8))

    # Define the grid for subplots
    ax_flux = plt.subplot2grid((2, 2), (0, 0), colspan=2)
    ax_time = plt.subplot2grid((2, 2), (1, 0))
    ax_freq = plt.subplot2grid((2, 2), (1, 1))

    # Plot 1: Input Magnetic Flux Stimulus (Top, full width)
    flux_values = [flux_func(t, **t_args) for t in t_eval]
    ax_flux.plot(t_eval, flux_values, color='green', label='Input Flux')
    ax_flux.set_title(f'Input Magnetic Flux Stimulus ({title})')
    ax_flux.set_xlabel('Time (ns)')
    ax_flux.set_ylabel('Magnetic Flux (Φ₀)')
    ax_flux.grid(True)
    ax_flux.legend()

    # Plot 2: Time-Domain Signal (Bottom-Left)
    ax_time.plot(sol.t, np.real(a_t))
    ax_time.set_title('Resulting Time-Domain Signal')
    ax_time.set_xlabel('Time (ns)')
    ax_time.set_ylabel('Amplitude (Re[<â>])')
    ax_time.grid(True)

    # Plot 3: Frequency-Domain Signal (Bottom-Right)
    fft_result = np.fft.fft(a_t)
    freqs = np.fft.fftfreq(len(sol.t), d=(t_eval[1] - t_eval[0]))
    ax_freq.plot(np.fft.fftshift(freqs), np.fft.fftshift(np.abs(fft_result)))
    ax_freq.set_title('Resulting Frequency-Domain Signal')
    ax_freq.set_xlabel('Frequency (GHz)')
    ax_freq.set_ylabel('Magnitude')
    ax_freq.grid(True)
    ax_freq.set_xlim(0, 15) # Zoom in on relevant frequencies

    plt.tight_layout()
    plt.show()

# --- 5. Main Execution ---
if __name__ == "__main__":
    # --- Simulation for Magnetic Flux Step ---
    t_span_step = [0, 10]
    t_eval_step = np.linspace(t_span_step[0], t_span_step[1], 1000)
    y0 = np.zeros(6, dtype=complex) # Initial conditions (A=B=C=D=E=F=0)
    # The flux step must cross the 0.5 Phi0 point to generate a pulse
    step_args = {'t_rise': 1.0, 'phi_initial': 0.6 * Phi0, 'phi_final': 0.0}
    run_and_plot(flux_step, step_args, t_span_step, t_eval_step, y0, "Magnetic Flux Step")

    # --- Simulation for Magnetic Flux Impulse Train ---
    t_span_train = [0, 50]
    t_eval_train = np.linspace(t_span_train[0], t_span_train[1], 5000)
    train_args = {'t_impulse': 0.5, 't_interval': 5, 'num_impulses': 8, 'phi_base': 0.0, 'phi_peak': 0.6 * Phi0}
    run_and_plot(flux_impulse_train, train_args, t_span_train, t_eval_train, y0, "Magnetic Flux Impulse Train")