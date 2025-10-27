import numpy as np
import matplotlib.pyplot as plt

def plot_signal_and_spectrum(t_dur_us=10.0):
    """
    Generates and plots a periodic signal and its frequency spectrum.

    Args:
        t_dur_us (float): The duration of the signal's ramp-down phase in microseconds.
    """
    if t_dur_us <= 0:
        print("Duration must be positive.")
        return

    # --- 1. Define Signal Parameters ---
    t_dur = t_dur_us * 1e-6  # Convert ramp duration to seconds
    T = 3 * t_dur            # Total period
    N = 2048                 # Number of samples for FFT
    dt = T / N               # Time resolution
    t = np.linspace(0, T - dt, N) # Time vector

    # --- 2. Generate the Signal in the Time Domain ---
    y = np.zeros_like(t)
    t1 = t_dur
    t2 = 2 * t_dur
    y[t < t1] = 1.0
    ramp_indices = (t >= t1) & (t < t2)
    y[ramp_indices] = (t2 - t[ramp_indices]) / t_dur
    y[t >= t2] = 1.0

    # --- 3. Calculate the Frequency Spectrum using FFT ---
    Y_fft = np.fft.fft(y)
    freqs = np.fft.fftfreq(N, d=dt)
    
    # Process for single-sided spectrum
    positive_freq_mask = freqs >= 0
    freqs_pos = freqs[positive_freq_mask]
    Y_fft_pos = Y_fft[positive_freq_mask]
    
    magnitude_spectrum = np.abs(Y_fft_pos) / N
    if N > 1:
        magnitude_spectrum[1:] *= 2

    # --- 4. Plot the Results ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Plot 1: Signal in the Time Domain
    ax1.plot(t * 1e6, y, 'b-')
    ax1.set_title(f'Periodic Signal (T = {T * 1e6:.1f} μs)')
    ax1.set_xlabel('Time (μs)')
    ax1.set_ylabel('Amplitude')
    ax1.set_xlim(0, T * 1e6)
    ax1.set_ylim(-0.1, 1.2)
    ax1.grid(True)

    # Plot 2: Frequency Spectrum (Magnitude)
    markerline, stemlines, baseline = ax2.stem(
        freqs_pos / 1e6, magnitude_spectrum, basefmt=" "
    )
    plt.setp(markerline, 'markersize', 4)
    ax2.set_title('Frequency Spectrum')
    ax2.set_xlabel('Frequency (MHz)')
    ax2.set_ylabel('Magnitude')
    
    fundamental_freq_Mhz = (1 / T) / 1e6 if T > 0 else 0
    ax2.set_xlim(0, 25 * fundamental_freq_Mhz)
    if len(magnitude_spectrum) > 1:
        max_harmonic_mag = np.max(magnitude_spectrum[1:])
        ax2.set_ylim(0, max_harmonic_mag * 1.25)
    
    ax2.grid(True)

    plt.tight_layout()
    plt.show()

# --- Call the function with the specified duration ---
plot_signal_and_spectrum(t_dur_us=10.0)