import numpy as np
import matplotlib.pyplot as plt

# --- Parameters to Adjust ---

# Microwave frequency in Hz
f_mw = 13e9  # 13 GHz

# Envelope width in seconds
width = 100e-9  # 100 ns

# Signal amplitude in Volts
amplitude = 0.050  # 50 mV


# --- Time Axis Generation ---

# We need a time axis that is fine enough to resolve the microwave frequency.
# For a smooth plot and accurate FFT, a high sampling rate is good. Let's use 20x the frequency.
sampling_rate = 20 * f_mw
time_span = 8 * width  # Extend time span for better frequency resolution

# Generate time points
num_points = int(time_span * sampling_rate)
if num_points % 2 == 0:
    num_points += 1
t = np.linspace(-time_span / 2, time_span / 2, num_points)
dt = t[1] - t[0] # time step


# --- Signal Generation Functions ---

def generate_rectangular_pulse(time_array, freq, pulse_width, amp):
    """Generates a sinusoidal signal with a rectangular envelope."""
    carrier = np.sin(2 * np.pi * freq * time_array)
    half_width = pulse_width / 2
    envelope = np.where((time_array >= -half_width) & (time_array <= half_width), 1, 0)
    signal = amp * carrier * envelope
    return signal, envelope * amp

def generate_gaussian_pulse(time_array, freq, fwhm_width, amp):
    """Generates a sinusoidal signal with a Gaussian envelope."""
    carrier = np.sin(2 * np.pi * freq * time_array)
    sigma = fwhm_width / (2 * np.sqrt(2 * np.log(2)))
    envelope = np.exp(-time_array**2 / (2 * sigma**2))
    signal = amp * carrier * envelope
    return signal, envelope * amp


# --- FFT Calculation Function ---

def calculate_fft(signal, time_step):
    """Calculates the FFT of a signal and returns the centered frequency axis and spectrum."""
    # Perform FFT
    fft_raw = np.fft.fft(signal)
    
    # Generate frequency axis and center it
    freq_axis = np.fft.fftfreq(len(signal), d=time_step)
    
    # Shift both the FFT output and the frequency axis to be centered around 0 Hz
    fft_shifted = np.fft.fftshift(fft_raw)
    freq_shifted = np.fft.fftshift(freq_axis)
    
    return freq_shifted, np.abs(fft_shifted)


# --- Generate Signals and their FFTs ---

# Generate signals
rect_signal, rect_envelope = generate_rectangular_pulse(t, f_mw, width, amplitude)
gauss_signal, gauss_envelope = generate_gaussian_pulse(t, f_mw, width, amplitude)

# Calculate FFTs
freq_axis, rect_fft = calculate_fft(rect_signal, dt)
_, gauss_fft = calculate_fft(gauss_signal, dt)


# --- Plotting ---

fig, axes = plt.subplots(2, 2, figsize=(15, 10))
((ax1, ax2), (ax3, ax4)) = axes

# Plot 1: Rectangular Pulse (Time Domain)
ax1.set_title(f'Rectangular Pulse: {amplitude*1000} mV, {width*1e9} ns width')
ax1.plot(t * 1e9, rect_signal * 1000, label='Signal', color='blue')
ax1.plot(t * 1e9, rect_envelope * 1000, label='Envelope', color='red', linestyle='--')
ax1.plot(t * 1e9, -rect_envelope * 1000, color='red', linestyle='--')
ax1.set_ylabel('Amplitude (mV)')
ax1.set_xlabel('Time (ns)')
ax1.legend()
ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

# Plot 2: Rectangular Pulse (Frequency Domain)
ax2.set_title('Frequency Spectrum')
ax2.plot(freq_axis / 1e9, rect_fft, color='blue')
ax2.set_ylabel('Magnitude (a.u.)')
ax2.set_xlabel('Frequency (GHz)')
ax2.set_xlim((f_mw - 0.2*f_mw)/1e9, (f_mw + 0.2*f_mw)/1e9) # Zoom into the carrier frequency
ax2.grid(True, which='both', linestyle='--', linewidth=0.5)

# Plot 3: Gaussian Pulse (Time Domain)
ax3.set_title(f'Gaussian Pulse: {amplitude*1000} mV, {width*1e9} ns FWHM')
ax3.plot(t * 1e9, gauss_signal * 1000, label='Signal', color='green')
ax3.plot(t * 1e9, gauss_envelope * 1000, label='Envelope', color='red', linestyle='--')
ax3.plot(t * 1e9, -gauss_envelope * 1000, color='red', linestyle='--')
ax3.set_ylabel('Amplitude (mV)')
ax3.set_xlabel('Time (ns)')
ax3.legend()
ax3.grid(True, which='both', linestyle='--', linewidth=0.5)

# Plot 4: Gaussian Pulse (Frequency Domain)
ax4.set_title('Frequency Spectrum')
ax4.plot(freq_axis / 1e9, gauss_fft, color='green')
ax4.set_ylabel('Magnitude (a.u.)')
ax4.set_xlabel('Frequency (GHz)')
ax4.set_xlim((f_mw - 0.2*f_mw)/1e9, (f_mw + 0.2*f_mw)/1e9) # Zoom into the carrier frequency
ax4.grid(True, which='both', linestyle='--', linewidth=0.5)


plt.tight_layout()
plt.show()
