import numpy as np
import matplotlib.pyplot as plt

def gaussian_pulse(t, amplitude, center, std_dev):
    """
    Generates a Gaussian pulse.

    Args:
        t (np.array): Time array.
        amplitude (float): The peak amplitude of the pulse.
        center (float): The time at which the pulse is centered.
        std_dev (float): The standard deviation of the pulse.

    Returns:
        np.array: The Gaussian pulse signal.
    """
    return amplitude * np.exp(-((t - center)**2) / (2 * std_dev**2))

# --- Signal Parameters (inspired by the provided figure) ---
peak_amplitude_mv = 220  # Peak amplitude in mV
pulse_center_ns = 50     # Center of the pulse in ns
std_deviation_ns = 15    # Standard deviation in ns, controlling the width

# --- Attenuation Levels ---
attenuation_db1 = -23
attenuation_db2 = -40

# --- Time Array ---
time_ns = np.linspace(0, 100, 1000)  # Time from 0 to 100 ns

# --- Generate the Input Signal ---
input_signal = gaussian_pulse(time_ns, peak_amplitude_mv, pulse_center_ns, std_deviation_ns)

# --- Calculate Attenuation Factors ---
scaling_factor1 = 10**(attenuation_db1 / 20)
scaling_factor2 = 10**(attenuation_db2 / 20)

# --- Generate the Attenuated Output Signals ---
output_signal1 = input_signal * scaling_factor1
output_signal2 = input_signal * scaling_factor2

# --- Plot the Signals ---
plt.figure(figsize=(10, 6))

plt.plot(time_ns, input_signal, label='Input Signal (Gaussian Pulse)', linewidth=2)
plt.plot(time_ns, output_signal1, label=f'Output Signal ({attenuation_db1} dB Attenuator)', linestyle='--')
plt.plot(time_ns, output_signal2, label=f'Output Signal ({attenuation_db2} dB Attenuator)', linestyle=':')

# --- Formatting the Plot ---
plt.title('Effect of Attenuators on a Gaussian Pulse')
plt.xlabel('Time (ns)')
plt.ylabel('Amplitude (mV)')
plt.legend()
plt.grid(True)
plt.show()