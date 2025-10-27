import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import kaiser
from scipy.fft import fft, fftshift, fftfreq

# --- Parameters ---
N = 512  # Number of points
T = 1.0  # Symbol period (for Raised Cosine)
t = np.linspace(-4*T, 4*T, N)

# --- 1. Generate Raised-Cosine Pulses ---
def raised_cosine(t, T, alpha):
    # Handle the t=0 and t = +/- T/(2*alpha) special cases to avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        sinc_part = np.sinc(t/T)
        cos_part = np.cos(np.pi * alpha * t / T)
        denom = 1 - (2 * alpha * t / T)**2
    
    # The value at t=0 is 1. np.sinc handles this. 
    # The value at the denominator singularity is pi/4 * sinc(1/(2*alpha))
    denom_zero_indices = np.where(np.abs(denom) < 1e-8)[0]
    for i in denom_zero_indices:
        cos_part[i] = (np.pi/4) * np.sinc(1/(2*alpha))
        denom[i] = 1 # Avoid division by zero by setting to 1, as we manually compute the value
        
    return sinc_part * cos_part / denom

rc_01 = raised_cosine(t, T, alpha=0.1)
rc_05 = raised_cosine(t, T, alpha=0.5)
rc_09 = raised_cosine(t, T, alpha=0.9)

# --- 2. Generate Kaiser Windows ---
# For a fair comparison, we'll use a subsection of the points
win_len = N // 2 
kaiser_2 = kaiser(win_len, beta=2)
kaiser_5 = kaiser(win_len, beta=5)
kaiser_8 = kaiser(win_len, beta=8)
# Pad with zeros to match the length of the RC pulses for FFT
kaiser_2_padded = np.pad(kaiser_2, (N - win_len) // 2)
kaiser_5_padded = np.pad(kaiser_5, (N - win_len) // 2)
kaiser_8_padded = np.pad(kaiser_8, (N - win_len) // 2)


# --- 3. Plotting ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Raised-Cosine Filter vs. Kaiser Window', fontsize=16)

# Time Domain Plots
axes[0, 0].set_title('Raised-Cosine Pulses (Time Domain)')
axes[0, 0].plot(t, rc_01, label=r'$\alpha=0.1$')
axes[0, 0].plot(t, rc_05, label=r'$\alpha=0.5$')
axes[0, 0].plot(t, rc_09, label=r'$\alpha=0.9$')
axes[0, 0].axvline(x=T, color='k', linestyle='--', alpha=0.5)
axes[0, 0].axvline(x=-T, color='k', linestyle='--', alpha=0.5)
axes[0, 0].text(T, 0.1, 'Zero ISI Point', ha='center')
axes[0, 0].grid(True), axes[0, 0].legend()
axes[0, 0].set_xlabel('Time (normalized to T)')

axes[1, 0].set_title('Kaiser Windows (Time Domain)')
axes[1, 0].plot(kaiser_2_padded, label=r'$\beta=2$')
axes[1, 0].plot(kaiser_5_padded, label=r'$\beta=5$')
axes[1, 0].plot(kaiser_8_padded, label=r'$\beta=8$')
axes[1, 0].grid(True), axes[1, 0].legend()
axes[1, 0].set_xlabel('Sample Number')


# Frequency Domain Plots
def get_db_spectrum(signal):
    Y = fftshift(fft(signal))
    Y_mag = np.abs(Y) / len(signal)
    Y_mag[Y_mag < 1e-15] = 1e-15
    return 20 * np.log10(Y_mag)

freqs = fftshift(fftfreq(N, d=(t[1]-t[0])))

axes[0, 1].set_title('Raised-Cosine Spectrum (Frequency Domain)')
axes[0, 1].plot(freqs, get_db_spectrum(rc_01), label=r'$\alpha=0.1$')
axes[0, 1].plot(freqs, get_db_spectrum(rc_05), label=r'$\alpha=0.5$')
axes[0, 1].plot(freqs, get_db_spectrum(rc_09), label=r'$\alpha=0.9$')
axes[0, 1].grid(True), axes[0, 1].legend()
axes[0, 1].set_ylim(-100, 5)
axes[0, 1].set_xlim(-2/T, 2/T)
axes[0, 1].set_xlabel('Frequency (normalized to 1/T)')


axes[1, 1].set_title('Kaiser Window Spectrum (Frequency Domain)')
axes[1, 1].plot(freqs, get_db_spectrum(kaiser_2_padded), label=r'$\beta=2$')
axes[1, 1].plot(freqs, get_db_spectrum(kaiser_5_padded), label=r'$\beta=5$')
axes[1, 1].plot(freqs, get_db_spectrum(kaiser_8_padded), label=r'$\beta=8$')
axes[1, 1].grid(True), axes[1, 1].legend()
axes[1, 1].set_ylim(-100, 5)
axes[1, 1].set_xlim(-2/T, 2/T)
axes[1, 1].set_xlabel('Frequency (normalized to 1/T)')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()