import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftshift, fftfreq
from matplotlib.widgets import Slider, Button

# --- 1. Raised-Cosine Pulse Generation Function ---
def raised_cosine(t, T, alpha, amplitude):
    """
    Generates a raised-cosine pulse.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        sinc_part = np.sinc(t / T)
        cos_part = np.cos(np.pi * alpha * t / T)
        denom = 1 - (2 * alpha * t / T)**2
    
    denom_zero_indices = np.where(np.abs(denom) < 1e-8)[0]
    if denom_zero_indices.size > 0:
        val_at_singularity = (np.pi / 4) * np.sinc(1 / (2 * alpha))
        cos_part[denom_zero_indices] = val_at_singularity
        denom[denom_zero_indices] = 1.0
        
    return amplitude * sinc_part * cos_part / denom

# --- 2. Setup the Interactive Plot ---
fig, (ax_time, ax_freq) = plt.subplots(1, 2, figsize=(15, 7))
plt.subplots_adjust(bottom=0.3)
fig.suptitle('Superposition of Two Raised-Cosine Pulses: y(t) = x(t) + x(t+t_d)', fontsize=16)

# --- Create empty line objects for plots ---
line_x1, = ax_time.plot([], [], 'k--', alpha=0.5, label='x(t)')
line_x2, = ax_time.plot([], [], 'gray', linestyle='--', alpha=0.5, label='x(t+t_d)')
line_y, = ax_time.plot([], [], 'r-', lw=2, label='y(t)')
line_freq_x, = ax_freq.plot([], [], color='gray', linestyle='--', label='Spectrum of x(t)')
line_freq_y, = ax_freq.plot([], [], 'b-', label='Spectrum of y(t)', lw=2)

# --- Configure subplots ---
ax_time.set_title('Time Domain')
ax_time.set_xlabel('Time (ns)'), ax_time.set_ylabel('Amplitude (mV)')
ax_time.grid(True), ax_time.legend()
ax_freq.set_title('Frequency Domain')
ax_freq.set_xlabel('Frequency (MHz)'), ax_freq.set_ylabel('Magnitude (dB)')
ax_freq.grid(True, which='both'), ax_freq.legend()
ax_freq.set_xlim(-1000, 1000)

# --- 3. Create Slider Widgets ---
ax_amp = plt.axes([0.25, 0.15, 0.65, 0.03])
ax_dur = plt.axes([0.25, 0.1, 0.65, 0.03])
ax_td = plt.axes([0.25, 0.05, 0.65, 0.03])

amp_slider = Slider(ax=ax_amp, label='Amplitude (mV)', valmin=10, valmax=200, valinit=100)
dur_slider = Slider(ax=ax_dur, label='Duration T (ns)', valmin=5, valmax=100, valinit=20)
td_slider = Slider(ax=ax_td, label='Shift t_d (ns)', valmin=1, valmax=100, valinit=10)

# --- 4. Define the Update Function ---
def update(val):
    amplitude_mv = amp_slider.val
    duration_ns = dur_slider.val
    td_ns = td_slider.val
    
    amplitude_v = amplitude_mv / 1000.0
    duration_s = duration_ns / 1e9
    td_s = td_ns / 1e9
    
    # --- Time Domain Calculation ---
    t_span = max(6 * duration_s + td_s, 500e-9) 
    t_center = -td_s / 2
    fs = 200e9
    n_points = int(t_span * fs)
    if n_points % 2 != 0: n_points += 1
    
    t = np.linspace(t_center - t_span/2, t_center + t_span/2, n_points)
    
    alpha = 0.5
    pulse_x1 = raised_cosine(t, duration_s, alpha, amplitude_v)
    pulse_x2 = raised_cosine(t + td_s, duration_s, alpha, amplitude_v)
    pulse_y = pulse_x1 + pulse_x2
    
    line_x1.set_data(t * 1e9, pulse_x1 * 1000)
    line_x2.set_data(t * 1e9, pulse_x2 * 1000)
    line_y.set_data(t * 1e9, pulse_y * 1000)

    # --- Frequency Domain Calculation ---
    freqs = fftshift(fftfreq(n_points, d=1/fs))
    def get_db_spectrum(time_signal):
        Y = fftshift(fft(time_signal))
        magnitude = np.abs(Y) / n_points
        magnitude[magnitude < 1e-15] = 1e-15
        return 20 * np.log10(magnitude)
        
    spectrum_x = get_db_spectrum(pulse_x1)
    spectrum_y = get_db_spectrum(pulse_y)
    
    line_freq_x.set_data(freqs / 1e6, spectrum_x)
    line_freq_y.set_data(freqs / 1e6, spectrum_y)

    # --- Rescale Axes ---
    # Autoscale the time-domain plot, then set a stable y-limit
    ax_time.relim()
    ax_time.autoscale_view()
    ax_time.set_ylim(-0.1 * amplitude_mv, 1.5 * amplitude_mv)

    # **CORRECTED SECTION**
    # The erroneous autoscale_view call is removed.
    # The set_ylim command correctly handles the dynamic y-axis for the frequency plot.
    ax_freq.set_ylim(-100, max(spectrum_y) + 5)
    
    fig.canvas.draw_idle()

# --- 5. Register callbacks and initialize ---
amp_slider.on_changed(update)
dur_slider.on_changed(update)
td_slider.on_changed(update)

ax_reset = plt.axes([0.8, 0.0, 0.1, 0.04])
button = Button(ax_reset, 'Reset', hovercolor='0.975')
def reset(event):
    amp_slider.reset(), dur_slider.reset(), td_slider.reset()
button.on_clicked(reset)

update(None)
plt.show()