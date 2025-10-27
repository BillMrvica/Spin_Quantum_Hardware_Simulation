import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fft, fftfreq, fftshift
from matplotlib.widgets import Slider, Button

# --- Setup the Interactive Plot ---
fig, axes = plt.subplots(3, 2, figsize=(15, 11), gridspec_kw={'width_ratios': [1.5, 1]})
plt.subplots_adjust(bottom=0.3, hspace=0.5, wspace=0.25)
fig.suptitle('Time and Frequency Domain Comparison of Pulse Shaping Filters', fontsize=16)

# Assign axes for clarity
(ax_time_rc, ax_freq_rc) = axes[0]
(ax_time_bessel, ax_freq_bessel) = axes[1]
(ax_time_gauss, ax_freq_gauss) = axes[2]

# --- Create empty line objects ---
line_in_rc, = ax_time_rc.plot([], [], 'r--', label='Input')
line_out_rc, = ax_time_rc.plot([], [], 'b-', label='RC Output', lw=2)
line_in_bessel, = ax_time_bessel.plot([], [], 'r--', label='Input')
line_out_bessel, = ax_time_bessel.plot([], [], 'g-', label='Bessel Output', lw=2)
line_in_gauss, = ax_time_gauss.plot([], [], 'r--', label='Input')
line_out_gauss, = ax_time_gauss.plot([], [], 'm-', label='Gaussian Output', lw=2)
line_freq_in_rc, = ax_freq_rc.plot([], [], color='gray', linestyle='--', label='Input Spectrum')
line_freq_out_rc, = ax_freq_rc.plot([], [], 'b-', label='RC Spectrum', lw=2)
line_freq_in_bessel, = ax_freq_bessel.plot([], [], color='gray', linestyle='--', label='Input Spectrum')
line_freq_out_bessel, = ax_freq_bessel.plot([], [], 'g-', label='Bessel Spectrum', lw=2)
line_freq_in_gauss, = ax_freq_gauss.plot([], [], color='gray', linestyle='--', label='Input Spectrum')
line_freq_out_gauss, = ax_freq_gauss.plot([], [], 'm-', label='Gaussian Spectrum', lw=2)

# --- Configure subplots with ROW TITLES ---
ax_time_rc.set_title('3rd-Order Cascaded RC Filter', fontsize=14, loc='left')
ax_time_rc.set_ylabel('Amplitude (mV)'), ax_time_rc.grid(True), ax_time_rc.legend()
ax_freq_rc.set_title('Frequency Spectrum', fontsize=12)
ax_freq_rc.set_ylabel('Magnitude (dB)'), ax_freq_rc.grid(True, which='both'), ax_freq_rc.legend()
ax_time_bessel.set_title('5th-Order LC Bessel Filter', fontsize=14, loc='left')
ax_time_bessel.set_ylabel('Amplitude (mV)'), ax_time_bessel.grid(True), ax_time_bessel.legend()
ax_freq_bessel.set_title('Frequency Spectrum', fontsize=12)
ax_freq_bessel.set_ylabel('Magnitude (dB)'), ax_freq_bessel.grid(True, which='both'), ax_freq_bessel.legend()
ax_time_gauss.set_title('5th-Order LC Gaussian Filter', fontsize=14, loc='left')
ax_time_gauss.set_xlabel('Time (ns)'), ax_time_gauss.set_ylabel('Amplitude (mV)')
ax_time_gauss.grid(True), ax_time_gauss.legend()
ax_freq_gauss.set_title('Frequency Spectrum', fontsize=12)
ax_freq_gauss.set_xlabel('Frequency (MHz)'), ax_freq_gauss.set_ylabel('Magnitude (dB)')
ax_freq_gauss.grid(True, which='both'), ax_freq_gauss.legend()

# --- Create Slider Widgets ---
ax_amp = plt.axes([0.25, 0.15, 0.65, 0.03])
ax_dur = plt.axes([0.1, 0.1, 0.8, 0.03])
ax_fc = plt.axes([0.25, 0.05, 0.65, 0.03])

amp_slider = Slider(ax=ax_amp, label='Pulse Amp (mV)', valmin=0, valmax=300, valinit=100, valstep=5)
dur_slider = Slider(ax=ax_dur, label='Pulse Duration (ns)', valmin=1, valmax=100, valinit=50, valstep=1)
fc_slider = Slider(ax=ax_fc, label='Filter Cutoff (MHz)', valmin=1, valmax=50, valinit=7.5, valstep=0.5)

# --- Define the Update Function ---
def update(val):
    amplitude_mv, duration_ns, cutoff_mhz = amp_slider.val, dur_slider.val, fc_slider.val
    
    f_c_hz, w_c = cutoff_mhz * 1e6, 2 * np.pi * cutoff_mhz * 1e6
    R, C = 50.0, 1 / (w_c * 50.0)
    num1, den1 = [1], [R*C, 1]
    system_rc = signal.lti([1], np.convolve(np.convolve(den1, den1), den1))
    system_bessel = signal.lti(*signal.bessel(5, w_c, 'low', analog=True))
    poles_norm = [-0.70494+1.20012j, -0.70494-1.20012j, -1.07128+0.71833j, -1.07128-0.71833j, -1.18737]
    a_gauss = np.poly(np.array(poles_norm) * w_c)
    system_gauss = signal.lti([a_gauss[-1]], a_gauss)

    amplitude_v, duration_s = amplitude_mv / 1000.0, duration_ns / 1e9
    t_duration = max(5 * duration_s, 400e-9)
    fs = 50e9
    t = np.arange(0, t_duration, 1/fs)
    pulse_start, pulse_end = 50e-9, 50e-9 + duration_s
    square_wave_sim = np.zeros_like(t)
    square_wave_sim[(t >= pulse_start) & (t < pulse_end)] = amplitude_v

    _, y_out_rc, _ = signal.lsim(system_rc, square_wave_sim, t)
    _, y_out_bessel, _ = signal.lsim(system_bessel, square_wave_sim, t)
    _, y_out_gauss, _ = signal.lsim(system_gauss, square_wave_sim, t)

    t_plot_sq = np.array([0, pulse_start*1e9, pulse_start*1e9, pulse_end*1e9, pulse_end*1e9, t_duration*1e9])
    y_plot_sq = np.array([0, 0, amplitude_mv, amplitude_mv, 0, 0])
    
    for line in [line_in_rc, line_in_bessel, line_in_gauss]: line.set_data(t_plot_sq, y_plot_sq)
    line_out_rc.set_data(t * 1e9, y_out_rc * 1000)
    line_out_bessel.set_data(t * 1e9, y_out_bessel * 1000)
    line_out_gauss.set_data(t * 1e9, y_out_gauss * 1000)

    # --- CORRECTED SECTION ---
    # Define n FIRST, then use it to define freqs
    n = len(t)
    freqs = fftshift(fftfreq(n, d=1/fs))
    
    def get_db_spectrum(time_signal):
        Y = fftshift(fft(time_signal))
        magnitude = np.abs(Y) / n
        magnitude[magnitude < 1e-15] = 1e-15
        return 20 * np.log10(magnitude)

    Y_in_db, Y_rc_db = get_db_spectrum(square_wave_sim), get_db_spectrum(y_out_rc)
    Y_bessel_db, Y_gauss_db = get_db_spectrum(y_out_bessel), get_db_spectrum(y_out_gauss)

    for line in [line_freq_in_rc, line_freq_in_bessel, line_freq_in_gauss]: line.set_data(freqs / 1e6, Y_in_db)
    line_freq_out_rc.set_data(freqs / 1e6, Y_rc_db)
    line_freq_out_bessel.set_data(freqs / 1e6, Y_bessel_db)
    line_freq_out_gauss.set_data(freqs / 1e6, Y_gauss_db)

    for ax in [ax_time_rc, ax_time_bessel, ax_time_gauss]:
        ax.relim(), ax.autoscale_view()
        ax.set_ylim(min(-10, -0.1*amplitude_mv), max(310, 1.1*amplitude_mv))
    
    for ax in [ax_freq_rc, ax_freq_bessel, ax_freq_gauss]:
        ax.relim(), ax.autoscale_view()
        ax.set_xlim(-5 * cutoff_mhz, 5 * cutoff_mhz)
        ax.set_ylim(-120, max(Y_in_db) + 5)

    fig.canvas.draw_idle()

# --- Register callbacks and initialize ---
amp_slider.on_changed(update)
dur_slider.on_changed(update)
fc_slider.on_changed(update)

ax_reset = plt.axes([0.8, 0.0, 0.1, 0.04])
button = Button(ax_reset, 'Reset', hovercolor='0.975')
def reset(event): amp_slider.reset(), dur_slider.reset(), fc_slider.reset()
button.on_clicked(reset)

update(None)
plt.show()