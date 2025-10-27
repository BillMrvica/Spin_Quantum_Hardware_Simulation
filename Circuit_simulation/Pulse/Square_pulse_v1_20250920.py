import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

# --- 1. Define the Analytical Fourier Transform of the Square Pulse ---
def fourier_transform_square_pulse(f, amplitude, duration):
    """
    Calculates the analytical Fourier transform of a square pulse starting at t=0.
    """
    # To avoid division by zero at f=0 for the sinc function, we handle it as a special case.
    with np.errstate(divide='ignore', invalid='ignore'):
        sinc_part = np.sin(np.pi * f * duration) / (np.pi * f * duration)
    sinc_part[f == 0] = 1.0  # Manually set the limit for f=0 where sinc(0)=1
    
    magnitude_part = amplitude * duration * sinc_part
    phase_part = np.exp(-1j * np.pi * f * duration)
    
    return magnitude_part * phase_part

# --- 2. Setup the Main Plot ---
# Create 3 subplots stacked vertically
fig, (ax_time, ax_mag, ax_phase) = plt.subplots(3, 1, figsize=(10, 9))
# Adjust layout to make room for sliders and prevent title overlap
plt.subplots_adjust(bottom=0.25, hspace=0.45)

# Initial parameters
initial_amplitude_mv = 100.0
initial_duration_ns = 20.0

# --- 3. Create the initial plot content ---
# Create empty line objects for each plot that will be updated
line_time, = ax_time.plot([], [], 'r-', lw=2)
line_mag, = ax_mag.plot([], [], 'b-')
line_phase, = ax_phase.plot([], [], 'g-')

# Configure Time-Domain Plot
ax_time.set_title('Time-Domain Signal', fontsize=14)
ax_time.set_xlabel('Time (ns)', fontsize=12)
ax_time.set_ylabel('Amplitude (mV)', fontsize=12)
ax_time.grid(True)

# Configure Magnitude Spectrum Plot
ax_mag.set_title('Double-Sided Magnitude Spectrum', fontsize=14)
ax_mag.set_ylabel('Spectral Density (µV / MHz)', fontsize=12)
ax_mag.grid(True)

# Configure Phase Spectrum Plot
ax_phase.set_title('Phase Spectrum', fontsize=14)
ax_phase.set_ylabel('Phase (degrees)', fontsize=12)
ax_phase.set_xlabel('Frequency (GHz)', fontsize=12)
ax_phase.grid(True)

# Add a text box for the DC value (will be updated)
dc_text = ax_mag.text(0.05, 0.95, '', transform=ax_mag.transAxes, verticalalignment='top',
                      bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

# --- 4. Create Slider Axes and Slider Widgets ---
ax_amp = plt.axes([0.25, 0.1, 0.65, 0.03])
ax_dur = plt.axes([0.25, 0.05, 0.65, 0.03])

amp_slider = Slider(
    ax=ax_amp, label='Amplitude (mV)', valmin=10, valmax=300,
    valinit=initial_amplitude_mv, valstep=5
)
dur_slider = Slider(
    ax=ax_dur, label='Duration (ns)', valmin=1, valmax=100,
    valinit=initial_duration_ns, valstep=1
)

# --- 5. Define the Update Function ---
def update(val):
    # Get current slider values
    amplitude_mv = amp_slider.val
    duration_ns = dur_slider.val
    
    # --- Update Time-Domain Plot ---
    # Create points for a clean square pulse shape for plotting
    t_pulse = np.array([-0.1 * duration_ns, 0, 0, duration_ns, duration_ns, 1.1 * duration_ns])
    y_pulse = np.array([0, 0, amplitude_mv, amplitude_mv, 0, 0])
    line_time.set_data(t_pulse, y_pulse)
    
    # --- Update Frequency-Domain Plots ---
    amplitude_v = amplitude_mv / 1000.0
    duration_s = duration_ns / 1e9
    
    # Recalculate frequency range and spectrum
    max_freq = 5 / duration_s if duration_s > 0 else 5e9
    f_new = np.linspace(-max_freq, max_freq, 4001)
    X_f_new = fourier_transform_square_pulse(f_new, amplitude_v, duration_s)
    
    # Update magnitude plot data
    line_mag.set_xdata(f_new / 1e9)
    line_mag.set_ydata(np.abs(X_f_new) * 1e-6) # in µV/MHz
    
    # Update phase plot data
    line_phase.set_xdata(f_new / 1e9)
    line_phase.set_ydata(np.unwrap(np.angle(X_f_new, deg=True)))
    
    # Update DC text
    dc_val = amplitude_v * duration_s * 1e-6
    dc_text.set_text(f'DC (f=0): {dc_val:.2f} µV/MHz')
    
    # Rescale all axes
    for ax in [ax_time, ax_mag, ax_phase]:
        ax.relim()
        ax.autoscale_view()
    
    # Redraw the figure
    fig.canvas.draw_idle()

# --- 6. Register the update function with the sliders ---
amp_slider.on_changed(update)
dur_slider.on_changed(update)

# --- 7. (Optional) Add a Reset Button ---
ax_reset = plt.axes([0.8, 0.0, 0.1, 0.04])
button = Button(ax_reset, 'Reset', hovercolor='0.975')

def reset(event):
    amp_slider.reset()
    dur_slider.reset()
button.on_clicked(reset)

# Call update once to draw the initial state
update(None)

# Show the plot
plt.show()