import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.constants import hbar, e
from matplotlib.widgets import Slider, Button

# --- 1. Define Physical Constants and Circuit Parameters ---
PHI_0 = np.pi * hbar / e  # More convenient form for the equations (h/2e = 2*pi*hbar/2e)

# Circuit Parameters
I_c = 100e-6   # Critical current: 100 µA
C_j = 0.5e-12  # Junction capacitance: 0.5 pF
R_j = 5.0      # Shunt resistance: 5 Ω
R_s = 50.0     # Source resistance: 50 Ω

# --- 2. Define the ODE for the Voltage-Driven RCSJ Model ---
def rcsj_voltage_driven_ode(t, y, I_c, C_j, R_j, R_s, V_in_func):
    """
    Represents the voltage-driven RCSJ model.
    y[0] = delta (phase)
    y[1] = d(delta)/dt (derivative of phase)
    """
    delta, delta_dot = y
    
    # Get the input voltage at the current time t
    V_in = V_in_func(t)
    
    # Constants for clarity
    V_j_per_delta_dot = PHI_0 / (2 * np.pi)
    
    # Rearranged second-order ODE derived from KCL
    delta_ddot = (1 / (C_j * V_j_per_delta_dot)) * (
        (V_in / R_s) - 
        (V_j_per_delta_dot * delta_dot * (1/R_s + 1/R_j)) - 
        I_c * np.sin(delta)
    )
    
    return [delta_dot, delta_ddot]

# --- 3. Setup the Interactive Plot ---
fig, (ax_v, ax_i) = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
plt.subplots_adjust(bottom=0.25)
fig.suptitle('Response of a Josephson Junction to a Voltage Pulse', fontsize=16)

# Create empty line objects
line_vin, = ax_v.plot([], [], 'r--', label='$V_{in}(t)$')
line_vout, = ax_v.plot([], [], 'b-', lw=2, label='$V_{out}(t)$ (across JJ)')
line_iin, = ax_i.plot([], [], 'g-', lw=2, label='$I_{in}(t)$ (Total Current)')

# Configure subplots
ax_v.set_title('Voltage Response')
ax_v.set_ylabel('Voltage (mV)'), ax_v.grid(True), ax_v.legend()
ax_i.set_title('Current Response')
ax_i.set_xlabel('Time (ns)'), ax_i.set_ylabel('Current (µA)'), ax_i.grid(True), ax_i.legend()

# --- 4. Create Slider Widgets ---
ax_amp = plt.axes([0.25, 0.1, 0.65, 0.03])
ax_dur = plt.axes([0.25, 0.05, 0.65, 0.03])

amp_slider = Slider(ax=ax_amp, label='Input Amp (mV)', valmin=10, valmax=200, valinit=50)
dur_slider = Slider(ax=ax_dur, label='Input Duration (ns)', valmin=1, valmax=100, valinit=20)

# --- 5. Define the Update Function ---
def update(val):
    amplitude_mv = amp_slider.val
    duration_ns = dur_slider.val
    
    amplitude_v = amplitude_mv / 1000.0
    duration_s = duration_ns / 1e9
    
    # Define the input voltage pulse (with smooth edges for numerical stability)
    t_rise = 0.1e-12 # 0.1 ps rise time
    def V_input_func(t):
        start_edge = 0.5 * (1 + np.tanh(t / t_rise))
        end_edge = 0.5 * (1 - np.tanh((t - duration_s) / t_rise))
        return amplitude_v * start_edge * end_edge

    # Simulation time span
    t_span = [-5 * duration_s, 5 * duration_s]
    t_eval = np.linspace(t_span[0], t_span[1], 2000)

    # Solve the ODE
    y0 = [0, 0] # Initial phase and voltage are zero
    sol = solve_ivp(
        fun=rcsj_voltage_driven_ode, t_span=t_span, y0=y0, t_eval=t_eval,
        args=(I_c, C_j, R_j, R_s, V_input_func), method='RK45'
    )
    
    # --- Post-process results ---
    delta_dot_t = sol.y[1]
    V_out_t = (PHI_0 / (2 * np.pi)) * delta_dot_t
    V_in_t = V_input_func(t_eval)
    I_in_t = (V_in_t - V_out_t) / R_s
    
    # --- Update plot data ---
    line_vin.set_data(t_eval * 1e9, V_in_t * 1000)
    line_vout.set_data(t_eval * 1e9, V_out_t * 1000)
    line_iin.set_data(t_eval * 1e9, I_in_t * 1e6)
    
    # Rescale axes
    for ax in [ax_v, ax_i]:
        ax.relim(), ax.autoscale_view()
    ax_v.set_ylim(-0.1*amplitude_mv, 1.1*amplitude_mv) # Keep y-axis stable
    
    fig.canvas.draw_idle()

# --- 6. Register callbacks and initialize ---
amp_slider.on_changed(update)
dur_slider.on_changed(update)

ax_reset = plt.axes([0.8, 0.0, 0.1, 0.04])
button = Button(ax_reset, 'Reset', hovercolor='0.975')
def reset(event):
    amp_slider.reset(), dur_slider.reset()
button.on_clicked(reset)

update(None)
plt.show()