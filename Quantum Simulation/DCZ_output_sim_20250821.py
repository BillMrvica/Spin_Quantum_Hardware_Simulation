import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider

# --- PART 1: FUNCTIONS TO DRAW THE QUANTUM CIRCUIT ---

def draw_gate(ax, x, y, width, height, label, color='skyblue'):
    """Draws a generic quantum gate box."""
    box = patches.Rectangle((x - width / 2, y - height / 2), width, height,
                            facecolor=color, edgecolor='black', linewidth=1.5, zorder=3)
    ax.add_patch(box)
    ax.text(x, y, label, ha='center', va='center', fontsize=12, zorder=4)

def draw_cphase_gate(ax, x, y1, y2, width, label='CPhase\n(t)'):
    """Draws the two-qubit CPhase gate block."""
    height = abs(y1 - y2) + 0.5
    center_y = (y1 + y2) / 2
    box = patches.Rectangle((x - width / 2, center_y - height / 2), width, height,
                            facecolor='lightcoral', edgecolor='black', linewidth=1.5, zorder=2)
    ax.add_patch(box)
    ax.text(x, center_y, label, ha='center', va='center', fontsize=12, zorder=4)

def draw_measurement(ax, x, y, width, height):
    """Draws the measurement symbol."""
    box = patches.Rectangle((x - width / 2, y - height / 2), width, height,
                            facecolor='lightgray', edgecolor='black', linewidth=1.5, zorder=3)
    ax.add_patch(box)
    arc = patches.Arc((x, y - height * 0.1), width * 0.6, height * 0.6,
                      theta1=180, theta2=0, color='black', linewidth=1.5, zorder=4)
    ax.add_patch(arc)
    arrow_angle = np.deg2rad(45)
    ax.plot([x, x + width * 0.3 * np.cos(arrow_angle)],
            [y - height * 0.1, y - height * 0.1 + height * 0.3 * np.sin(arrow_angle)],
            color='black', lw=1.5, zorder=4)

# --- Main Figure and Subplot Setup ---
fig = plt.figure(figsize=(10, 15))
ax_circuit = fig.add_subplot(3, 1, 1)
ax_sphere = fig.add_subplot(3, 1, 2, projection='3d')
ax_prob = fig.add_subplot(3, 1, 3)
plt.subplots_adjust(bottom=0.18, hspace=0.45) # Make room for plots and sliders

# --- PART 2: DRAW THE CIRCUIT DIAGRAM (ax_circuit) ---
y_q0, y_q1 = 2, 1
ax_circuit.plot([0, 9], [y_q0, y_q0], color='black', zorder=1)
ax_circuit.plot([0, 9], [y_q1, y_q1], color='black', zorder=1)
ax_circuit.text(-0.5, y_q0, r'$|q_0\rangle$', ha='center', va='center', fontsize=16)
ax_circuit.text(-0.5, y_q1, r'$|q_1\rangle$', ha='center', va='center', fontsize=16)

gate_w, gate_h = 0.8, 0.8
x_coords = {'start_x2': 1, 'cphase1': 2.5, 'x_pi': 4, 'cphase2': 5.5, 'end_x2': 7, 'measure': 8.2}

draw_gate(ax_circuit, x_coords['start_x2'], y_q0, gate_w, gate_h, r'$X(\frac{\pi}{2})$')
draw_cphase_gate(ax_circuit, x_coords['cphase1'], y_q0, y_q1, width=1.5)
ax_circuit.annotate('', xy=(1.75, 2.7), xytext=(3.25, 2.7), arrowprops=dict(arrowstyle='<->', lw=1.5))
ax_circuit.text(x_coords['cphase1'], 2.9, r'$t_0/2$', ha='center', fontsize=12)
draw_gate(ax_circuit, x_coords['x_pi'], y_q0, gate_w, gate_h, r'$X(\pi)$', color='lightgreen')
draw_gate(ax_circuit, x_coords['x_pi'], y_q1, gate_w, gate_h, r'$X(\pi)$', color='lightgreen')
draw_cphase_gate(ax_circuit, x_coords['cphase2'], y_q0, y_q1, width=1.5)
ax_circuit.annotate('', xy=(4.75, 2.7), xytext=(6.25, 2.7), arrowprops=dict(arrowstyle='<->', lw=1.5))
ax_circuit.text(x_coords['cphase2'], 2.9, r'$t_0/2$', ha='center', fontsize=12)
draw_gate(ax_circuit, x_coords['end_x2'], y_q0, gate_w, gate_h, r'$X(\frac{\pi}{2})$')
draw_measurement(ax_circuit, x_coords['measure'], y_q0, gate_w, gate_h)

ax_circuit.set_xlim(0, 9)
ax_circuit.set_ylim(0, 3.5)
ax_circuit.set_title('1. Quantum Circuit for Decoupled CPhase Gate', fontsize=16)
ax_circuit.set_aspect('equal')
ax_circuit.axis('off')

# --- PART 3: SETUP BLOCH SPHERE and PROBABILITY PLOTS (ax_sphere, ax_prob) ---

# --- Bloch Sphere (Static Elements) ---
ax_sphere.set_box_aspect([1, 1, 1])
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x_sphere = np.outer(np.cos(u), np.sin(v))
y_sphere = np.outer(np.sin(u), np.sin(v))
z_sphere = np.outer(np.ones(np.size(u)), np.cos(v))
ax_sphere.plot_wireframe(x_sphere, y_sphere, z_sphere, color='gray', alpha=0.2, rstride=10, cstride=10)
ax_sphere.quiver(0, 0, 0, 1.3, 0, 0, color='black', arrow_length_ratio=0.08)
ax_sphere.quiver(0, 0, 0, 0, 1.3, 0, color='black', arrow_length_ratio=0.08)
ax_sphere.quiver(0, 0, 0, 0, 0, 1.3, color='black', arrow_length_ratio=0.08)
ax_sphere.text(1.4, 0, 0, 'X', fontsize=12)
ax_sphere.text(0, 1.4, 0, 'Y', fontsize=12)
ax_sphere.text(0, 0, 1.4, r'$|0\rangle$ (Z)', fontsize=14)
ax_sphere.set_axis_off()
ax_sphere.set_title('2. Interactive Qubit State Visualization', fontsize=16)

# --- Probability Plot (Static Elements) ---
T2_star = 50.0  # Dephasing time in μs
t_max = 150     # Max time in μs to show the decay effect
t_space = np.linspace(0, t_max, 500)
J1 = 20; J2 = 5 # kHz

# Define the decaying probability function
def decaying_prob(t, J_khz):
    J_mhz = J_khz / 1000.0
    oscillation = 0.5 * np.cos(np.pi * J_mhz * t)
    decay = np.exp(-(t / T2_star)**2)
    return 0.5 + oscillation * decay

ax_prob.plot(t_space, decaying_prob(t_space, J1), label=f'J = {J1} kHz', color='cyan')
ax_prob.plot(t_space, decaying_prob(t_space, J2), label=f'J = {J2} kHz', color='orange')
ax_prob.set_xlabel('Evolution Time $t_0$ (μs)', fontsize=12)
ax_prob.set_ylabel('Spin-up Probability P($|0\\rangle$)', fontsize=12)
ax_prob.set_title('3. Probability Evolution with Dephasing', fontsize=16)
ax_prob.set_ylim(-0.1, 1.1)
ax_prob.set_xlim(0, t_max)
ax_prob.grid(True)

# --- PART 4: INTERACTIVE ELEMENTS AND UPDATE LOGIC ---
J_init = 25.0; t0_init = 0.0

def calculate_coords_and_prob(J_khz, t0_us):
    J_mhz = J_khz / 1000.0
    theta = np.pi * J_mhz * t0_us
    x = np.sin(theta); y = 0; z = np.cos(theta)
    prob = decaying_prob(t0_us, J_khz) # Use the new decaying function
    return x, y, z, prob

x_init, y_init, z_init, _ = calculate_coords_and_prob(J_init, t0_init)
quiver_handle = [ax_sphere.quiver(0, 0, 0, x_init, y_init, z_init, color='red', lw=3, arrow_length_ratio=0.1)]

vline = ax_prob.axvline(t0_init, color='red', linestyle='--', lw=2)
current_J_curve, = ax_prob.plot(t_space, decaying_prob(t_space, J_init), color='red', lw=2, label=f'Current J = {J_init:.1f} kHz')
point, = ax_prob.plot(t0_init, decaying_prob(t0_init, J_init), 'ro', markersize=8)
ax_prob.legend()

# --- Sliders ---
ax_j = plt.axes([0.25, 0.08, 0.65, 0.02])
ax_t0 = plt.axes([0.25, 0.04, 0.65, 0.02])
j_slider = Slider(ax=ax_j, label='$J_{res}$ (kHz)', valmin=0, valmax=50, valinit=J_init)
t0_slider = Slider(ax=ax_t0, label='$t_0$ (μs)', valmin=0, valmax=t_max, valinit=t0_init)

# --- Update Function ---
def update(val):
    J = j_slider.val
    t0 = t0_slider.val
    
    # Update Bloch Sphere
    quiver_handle[0].remove()
    x_new, y_new, z_new, prob_new = calculate_coords_and_prob(J, t0)
    quiver_handle[0] = ax_sphere.quiver(0, 0, 0, x_new, y_new, z_new, color='red', lw=3, arrow_length_ratio=0.1)
    
    # Update Probability Plot
    vline.set_xdata([t0, t0])
    current_J_curve.set_ydata(decaying_prob(t_space, J))
    current_J_curve.set_label(f'Current J = {J:.1f} kHz')
    ax_prob.legend()
    point.set_data(t0, prob_new)
    
    fig.canvas.draw_idle()

j_slider.on_changed(update)
t0_slider.on_changed(update)

plt.show()
