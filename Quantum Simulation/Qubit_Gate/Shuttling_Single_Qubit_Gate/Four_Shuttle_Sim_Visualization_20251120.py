import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg
import os

# --- Constants ---
DEG = np.pi/180
s0 = np.eye(2)
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]])

# --- Initial Parameters ---
params = {
    't_q20': -2.28847019,
    't_q30': -1.54215921,
    't_q2_res': 1.1601954,
    'theta_in': 44.71842756,
    'A': 0.84173585,
    'B': 0.08989211,
    'fq2': 0.070926,
    'fq3': 0.06203767
}

# ==========================================
# Part 1: Physics Functions
# ==========================================

def get_hamiltonians_and_axis(p):
    theta, phi = DEG * p['theta_in'], 0
    qax = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    qx, qy, qz = qax
    
    H_q3 = np.array([[qz, qx-1j*qy], [qx+1j*qy, -qz]]) * p['fq3'] * 0.5
    H_q2 = np.array([[1, 0], [0, -1]]) * p['fq2'] * 0.5
    return H_q2, H_q3, qax

def func_U_matrix(t_q2, t_q3, p):
    H_q2, H_q3, _ = get_hamiltonians_and_axis(p)
    d_res = p['t_q2_res']
    d_3 = t_q3 - p['t_q30']
    d_2 = t_q2 - p['t_q20']
    
    U_res = scipy.linalg.expm(-2*np.pi*1j * H_q2 * d_res)
    U_3   = scipy.linalg.expm(-2*np.pi*1j * H_q3 * d_3)
    U_2   = scipy.linalg.expm(-2*np.pi*1j * H_q2 * d_2)
    
    U_single = U_res @ U_3 @ U_2 @ U_3 @ U_res
    return U_single

def func2_SU2decompose(mat):
    u0 = np.real(np.trace(mat @ s0)/2)
    ux = np.imag(np.trace(mat @ sx)/2)
    uy = np.imag(np.trace(mat @ sy)/2)
    uz = np.imag(np.trace(mat @ sz)/2)
    uvec = np.array([ux, uy, uz])
    norm = np.sqrt(np.sum(uvec**2))
    if norm <= 1e-14: uvec = np.array([0,0,1])
    else: uvec = uvec/norm
    rot_angle = 2*np.arccos(u0)
    gate_polar_angle = np.arctan2(np.sqrt(uvec[0]**2 + uvec[1]**2), uvec[2]) 
    return gate_polar_angle, rot_angle

def calculate_trajectory_segments(t2_in, t3_in, lab_time_limit, p):
    H_q2, H_q3, _ = get_hamiltonians_and_axis(p)
    d_res = max(0, p['t_q2_res'])
    d_3   = max(0, t3_in - p['t_q30'])
    d_2   = max(0, t2_in - p['t_q20'])
    d_tau = 1.0 / p['fq2'] if p['fq2'] > 1e-9 else 0
    
    # Structure: (Duration, Hamiltonian, Dot_ID)
    # Dot_ID: 0 = Dot 2 (Idle/Precession), 1 = Dot 3 (Pulse)
    pulse_structure = [
        (d_res, H_q2, 0), (d_3, H_q3, 1), (d_2, H_q2, 0), (d_3, H_q3, 1), (d_res, H_q2, 0)
    ]
    sequence = []
    sequence.extend(pulse_structure)
    sequence.append((d_tau, H_q2, 0))
    sequence.extend(pulse_structure)
    sequence.append((999.0, H_q2, 0)) 
    
    psi = np.array([0, 1], dtype=complex)
    plot_segments = []
    time_remaining = lab_time_limit
    
    # Colors
    c_dot2 = '#1f77b4' # Blue
    c_dot3 = '#d62728' # Red
    
    for dur, H, dot_id in sequence:
        if time_remaining <= 0: break
        
        # Logic for Sim Time
        if time_remaining >= dur:
            sim_time = dur
            is_completed = True # This segment is fully finished in the past
        else:
            sim_time = time_remaining
            is_completed = False # This segment is currently happening (active)
            
        steps = int(max(5, 40 * (sim_time/5.0)))
        dt_arr = np.linspace(0, sim_time, steps+1)
        xs, ys, zs = [], [], []
        start_psi = psi.copy()
        for dt in dt_arr:
            U = scipy.linalg.expm(-1j * 2 * np.pi * H * dt)
            psi_t = U @ start_psi
            rho = np.outer(psi_t, psi_t.conj())
            xs.append(np.real(np.trace(rho @ sx)))
            ys.append(np.real(np.trace(rho @ sy)))
            zs.append(np.real(np.trace(rho @ sz)))
        
        # Determine Line Style
        # Requirement: 
        # 1. Precession (Dot 2, id=0) -> ALWAYS Dashed ('--')
        # 2. Rotation (Dot 3, id=1):
        #    - If completed -> Dashed ('--')
        #    - If active (currently drawing) -> Solid ('-')
        
        if dot_id == 0:
            linestyle = '--' 
        else:
            if is_completed:
                linestyle = '--'
            else:
                linestyle = '-'
        
        plot_segments.append({
            'x': np.array(xs), 'y': np.array(ys), 'z': np.array(zs),
            'c': c_dot3 if dot_id == 1 else c_dot2,
            'ls': linestyle
        })
        psi = scipy.linalg.expm(-1j * 2 * np.pi * H * sim_time) @ start_psi
        time_remaining -= sim_time
        
    last_pt = (plot_segments[-1]['x'][-1], plot_segments[-1]['y'][-1], plot_segments[-1]['z'][-1]) if plot_segments else (0,0,-1)
    return plot_segments, last_pt

def get_angle_plot_data(lx, ly, lz):
    """Calculates arcs for theta and phi visualization."""
    # Radius for arcs
    r_theta = 0.4
    r_phi = 0.4
    
    # 1. Phi Arc (XY Plane)
    phi = np.arctan2(ly, lx)
    t_phi = np.linspace(0, phi, 40)
    px = r_phi * np.cos(t_phi)
    py = r_phi * np.sin(t_phi)
    pz = np.zeros_like(px)
    
    # 2. Theta Arc (Vertical Plane defined by phi)
    # Angle from Z axis (0) to current theta
    theta = np.arccos(np.clip(lz, -1, 1))
    t_theta = np.linspace(0, theta, 40)
    # Parametric equation for circle in vertical plane rotated by phi
    tx = r_theta * np.sin(t_theta) * np.cos(phi)
    ty = r_theta * np.sin(t_theta) * np.sin(phi)
    tz = r_theta * np.cos(t_theta)
    
    # Label positions (at the midpoint of the arcs)
    phi_mid = phi / 2
    label_phi_pos = (r_phi*1.2 * np.cos(phi_mid), r_phi*1.2 * np.sin(phi_mid), 0)
    
    theta_mid = theta / 2
    label_theta_pos = (r_theta*1.2 * np.sin(theta_mid) * np.cos(phi), 
                       r_theta*1.2 * np.sin(theta_mid) * np.sin(phi), 
                       r_theta*1.2 * np.cos(theta_mid))
    
    return (px, py, pz), (tx, ty, tz), label_phi_pos, label_theta_pos


# ==========================================
# Part 2: Visualization Layout
# ==========================================

fig = plt.figure(figsize=(17, 9))
gs = fig.add_gridspec(2, 3, height_ratios=[0.7, 1.3], width_ratios=[0.25, 1, 1], bottom=0.22, wspace=0.3)

ax_inputs_dummy = fig.add_subplot(gs[:, 0]) 
ax_inputs_dummy.axis('off')
ax_img   = fig.add_subplot(gs[0, 1]) 
ax_map   = fig.add_subplot(gs[1, 1]) 
ax_bloch = fig.add_subplot(gs[:, 2], projection='3d') 

# --- 1. Image Panel ---
ax_img.axis('off')
image_filename = "C:\\Users\\Administrator\\Desktop\\Spin_Quantum_Hardware_Simulation\\Quantum Simulation\\Qubit_Gate\\Shuttling_Single_Qubit_Gate\\schematic.png" 
if os.path.exists(image_filename):
    img = mpimg.imread(image_filename)
    ax_img.imshow(img)
    ax_img.set_title("Pulse Sequence Schematic")
else:
    ax_img.text(0.5, 0.5, "Image not found", ha='center', va='center', color='red')

# --- 2. Heatmap Setup ---
map_res = 150 
t_q2s = np.linspace(0, 50, map_res)
t_q3s = np.linspace(0, 50, map_res)

im_obj = None
contour_objs = []
crosshair, = ax_map.plot([], [], 'r+', markersize=20, markeredgewidth=3, zorder=10)

def recalculate_heatmap(event=None):
    global im_obj, contour_objs
    print("Recalculating map...")
    
    for c in contour_objs:
        if hasattr(c, 'collections'):
            for coll in c.collections: 
                try: coll.remove()
                except ValueError: pass 
        else:
            try: c.remove()
            except ValueError: pass
    contour_objs = []
    
    Pups = np.zeros((map_res, map_res))
    GatePol = np.zeros((map_res, map_res))
    RotAng = np.zeros((map_res, map_res))
    
    for i in range(map_res):
        for j in range(map_res):
            t2 = t_q2s[i]
            t3 = t_q3s[j]
            U_single = func_U_matrix(t2, t3, params)
            gp, ra = func2_SU2decompose(U_single)
            stateN = U_single @ U_single @ np.array([0,1])
            p_val = (np.real(np.abs(stateN[0]))**2 ) * params['A'] + params['B']
            Pups[j, i] = p_val 
            GatePol[j, i] = gp
            RotAng[j, i] = ra

    if im_obj is None:
        im_obj = ax_map.pcolormesh(t_q2s, t_q3s, Pups, shading='auto', cmap='viridis')
        cbar = fig.colorbar(im_obj, ax=ax_map)
        cbar.set_label('Spin-Up Probability')
    else:
        im_obj.set_array(Pups.ravel())
        im_obj.set_clim(np.min(Pups), np.max(Pups))
    
    c1 = ax_map.contour(t_q2s, t_q3s, GatePol/DEG, levels=[90], colors='C9', linewidths=1.5)
    c2 = ax_map.contour(t_q2s, t_q3s, RotAng/DEG, levels=[90], colors='C1', linewidths=1.5)
    c3 = ax_map.contour(t_q2s, t_q3s, RotAng/DEG, levels=[270], colors='C0', linewidths=1.5)
    contour_objs.extend([c1, c2, c3])
    
    ax_map.set_title(f"Map ($f_{{q2}}={params['fq2']:.3f}, \\theta={params['theta_in']:.1f}^\circ$)")
    update_all(0)
    fig.canvas.draw_idle()
    print("Calculation done.")

ax_map.set_xlabel(r'$t_2$ (ns)')
ax_map.set_ylabel(r'$t_3$ (ns)')
ax_map.legend(handles=[
    Line2D([0], [0], color='C9', lw=1.5, label=r'Polar $90^\circ$'),
    Line2D([0], [0], color='C1', lw=1.5, label=r'Rot $90^\circ$'),
    Line2D([0], [0], color='C0', lw=1.5, label=r'Rot $270^\circ$'),
    Line2D([0], [0], marker='+', color='r', lw=0, markersize=10, label='Point')
], loc='upper right', fontsize='x-small', framealpha=0.8)

# --- 3. Bloch Sphere Setup ---
u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:15j]
x_sphere, y_sphere, z_sphere = np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)
ax_bloch.plot_wireframe(x_sphere, y_sphere, z_sphere, color="gray", alpha=0.1, linewidth=0.5)

# Colored Axes
ax_bloch.plot([-1.2, 1.2], [0, 0], [0, 0], color='r', linestyle='-', lw=1.5) # X
ax_bloch.text(1.3, 0, 0, "x", color='r', fontsize=12, fontweight='bold')

ax_bloch.plot([0, 0], [-1.2, 1.2], [0, 0], color='g', linestyle='-', lw=1.5) # Y
ax_bloch.text(0, 1.3, 0, "y", color='g', fontsize=12, fontweight='bold')

ax_bloch.plot([0, 0], [0, 0], [-1.2, 1.2], color='b', linestyle='-', lw=1.5) # Z
ax_bloch.text(0, 0, 1.3, "z $|0\\rangle$", color='b', fontsize=12, fontweight='bold')
ax_bloch.text(0, 0, -1.4, "$|1\\rangle$", color='b', fontsize=12)

ax_bloch.set_axis_off()
ax_bloch.set_title('Interactive Bloch Sphere')
ax_bloch.set_aspect('equal')

# Plot Objects
q_axis_arrow = ax_bloch.quiver(0, 0, 0, 0, 0, 1, color='purple', linestyle='-.') 
lines_3d = [ax_bloch.plot([], [], [], linewidth=2, alpha=0.9)[0] for _ in range(25)]
state_arrow = None
proj_line, = ax_bloch.plot([], [], [], 'k:', lw=1, alpha=0.6) # Projection line for phi
arc_phi, = ax_bloch.plot([], [], [], color='k', lw=1)
arc_theta, = ax_bloch.plot([], [], [], color='k', lw=1)
txt_theta_lbl = ax_bloch.text(0,0,0, r"$\theta$", color='k', fontsize=10)
txt_phi_lbl = ax_bloch.text(0,0,0, r"$\phi$", color='k', fontsize=10)

# Info Box
angle_text = ax_bloch.text2D(0.02, 0.98, "", transform=ax_bloch.transAxes, fontsize=10, 
                             va='top', bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray', boxstyle='round'))

ax_bloch.legend(handles=[
    Line2D([0], [0], color='#1f77b4', lw=2, ls='--', label='Idle/Past (Dashed)'),
    Line2D([0], [0], color='#d62728', lw=2, ls='-', label='Active R(n,$\\theta$) (Solid)'),
    Line2D([0], [0], color='purple', linestyle='-.', lw=1, label='Quant. Axis'),
    Line2D([0], [0], color='k', marker='<', markersize=8, linestyle='-', lw=2, label='State')
], loc='lower left', fontsize='x-small')

# ==========================================
# Part 4: Inputs & Compute
# ==========================================
input_left, input_width, input_height, start_y, gap_y = 0.05, 0.12, 0.04, 0.80, 0.06
inputs = {}
input_labels = [
    ('fq2', r'$f_{q2}$ (GHz)'), ('fq3', r'$f_{q3}$ (GHz)'), ('theta_in', r'$\theta_{in}$ (deg)'),
    ('t_q2_res', r'$t_{res}$ (ns)'), ('t_q20', r'$t_{q20}$ (ns)'), ('t_q30', r'$t_{q30}$ (ns) ')
]

def make_update_param_callback(key):
    def update_param(text):
        try:
            val = float(text)
            params[key] = val
        except ValueError: print(f"Invalid input for {key}")
    return update_param

for i, (key, label) in enumerate(input_labels):
    y_pos = start_y - i * gap_y
    ax_box = plt.axes([input_left + 0.06, y_pos, input_width, input_height]) 
    fig.text(input_left + 0.05, y_pos + 0.01, label + ":", ha='right', va='bottom', fontsize=10)
    text_box = TextBox(ax_box, "", initial=str(params[key]))
    text_box.on_submit(make_update_param_callback(key))
    inputs[key] = text_box

y_compute = start_y - len(input_labels) * gap_y - 0.02
ax_compute = plt.axes([input_left + 0.06, y_compute, input_width, 0.05])
btn_compute = Button(ax_compute, "Compute Heatmap", color='lightblue', hovercolor='0.9')
btn_compute.on_clicked(recalculate_heatmap)
fig.text(input_left, start_y + 0.06, "Parameters", fontsize=12, fontweight='bold')

# ==========================================
# Part 5: Sliders & Buttons
# ==========================================
axcolor = 'lightgoldenrodyellow'
x_slider_start, w_slider, w_btn, h_ctrl = 0.30, 0.45, 0.03, 0.03
y_t2, y_t3, y_tm = 0.12, 0.08, 0.04

ax_t2 = plt.axes([x_slider_start, y_t2, w_slider, h_ctrl], facecolor=axcolor)
ax_t3 = plt.axes([x_slider_start, y_t3, w_slider, h_ctrl], facecolor=axcolor)
ax_tm = plt.axes([x_slider_start, y_tm, w_slider, h_ctrl], facecolor=axcolor)

x_btn_l = x_slider_start + w_slider + 0.01
x_btn_r = x_btn_l + w_btn + 0.01
x_val_text = x_btn_r + w_btn + 0.01

ax_b_t2_l, ax_b_t2_r = plt.axes([x_btn_l, y_t2, w_btn, h_ctrl]), plt.axes([x_btn_r, y_t2, w_btn, h_ctrl])
ax_b_t3_l, ax_b_t3_r = plt.axes([x_btn_l, y_t3, w_btn, h_ctrl]), plt.axes([x_btn_r, y_t3, w_btn, h_ctrl])
ax_b_tm_l, ax_b_tm_r = plt.axes([x_btn_l, y_tm, w_btn, h_ctrl]), plt.axes([x_btn_r, y_tm, w_btn, h_ctrl])

s_t2 = Slider(ax_t2, r'$t_2$ (ns) ', 0.0, 50.0, valinit=10.0)
s_t3 = Slider(ax_t3, r'$t_3$ (ns) ', 0.0, 50.0, valinit=10.0)
s_time = Slider(ax_tm, 'Lab Time (ns) ', 0.0, 400.0, valinit=60.0)
for s in [s_t2, s_t3, s_time]: s.valtext.set_visible(False)

txt_t2 = fig.text(x_val_text, y_t2 + h_ctrl/4, "10.00 ns", fontsize=11)
txt_t3 = fig.text(x_val_text, y_t3 + h_ctrl/4, "10.00 ns", fontsize=11)
txt_tm = fig.text(x_val_text, y_tm + h_ctrl/4, "60.00 ns", fontsize=11)

btn_props = dict(color=axcolor, hovercolor='0.975')
btn_t2_l, btn_t2_r = Button(ax_b_t2_l, '<', **btn_props), Button(ax_b_t2_r, '>', **btn_props)
btn_t3_l, btn_t3_r = Button(ax_b_t3_l, '<', **btn_props), Button(ax_b_t3_r, '>', **btn_props)
btn_tm_l, btn_tm_r = Button(ax_b_tm_l, '<', **btn_props), Button(ax_b_tm_r, '>', **btn_props)

def create_stepper(slider, step_val):
    def step_func(event):
        new_val = np.clip(slider.val + step_val, slider.valmin, slider.valmax)
        slider.set_val(new_val)
    return step_func

btn_t2_l.on_clicked(create_stepper(s_t2, -0.1))
btn_t2_r.on_clicked(create_stepper(s_t2, 0.1))
btn_t3_l.on_clicked(create_stepper(s_t3, -0.1))
btn_t3_r.on_clicked(create_stepper(s_t3, 0.1))
btn_tm_l.on_clicked(create_stepper(s_time, -0.5))
btn_tm_r.on_clicked(create_stepper(s_time, 0.5))

# ==========================================
# Part 6: Update Logic
# ==========================================

def update_all(val):
    global state_arrow, q_axis_arrow
    
    t2_val = s_t2.val
    t3_val = s_t3.val
    lab_time = s_time.val
    
    txt_t2.set_text(f"{t2_val:5.2f} ns")
    txt_t3.set_text(f"{t3_val:5.2f} ns")
    txt_tm.set_text(f"{lab_time:5.1f} ns")
    
    crosshair.set_data([t2_val], [t3_val])
    
    # 1. Trajectory with Styles
    segments_data, last_pt = calculate_trajectory_segments(t2_val, t3_val, lab_time, params)
    for line in lines_3d:
        line.set_data([], [])
        line.set_3d_properties([])
    for i, dat in enumerate(segments_data):
        if i < len(lines_3d):
            lines_3d[i].set_data(dat['x'], dat['y'])
            lines_3d[i].set_3d_properties(dat['z'])
            lines_3d[i].set_color(dat['c'])
            lines_3d[i].set_linestyle(dat['ls']) # Apply dashed/solid
            
    # 2. Arrows
    lx, ly, lz = last_pt
    if state_arrow: state_arrow.remove()
    state_arrow = ax_bloch.quiver(0, 0, 0, lx, ly, lz, length=1.0, color='k', pivot='tail', arrow_length_ratio=0.1, linewidth=2.0)
    
    if q_axis_arrow: q_axis_arrow.remove()
    _, _, qax = get_hamiltonians_and_axis(params)
    qv = qax * 1.2
    q_axis_arrow = ax_bloch.quiver(0, 0, 0, qv[0], qv[1], qv[2], color='purple', arrow_length_ratio=0.1, linestyle='-.')
    
    # 3. Angle Arcs & Labels
    phi_arc_dat, theta_arc_dat, l_phi, l_theta = get_angle_plot_data(lx, ly, lz)
    arc_phi.set_data(phi_arc_dat[0], phi_arc_dat[1])
    arc_phi.set_3d_properties(phi_arc_dat[2])
    
    arc_theta.set_data(theta_arc_dat[0], theta_arc_dat[1])
    arc_theta.set_3d_properties(theta_arc_dat[2])
    
    proj_line.set_data([lx, lx, 0], [ly, ly, 0])
    proj_line.set_3d_properties([lz, 0, 0])
    
    txt_phi_lbl.set_position((l_phi[0], l_phi[1]))
    txt_phi_lbl.set_3d_properties(l_phi[2])
    
    txt_theta_lbl.set_position((l_theta[0], l_theta[1]))
    txt_theta_lbl.set_3d_properties(l_theta[2])

    # 4. Info Text
    theta_deg = np.degrees(np.arccos(np.clip(lz, -1, 1)))
    phi_deg = np.degrees(np.arctan2(ly, lx))
    
    U_single = func_U_matrix(t2_val, t3_val, params)
    U_double = U_single @ U_single 
    _, rot_rad_single = func2_SU2decompose(U_single)
    rot_deg_single = np.degrees(rot_rad_single)
    
    psi_final = U_double @ np.array([0, 1])
    fidelity = np.abs(psi_final[0])**2
    infidelity = 1.0 - fidelity
    
    info_str = (f"Current Position:\n"
                f" $\\theta$ = {theta_deg:.1f}°\n"
                f" $\\phi$ = {phi_deg:.1f}°\n"
                f"-----------------\n"
                f"Single $R(\\hat{{n}}, \\theta)$ Angle: {rot_deg_single:.2f}°\n"
                f"Target Single Angle: 90°\n"
                f"-----------------\n"
                f"Seq Fidelity ($P_{{\\uparrow}}$): {fidelity:.5f}\n"
                f"Infidelity: {infidelity:.4e}")
    angle_text.set_text(info_str)

s_t2.on_changed(update_all)
s_t3.on_changed(update_all)
s_time.on_changed(update_all)

print("Initializing...")
recalculate_heatmap()
update_all(0)
plt.show()  