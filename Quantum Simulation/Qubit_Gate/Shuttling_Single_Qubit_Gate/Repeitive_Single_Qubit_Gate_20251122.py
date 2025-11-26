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
    'fq3': 0.06203767,
    't_wait1': 0.0,  
    't_wait2': 0.0,  
    'N_rep': 1      
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
    d_w1 = p['t_wait1']
    d_w2 = p['t_wait2']
    
    U_res = scipy.linalg.expm(-2*np.pi*1j * H_q2 * d_res)
    U_3   = scipy.linalg.expm(-2*np.pi*1j * H_q3 * d_3)
    U_2   = scipy.linalg.expm(-2*np.pi*1j * H_q2 * d_2)
    U_z1  = scipy.linalg.expm(-2*np.pi*1j * H_q2 * d_w1)
    U_z2  = scipy.linalg.expm(-2*np.pi*1j * H_q2 * d_w2)
    
    U_single = U_res @ U_z2 @ U_3 @ U_2 @ U_3 @ U_z1 @ U_res
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
    d_w1  = max(0, p['t_wait1'])
    d_w2  = max(0, p['t_wait2'])
    d_gap = 1.0 / p['fq2'] if p['fq2'] > 1e-9 else 0
    
    single_gate_structure = [
        (d_res, H_q2, 0), (d_w1,  H_q2, 0), (d_3,   H_q3, 1), 
        (d_2,   H_q2, 0), (d_3,   H_q3, 1), (d_w2,  H_q2, 0), (d_res, H_q2, 0)
    ]
    
    sequence = []
    gate_end_indices = set()
    n_reps = int(p.get('N_rep', 1))
    current_seg_idx = 0
    
    for i in range(n_reps):
        for seg in single_gate_structure:
            sequence.append(seg)
            current_seg_idx += 1
        gate_end_indices.add(current_seg_idx - 1)
        if i < n_reps - 1:
            sequence.append((d_gap, H_q2, 2)) 
            current_seg_idx += 1
            
    sequence.append((999.0, H_q2, 0))
    
    psi = np.array([0, 1], dtype=complex)
    plot_segments = []
    waypoints = []
    time_remaining = lab_time_limit
    
    colors = {0: '#1f77b4', 1: '#d62728', 2: '#2ca02c'}
    linestyles = {0: '--', 1: '-', 2: '-.'}
    
    gate_counter = 0
    
    for idx, (dur, H, dot_id) in enumerate(sequence):
        if time_remaining <= 0: break
        
        sim_time = dur if time_remaining >= dur else time_remaining
        is_completed = (time_remaining >= dur)
        
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
        
        ls = linestyles[dot_id]
        if dot_id == 1 and is_completed: ls = '--'
        
        plot_segments.append({
            'x': np.array(xs), 'y': np.array(ys), 'z': np.array(zs),
            'c': colors[dot_id], 'ls': ls
        })
        
        psi = scipy.linalg.expm(-1j * 2 * np.pi * H * sim_time) @ start_psi
        time_remaining -= sim_time
        
        if idx in gate_end_indices and is_completed:
            gate_counter += 1
            waypoints.append({'pos': (xs[-1], ys[-1], zs[-1]), 'n': gate_counter})
            
    last_pt = (plot_segments[-1]['x'][-1], plot_segments[-1]['y'][-1], plot_segments[-1]['z'][-1]) if plot_segments else (0,0,-1)
    return plot_segments, last_pt, waypoints

def get_angle_plot_data(lx, ly, lz):
    r = 0.4
    phi = np.arctan2(ly, lx)
    theta = np.arccos(np.clip(lz, -1, 1))
    
    t_phi = np.linspace(0, phi, 40)
    px, py, pz = r * np.cos(t_phi), r * np.sin(t_phi), np.zeros_like(t_phi)
    
    t_theta = np.linspace(0, theta, 40)
    tx = r * np.sin(t_theta) * np.cos(phi)
    ty = r * np.sin(t_theta) * np.sin(phi)
    tz = r * np.cos(t_theta)
    
    l_phi = (r*1.2 * np.cos(phi/2), r*1.2 * np.sin(phi/2), 0)
    l_theta = (r*1.2 * np.sin(theta/2) * np.cos(phi), r*1.2 * np.sin(theta/2) * np.sin(phi), r*1.2 * np.cos(theta/2))
    return (px, py, pz), (tx, ty, tz), l_phi, l_theta

# ==========================================
# Part 2: Visualization Layout
# ==========================================

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, height_ratios=[0.7, 1.1], width_ratios=[0.25, 1, 1], bottom=0.32, wspace=0.3)

ax_inputs_dummy = fig.add_subplot(gs[:, 0]); ax_inputs_dummy.axis('off')
ax_img   = fig.add_subplot(gs[0, 1]); ax_img.axis('off')
ax_map   = fig.add_subplot(gs[1, 1]) 
ax_bloch = fig.add_subplot(gs[:, 2], projection='3d') 

# --- Image ---
image_filename = "C:\\Users\\Administrator\\Desktop\\Spin_Quantum_Hardware_Simulation\\Quantum Simulation\\Qubit_Gate\\Shuttling_Single_Qubit_Gate\\schematic.png" 
if os.path.exists(image_filename):
    img = mpimg.imread(image_filename)
    ax_img.imshow(img); ax_img.set_title("Pulse Sequence Schematic")
else:
    ax_img.text(0.5, 0.5, "Image not found", ha='center', va='center', color='red')

# --- Heatmap ---
map_res = 60
t_q2s = np.linspace(0, 50, map_res)
t_q3s = np.linspace(0, 50, map_res)
im_obj = None; contour_objs = []
crosshair, = ax_map.plot([], [], 'r+', markersize=20, markeredgewidth=3, zorder=10)

def recalculate_heatmap(event=None):
    global im_obj, contour_objs
    print("Recalculating map...")
    for c in contour_objs:
        try: c.remove()
        except: 
            for coll in c.collections: coll.remove()
    contour_objs = []
    
    Pups = np.zeros((map_res, map_res))
    GatePol = np.zeros((map_res, map_res))
    RotAng = np.zeros((map_res, map_res))
    
    for i in range(map_res):
        for j in range(map_res):
            t2, t3 = t_q2s[i], t_q3s[j]
            U_single = func_U_matrix(t2, t3, params)
            gp, ra = func2_SU2decompose(U_single)
            stateN = U_single @ U_single @ np.array([0,1])
            Pups[j, i] = (np.real(np.abs(stateN[0]))**2 ) * params['A'] + params['B']
            RotAng[j, i] = ra
            GatePol[j, i] = gp

    if im_obj is None:
        im_obj = ax_map.pcolormesh(t_q2s, t_q3s, Pups, shading='auto', cmap='viridis')
        fig.colorbar(im_obj, ax=ax_map).set_label('Spin-Up Prob (2 Gates)')
    else:
        im_obj.set_array(Pups.ravel()); im_obj.set_clim(np.min(Pups), np.max(Pups))
    
    c1 = ax_map.contour(t_q2s, t_q3s, GatePol/DEG, levels=[90], colors='C9', linewidths=1.5)
    c2 = ax_map.contour(t_q2s, t_q3s, RotAng/DEG, levels=[90], colors='C1', linewidths=1.5)
    c3 = ax_map.contour(t_q2s, t_q3s, RotAng/DEG, levels=[270], colors='C0', linewidths=1.5)
    contour_objs.extend([c1, c2, c3])
    ax_map.set_title(f"Map (w1={params['t_wait1']}, w2={params['t_wait2']})")
    update_all(0); fig.canvas.draw_idle()

ax_map.set_xlabel(r'$t_2$ (ns)'); ax_map.set_ylabel(r'$t_3$ (ns)')
ax_map.legend(handles=[
    Line2D([0],[0], color='C9', lw=1.5, label=r'Axis Equator ($90^\circ$)'),
    Line2D([0],[0], color='C1', lw=1.5, label=r'Rot $90^\circ$'),
    Line2D([0],[0], color='C0', lw=1.5, label=r'Rot $270^\circ$')
], loc='upper right', fontsize='x-small', framealpha=0.8)

# --- Bloch Sphere ---
u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:15j]
x_sphere, y_sphere, z_sphere = np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)
ax_bloch.plot_wireframe(x_sphere, y_sphere, z_sphere, color="gray", alpha=0.1, linewidth=0.5)
ax_bloch.plot([-1.2, 1.2], [0, 0], [0, 0], 'r'); ax_bloch.text(1.3, 0, 0, "x", color='r')
ax_bloch.plot([0, 0], [-1.2, 1.2], [0, 0], 'g'); ax_bloch.text(0, 1.3, 0, "y", color='g')
ax_bloch.plot([0, 0], [0, 0], [-1.2, 1.2], 'b'); ax_bloch.text(0, 0, 1.3, "z", color='b')
ax_bloch.set_axis_off(); ax_bloch.set_aspect('equal')

lines_3d = [ax_bloch.plot([], [], [], linewidth=2, alpha=0.9)[0] for _ in range(100)]
waypoint_scatter = ax_bloch.scatter([], [], [], s=50, c='orange', marker='o', depthshade=False)
waypoint_texts = []
state_arrow = None
q_axis_arrow = ax_bloch.quiver(0, 0, 0, 0, 0, 1, color='purple', linestyle='-.') 
proj_line, = ax_bloch.plot([], [], [], 'k:', lw=1, alpha=0.6)
arc_phi, = ax_bloch.plot([], [], [], 'k', lw=1)
arc_theta, = ax_bloch.plot([], [], [], 'k', lw=1)
txt_lbls = [ax_bloch.text(0,0,0,t,color='k',fontsize=10) for t in [r"$\theta$", r"$\phi$"]]
angle_text = ax_bloch.text2D(0.02, 0.98, "", transform=ax_bloch.transAxes, fontsize=9, va='top', 
                             bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray', boxstyle='round'))

ax_bloch.legend(handles=[
    Line2D([0],[0], color='#1f77b4', lw=2, ls='--', label='Idle'),
    Line2D([0],[0], color='#d62728', lw=2, ls='-', label='Pulse'),
    Line2D([0],[0], color='#2ca02c', lw=2, ls='-.', label='Gap'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='orange', markersize=8, label='Gate End'),
    Line2D([0],[0], color='purple', linestyle='-.', lw=1, label='Quant. Axis')
], loc='lower left', fontsize='x-small')

# ==========================================
# Part 4: Inputs
# ==========================================
input_left, input_width, start_y, gap_y = 0.05, 0.12, 0.80, 0.06
inputs = {}
labels = [('fq2', '$f_{q2}$'), ('fq3', '$f_{q3}$'), ('theta_in', r'$\theta_{in}$'),
          ('t_q2_res', '$t_{res}$'), ('t_q20', '$t_{q20}$'), ('t_q30', '$t_{q30}$'),
          ('t_wait1', '$t_{w1}$'), ('t_wait2', '$t_{w2}$')]

def make_cb(key):
    def update(text):
        try:
            val = float(text)
            params[key] = val
            if key == 't_wait1': s_w1.set_val(val)
            elif key == 't_wait2': s_w2.set_val(val)
            update_all(0)
            fig.canvas.draw_idle()
        except ValueError: pass
    return update

for i, (k, l) in enumerate(labels):
    y = start_y - i * gap_y
    fig.text(input_left+0.05, y+0.01, l+":", ha='right', fontsize=9)
    tb = TextBox(plt.axes([input_left+0.06, y, input_width, 0.04]), "", initial=str(params[k]))
    tb.on_submit(make_cb(k)); inputs[k] = tb

btn_comp = Button(plt.axes([input_left+0.06, start_y - len(labels)*gap_y - 0.03, input_width, 0.05]), "Map", color='lightblue')
btn_comp.on_clicked(recalculate_heatmap)

# ==========================================
# Part 5: Controls & Button Fix
# ==========================================
x_s, w_s, h_c = 0.30, 0.45, 0.03
y_b, y_st = 0.04, 0.045
sliders = {}

def create_slider(name, label, vmin, vmax, vinit, y_idx, fmt=None, step=None):
    y = y_b + y_idx * y_st
    ax = plt.axes([x_s, y, w_s, h_c], facecolor='lightgoldenrodyellow')
    s = Slider(ax, label, vmin, vmax, valinit=vinit, valfmt=fmt, valstep=step)
    s.valtext.set_visible(False)
    
    # --- Buttons ---
    ax_l = plt.axes([x_s + w_s + 0.01, y, 0.03, h_c])
    ax_r = plt.axes([x_s + w_s + 0.05, y, 0.03, h_c])
    
    # 创建按钮对象
    bl = Button(ax_l, '<', color='0.85', hovercolor='0.95')
    br = Button(ax_r, '>', color='0.85', hovercolor='0.95')
    
    val_txt = fig.text(x_s + w_s + 0.09, y + h_c/4, "", fontsize=10)
    
    delta = step if step else (vmax-vmin)/100
    
    # 绑定事件
    # 注意：lambda 需要绑定 s=s, delta=delta 避免闭包问题（虽然这里是立即执行的，但好习惯）
    bl.on_clicked(lambda e, s=s, d=delta: s.set_val(s.val - d))
    br.on_clicked(lambda e, s=s, d=delta: s.set_val(s.val + d))
    
    # --- 关键修复 ---
    # 必须保存 bl 和 br 对象的引用，否则它们会被垃圾回收，导致点击无反应
    sliders[name] = {
        's': s, 
        'txt': val_txt, 
        'bl': bl, # 保存左按钮引用
        'br': br  # 保存右按钮引用
    }
    return s

s_t2  = create_slider('t2', '$t_2$', 0, 50, 10.0, 5)
s_t3  = create_slider('t3', '$t_3$', 0, 50, 10.0, 4)
s_w1  = create_slider('w1', '$t_{w1}$', -3, 30, 0.0, 3)
s_w2  = create_slider('w2', '$t_{w2}$', -3, 30, 0.0, 2)
s_rep = create_slider('rep', 'Rep', 1, 20, 1, 1, '%0.0f', 1)
s_tm  = create_slider('tm', 'Time', 0, 400, 60.0, 0, None, 1.0)

# ==========================================
# Part 6: Update
# ==========================================
def update_all(val):
    global state_arrow, q_axis_arrow, waypoint_texts
    params['t_wait1'] = s_w1.val; params['t_wait2'] = s_w2.val; params['N_rep'] = int(s_rep.val)
    
    for k, v in sliders.items(): 
        fmt = "{:.0f}" if k == 'rep' else "{:.2f}ns"
        if k == 'tm': fmt = "{:.1f}ns"
        v['txt'].set_text(fmt.format(v['s'].val))
    
    crosshair.set_data([s_t2.val], [s_t3.val])
    segments, last_pt, waypoints = calculate_trajectory_segments(s_t2.val, s_t3.val, s_tm.val, params)
    
    for line in lines_3d: line.set_data([],[]); line.set_3d_properties([])
    for i, dat in enumerate(segments):
        if i < len(lines_3d):
            lines_3d[i].set_data(dat['x'], dat['y'])
            lines_3d[i].set_3d_properties(dat['z'])
            lines_3d[i].set_color(dat['c'])
            lines_3d[i].set_linestyle(dat['ls'])
            
    for t in waypoint_texts: t.remove()
    waypoint_texts = []
    if waypoints:
        wx, wy, wz = zip(*[w['pos'] for w in waypoints])
        waypoint_scatter._offsets3d = (wx, wy, wz)
        for w in waypoints:
            px, py, pz = w['pos']
            theta_deg = np.degrees(np.arccos(np.clip(pz, -1, 1)))
            phi_deg = np.degrees(np.arctan2(py, px))
            t_obj = ax_bloch.text(px, py, pz+0.15, f"n={w['n']}\n$\\theta$={theta_deg:.0f}°\n$\\phi$={phi_deg:.0f}°", 
                                  fontsize=8, color='darkorange', ha='center', fontweight='bold')
            waypoint_texts.append(t_obj)
    else: waypoint_scatter._offsets3d = ([], [], [])

    lx, ly, lz = last_pt
    if state_arrow: state_arrow.remove()
    state_arrow = ax_bloch.quiver(0,0,0, lx, ly, lz, color='k', pivot='tail', length=1.0)
    
    if q_axis_arrow: q_axis_arrow.remove()
    _, _, qax = get_hamiltonians_and_axis(params)
    qv = qax * 1.2
    q_axis_arrow = ax_bloch.quiver(0, 0, 0, qv[0], qv[1], qv[2], color='purple', arrow_length_ratio=0.1, linestyle='-.')
    
    p_dat, t_dat, l_phi, l_theta = get_angle_plot_data(lx, ly, lz)
    arc_phi.set_data(p_dat[0], p_dat[1]); arc_phi.set_3d_properties(p_dat[2])
    arc_theta.set_data(t_dat[0], t_dat[1]); arc_theta.set_3d_properties(t_dat[2])
    proj_line.set_data([lx, lx, 0], [ly, ly, 0]); proj_line.set_3d_properties([lz, 0, 0])
    txt_lbls[1].set_position(l_phi[:2]); txt_lbls[1].set_3d_properties(l_phi[2])
    txt_lbls[0].set_position(l_theta[:2]); txt_lbls[0].set_3d_properties(l_theta[2])

    U_single = func_U_matrix(s_t2.val, s_t3.val, params)
    U_tot = np.eye(2, dtype=complex)
    for _ in range(int(s_rep.val)): U_tot = U_single @ U_tot
    fid = np.abs((U_tot @ [0,1])[0])**2
    angle_text.set_text(f"Current (t={s_tm.val:.1f}):\n $\\theta$={np.degrees(np.arccos(np.clip(lz,-1,1))):.1f}°\n $\\phi$={np.degrees(np.arctan2(ly,lx)):.1f}°\n Seq Fid ($P_{{\\uparrow}}$): {fid:.4f}")

for s in sliders.values(): s['s'].on_changed(update_all)
recalculate_heatmap()
plt.show()