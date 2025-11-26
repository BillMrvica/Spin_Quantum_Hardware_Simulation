import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
import Baseband_Control_Package as phys 

# ==========================================
# Helper: SU(2) Decomposition
# ==========================================
def decompose_unitary(U):
    tr = np.trace(U)
    cos_half_th = np.real(tr) / 2.0
    cos_half_th = np.clip(cos_half_th, -1.0, 1.0)
    rot_ang = 2 * np.arccos(cos_half_th)
    sin_half_th = np.sin(rot_ang / 2.0)
    
    if abs(sin_half_th) < 1e-6: return 0.0, 0.0
        
    sz_mat = np.array([[1, 0], [0, -1]], dtype=complex)
    nz = np.imag(np.trace(U @ sz_mat)) / (-2 * sin_half_th)
    axis_polar = np.arccos(np.clip(nz, -1.0, 1.0))
    
    return np.degrees(rot_ang), np.degrees(axis_polar)

# ==========================================
# Shared Data Container
# ==========================================
class SystemState:
    def __init__(self):
        # Physics Parameters
        self.pA = {'fq2': 0.0709, 'fq3': 0.062, 'th': 44.7, 'res': 1.16, 'tq20': -2.288, 'tq30': -1.542, 'A': 0.8417, 'B': 0.0899}
        self.pB = {'fq2': 0.0715, 'fq3': 0.033, 'th': 41.5, 'res': 1.16, 'tq20': -2.288, 'tq30': -1.542, 'A': 0.8417, 'B': 0.0899}
        self.tune_t2_A = 10.0; self.tune_t3_A = 10.0; self.tune_z_A  = 0.0
        self.tune_t2_B = 10.0; self.tune_t3_B = 10.0; self.tune_z_B  = 0.0
        self.cz_duration = 5.0
        self.circuit = [] 
        self.engA = phys.ShuttlingEngine(*self.unpack(self.pA))
        self.engB = phys.ShuttlingEngine(*self.unpack(self.pB))

    def unpack(self, p):
        return p['fq2'], p['fq3'], p['th'], p['res'], p['tq20'], p['tq30']

    def update_engines(self):
        self.engA.update(*self.unpack(self.pA))
        self.engB.update(*self.unpack(self.pB))

state = SystemState()

# ==========================================
# GUI Manager
# ==========================================
class QuantumStudio:
    def __init__(self):
        self.fig = plt.figure(figsize=(18, 10))
        self.fig.canvas.manager.set_window_title('Quantum Spin Studio')
        self.current_tab = 1
        
        # Objects storage to prevent GC
        self.cross_A = None; self.cross_B = None
        self.lines_A = []; self.lines_B = []
        self.cbA = None; self.cbB = None
        self.sliders_cache = []; self.buttons_cache = []
        
        self.show_tab_1() 

    def clear_ui(self):
        self.fig.clf()
        self.cbA = None; self.cbB = None
        self.sliders_cache = []; self.buttons_cache = []
        
        ax_nav = self.fig.add_axes([0.0, 0.95, 1.0, 0.05], facecolor='#e0e0e0')
        ax_nav.axis('off')
        
        active = '#87cefa'
        c = [active if self.current_tab == i else 'white' for i in range(1,4)]

        self.btn_t1 = Button(plt.axes([0.05, 0.955, 0.25, 0.04]), '1. Single Qubit Tuning', color=c[0])
        self.btn_t2 = Button(plt.axes([0.35, 0.955, 0.25, 0.04]), '2. Two-Qubit (CZ) Config', color=c[1])
        self.btn_t3 = Button(plt.axes([0.65, 0.955, 0.25, 0.04]), '3. Circuit Composer', color=c[2])
        
        self.btn_t1.on_clicked(lambda x: self.switch_tab(1))
        self.btn_t2.on_clicked(lambda x: self.switch_tab(2))
        self.btn_t3.on_clicked(lambda x: self.switch_tab(3))

    def switch_tab(self, tab_id):
        self.current_tab = tab_id
        if tab_id == 1: self.show_tab_1()
        elif tab_id == 2: self.show_tab_2()
        elif tab_id == 3: self.show_tab_3()

    # ==========================================
    # Tab 1: Single Qubit Tuning
    # ==========================================
    def show_tab_1(self):
        self.clear_ui()
        
        gs = self.fig.add_gridspec(2, 3, 
                                   width_ratios=[0.15, 0.35, 0.5], 
                                   height_ratios=[1, 1], 
                                   top=0.92, bottom=0.08, 
                                   wspace=0.25, hspace=0.3)

        # --- Column 1: Parameters ---
        ax_col1_top = self.fig.add_subplot(gs[0, 0]); ax_col1_top.axis('off')
        self.tbs = [] 

        def add_box(fig_obj, label, key, p_dict, x, y):
            fig_obj.text(x, y, label, fontsize=9)
            ax = plt.axes([x+0.05, y-0.01, 0.07, 0.03])
            tb = TextBox(ax, "", initial=str(p_dict[key]))
            def submit(val): 
                try: p_dict[key] = float(val); state.update_engines()
                except: pass
            tb.on_submit(submit)
            self.tbs.append(tb)

        self.fig.text(0.02, 0.88, "Qubit A (Blue)", fontweight='bold', color='#1f77b4')
        add_box(self.fig, "fq2", 'fq2', state.pA, 0.02, 0.84)
        add_box(self.fig, "fq3", 'fq3', state.pA, 0.02, 0.80)
        add_box(self.fig, "theta", 'th', state.pA, 0.02, 0.76)

        self.fig.text(0.02, 0.45, "Qubit B (Orange)", fontweight='bold', color='#ff7f0e')
        add_box(self.fig, "fq2", 'fq2', state.pB, 0.02, 0.41)
        add_box(self.fig, "fq3", 'fq3', state.pB, 0.02, 0.37)
        add_box(self.fig, "theta", 'th', state.pB, 0.02, 0.33)

        ax_calc = plt.axes([0.02, 0.15, 0.12, 0.05])
        self.btn_calc = Button(ax_calc, "Compute Heatmaps", color='lightgreen')

        # --- Column 2: Heatmaps ---
        ax_hm_A = self.fig.add_subplot(gs[0, 1])
        ax_hm_B = self.fig.add_subplot(gs[1, 1])
        
        def calc_prob_grid(engine, p_dict, N, t_max):
            t2s = np.linspace(0, t_max, N); t3s = np.linspace(0, t_max, N)
            grid_prob = np.zeros((N, N))
            grid_rot = np.zeros((N, N))
            grid_ax = np.zeros((N, N))
            
            init_state = phys.basis_1 # |1> Spin Down
            
            for i, t2 in enumerate(t2s):
                for j, t3 in enumerate(t3s):
                    U = engine.get_unitary(t2, t3)
                    # Double Pulse Sequence
                    psi_final = U @ U @ init_state
                    p_up_raw = np.abs(psi_final[0])**2
                    p_final = p_up_raw * p_dict['A'] + p_dict['B']
                    
                    grid_prob[j, i] = p_final
                    
                    r, a = decompose_unitary(U)
                    grid_rot[j, i] = r
                    grid_ax[j, i] = a
            return t2s, t3s, grid_prob, grid_rot, grid_ax

        def calc_maps(event):
            print("Computing Heatmaps...")
            N = 50; t_max = 50
            
            t2s, t3s, pA, rA, axA = calc_prob_grid(state.engA, state.pA, N, t_max)
            t2s, t3s, pB, rB, axB = calc_prob_grid(state.engB, state.pB, N, t_max)
            
            def plot_hm(ax, Z, R, Ax, title, is_A):
                ax.clear()
                im = ax.pcolormesh(t2s, t3s, Z, shading='auto', cmap='viridis', vmin=0, vmax=1)
                
                ax.contour(t2s, t3s, R, levels=[90], colors='C1', linewidths=1.5)
                ax.contour(t2s, t3s, R, levels=[270], colors='C0', linewidths=1.5)
                ax.contour(t2s, t3s, Ax, levels=[90], colors='cyan', linewidths=1.5)
                
                ax.set_title(f"{title} Spin-Up Prob (Double Pulse)")
                ax.set_ylabel("t3 (ns)")
                ax.set_aspect('equal')
                ax.set_xlim(0, t_max); ax.set_ylim(0, t_max)
                
                if is_A:
                    if self.cbA: self.cbA.remove()
                    self.cbA = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="P_up")
                else:
                    if self.cbB: self.cbB.remove()
                    self.cbB = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="P_up")

            plot_hm(ax_hm_A, pA, rA, axA, "Qubit A", True)
            plot_hm(ax_hm_B, pB, rB, axB, "Qubit B", False)
            ax_hm_B.set_xlabel("t2 (ns)")

            legend_elements = [
                Line2D([0], [0], color='C1', lw=1.5, label='Rot=90°'),
                Line2D([0], [0], color='C0', lw=1.5, label='Rot=270°'),
                Line2D([0], [0], color='cyan', lw=1.5, label='Axis=Equator'),
                Line2D([0], [0], marker='+', color='r', lw=0, markersize=10, label='Point')
            ]
            ax_hm_A.legend(handles=legend_elements, loc='upper right', fontsize=7, framealpha=0.8)

            self.cross_A, = ax_hm_A.plot(state.tune_t2_A, state.tune_t3_A, 'r+', ms=15, mew=2)
            self.cross_B, = ax_hm_B.plot(state.tune_t2_B, state.tune_t3_B, 'b+', ms=15, mew=2)
            self.fig.canvas.draw_idle()

        self.btn_calc.on_clicked(calc_maps)

        # --- Column 3: Bloch & Sliders ---
        ax_bloch_A = self.fig.add_subplot(gs[0, 2], projection='3d')
        pos_A = ax_bloch_A.get_position()
        ax_bloch_A.set_position([pos_A.x0, pos_A.y0 + 0.05, pos_A.width, pos_A.height * 0.85])

        ax_bloch_B = self.fig.add_subplot(gs[1, 2], projection='3d')
        pos_B = ax_bloch_B.get_position()
        ax_bloch_B.set_position([pos_B.x0, pos_B.y0 + 0.05, pos_B.width, pos_B.height * 0.85])
        
        def init_bloch(ax, title):
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x, y, z = np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)
            ax.plot_wireframe(x, y, z, color="gray", alpha=0.1)
            ax.plot([-1.2,1.2],[0,0],[0,0],'r-',lw=1); ax.text(1.3,0,0,'x',color='r')
            ax.plot([0,0],[-1.2,1.2],[0,0],'g-',lw=1); ax.text(0,1.3,0,'y',color='g')
            ax.plot([0,0],[0,0],[-1.2,1.2],'b-',lw=1); 
            ax.text(0,0,1.4,'z |0>',color='b')
            ax.text(0,0,-1.4,'|1> (Start)',color='b')
            
            ax.set_axis_off(); ax.set_aspect('equal'); ax.set_title(title, fontweight='bold')
            
            proj, = ax.plot([],[],[], 'k:', lw=1, alpha=0.6)
            arc_th, = ax.plot([],[],[], 'k', lw=1)
            arc_ph, = ax.plot([],[],[], 'k', lw=1)
            txt_th = ax.text(0,0,0, r"$\theta$", fontsize=9)
            txt_ph = ax.text(0,0,0, r"$\phi$", fontsize=9)
            
            return ax.quiver(0,0,0,0,0,1, color='k', lw=2.5), proj, arc_th, arc_ph, txt_th, txt_ph

        self.vecA, self.projA, self.arcA_th, self.arcA_ph, self.txtA_th, self.txtA_ph = init_bloch(ax_bloch_A, "Qubit A")
        self.vecB, self.projB, self.arcB_th, self.arcB_ph, self.txtB_th, self.txtB_ph = init_bloch(ax_bloch_B, "Qubit B")
        
        qA = state.engA.get_quantization_axis() * 1.2
        self.q_axis_A = ax_bloch_A.plot([0, qA[0]], [0, qA[1]], [0, qA[2]], color='purple', ls='--', lw=1.5)[0]
        qB = state.engB.get_quantization_axis() * 1.2
        self.q_axis_B = ax_bloch_B.plot([0, qB[0]], [0, qB[1]], [0, qB[2]], color='purple', ls='--', lw=1.5)[0]

        self.info_txt_A = ax_bloch_A.text2D(0.0, 0.9, "", transform=ax_bloch_A.transAxes, fontsize=9)
        self.info_txt_B = ax_bloch_B.text2D(0.0, 0.9, "", transform=ax_bloch_B.transAxes, fontsize=9)

        # --- Sliders with Buttons ---
        sx_base = 0.58; sw = 0.22; sh = 0.015; btn_w = 0.025; gap = 0.005
        y_A_start = 0.65; y_B_start = 0.20; dy = 0.035 

        def add_control(y_pos, label, v_init, v_max, step=0.1):
            ax_s = plt.axes([sx_base, y_pos, sw, sh])
            sl = Slider(ax_s, label, 0, v_max, valinit=v_init)
            sl.valtext.set_visible(False) 
            ax_l = plt.axes([sx_base + sw + gap, y_pos, btn_w, sh])
            ax_r = plt.axes([sx_base + sw + gap + btn_w + gap, y_pos, btn_w, sh])
            bl = Button(ax_l, '<', color='0.9', hovercolor='0.8')
            br = Button(ax_r, '>', color='0.9', hovercolor='0.8')
            def dec(e): sl.set_val(sl.val - step)
            def inc(e): sl.set_val(sl.val + step)
            bl.on_clicked(dec); br.on_clicked(inc)
            self.sliders_cache.append(sl); self.buttons_cache.extend([bl, br])
            return sl

        s_t2a = add_control(y_A_start, 'A t2', state.tune_t2_A, 50)
        s_t3a = add_control(y_A_start - dy, 'A t3', state.tune_t3_A, 50)
        s_za  = add_control(y_A_start - 2*dy, 'A Z', state.tune_z_A, 180, step=1.0); s_za.valmin = -180

        s_t2b = add_control(y_B_start, 'B t2', state.tune_t2_B, 50)
        s_t3b = add_control(y_B_start - dy, 'B t3', state.tune_t3_B, 50)
        s_zb  = add_control(y_B_start - 2*dy, 'B Z', state.tune_z_B, 180, step=1.0); s_zb.valmin = -180

        s_lab = add_control(0.03, 'Lab Time', 0, 300, step=1.0) 
        s_lab.ax.set_position([0.60, 0.03, sw, 0.02]) 

        txt_x = sx_base + sw + gap + btn_w + gap + btn_w + 0.01
        self.txt_t2a = self.fig.text(txt_x, y_A_start, "")
        self.txt_t3a = self.fig.text(txt_x, y_A_start - dy, "")
        self.txt_za  = self.fig.text(txt_x, y_A_start - 2*dy, "")
        self.txt_t2b = self.fig.text(txt_x, y_B_start, "")
        self.txt_t3b = self.fig.text(txt_x, y_B_start - dy, "")
        self.txt_zb  = self.fig.text(txt_x, y_B_start - 2*dy, "")
        self.txt_lab = self.fig.text(txt_x, 0.03, "")

        # --- Helper: Update Arcs ---
        def update_arcs(lx, ly, lz, proj, arc_th, arc_ph, txt_th, txt_ph):
            # Projection
            proj.set_data([lx, lx, 0], [ly, ly, 0])
            proj.set_3d_properties([lz, 0, 0])
            
            r = 0.4
            # Phi arc
            phi = np.arctan2(ly, lx)
            t_phi = np.linspace(0, phi, 20)
            arc_ph.set_data(r*np.cos(t_phi), r*np.sin(t_phi))
            arc_ph.set_3d_properties(np.zeros_like(t_phi))
            txt_ph.set_position((r*1.2*np.cos(phi/2), r*1.2*np.sin(phi/2)))
            txt_ph.set_3d_properties(0)
            
            # Theta arc
            theta = np.arccos(np.clip(lz, -1, 1))
            t_th = np.linspace(0, theta, 20)
            tx = r*np.sin(t_th)*np.cos(phi)
            ty = r*np.sin(t_th)*np.sin(phi)
            tz = r*np.cos(t_th)
            arc_th.set_data(tx, ty)
            arc_th.set_3d_properties(tz)
            txt_th.set_position((r*1.2*np.sin(theta/2)*np.cos(phi), r*1.2*np.sin(theta/2)*np.sin(phi)))
            txt_th.set_3d_properties(r*1.2*np.cos(theta/2))

        def plot_segmented_trajectory(ax, segments, line_storage, current_time):
            if isinstance(line_storage, list):
                for l in line_storage: l.remove()
                line_storage.clear()
            
            t_acc = 0
            for seg in segments:
                if t_acc > current_time: break
                pts, typ, dur = seg['points'], seg['type'], seg['duration']
                
                if typ=='q2': c,ls,lw='#1f77b4','--',1.5
                elif typ=='q3': c,ls,lw='#d62728','-',2.0
                elif typ=='z': c,ls,lw='#FF00FF',':',2.5 
                
                plot_pts = pts
                if dur > 0 and current_time < t_acc + dur:
                    idx = int(((current_time-t_acc)/dur) * len(pts))
                    plot_pts = pts[:idx+1]
                
                l, = ax.plot(plot_pts[:,0], plot_pts[:,1], plot_pts[:,2], color=c, ls=ls, lw=lw)
                line_storage.append(l)
                t_acc += dur

        def get_current_vector_and_info(engine, params_d, t2, t3, z, lab_t):
            init_vec = phys.basis_1 
            segs = engine.get_double_pulse_trajectory(t2, t3, z, init_state=init_vec)
            
            found_pt = [0,0,-1]
            t_acc = 0
            for seg in segs:
                d = seg['duration']; pts = seg['points']
                if lab_t < t_acc + d:
                    if d > 0:
                        idx = int(((lab_t - t_acc)/d) * len(pts))
                        found_pt = pts[np.clip(idx, 0, len(pts)-1)]
                    else: found_pt = pts[-1]
                    break
                t_acc += d
                found_pt = pts[-1]

            U = engine.get_unitary(t2, t3, z)
            final = U @ U @ init_vec
            pup = np.abs(final[0])**2 * params_d['A'] + params_d['B']
            fid = phys.get_fidelity(U)
            
            theta = np.degrees(np.arccos(np.clip(found_pt[2], -1, 1)))
            phi = np.degrees(np.arctan2(found_pt[1], found_pt[0]))
            return segs, found_pt, theta, phi, pup, fid

        self.lines_A = []; self.lines_B = []

        def update_bloch(val):
            state.tune_t2_A, state.tune_t3_A, state.tune_z_A = s_t2a.val, s_t3a.val, s_za.val
            state.tune_t2_B, state.tune_t3_B, state.tune_z_B = s_t2b.val, s_t3b.val, s_zb.val
            
            # Update Text
            self.txt_t2a.set_text(f"{state.tune_t2_A:.2f}")
            self.txt_t3a.set_text(f"{state.tune_t3_A:.2f}")
            self.txt_za.set_text(f"{state.tune_z_A:.1f}")
            self.txt_t2b.set_text(f"{state.tune_t2_B:.2f}")
            self.txt_t3b.set_text(f"{state.tune_t3_B:.2f}")
            self.txt_zb.set_text(f"{state.tune_z_B:.1f}")
            self.txt_lab.set_text(f"{s_lab.val:.1f}")

            if self.cross_A: self.cross_A.set_data([s_t2a.val], [s_t3a.val])
            if self.cross_B: self.cross_B.set_data([s_t2b.val], [s_t3b.val])

            # Update A
            segsA, vecA_pos, thA, phA, pA, fidA = get_current_vector_and_info(state.engA, state.pA, s_t2a.val, s_t3a.val, s_za.val, s_lab.val)
            plot_segmented_trajectory(ax_bloch_A, segsA, self.lines_A, s_lab.val)
            self.vecA.remove()
            self.vecA = ax_bloch_A.quiver(0,0,0, vecA_pos[0], vecA_pos[1], vecA_pos[2], color='k', lw=3)
            update_arcs(vecA_pos[0], vecA_pos[1], vecA_pos[2], self.projA, self.arcA_th, self.arcA_ph, self.txtA_th, self.txtA_ph)
            self.info_txt_A.set_text(f"$\\theta={thA:.1f}^\\circ, \\phi={phA:.1f}^\\circ$\nSeq Prob($\\uparrow$): {pA:.4f}\nFidelity(X90): {fidA:.4f}")

            # Update B
            segsB, vecB_pos, thB, phB, pB, fidB = get_current_vector_and_info(state.engB, state.pB, s_t2b.val, s_t3b.val, s_zb.val, s_lab.val)
            plot_segmented_trajectory(ax_bloch_B, segsB, self.lines_B, s_lab.val)
            self.vecB.remove()
            self.vecB = ax_bloch_B.quiver(0,0,0, vecB_pos[0], vecB_pos[1], vecB_pos[2], color='k', lw=3)
            update_arcs(vecB_pos[0], vecB_pos[1], vecB_pos[2], self.projB, self.arcB_th, self.arcB_ph, self.txtB_th, self.txtB_ph)
            self.info_txt_B.set_text(f"$\\theta={thB:.1f}^\\circ, \\phi={phB:.1f}^\\circ$\nSeq Prob($\\uparrow$): {pB:.4f}\nFidelity(X90): {fidB:.4f}")

            self.fig.canvas.draw_idle()

        for s in self.sliders_cache: s.on_changed(update_bloch)
        update_bloch(0)

    # ==========================================
    # Tab 2: CZ Config
    # ==========================================
    def show_tab_2(self):
        self.clear_ui()
        self.fig.text(0.5, 0.8, "Two-Qubit Gate Configuration", ha='center', fontsize=16)
        ax_box = plt.axes([0.45, 0.6, 0.2, 0.05])
        self.fig.text(0.3, 0.61, "CZ Duration (ns):", fontsize=12)
        tb = TextBox(ax_box, "", initial=str(state.cz_duration))
        def submit(val):
            try: state.cz_duration = float(val)
            except: pass
        tb.on_submit(submit)
        self.tb_cz = tb 
        self.fig.canvas.draw_idle()

    # ==========================================
    # Tab 3: Circuit Composer
    # ==========================================
    def show_tab_3(self):
        self.clear_ui()
        gs = self.fig.add_gridspec(3, 2, height_ratios=[0.4, 0.4, 0.2], hspace=0.4)
        ax_bloch_A = self.fig.add_subplot(gs[0, 0], projection='3d')
        ax_bloch_B = self.fig.add_subplot(gs[0, 1], projection='3d')
        ax_timeline = self.fig.add_subplot(gs[1, :])
        
        def init_bloch(ax, title):
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x, y, z = np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)
            ax.plot_wireframe(x, y, z, color="gray", alpha=0.1)
            ax.plot([-1,1],[0,0],[0,0],'r',alpha=0.3); ax.plot([0,0],[-1,1],[0,0],'g',alpha=0.3); ax.plot([0,0],[0,0],[-1,1],'b',alpha=0.3)
            ax.set_axis_off(); ax.set_aspect('equal'); ax.set_title(title, fontweight='bold')
            return ax.quiver(0,0,0,0,0,1, color='k', lw=3)

        self.qvA = init_bloch(ax_bloch_A, "Qubit A")
        self.qvB = init_bloch(ax_bloch_B, "Qubit B")
        
        curr_t = 0
        for gate in state.circuit:
            if gate['type'] == 'A':
                dur = state.engA.get_duration(gate['t2'], gate['t3'])
                ax_timeline.add_patch(plt.Rectangle((curr_t, 0.6), dur, 0.8, color='#1f77b4', alpha=0.7))
                curr_t += dur
            elif gate['type'] == 'B':
                dur = state.engB.get_duration(gate['t2'], gate['t3'])
                ax_timeline.add_patch(plt.Rectangle((curr_t, -0.4), dur, 0.8, color='#ff7f0e', alpha=0.7))
                curr_t += dur
            elif gate['type'] == 'CZ':
                dur = state.cz_duration
                ax_timeline.add_patch(plt.Rectangle((curr_t, 0.6), dur, 0.8, color='green', alpha=0.5))
                ax_timeline.add_patch(plt.Rectangle((curr_t, -0.4), dur, 0.8, color='green', alpha=0.5))
                ax_timeline.plot([curr_t+dur/2]*2, [0,1], 'k-', lw=2)
                curr_t += dur
        
        total_time = max(50, curr_t + 10)
        ax_timeline.set_xlim(0, total_time); ax_timeline.set_ylim(-1, 2)
        ax_timeline.set_yticks([0, 1]); ax_timeline.set_yticklabels(['Qubit B', 'Qubit A'])
        ax_timeline.set_xlabel("Time (ns)"); ax_timeline.grid(True, axis='x', alpha=0.3)
        line_cursor = ax_timeline.axvline(0, color='red', lw=2)

        ax_btn_A = plt.axes([0.1, 0.12, 0.15, 0.05]); ax_btn_B = plt.axes([0.26, 0.12, 0.15, 0.05])
        ax_btn_CZ = plt.axes([0.42, 0.12, 0.15, 0.05]); ax_btn_Clr = plt.axes([0.8, 0.12, 0.1, 0.05])
        btn_A = Button(ax_btn_A, f"Add A", color='#ADD8E6')
        btn_B = Button(ax_btn_B, f"Add B", color='#FFDAB9')
        btn_CZ = Button(ax_btn_CZ, "Add CZ", color='#90EE90')
        btn_Clr = Button(ax_btn_Clr, "Clear", color='#ffcccb')
        
        ax_scroll = plt.axes([0.1, 0.05, 0.8, 0.03])
        s_scroll = Slider(ax_scroll, "Lab Time", 0, total_time, valinit=0)

        def redraw_circuit(): self.show_tab_3()
        def add_A(e): state.circuit.append({'type': 'A', 't2': state.tune_t2_A, 't3': state.tune_t3_A, 'z': state.tune_z_A}); redraw_circuit()
        def add_B(e): state.circuit.append({'type': 'B', 't2': state.tune_t2_B, 't3': state.tune_t3_B, 'z': state.tune_z_B}); redraw_circuit()
        def add_CZ(e): state.circuit.append({'type': 'CZ'}); redraw_circuit()
        def clear(e): state.circuit = []; redraw_circuit()

        btn_A.on_clicked(add_A); btn_B.on_clicked(add_B)
        btn_CZ.on_clicked(add_CZ); btn_Clr.on_clicked(clear)

        def update_evolution(val):
            curr_lab_t = s_scroll.val
            line_cursor.set_xdata([curr_lab_t])
            psi = np.kron(phys.basis_0, phys.basis_0) 
            elapsed = 0
            
            for gate in state.circuit:
                if gate['type'] == 'A': dur = state.engA.get_duration(gate['t2'], gate['t3'])
                elif gate['type'] == 'B': dur = state.engB.get_duration(gate['t2'], gate['t3'])
                else: dur = state.cz_duration
                
                if curr_lab_t >= elapsed + dur:
                    if gate['type'] == 'A': psi = np.kron(state.engA.get_unitary(gate['t2'], gate['t3'], gate['z']), np.eye(2)) @ psi
                    elif gate['type'] == 'B': psi = np.kron(np.eye(2), state.engB.get_unitary(gate['t2'], gate['t3'], gate['z'])) @ psi
                    elif gate['type'] == 'CZ': psi = np.diag([1,1,1,-1]) @ psi
                    elapsed += dur
                else: break
            
            va, vb = phys.partial_trace(psi)
            self.qvA.remove(); self.qvB.remove()
            self.qvA = ax_bloch_A.quiver(0,0,0, va[0], va[1], va[2], color='#1f77b4', lw=3)
            self.qvB = ax_bloch_B.quiver(0,0,0, vb[0], vb[1], vb[2], color='#ff7f0e', lw=3)
            self.fig.canvas.draw_idle()

        s_scroll.on_changed(update_evolution)
        self.t3_widgets = [btn_A, btn_B, btn_CZ, btn_Clr, s_scroll]
        self.fig.canvas.draw_idle()

if __name__ == "__main__":
    app = QuantumStudio()
    plt.show()