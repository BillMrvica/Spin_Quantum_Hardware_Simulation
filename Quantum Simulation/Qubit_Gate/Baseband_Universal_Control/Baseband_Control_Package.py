import numpy as np
import scipy.linalg

# --- Constants ---
DEG = np.pi/180
s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

# Basis States
basis_0 = np.array([1, 0], dtype=complex) # Spin Up |0>
basis_1 = np.array([0, 1], dtype=complex) # Spin Down |1>
state_00 = np.kron(basis_0, basis_0)

class ShuttlingEngine:
    def __init__(self, fq2, fq3, theta, t_res, t_q20, t_q30):
        self.update(fq2, fq3, theta, t_res, t_q20, t_q30)

    def update(self, fq2, fq3, theta, t_res, t_q20, t_q30):
        self.fq2 = fq2
        self.fq3 = fq3
        self.theta = theta
        self.t_res = t_res
        self.t_q20 = t_q20
        self.t_q30 = t_q30
        
        # Hamiltonians
        th_rad = theta * DEG
        qx, qy, qz = np.sin(th_rad), 0, np.cos(th_rad)
        self.H_q3 = np.array([[qz, qx], [qx, -qz]]) * fq3 * 0.5
        self.H_q2 = np.array([[1, 0], [0, -1]]) * fq2 * 0.5

    def get_quantization_axis(self):
        """Returns the vector representing the H_q3 axis."""
        th_rad = self.theta * DEG
        return np.array([np.sin(th_rad), 0, np.cos(th_rad)])

    def get_unitary(self, t2, t3, z_deg=0):
        """
        Calculates unitary for ONE full pulse sequence.
        Order: U_res @ U_z @ U_3 @ U_2 @ U_3 @ U_res
        """
        d_res = max(0, self.t_res)
        d_3 = max(0, t3 - self.t_q30)
        d_2 = max(0, t2 - self.t_q20)
        
        U_res = scipy.linalg.expm(-2*np.pi*1j * self.H_q2 * d_res)
        U_3   = scipy.linalg.expm(-2*np.pi*1j * self.H_q3 * d_3)
        U_2   = scipy.linalg.expm(-2*np.pi*1j * self.H_q2 * d_2)
        U_z   = scipy.linalg.expm(-1j * (z_deg*DEG/2) * sz)
        
        # Modified Order: Z is applied after the core sequence (before the last Res if viewing from right)
        # Time: Res_start -> U3 -> U2 -> U3 -> Z -> Res_end
        return U_res @ U_z @ U_3 @ U_2 @ U_3 @ U_res

    def get_duration(self, t2, t3):
        d_res = max(0, self.t_res)
        d_3 = max(0, t3 - self.t_q30)
        d_2 = max(0, t2 - self.t_q20)
        return 2*d_res + 2*d_3 + d_2

    def get_free_evolution_segment(self, start_state, duration, steps=50):
        """Calculates trajectory for free precession (H_q2)."""
        if duration <= 0: return None
        t_steps = np.linspace(0, duration, steps)
        seg_points = []
        psi = start_state.copy()
        
        # Initial point
        rho = np.outer(psi, psi.conj())
        seg_points.append([np.real(np.trace(rho @ sx)), np.real(np.trace(rho @ sy)), np.real(np.trace(rho @ sz))])
        
        for dt in t_steps[1:]:
            U = scipy.linalg.expm(-2*np.pi*1j * self.H_q2 * dt)
            pt = U @ psi
            rho = np.outer(pt, pt.conj())
            seg_points.append([np.real(np.trace(rho @ sx)), np.real(np.trace(rho @ sy)), np.real(np.trace(rho @ sz))])
            
        return {'points': np.array(seg_points), 'type': 'q2', 'duration': duration}

    def get_trajectory_segments(self, t2, t3, z_deg, init_state=None, steps_per_seg=30):
        """
        Returns (segments, final_state).
        Trajectory matches the unitary order: Res -> U3 -> U2 -> U3 -> Z -> Res
        """
        d_res = max(0, self.t_res)
        d_3 = max(0, t3 - self.t_q30)
        d_2 = max(0, t2 - self.t_q20)
        
        # Definition: (Duration, Operator, TypeLabel)
        segments_def = [
            (d_res, self.H_q2, 'q2'),       # 1. Res
            (d_3, self.H_q3, 'q3'),         # 2. U3
            (d_2, self.H_q2, 'q2'),         # 3. U2
            (d_3, self.H_q3, 'q3'),         # 4. U3
            (0, 'Z', 'z'),                  # 5. Z (Instant)
            (d_res, self.H_q2, 'q2')        # 6. Res
        ]
        
        if init_state is None:
            psi = basis_1.copy() 
        else:
            psi = init_state.copy()
            
        segments_out = []
        
        for duration, Op, label in segments_def:
            seg_points = []
            
            # Handle Instantaneous Z jump
            if isinstance(Op, str) and Op == 'Z':
                rho_before = np.outer(psi, psi.conj())
                p_before = [np.real(np.trace(rho_before @ sx)), np.real(np.trace(rho_before @ sy)), np.real(np.trace(rho_before @ sz))]
                
                psi = scipy.linalg.expm(-1j * (z_deg*DEG/2) * sz) @ psi
                
                rho_after = np.outer(psi, psi.conj())
                p_after = [np.real(np.trace(rho_after @ sx)), np.real(np.trace(rho_after @ sy)), np.real(np.trace(rho_after @ sz))]
                
                segments_out.append({'points': np.array([p_before, p_after]), 'type': 'z', 'duration': 0.0})
                continue
                
            # Handle Time Evolution
            n_steps = int(max(2, steps_per_seg * (duration/5.0) + 2))
            t_steps = np.linspace(0, duration, n_steps)
            start_psi = psi.copy()
            
            rho = np.outer(start_psi, start_psi.conj())
            seg_points.append([np.real(np.trace(rho @ sx)), np.real(np.trace(rho @ sy)), np.real(np.trace(rho @ sz))])

            for dt in t_steps[1:]:
                U = scipy.linalg.expm(-2*np.pi*1j * Op * dt)
                pt = U @ start_psi
                rho = np.outer(pt, pt.conj())
                seg_points.append([np.real(np.trace(rho @ sx)), np.real(np.trace(rho @ sy)), np.real(np.trace(rho @ sz))])
            
            psi = scipy.linalg.expm(-2*np.pi*1j * Op * duration) @ start_psi
            
            segments_out.append({'points': np.array(seg_points), 'type': label, 'duration': duration})
            
        return segments_out, psi

    def get_double_pulse_trajectory(self, t2, t3, z_deg, init_state=None):
        """
        Returns trajectory for: Gate 1 -> Wait(tau) -> Gate 2.
        Wait time tau = 1 / fq2 (One spin procession period).
        """
        # 1. First Pulse (Using new return signature)
        segs1, psi_end1 = self.get_trajectory_segments(t2, t3, z_deg, init_state)
        
        # 2. Free Evolution (tau)
        if self.fq2 > 1e-9:
            tau = 1.0 / self.fq2
        else:
            tau = 0.0
            
        seg_wait = self.get_free_evolution_segment(psi_end1, tau, steps=40)
        
        # Update state after wait
        if seg_wait:
            U_wait = scipy.linalg.expm(-2*np.pi*1j * self.H_q2 * tau)
            psi_end_wait = U_wait @ psi_end1
        else:
            psi_end_wait = psi_end1
            
        # 3. Second Pulse
        segs2, psi_final = self.get_trajectory_segments(t2, t3, z_deg, init_state=psi_end_wait)
        
        # Combine lists
        all_segs = segs1 + ([seg_wait] if seg_wait else []) + segs2
        return all_segs

# Helper functions
def partial_trace(psi_4x1):
    rho = np.outer(psi_4x1, psi_4x1.conj())
    ra = np.zeros((2,2), dtype=complex)
    ra[0,0] = rho[0,0]+rho[1,1]; ra[0,1] = rho[0,2]+rho[1,3]
    ra[1,0] = rho[2,0]+rho[3,1]; ra[1,1] = rho[2,2]+rho[3,3]
    rb = np.zeros((2,2), dtype=complex)
    rb[0,0] = rho[0,0]+rho[2,2]; rb[0,1] = rho[0,1]+rho[2,3]
    rb[1,0] = rho[1,0]+rho[3,2]; rb[1,1] = rho[1,1]+rho[3,3]
    def vec(r): return [np.real(np.trace(r@sx)), np.real(np.trace(r@sy)), np.real(np.trace(r@sz))]
    return vec(ra), vec(rb)

def get_fidelity(U_actual):
    # Target: X90 (from |0> to |-i>)
    psi_target = np.array([1/np.sqrt(2), -1j/np.sqrt(2)]) 
    psi_actual = U_actual @ basis_0
    return np.abs(np.vdot(psi_target, psi_actual))**2