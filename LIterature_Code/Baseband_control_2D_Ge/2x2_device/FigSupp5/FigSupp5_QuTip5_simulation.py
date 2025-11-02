# this file is updated for QuTiP 5
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
# The qutip.interpolate module has been removed in QuTiP 5.
# We no longer need Cubic_Spline.
import time
from scipy.optimize import fsolve
import scipy
import pickle



def pulse3(x, xe0, xe1, xe2, xe3):
    x0, e0 = xe0
    x1, e1 = xe1
    x2, e2 = xe2
    x3, e3 = xe3
    if np.isscalar(x):
        if x<x0:
            return e0
        elif x<x1:
            return (-(x-x0)*e1 + (x-x1)*e0)/(x0-x1)
        elif x<x2:
            return (-(x-x1)*e2 + (x-x2)*e1)/(x1-x2)
        elif x<x3:
            return (-(x-x2)*e3 + (x-x3)*e2)/(x2-x3)
        else:
            return e3
    else:   
        y = np.zeros(x.shape)
        
        mask0 = x<x0
        mask01 = np.logical_and(x>=x0,x<x1)
        mask12 = np.logical_and(x>=x1,x<x2)
        mask23 = np.logical_and(x>=x2,x<x3)
        mask3 = x>=x3
        y[mask0] = e0
        y[mask01] = (-(x[mask01]-x0)*e1 + (x[mask01]-x1)*e0)/(x0-x1)
        y[mask12] = (-(x[mask12]-x1)*e2 + (x[mask12]-x2)*e1)/(x1-x2)
        y[mask23] = (-(x[mask23]-x2)*e3 + (x[mask23]-x3)*e2)/(x2-x3)
        y[mask3] = e3

        return y
    

sq2 = np.sqrt(2)
sq3 = np.sqrt(3)

ueV2GHz = (1e-6*1.6e-19/6.62e-34)/1e9
V0 = Qobj(np.diag([1,0]))
s_z = np.array([[1,0],[0,-1]])
s_x = np.array([[0,1],[1, 0]])
s_y = np.array([[0,-1j],[1j, 0]])
s_0 = np.eye(2)
s0, sx, sy, sz = s_0, s_x, s_y, s_z


def SU2decompose(mat):
    u0 = np.real( np.trace( mat @ s0 )/2  )
    ux = np.imag( np.trace( mat @ sx )/2  )
    uy = np.imag( np.trace( mat @ sy )/2  )
    uz = np.imag( np.trace( mat @ sz )/2  )
    
    U = s0* u0 +  (ux*sx + uy*sy + uz*sz) 
    uvec = np.array([ux, uy, uz])
    # if  i==0:
    #     print(  np.sqrt(np.sum(uvec**2))  )    
    if np.sqrt(np.sum(uvec**2)) <= 1e-14:
        uvec = np.array([0,0,1])
    else:
        uvec = uvec/np.sqrt(np.sum(uvec**2))
    direction = uvec
    rot_angle = 2*np.arccos( u0)
    return direction, rot_angle

# pair = '14'
pair = '23'

results = {}
for print_residual in [True, False]:

    
    
    if pair  not in ['14', '23']:
        raise ValueError
    
    if pair == '14':
        t_in = 27 
        phi_in = 65
    else:
        t_in = 20 
        phi_in = 51
    # qubit 1 LZ transition
    # np.exp(-2*np.pi*np.pi*26.9*26.9/((563.87/2)))
    # >> 9.945397242660604e-23
    # qubit 2 LZ transition
    # np.exp(-2*np.pi*np.pi*20*20/((377.026586102719/2)))
    # >> 6.456937929755734e-19
    
    
    b_0, b_1 = 0, 1
    M_pts = 5001   # accuracy for charge degree of freedom
    
    ramptimes = np.array([2])
    if print_residual:
        if pair == '14':
            waittimes = np.array([8.8]) #
        else:
            waittimes = np.array([6.0]) #
    else:
        waittimes = np.linspace(0,30,61)
    
    
    
    finalstate_q2_f = np.zeros((len(ramptimes), len(waittimes),4))
    finalstate_q2_reduction_f = np.zeros((len(ramptimes), len(waittimes),2))
    finalstate_q2_simple_f = np.zeros((len(ramptimes), len(waittimes),2))
    
    finalstate_q2_xyz = np.zeros((len(ramptimes), len(waittimes),3))
    finalstate_q2_reduction_xyz = np.zeros((len(ramptimes), len(waittimes),3))
    finalstate_q2_simple_xyz = np.zeros((len(ramptimes), len(waittimes),3))
    
    
    jt = 0
    for it in range(len(waittimes)):
        print(f'it = {it}')
        
        ramptime =  ramptimes[jt]  # unit: ns
        waittime =  waittimes[it]  # unit: ns
    
        ts = np.linspace(0, ramptime+waittime+ramptime, int((ramptime+waittime+ramptime)*M_pts)  +1) 
        p2 = pulse3(ts, (0, b_0), (ramptime, b_1), (ramptime+waittime, b_1), (ramptime+waittime+ramptime, b_0))   
        # vP3 = 39.6+p2*(65-39.6)    
        if pair == '14':    
            vP3 = 39.6+p2*(65-39.6) 
        else:
            vP3 = 27+p2*(46-27) 
        
        
        finalstate_q2 = np.zeros((len(ramptimes), len(waittimes),len(ts),4))
        finalstate_q2_reduction = np.zeros((len(ramptimes), len(waittimes),len(ts),2))
        finalstate_q2_simple = np.zeros((len(ramptimes), len(waittimes),len(ts),2))
        
        
        t = t_in
        phi0 = phi_in*np.pi/180
        
        if pair == '14':
            e_offset = 54.8  
            leverarm = 0.078*1 + 0.094*28.5/25.4
        else:
            e_offset = 36.6 
            leverarm = (0.084 + 0.08)
        detuning = vP3 
        
        es =   0.5* (detuning - e_offset) * leverarm * (1e-3*1.602e-19 / (1e9*6.62e-34))  
        evs = np.zeros((len(es),4))
        
        
        if pair == '14':
            q3 = 0.0715   
            q2 = 0.0338    
        else:
            q3 = 0.0620   
            q2 = 0.0709        

        # QuTiP 4 (old): Create a spline object
        # S_e = Cubic_Spline(ts[0], ts[-1], es)

        # QuTiP 5 (new): The NumPy array 'es' is used directly as the coefficient.
    
        Hsxs0 = Qobj(np.kron(s_x, s_0))
        Hszs0 = Qobj(np.kron(s_z, s_0))
        H_Ez = Qobj(np.array([[ q2,   0,              0,                 0],
                            [  0,  -q2,             0,                 0],
                            [  0,   0, q3*np.cos(phi0),  q3*np.sin(phi0)],
                            [  0,   0, q3*np.sin(phi0), -q3*np.cos(phi0)]]) /2)
    
        # QuTiP 4 (old): H = [..., [..., S_e]]
        # H = [2*np.pi*t*Hsxs0 + 2*np.pi*H_Ez, 
        #      [2*np.pi*Hszs0, S_e],   
        #       ]
        
        # QuTiP 5 (new): The coefficient is the NumPy array 'es' itself.
        # The solver will map it to the time points 'ts'.
        H = [2*np.pi*t*Hsxs0 + 2*np.pi*H_Ez, 
             [2*np.pi*Hszs0, es],   
              ]
        
        
        # QuTiP 4 (old): H[1][1] was the spline function, so you could call it: H[1][1](ts[0])
        # evals, ekets = (H[0] + H[1][0]*H[1][1](ts[0])  ).eigenstates()
        
        # QuTiP 5 (new): H[1][1] is now the NumPy array 'es'. We access its value by indexing.
        evals, ekets = (H[0] + H[1][0] * es[0]).eigenstates()
        psi0 = ekets[0] 
         
        # output = mesolve(H, psi0, ts)
        # Uprop = propagator(H, ts) 
        # NEW, CORRECTED CODE
        H_evo = QobjEvo(H, tlist=ts)  # Manually create the time-dependent object
        output = mesolve(H_evo, psi0, ts)
        Uprop = propagator(H_evo, ts)
        if print_residual:
            r = range(len(ts))
        else: 
            r = [len(ts)-1]
        for idx in r:        
            # QuTiP 5 (new): Access the coefficient value by indexing the 'es' array at 'idx'.
            evals, ekets = (H[0] + H[1][0] * es[idx]).eigenstates()
            finalstate_q2[jt,it,idx, 0] = abs(  output.states[idx].overlap( ekets[0])   )**2
            finalstate_q2[jt,it,idx,1] = abs(  output.states[idx].overlap( ekets[1])   )**2 
            finalstate_q2[jt,it,idx,2] = abs(  output.states[idx].overlap( ekets[2])   )**2 
            finalstate_q2[jt,it,idx,3] = abs(  output.states[idx].overlap( ekets[3])   )**2 
    
        finalstate_q2_f[jt,it,:] = np.array(finalstate_q2[jt,it,-1,:])
        #########################################################
        #########################################################
        
        
        sinetheta2 = (1-es/np.sqrt(t**2 + es**2))/2
        # QuTiP 4 (old):
        # S_sinetheta2 = Cubic_Spline(ts[0], ts[-1], sinetheta2)

        # QuTiP 5 (new): The NumPy array 'sinetheta2' is used directly.

        H_q2 = Qobj( q2* np.array([[ 1,   0,      ],
                            [  0,  -1,       ]])/2  )
        H_q3 = Qobj( q3 *   np.array([[  np.cos(phi0),  np.sin(phi0)],
                                      [  np.sin(phi0), -np.cos(phi0)]])/2 )
        
        
        # QuTiP 5 (new): Use the 'sinetheta2' array directly as the coefficient.
        H_r = [2*np.pi*H_q3 , 
             [ 2*np.pi*(-H_q3 + H_q2), sinetheta2],   
              ]
        
        
        # QuTiP 5 (new): Access the coefficient value by indexing at [0].
        evals, ekets = (H_r[0] + H_r[1][0] * sinetheta2[0]).eigenstates()
        psi0 = ekets[0] 
            
        # output_r = mesolve(H_r, psi0, ts)
        # Uprop_r  = propagator(H_r, ts) 
        # NEW, CORRECTED CODE
        H_r_evo = QobjEvo(H_r, tlist=ts) # Manually create the time-dependent object
        output_r = mesolve(H_r_evo, psi0, ts)
        Uprop_r  = propagator(H_r_evo, ts)
        
        if print_residual:
            r = range(len(ts))
        else: 
            r = [len(ts)-1]
        for idx in r:
            # QuTiP 5 (new): Access the coefficient value by indexing at [idx].
            evals, ekets = (H_r[0] + H_r[1][0] * sinetheta2[idx]).eigenstates()
            finalstate_q2_reduction[jt,it,idx,0] = abs(  output_r.states[idx].overlap( ekets[0])   )**2
            finalstate_q2_reduction[jt,it,idx,1] = abs(  output_r.states[idx].overlap( ekets[1])   )**2  
        finalstate_q2_reduction_f[jt,it,:] = finalstate_q2_reduction[jt,it,-1,:]
            
        uvecs = np.zeros((len(Uprop_r),3), )
        thetas = np.zeros((len(Uprop_r),), )
        
        
        if print_residual:
            r = range(len(ts))
        else: 
            r = [len(ts)-1]
        for i in r:    
            temp = Uprop_r[i].full()
            direction, rot_angle = SU2decompose(temp)
            uvecs[i, :] = direction
            thetas[i] = rot_angle
    
        #########################################################
        #########################################################
        
    
        e_temp = vP3 - e_offset
        i_23AC, i_32AC = np.where( e_temp[:-1] * e_temp[1:] <0 )[0]
        
        t_23AC = ( ts[i_23AC+1]*abs(e_temp[i_23AC]) + ts[i_23AC]*abs(e_temp[i_23AC+1]) )/( abs(e_temp[i_23AC]) + abs(e_temp[i_23AC+1]) )
        t_32AC = ( ts[i_32AC+1]*abs(e_temp[i_32AC]) + ts[i_32AC]*abs(e_temp[i_32AC+1]) )/( abs(e_temp[i_32AC]) + abs(e_temp[i_32AC+1]) )
        
        # QuTiP 5 (new): Access coefficients by indexing the 'sinetheta2' array.
        H_q2_simple = (H_r[0] + H_r[1][0] * sinetheta2[0]) / (2*np.pi)
        H_q3_simple = (H_r[0] + H_r[1][0] * sinetheta2[ramptime * M_pts]) / (2*np.pi)
        U_0_q2 = scipy.linalg.expm(  -2*np.pi*1j* H_q2_simple.full() * (t_23AC-ts[0])  )  
        U_1_q3 = scipy.linalg.expm(  -2*np.pi*1j* H_q3_simple.full() * (t_32AC-t_23AC)  ) 
        U_2_q2 = scipy.linalg.expm(  -2*np.pi*1j* H_q2_simple.full() * (ts[-1]-t_32AC)  ) 
        
        
        U_ramp_q2q3 = scipy.linalg.expm(  -2*np.pi*1j* H_q3.full() * (ramptime-t_23AC)  )     @     scipy.linalg.expm(  -2*np.pi*1j* H_q2.full() * (t_23AC-ts[0])  )  
    
        # QuTiP 5 (new): Access coefficient by indexing the 'sinetheta2' array.
        evals, ekets = (H_r[0] + H_r[1][0] * sinetheta2[0]).eigenstates()
        psi0 = ekets[0].full()
    
        
        if print_residual:
            r = range(len(ts))
        else: 
            r = [len(ts)-1]
        for idx in r:
            if idx <= i_23AC:
                state = scipy.linalg.expm(  -2*np.pi*1j* H_q2_simple.full() * (ts[idx]-ts[0])  )  @ psi0
            elif idx <= i_32AC:
                state = scipy.linalg.expm(  -2*np.pi*1j* H_q3_simple.full() * (ts[idx]-t_23AC)  )  @ U_0_q2 @ psi0
            else:
                state = scipy.linalg.expm(  -2*np.pi*1j* H_q2_simple.full() * (ts[idx]-t_32AC)  ) @ U_1_q3  @ U_0_q2 @ psi0
    
            # QuTiP 5 (new): Access coefficient by indexing the 'sinetheta2' array.
            evals, ekets = (H_r[0] + H_r[1][0] * sinetheta2[idx]).eigenstates()
            finalstate_q2_simple[jt,it,idx,0] = abs(  np.dot(state.T , ekets[0].full())   )**2
            finalstate_q2_simple[jt,it,idx,1] = abs(  np.dot(state.T , ekets[1].full())    )**2  
        
        finalstate_q2_simple_f[jt,it,:] = finalstate_q2_simple[jt,it,-1,:]    
    if print_residual == False:
        results['print_residual_False'] = dict(finalstate_q2_f=finalstate_q2_f,
                        finalstate_q2_reduction_f=finalstate_q2_reduction_f,
                         finalstate_q2_simple_f=finalstate_q2_simple_f)  
        results['print_residual_False']['waittimes'] = waittimes        
    else:
        results['print_residual_True'] = dict(finalstate_q2=finalstate_q2,
                        finalstate_q2_reduction=finalstate_q2_reduction,
                         finalstate_q2_simple=finalstate_q2_simple) 
        results['print_residual_True']['ts'] = ts     

# %%

results['pair'] = pair

import pickle
# Using a more standard filename format
pickle.dump(results , open('FigSupp5_qutip_simulation_results_QuTip5_23.pickle', 'wb'))

print("Simulation finished and results saved to FigSupp5_qutip_simulation_results_QuTip5.pickle")