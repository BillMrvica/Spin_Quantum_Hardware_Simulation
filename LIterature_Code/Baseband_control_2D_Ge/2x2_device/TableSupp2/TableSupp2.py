# this file is modified from Supp/unequal_time/q1x90_data_sim.py


#%% Import for plotting

import matplotlib.pyplot as plt
import numpy as np
import scipy
sx = np.array([[0,1],[1,0]],dtype=complex)
sy = np.array([[0,-1j],[1j,0]],dtype=complex)
sz = np.array([[1,0],[0,-1]],dtype=complex)
sm = np.array([[0,0],[1,0]],dtype=complex)
s0 = np.eye(2)


#%%
import numpy as np
import scipy.linalg
dtq2 = 2.28847028
dtq3 = 1.54215906
DEG = np.pi/180
fq2 = 0.070926 *5/4 
fq3     = 0.06203767*5/4   # qubit angular frequency
theta = 44.718 * DEG       # q3 angle from sigma_z axis (toward sigma_x axis)
phi = 0*DEG

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
    rot_angle = 2*np.arccos( u0)
    gate_polar_angle = np.arctan2( np.sqrt(uvec[0]**2 + uvec[1]**2), uvec[2]   ) 
    gate_azu_angle = np.arctan2( uvec[1], uvec[0] )
    # print( rot_angle/DEG,     gate_polar_angle /DEG , gate_azu_angle/DEG  )
    return np.array( [rot_angle/DEG,     gate_polar_angle /DEG , gate_azu_angle/DEG])


def Q2Xgate_uneq(fq2, fq3, dt_es, theta=44.718 * DEG):
    H2_0 = (1/2) *  sz
    H3_0 = (1/2) * (np.cos(theta)*sz + np.sin(theta)*sx   )
    Tq2 = 11.279361588134112
    output = scipy.linalg.expm(  -2*np.pi*1j* H2_0 *fq2* ( Tq2 + dt_es)  )
    output = scipy.linalg.expm(  -2*np.pi*1j* H3_0 *fq3* (wait_q3+dtq3  - 2*dt_es)     )  @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H2_0 *fq2* (wait_q2+dtq2 + Tq2 +dt_es )  ) @ output
    # U = scipy.linalg.expm(  -1j* 0.5*sz* ( )*np.pi/180    )
    # return U.T.conj() @  output @ U
    
    return output

def Q2Xgate_eq(fq2, fq3, dt_es, theta=44.718 * DEG):
    fq2, fq3 = 0.8*fq2, 0.8*fq3
    H2_0 = (1/2) *  sz
    H3_0 = (1/2) * (np.cos(theta)*sz + np.sin(theta)*sx   )
    Tq2 = 11.279361588134112 / 0.8
    output = scipy.linalg.expm(  -2*np.pi*1j* H2_0 *fq2 *  (Tq2+ dtq2/2 + dt_es)  ) 
    output = scipy.linalg.expm(  -2*np.pi*1j* H3_0 *fq3 *  (  4.94  +dtq3 - 2*dt_es) ) @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H2_0 *fq2 *  (  10.125  +dtq2 + 2*dt_es)  ) @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H3_0 *fq3 *  (  4.94  +dtq3 - 2*dt_es)  ) @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H2_0 *fq2 *  (Tq2 + dtq2/2 + dt_es)  ) @ output
    # U = scipy.linalg.expm(  -1j* 0.5*sz* ( )*np.pi/180    )
    # return U.T.conj() @  output @ U
    
    return output



#%%
dfq2, dfq3 = 0.1e-3, 0.1e-3    # unit in GHz

wait_q3 = 4.86 
wait_q2 = 3.42

Q2Xgate = Q2Xgate_uneq
print('\n' + 'qubit A unequal time ')
print(  *SU2decompose(  Q2Xgate(fq2, fq3, 0)  )   )
print('Dot1 Larmor frequency fluctuation ')
print(  (SU2decompose(  Q2Xgate(fq2+dfq2, fq3, 0)  ) - SU2decompose(  Q2Xgate(fq2-dfq2, fq3, 0)  ) )/2   )
print('Dot4 Larmor frequency fluctuation ')
print(  (SU2decompose(  Q2Xgate(fq2, fq3+dfq3, 0)  ) - SU2decompose(  Q2Xgate(fq2, fq3-dfq3, 0)  ) )/2   )

dtheta = 0.1*DEG  
print('quantization axis angle theta_14 fluctuation')
print(  (SU2decompose(  Q2Xgate(fq2, fq3, 0, theta=theta+dtheta)  ) - SU2decompose(  Q2Xgate(fq2, fq3, 0, theta=theta-dtheta )  ) )/2   )

# detuning noise 1ueV
# es[0] = -190.49 GHz
# es[600] = 186.52 GHz
# des = 0.24169 GHz
# -> dt = tramp* 0.24169/(190.49+186.52) = 0.00128214 ns (in the case of tramp = 2ns)
# dt_es = 0.00128214
dt_es = 0.0128214 # detuning noise 10ueV
print('detuning noise epsilon_14')
print(  (SU2decompose(  Q2Xgate(fq2, fq3, dt_es)  ) - SU2decompose(  Q2Xgate(fq2, fq3, -dt_es)  ) )/2   )


#%%
dfq2, dfq3 = 0.1e-3, 0.1e-3    # unit in GHz

Q2Xgate = Q2Xgate_eq
print('\n' + 'qubit A equal time ')
print(  *SU2decompose(  Q2Xgate(fq2, fq3, 0)  )   )
print('Dot2 Larmor frequency fluctuation ')
print(  (5/4)*(SU2decompose(  Q2Xgate(fq2+dfq2, fq3, 0)  ) - SU2decompose(  Q2Xgate(fq2-dfq2, fq3, 0)  ) )/2   )
print('Dot3 Larmor frequency fluctuation ')
print(  (5/4)*(SU2decompose(  Q2Xgate(fq2, fq3+dfq3, 0)  ) - SU2decompose(  Q2Xgate(fq2, fq3-dfq3, 0)  ) )/2   )


dtheta = 0.1*DEG 
print('quantization axis angle theta_23 fluctuation')
print(  (SU2decompose(  Q2Xgate(fq2, fq3, 0, theta=theta+dtheta)  ) - SU2decompose(  Q2Xgate(fq2, fq3, 0, theta=theta-dtheta )  ) )/2   )


# detuning noise 1ueV
# es[0] = -190.49 GHz
# es[600] = 186.52 GHz
# des = 0.24169 GHz
# -> dt = tramp* 0.24169/(190.49+186.52) = 0.00128214 ns (in the case of tramp = 2ns)
# dt_es = 0.00128214
dt_es = 0.0128214 # detuning noise 10ueV
print('detuning noise epsilon_23')
print(  (5/4)*(SU2decompose(  Q2Xgate(fq2, fq3, dt_es)  ) - SU2decompose(  Q2Xgate(fq2, fq3, -dt_es)  ) )/2   )













