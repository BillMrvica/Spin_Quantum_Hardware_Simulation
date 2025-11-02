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
dtq1 = 1.9489
dtq4 = 1.9411
DEG = np.pi/180
fq1 = 3.38324710e-02*5/4  
fq4     = 7.15275643e-02*5/4   # qubit angular frequency, unit in GHz
theta = 41.517 * DEG       # q3 angle from sigma_z axis (toward sigma_x axis)
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


def Q1Xgate_uneq(fq1, fq4, dt_es, theta=41.517 * DEG):
    H1_0 = (1/2) *  sz
    H4_0 = (1/2) * (np.cos(theta)*sz + np.sin(theta)*sx   )
    Tq1 = 23.645922876871747
    output = scipy.linalg.expm(  -2*np.pi*1j* H1_0 *fq1* ( Tq1 +dt_es )  )
    output = scipy.linalg.expm(  -2*np.pi*1j* H4_0 *fq4* (wait_q4_1+dtq4  - 2*dt_es)     )  @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H1_0 *fq1* (wait_q1+dtq1 +2*dt_es )  ) @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H4_0 *fq4* (wait_q4_2+dtq4 -2*dt_es )  ) @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H1_0 *fq1* (add_wait_q1+dtq1 + Tq1 +dt_es )  ) @ output
    # U = scipy.linalg.expm(  -1j* 0.5*sz* (-109.87825962836682)*np.pi/180    )
    # return U.T.conj() @  output @ U
    
   
    return output


def Q1Xgate_eq(fq1, fq4, dt_es, theta=41.517 * DEG):
    fq1, fq4 = 0.8*fq1, 0.8*fq4
    H1_0 = (1/2) *  sz
    H4_0 = (1/2) * (np.cos(theta)*sz + np.sin(theta)*sx   )
    Tq1 = 23.645922876871747/0.8
    output = scipy.linalg.expm(  -2*np.pi*1j* H1_0 *fq1* ( Tq1 + dtq1/2 +dt_es )  )
    output = scipy.linalg.expm(  -2*np.pi*1j* H4_0 *fq4* (9.16+dtq4  - 2*dt_es)     )  @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H1_0 *fq1* (4.744+dtq1 +2*dt_es )  ) @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H4_0 *fq4* (9.16+dtq4 -2*dt_es )  ) @ output
    output = scipy.linalg.expm(  -2*np.pi*1j* H1_0 *fq1* ( Tq1 +dtq1/2 +dt_es )  ) @ output
    # U = scipy.linalg.expm(  -1j* 0.5*sz* (-109.87825962836682)*np.pi/180    )
    # return U.T.conj() @  output @ U
     
    return output

#%%
dfq1, dfq4 = 0.1e-3, 0.1e-3    # unit in GHz
 
wait_q4_1=3.747 
wait_q1 = 19.33 
wait_q4_2 = 10.17  
add_wait_q1 =  9.4 

Q1Xgate = Q1Xgate_uneq
print('\n' + 'qubit A unequal time ')
print(  *SU2decompose(  Q1Xgate(fq1, fq4, 0)  )   )
print('Dot1 Larmor frequency fluctuation ')
print(  (SU2decompose(  Q1Xgate(fq1+dfq1, fq4, 0)  ) - SU2decompose(  Q1Xgate(fq1-dfq1, fq4, 0)  ) )/2   )
print('Dot4 Larmor frequency fluctuation ')
print(  (SU2decompose(  Q1Xgate(fq1, fq4+dfq4, 0)  ) - SU2decompose(  Q1Xgate(fq1, fq4-dfq4, 0)  ) )/2   )

dtheta = 0.1*DEG  
print('quantization axis angle theta_14 fluctuation')
print(  (SU2decompose(  Q1Xgate(fq1, fq4, 0, theta=theta+dtheta)  ) - SU2decompose(  Q1Xgate(fq1, fq4, 0, theta=theta-dtheta )  ) )/2   )

# detuning noise 1ueV
# es[0] = -337.43 GHz
# es[600] = 226.43 GHz
# des = 0.24169 GHz
# -> dt = tramp* 0.24169/(337.43+226.43) = 0.000857269 ns (in the case of tramp = 2ns)
# dt_es = 0.000857269 
dt_es = 0.00857269 # detuning noise 10ueV
print('detuning noise epsilon_14')
print(  (SU2decompose(  Q1Xgate(fq1, fq4, dt_es)  ) - SU2decompose(  Q1Xgate(fq1, fq4, -dt_es)  ) )/2   )


#%%
dfq1, dfq4 = 0.1e-3, 0.1e-3    # unit in GHz

Q1Xgate = Q1Xgate_eq
print('\n' + 'qubit A equal time ')
print(  *SU2decompose(  Q1Xgate(fq1, fq4, 0)  )   )
print('Dot1 Larmor frequency fluctuation ')
print(  (5/4)*(SU2decompose(  Q1Xgate(fq1+dfq1, fq4, 0)  ) - SU2decompose(  Q1Xgate(fq1-dfq1, fq4, 0)  ) )/2   )
print('Dot4 Larmor frequency fluctuation ')
print(  (5/4)*(SU2decompose(  Q1Xgate(fq1, fq4+dfq4, 0)  ) - SU2decompose(  Q1Xgate(fq1, fq4-dfq4, 0)  ) )/2   )


dtheta = 0.1*DEG 
print('quantization axis angle theta_14 fluctuation')
print(  (SU2decompose(  Q1Xgate(fq1, fq4, 0, theta=theta+dtheta)  ) - SU2decompose(  Q1Xgate(fq1, fq4, 0, theta=theta-dtheta )  ) )/2   )


# detuning noise 1ueV
# es[0] = -337.43 GHz
# es[600] = 226.43 GHz
# des = 0.24169 GHz
# -> dt = tramp* 0.24169/(337.43+226.43) = 0.000857269 ns (in the case of tramp = 2ns)
# dt_es = 0.000857269
dt_es = 0.00857269 # detuning noise 10ueV
print('detuning noise epsilon_14')
print(  (5/4)*(SU2decompose(  Q1Xgate(fq1, fq4, dt_es)  ) - SU2decompose(  Q1Xgate(fq1, fq4, -dt_es)  ) )/2   )





