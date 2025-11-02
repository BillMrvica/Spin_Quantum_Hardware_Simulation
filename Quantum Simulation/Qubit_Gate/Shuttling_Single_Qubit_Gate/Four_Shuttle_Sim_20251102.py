import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import Divider, Size
import scipy.linalg

DEG = np.pi/180
s0 = np.eye(2)
sx = np.array([[0, 1],
               [1, 0]])
sy = np.array([[0, -1j],
               [1j, 0]])
sz = np.array([[1, 0],
               [0,-1]])

# 
def func_U(t_q2, t_q3, t_q20, t_q30, t_q2_res, theta_in, A, B, fq2, fq3 ):
    theta, phi = DEG*theta_in, DEG*0 

    qax = np.array([np.sin(theta)*np.cos(phi), 
                    np.sin(theta)*np.sin(phi),
                    np.cos(theta)])
    
    qx, qy, qz = qax
    dHq3 = np.array([[      qz, qx-1j*qy ],
                   [qx+1j*qy,      -qz]])  *  fq3 * 0.5 
    U3 = scipy.linalg.expm(  -2*np.pi*1j* dHq3 * (t_q3-t_q30)  )

    dHq2 = np.array([[      1,       0 ],
                   [      0,      -1 ]])  *  fq2 * 0.5 
    U2 = scipy.linalg.expm(  -2*np.pi*1j* dHq2 * (t_q2-t_q20)  )  
    U2_residual = scipy.linalg.expm(  -2*np.pi*1j* dHq2 * t_q2_res  )  
    U323 = U2_residual @ U3 @ U2 @ U3  @ U2_residual
    return U323

def func2(t_q2, t_q3, t_q20, t_q30, t_q2_res, theta_in, A, B, fq2, fq3 ):      
    init_state = np.array([0, 1]) # Initial state |1> (spin-down)
    U323 = func_U(t_q2, t_q3, t_q20, t_q30, t_q2_res, theta_in, A, B, fq2, fq3 )
    stateN = U323 @ U323 @ init_state # Apply the pi/2 single-qubit gate twice
    # Probability of being in the |0> state (spin-up)
    Pup = (np.real(np.abs(stateN[0]))**2 ) * A + B
    return Pup

def func2_wrap(M, *args):
    x, y = M
    arr = np.zeros(x.shape)
    for i in range(len(x)):
        arr[i] = func2(x[i], y[i], *args)
    return arr

def func2_SU2decompose(mat):
    u0 = np.real( np.trace( mat @ s0 )/2  )
    ux = np.imag( np.trace( mat @ sx )/2  )
    uy = np.imag( np.trace( mat @ sy )/2  )
    uz = np.imag( np.trace( mat @ sz )/2  )
    
    uvec = np.array([ux, uy, uz])
    if np.sqrt(np.sum(uvec**2)) <= 1e-14:
        uvec = np.array([0,0,1])
    else:
        uvec = uvec/np.sqrt(np.sum(uvec**2))
    rot_angle = 2*np.arccos( u0)
    gate_polar_angle = np.arctan2( np.sqrt(uvec[0]**2 + uvec[1]**2), uvec[2]   ) 
    return gate_polar_angle, rot_angle

# Fitted parameters from Chien-An
t_q20 = -2.28847019
t_q30 = -1.54215921
t_q2_res = 1.1601954
theta_in = 44.71842756
A = 0.84173585
B = 0.08989211
fq2 = 0.070926
fq3 = 0.06203767
popt = [t_q20, t_q30,  t_q2_res,  theta_in,  A,  B,  fq2,  fq3 ]

# Define the range for t_q2 and t_q3
t_q2s = np.linspace(0, 30, 120)
t_q3s = np.linspace(0, 30, 120)

# Initialize arrays to store results
gate_polar_angles = np.zeros( (len(t_q2s), len(t_q3s)) )
rot_angles = np.zeros( (len(t_q2s), len(t_q3s)) )
Pups = np.zeros( (len(t_q2s), len(t_q3s)) )

# Calculate Pups and gate parameters for each combination of t_q2 and t_q3
print("Calculating Pups and gate parameters...")
for i2 in range(len(t_q2s)):
    for i3 in range(len(t_q3s)):
        U_matrix = func_U(t_q2s[i2], t_q3s[i3], *popt)
        gate_polar_angle, rot_angle = func2_SU2decompose(U_matrix)
        Pup = func2(t_q2s[i2], t_q3s[i3],  *popt)
        gate_polar_angles[i2, i3] = gate_polar_angle
        rot_angles[i2, i3] = rot_angle
        Pups[i2, i3] = Pup
print("Calculation complete.")

# --- Plotting Pups with contours ---
fig, ax = plt.subplots(figsize=(8, 6))

# Use pcolormesh for a 2D plot of Pups
im = ax.pcolormesh(t_q2s, t_q3s, Pups.T, shading='auto', cmap='viridis') # Transpose Pups for correct orientation

# Add a color bar
cbar = fig.colorbar(im, ax=ax, orientation='vertical')
cbar.set_label('Spin-Up Probability (Pups)')

# Add contours
cs = ax.contour( t_q2s, t_q3s, gate_polar_angles.T/DEG, levels=[90,] , colors='C9' , linewidths=[0.5,]  )
cs = ax.contour( t_q2s, t_q3s, rot_angles.T/DEG, levels=[90,]  , colors='C1' , linewidths=[0.5,] )    
cs = ax.contour( t_q2s, t_q3s, rot_angles.T/DEG, levels=[270,]  , colors='C0' , linewidths=[0.5,] ) 

# Set labels and title
ax.set_xlabel(r'$\rm t_2 (ns)$')
ax.set_ylabel(r'$\rm t_3 (ns)$')
ax.set_title('Simulated Spin-Up Probability (Pups) with Contours')

# Set ticks (optional, as pcolormesh handles it well)
ax.set_xticks(np.arange(0, 31, 5))
ax.set_yticks(np.arange(0, 31, 5))

# Ensure tight layout
plt.tight_layout()
plt.show()

# The following lines from the original code are for debugging/printing specific values.
print(func_U(10, 10, *popt))
print(func2(10, 10, *popt))
print(func2_SU2decompose(func_U(10, 10, *popt))[0]*180/np.pi, func2_SU2decompose(func_U(10, 10, *popt))[1]*180/np.pi)