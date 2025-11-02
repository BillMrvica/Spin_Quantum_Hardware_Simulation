# this file is modified from Supp/fit_tc_angle/H4x4_example.py
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import time
from scipy.optimize import fsolve
# from tqdm import tqdm 
import scipy
import pickle
# file = open('extracted_Jx_Jy.pickle','rb')
# file = open('fit_params.pickle', 'rb')
# d = pickle.load(file)


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

print_residual = True
# print_residual = False


if pair  not in ['14', '23']:
    raise ValueError

if pair == '14':
    t_in = 27 #54 
    phi_in = 90 #65
else:
    t_in = 20 #40  
    # phi_in = 51
    # phi_ins = [0, 30, 60, 90, ]
    phi_ins = [0, 45, 90, 135]



p2 = np.linspace(0, 1, 1001)
# vP3 = 39.6+p2*(65-39.6)    
if pair == '14':    
    vP3 = 39.6+p2*(65-39.6) 
else:
    vP3 = 27+p2*(46-27) 
    
evs = np.zeros((len(phi_ins), len(p2),4))
for ii in range(len(phi_ins)):
    phi_in = phi_ins[ii]
    t = t_in
    phi0 = phi_in*np.pi/180
    
    if pair == '14':
        e_offset = 54.8  
        leverarm = 0.078*1 + 0.094*28.5/25.4
    else:
        e_offset = 36.6 
        leverarm = (0.084 + 0.08)
    detuning = vP3 #np.linspace(10, 50, 4001)
    
    es =   -0.5* (detuning - e_offset) * leverarm * (1e-3*1.602e-19 / (1e9*6.62e-34))  
    # evs = np.zeros((len(es),4))
    
    
    if pair == '14':
        q3 = 0.0715   
        q2 = 0.0338    
    else:
        q3 = 0.0620   *5/4
        q2 = 0.0709   *5/4     
        
    
    
    
    Hsxs0 = np.kron(s_x, s_0)
    Hszs0 = np.kron(s_z, s_0)
    H_Ez = np.array([[ q2,   0,              0,                 0],
                        [  0,  -q2,             0,                 0],
                        [  0,   0, q3*np.cos(phi0),  q3*np.sin(phi0)],
                        [  0,   0, q3*np.sin(phi0), -q3*np.cos(phi0)]]) /2
    
    
    for i in range(len(es)):
        H4x4 = t*Hsxs0 + H_Ez + Hszs0*es[i]
        ev, evec = np.linalg.eigh(H4x4)
        evs[ii,i,:] = ev.copy()




from mpl_toolkits.axes_grid1 import Divider, Size





fig = plt.figure( figsize=(10,4))
HSPACE, VSPACE = 0.5, 0.5
h = [ s for _ in range(30) for s in [Size.Fixed(HSPACE), Size.Fixed(HSPACE)]  ]
v = [ s for _ in range(8) for s in [Size.Fixed(VSPACE), Size.Fixed(VSPACE)]  ]

divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1, nx1=4, ny1=4))

ax31_up = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=6, ny=3, nx1=7, ny1=4))

ax31_dn = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=6, ny=1, nx1=7, ny1=2))

ax51 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=9, ny=1, nx1=12, ny1=4))

factor = 0.5* leverarm * (1e-3*1.602e-19 / (1e9*6.62e-34))

ax = ax11
ax.plot( es, evs[0,:,0:2], 'k-', linewidth=0.75)
ax.plot( es, evs[0,:,2:4], 'k-', linewidth=0.75)
ax.set_xlim(-200, 200)
ax.set_xticks([-200, 0, 200])


# ax = ax31_dn
# ax.plot( vP3, evs[0,:,0:2], 'k-', linewidth=0.75)
# ax.set_xlim(46-20/2000, 46+ 20/2000)
# ax.set_ylim(-187.64239479 -0.05, -187.55380149 + 0.15)
# # ax.set_xlim(46-20/1000, 46)
# # ax.set_xticks([27, 36.5, 46])

# ax = ax31_up
# ax.plot( vP3, evs[0,:,2:4], 'k-', linewidth=0.75)
# ax.set_xlim(46-20/2000, 46+ 20/2000)
# ax.set_ylim(187.55933229 -0.2, 187.63686399 + 0.2)
# ax.set_xticks([27, 36.5, 46])



ax = ax31_dn
ax.plot( es, evs[3,:,0:2], 'k-', linewidth=0.75)
ax.set_xlim( -2,  2)
ax.set_ylim( -20.18, -19.93)
# ax.set_xlim(46-20/1000, 46)
# ax.set_xticks([27, 36.5, 46])

ax = ax31_up
ax.plot( es, evs[3,:,2:4], 'k-', linewidth=0.75)
ax.set_xlim( -2,  2)
ax.set_ylim( 19.93, 20.18)


ax = ax51
for i in [3,2,1,0]:
    ax.plot( es, (evs[i,:,1]-evs[i,:,0]), linewidth=0.75 )

ax.set_xlim(-200, 200)
ax.set_xticks([-200, 0, 200])
# ax.set_ylim(0.057, 0.093)

plt.show()

# filename = 'FigureSupp3' # 'H4x4_example'
# plt.savefig(filename+'.pdf')
