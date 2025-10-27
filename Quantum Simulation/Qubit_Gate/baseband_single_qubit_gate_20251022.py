
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import Divider, Size
import pickle
import scipy.linalg
import numpy as np
import matplotlib.pyplot as plt
DEG = np.pi/180
s0 = np.eye(2)
sx = np.array([[0, 1],
               [1, 0]])
sy = np.array([[0, -1j],
               [1j, 0]])
sz = np.array([[1, 0],
               [0,-1]])

DEG = np.pi/180

#%%
d_exp = dict(xdata=xdata, ydata=ydata, zdatas=zdatas)
t_q2s = xdata
t_q3s = ydata

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
    init_state = np.array([0, 1])
    U323 = func_U(t_q2, t_q3, t_q20, t_q30, t_q2_res, theta_in, A, B,fq2, fq3 )
    stateN = U323 @ U323 @ init_state
    return (np.real(np.abs(stateN[0]))**2   ) * A + B


def func2_wrap(M, *args):
    
    x, y = M
    arr = np.zeros(x.shape)
    for i in range(len(x)):
        arr[i] = func2(x[i], y[i], *args)
    return arr
#%%
from scipy.optimize import curve_fit


###
t_q2sL, t_q2sR = 0, 30   
maskx = np.logical_and( t_q2sL <= t_q2s, t_q2s <= t_q2sR)
i_t_q2sL, i_t_q2sR = np.min(np.where(maskx)[0]), np.max(np.where(maskx)[0])
xdata = t_q2s[ i_t_q2sL: i_t_q2sR+1  ]

t_q3sL, t_q3sR = 0, 30  
masky = np.logical_and( t_q3sL <= t_q3s, t_q3s <= t_q3sR)
i_t_q3sL, i_t_q3sR = np.min(np.where(masky)[0]), np.max(np.where(masky)[0])
ydata = t_q3s[ i_t_q3sL: i_t_q3sR+1    ]
zdatas = d_exp['zdatas']
X, Y = np.meshgrid(xdata, ydata)
Z = zdatas[  i_t_q3sL: i_t_q3sR+1 ,   i_t_q2sL: i_t_q2sR+1  ]
###


xxdata = np.vstack((X.ravel(), Y.ravel()))
p0 = np.array([-2.28847028, -1.54215906,  1.16019542, 44.71842728,  0.84173585, 0.08989211,  0.070926  ,  0.06203767])  
bounds = (  (-np.inf, -np.inf, -np.inf, -np.inf, 0.8, 0.01, 0.064, 0.056), (np.inf, np.inf, np.inf, np.inf, 0.99, 0.15, 0.078, 0.068)   )
popt, pcov = curve_fit(func2_wrap, xxdata, Z.ravel(), p0, bounds = bounds )
pstd = np.sqrt(np.diag(pcov))
print(popt)
print(pstd)
# popt = [2.28847019 -1.54215921  1.1601954  44.71842756  0.84173585  0.08989211
#   0.070926    0.06203767]
# pstd = [4.56254326e-03 3.40969332e-03 2.44957553e-03 2.36251214e-02
#  9.58524976e-04 3.69579981e-04 1.38538549e-05 7.98846076e-06]


p0 = np.array([-2.28847028, -1.54215906,  1.16019542, 44.71842728,  0.84173585, 0.08989211,  0.070926  ,  0.06203767])  
popt, pcov = curve_fit(func2_wrap, xxdata, Z.ravel(), p0, bounds = bounds, sigma=[0.02]*len(Z.ravel()), absolute_sigma=True )
pstd = np.sqrt(np.diag(pcov))
print(popt)
print(pstd)
# popt = [-2.28847019 -1.54215921  1.1601954  44.71842756  0.84173585  0.08989211
#   0.070926    0.06203767]
# pstd = [2.86451636e-03 2.14071895e-03 1.53792502e-03 1.48326368e-02
#  6.01793848e-04 2.32034607e-04 8.69791015e-06 5.01542143e-06]


# Z_fit = func2_wrap(xxdata, *popt).reshape(len(ydata),len(xdata))

def func2_SU2decompose(mat):
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
    return gate_polar_angle, rot_angle


gate_polar_angles = np.zeros( (len(t_q2s), len(t_q3s), ) )
rot_angles = np.zeros( (len(t_q2s), len(t_q3s), ) )
Pups = np.zeros( (len(t_q2s), len(t_q3s), ) )
for i2 in range(len(t_q2s)):
    for i3 in range(len(t_q3s)):
        # U323 = func_U(t_q2s[i2], t_q3s[i3], *popt)
        gate_polar_angle, rot_angle = func2_SU2decompose(func_U(t_q2s[i2], t_q3s[i3], *popt))
        Pup = func2(t_q2s[i2], t_q3s[i3],  *popt)
        gate_polar_angles[i2, i3] = gate_polar_angle
        rot_angles[i2, i3] = rot_angle
        Pups[i2, i3] = Pup
  

fig = plt.figure( figsize=(7,3.5))
h = [Size.Fixed(1), Size.Fixed(1.3), Size.Fixed(1), Size.Fixed(1.3)]
v = [Size.Fixed(0.5), Size.Fixed(1.3)]

divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1))
ax31 = fig.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=3, ny=1))

ax = ax11
img = ax.pcolormesh( xdata, ydata, Z, shading='auto', rasterized=True, vmin=0, vmax=1)

ax.set_xticks([0,  15, 30])
ax.set_yticks([0,  15, 30])
ax.set_xlabel(r'$\rm t_2 (ns)$')
ax.set_ylabel(r'$\rm t_3 (ns)$')


ax = ax31
img = ax.pcolormesh( xdata, ydata, Pups.T, shading='auto', rasterized=True, vmin=0, vmax=1)


cb=fig.colorbar(img, location='right',ticks=[0.2,0.8], shrink=0.12, aspect=5, )
cb.ax.tick_params(labelsize=12) 
ax.set_xticks([0,  15, 30])
ax.set_yticks([0,  15, 30])
ax.set_xlabel(r'$\rm t_2 (ns)$')
ax.set_ylabel(r'$\rm t_3 (ns)$')


cs = plt.contour( t_q2s, t_q3s, gate_polar_angles.T/DEG, levels=[90,] , colors='C9' , linewidths=[0.5,]  )
cs = plt.contour( t_q2s, t_q3s, rot_angles.T/DEG, levels=[90,]  , colors='C1' , linewidths=[0.5,] )    
cs = plt.contour( t_q2s, t_q3s, rot_angles.T/DEG, levels=[270,]  , colors='C0' , linewidths=[0.5,] ) 

plt.show()

# filename = 'FigureSupp6c' # 'Supp_fit2D_contour_q2x90'
# plt.savefig(filename+'.pdf')





