# this file is modiefied from Supp/T2/T2driven_fig1_v1.py
import matplotlib.pyplot as plt
import numpy as np

import random

from scipy.optimize import curve_fit


# values obtained from FigSupp10
T2_driven_q1_err = np.array([  
[1, 2.34, 2.62, 0.00472, 0.305],
[2, 3.72, 13.1, 0.00142, 0.545],
[3, 5.39, 19.1, 0.000531, 1.074],
[5, 8.75, 18.0, 0.000574, 0.889],
[10, 17.2, 13.3, 0.000791, 0.419],
[25, 42.6, 6.96, 0.00122, 0.177],
[40, 68.0, 4.56, 0.00195, 0.127]
])

T2_driven_q2_err = np.array([  
[1, 4.46, 4.85, 0.00409,  0.662],
[2, 7.72, 13.5, 0.000822,  0.832],
[3, 11.3, 15.4, 0.000805,  0.571],
[5, 18.4, 12.3, 0.00126,  0.606],
[10, 36.2, 8.98, 0.00153,  0.317],
[25, 89.5, 4.22, 0.00200,  0.154],
[40, 143, 2.72, 0.00291,  0.081]
])


#%%
from scipy.optimize import curve_fit
def B2freq(B,B0,B2f,f_orth):
    return np.sqrt( (B2f*(B+B0))**2 + f_orth**2)
    

plt.figure( figsize=(5,5))

x = T2_driven_q1_err[:,0]
y = T2_driven_q1_err[:,1]
y_err = T2_driven_q1_err[:,3]
plt.plot(x, y, 'o')
coeffs = np.polyfit(x, y, deg=1)
p0 = [0.15, 1.68, 0.4]
popt, pcov = curve_fit(B2freq, x, y, p0=p0 )
# popt, pcov = curve_fit(B2freq, x, y, p0=p0, sigma=y_err, absolute_sigma=True)

xxx = np.linspace(0, 40, 401)
plt.plot(xxx, np.poly1d( coeffs)(xxx) )
plt.plot(xxx, B2freq(xxx,*popt) )
print( y - np.poly1d( coeffs)(x)   )
print( y - B2freq(x,*popt)   )
print(popt ,  ' +- ' , np.sqrt(np.diag(pcov))   )
# [0.08367998 1.69660382 1.3762088 ]  +-  [0.02992508 0.00191272 0.12206411]
q1_param = popt



x = T2_driven_q2_err[:,0]
y = T2_driven_q2_err[:,1]
y_err = T2_driven_q2_err[:,3]
plt.plot(x, y, 'o')
coeffs = np.polyfit(x, y, deg=1)
p0 = [0.15, 1.68, 0.4]
popt, pcov = curve_fit(B2freq, x, y, p0=p0)
# popt, pcov = curve_fit(B2freq, x, y, p0=p0, sigma=y_err, absolute_sigma=True)
q2_param = popt

xxx = np.linspace(0, 40, 401)
plt.plot(xxx, np.poly1d( coeffs)(xxx) )
plt.plot(xxx, B2freq(xxx,*popt) )
print( y - np.poly1d( coeffs)(x)   )
print( y - B2freq(x,*popt)   )
print(popt ,  ' +- ' , np.sqrt(np.diag(pcov))   )
# [0.13183808 3.56253521 1.81241179]  +-  [0.01799204 0.00242655 0.24078626]
plt.xlim(0,5.5)
plt.ylim(0, 20)

#%%

CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
MRSIZE = 2

from mpl_toolkits.axes_grid1 import Divider, Size
from functools import partial

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.2), Size.Fixed(0.7), Size.Fixed(1.2), Size.Fixed(0.7), Size.Fixed(1.5)]
v = [Size.Fixed(0.5), Size.Fixed(1.0), Size.Fixed(0.7), Size.Fixed(1.0)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax11 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=1))


def q1T2_v1(B, dg, fn, fn_par, params):
    B2f, B0, f_orth = params['B2f'], params['B0'], params['f_orth']
    # B2f = 1.6966
    # B0 = 0.08368
    # f_orth =  1.376 
    dfn = np.sqrt( (B2f*(B+B0))**2 + (fn+f_orth)**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    dfso = np.sqrt( (B2f*(B+B0)*(1+dg))**2 + f_orth**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    dfn_par = np.sqrt( (B2f*(B+B0)+fn_par)**2 + f_orth**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    df = np.sqrt(dfso**2 + dfn**2 + dfn_par**2)
    return     np.sqrt(2)/(2*np.pi*df)

def q2T2_v1(B, dg, fn, fn_par, params):
    B2f, B0, f_orth = params['B2f'], params['B0'], params['f_orth']
    # B2f = 3.562
    # B0 = 0.1318
    # f_orth =  1.812 
    dfn = np.sqrt( (B2f*(B+B0))**2 + (fn+f_orth)**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    dfso = np.sqrt( (B2f*(B+B0)*(1+dg))**2 + f_orth**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    dfn_par = np.sqrt( (B2f*(B+B0)+fn_par)**2 + f_orth**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    df = np.sqrt(dfso**2 + dfn**2 + dfn_par**2)
    return     np.sqrt(2)/(2*np.pi*df)


q1T2_v3 = partial(q1T2_v1, fn_par=0, params=dict(B0= q1_param[0], B2f=q1_param[1], f_orth=q1_param[2] ))
q2T2_v3 = partial(q2T2_v1, fn_par=0, params=dict(B0= q2_param[0], B2f=q2_param[1], f_orth=q2_param[2] ))




x = T2_driven_q1_err[:,0]
y = T2_driven_q1_err[:,2]

yerr = T2_driven_q1_err[:,4]
for i in range(len(x)):
    plt.errorbar(x[i], y[i],
                 yerr = yerr[i], c=CQ1,
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=0.75)

    
    




p0 = [0.0008, 0.34,]
popt, pcov = curve_fit( q1T2_v3, x, y, p0=p0)


Bs = np.linspace(1,50,101)  
plt.plot(Bs, q1T2_v3(Bs, *popt),  c=CQ1)
print('q1', popt)
print('q1', np.sqrt(np.diag(pcov)))



x = T2_driven_q2_err[:,0]
y = T2_driven_q2_err[:,2]

yerr = T2_driven_q2_err[:,4]
for i in range(len(x)):
    plt.errorbar(x[i], y[i],
                 yerr = yerr[i], c=CQ2,
                 fmt ='o', markersize=MRSIZE, linestyle='', linewidth=0.75)



p0 = [0.0008, 0.34, ]
popt, pcov = curve_fit(q2T2_v3, x, y, p0=p0)


Bs = np.linspace(1,50,101)  

plt.plot(Bs, q2T2_v3(Bs, *popt),  c=CQ2)
print('q2', popt)
print('q2', np.sqrt(np.diag(pcov)))


ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('field (mT)')
ax.set_ylabel('T2 (us)')

ax.set_xlim(0.6, 60)
ax.set_ylim(1.5, 28)

# filename = 'Figure1g'   # 'T2driven_fig1_v1'
# plt.savefig(filename+'.pdf') 