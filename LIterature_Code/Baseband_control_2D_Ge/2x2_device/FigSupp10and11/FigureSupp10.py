# this file is modified from Supp/T2/T2_Ge2x2_1mT_40mT.py
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':15})

# last column: experiment time in minutes
T2_driven_q1 = np.array([  
[1, 2.34, 2.62, 0.984, 17.2, 3.48, 13.72],
[2, 3.72, 13.1, 2.73, 48.1, 1.35, 10.38],
[3, 5.39, 19.1, 1.31, 110, 0.91, 13.63],
[5, 8.75, 18.0, 1.45, 108, 1.29, 12.55],
[10, 17.2, 13.3, 2.02, 72.7, 1.54, 11.98],
[25, 42.6, 6.96, 2.02, 28.6, 1.28, 11.92],
[40, 68.0, 4.56, 1.93, 19.4, 1.47, 11.67]
])

T2_driven_q2 = np.array([  
[1, 4.46, 4.85, 1.11, 39.5, 2.09, 9.12],
[2, 7.72, 13.5, 1.36, 50.5, 1.81, 10.48],
[3, 11.3, 15.4, 1.99, 65.0, 0.90, 17.18],
[5, 18.4, 12.3, 1.92, 61.2, 0.97, 12.72],
[10, 36.2, 8.98, 2.29, 45.8, 1.05, 12.9],
[25, 89.5, 4.22, 1.52, 26.4, 1.48, 13.28],
[40, 143, 2.72, 1.67, 15.4, 1.22, 13.03]
])

# last column: fitting standard deviation of T2star
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

# last column: experiment time in minutes
T2star_persist_q1 = np.array([  
[5, 8.63, 24.1, 1.75, 12.35],
[10, 16.9, 17.1, 1.53, 13.82],
[20, 33.9, 9.97, 2.07, 17.8],
[25, 42.5, 7.00, 1.94, 15.67],
[40, 67.6, 4.94, 2.04, 16.77]
])

T2star_persist_q2 = np.array([  
[5, 18.1, 15.4, 1.54, 12.88],
[10, 35.4, 11.1, 1.66, 14.4],
[20, 70.9, 5.85, 1.61, 19.0],
[25, 89.4, 4.46, 1.68, 17.08],
[40, 141, 3.06, 1.79, 17.72]
])


T2Hahn_persist_q1 = np.array([  
[5, 8.63, 122, 2.09],
[10, 16.9, 83.7, 1.55],
[20, 33.9, 50.3, 1.32],
[25, 42.5, 32.1, 1.42],
])

T2Hahn_persist_q2 = np.array([  
[5, 18.1, 62.6, 1.02],
[10, 35.4, 50.0, 1.38],
[20, 70.9, 30.6, 1.31],
[25, 89.4, 23.8, 1.27],
])




T2cpmg_persist_q1 = np.array([
    [5, 8.63, 117, 205, 117, 105, 185, 358, 668, 1.2e+03, 1.99e+03, 3.29e+03, ],  # 5mT,  01-09-2023.pptx    
    [10, 16.9, 82.7, 100, 136, 117, 188, 331, 580, 1.03e+03, 1.71e+03, 2.79e+03],      # 10mT,  01-08-2023.pptx     
    [20, 33.9, 48.2, 71.4, 99.8, 108, 175, 306, 424, 770, 1.28e+03, 1.84e+03  ],        #  20mT, 22-07-2023.pptx    
    [25, 42.5, 31.9, 58.4, 84.1, 96.2, 148, 274, 473, 822, 1.25e+03, 2.10e+03],       # 25mT,  01-09-2023.pptx
    ])

T2cpmg_persist_q2 = np.array([
    [5, 18.1, 61.7, 68.8, 88.8, 130, 201, 364, 675, 1.21e+03, 2.17e+03, 3.65e+03 ] ,  # 5mT,  01-09-2023.pptx   
    [10, 35.4, 53.1, 64.1, 79, 123, 201, 328, 623, 1.07e+03, 1.91e+03, 3.02e+03,] ,    # 10mT,  01-08-2023.pptx   
    [20, 70.9, 29.2, 49.5, 70.2, 94.1, 166, 293, 410, 838, 1.27e+03, 2.12e+03] ,        #  20mT, 22-07-2023.pptx    
    [25, 89.4, 24.7, 37.4, 55.5, 83.3, 143, 256, 423, 708, 1.1e+03, 1.72e+03],       # 25mT,  01-09-2023.pptx
    ])






#%%
plt.rcParams.update({'lines.markersize':4})
plt.rcParams.update({'lines.linewidth': 0.5})
plt.rcParams.update({'font.size':15})

plt.figure( figsize=(14,8))
plt.subplots_adjust(hspace=0.3, wspace=0.3)



from scipy.optimize import curve_fit
def B2freq(B,B0,B2f,f_orth):
    return np.sqrt( (B2f*(B+B0))**2 + f_orth**2)

plt.subplot(231)
# x = T2_driven_q1[:,0]
# y = T2_driven_q1[:,1]
# plt.plot(x, y, 'bo')
# coeffs = np.polyfit(x, y, deg=1)
# print(coeffs)
# plt.plot(x, np.poly1d( coeffs)(x) , 'b')

# x = T2_driven_q2[:,0]
# y = T2_driven_q2[:,1]
# plt.plot(x, y, 'ro')
# coeffs = np.polyfit(x, y, deg=1)
# print(coeffs)
# plt.plot(x, np.poly1d( coeffs)(x) , 'r')


x = T2_driven_q1_err[:,0]
y = T2_driven_q1_err[:,1]
y_err = T2_driven_q1_err[:,3]
plt.plot(x, y, 'bo', markersize=5)
popt, pcov = curve_fit(B2freq, x, y, p0=[0.15, 1.68, 0.4] )
# popt, pcov = curve_fit(B2freq, x, y, p0=p0, sigma=y_err, absolute_sigma=True)

xxx = np.linspace(0, 40, 401)
plt.plot(xxx, B2freq(xxx,*popt) , 'b')
print(popt ,  ' +- ' , np.sqrt(np.diag(pcov))   )
# [0.08367998 1.69660382 1.3762088 ]  +-  [0.02992508 0.00191272 0.12206411]
q1_param = popt


x = T2_driven_q2_err[:,0]
y = T2_driven_q2_err[:,1]
y_err = T2_driven_q2_err[:,3]
plt.plot(x, y, 'ro', markersize=5)
popt, pcov = curve_fit(B2freq, x, y, p0=[0.15, 3.5, 0.4])
# popt, pcov = curve_fit(B2freq, x, y, p0=p0, sigma=y_err, absolute_sigma=True)

xxx = np.linspace(0, 40, 401)
plt.plot(xxx, B2freq(xxx,*popt) , 'r')
print(popt ,  ' +- ' , np.sqrt(np.diag(pcov))   )
# [0.13183808 3.56253521 1.81241179]  +-  [0.01799204 0.00242655 0.24078626]
q2_param = popt

plt.xlim(0,5.5)
plt.ylim(0, 20)



# plt.xscale('log')
# plt.yscale('log')
plt.xlabel('field (mT)')
plt.ylabel('frequency (MHz)')
plt.grid()





plt.subplot(232)
x = T2_driven_q1[:,0]
y = T2_driven_q1[:,1]
plt.plot(x, y, 'bo', markersize=5)
# coeffs = np.polyfit(x, y, deg=1)
# print(coeffs)
# plt.plot(x, np.poly1d( coeffs)(x) , 'b')
plt.plot(x, B2freq(x,*q1_param) , 'b')

x = T2_driven_q2[:,0]
y = T2_driven_q2[:,1]
plt.plot(x, y, 'ro', markersize=5)
# coeffs = np.polyfit(x, y, deg=1)
# print(coeffs)
# plt.plot(x, np.poly1d( coeffs)(x) , 'r')
plt.plot(x, B2freq(x,*q2_param) , 'r')

plt.xlim(0,)
plt.ylim(0,) 
plt.xlabel('field (mT)')
plt.ylabel('frequency (MHz)')
plt.grid()




################## copy from T2driven_fig1_v1.py ######################
plt.subplot(234)

from functools import partial
def q1T2_v1(B, dg, fn, fn_par, params):
    B2f, B0, f_orth = params['B2f'], params['B0'], params['f_orth']
    dfn = np.sqrt( (B2f*(B+B0))**2 + (fn+f_orth)**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    dfso = np.sqrt( (B2f*(B+B0)*(1+dg))**2 + f_orth**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    dfn_par = np.sqrt( (B2f*(B+B0)+fn_par)**2 + f_orth**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    df = np.sqrt(dfso**2 + dfn**2 + dfn_par**2)
    return     np.sqrt(2)/(2*np.pi*df)

def q2T2_v1(B, dg, fn, fn_par, params):
    B2f, B0, f_orth = params['B2f'], params['B0'], params['f_orth']
    dfn = np.sqrt( (B2f*(B+B0))**2 + (fn+f_orth)**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    dfso = np.sqrt( (B2f*(B+B0)*(1+dg))**2 + f_orth**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    dfn_par = np.sqrt( (B2f*(B+B0)+fn_par)**2 + f_orth**2 ) - np.sqrt( (B2f*(B+B0))**2 + f_orth**2  )
    df = np.sqrt(dfso**2 + dfn**2 + dfn_par**2)
    return     np.sqrt(2)/(2*np.pi*df)

q1T2_v3 = partial(q1T2_v1, fn_par=0, params=dict(B0= q1_param[0], B2f=q1_param[1], f_orth=q1_param[2] ))
q2T2_v3 = partial(q2T2_v1, fn_par=0, params=dict(B0= q2_param[0], B2f=q2_param[1], f_orth=q2_param[2] ))


x = T2_driven_q1_err[:,0]
y = T2_driven_q1_err[:,2]
plt.plot(x, y, 'o', c='b', markersize=4, )
yerr = T2_driven_q1_err[:,4]
# for i in range(len(x)):
#     plt.errorbar(x[i], y[i],
#                  yerr = yerr[i], c='b',
#                  fmt ='o', markersize=4, linestyle='', linewidth=0.75)

Bs = np.linspace(1,41,101)  
popt = np.array([ 0.00096053, -0.05152972])
plt.plot(Bs, q1T2_v3(Bs, *popt),  c='b')

x = T2_driven_q2_err[:,0]
y = T2_driven_q2_err[:,2]
plt.plot(x, y, 'o', c='r', markersize=4)
yerr = T2_driven_q2_err[:,4]
# for i in range(len(x)):
#     plt.errorbar(x[i], y[i],
#                  yerr = 0, c='r',
#                  fmt ='o', markersize=4, linestyle='', linewidth=0.75)

Bs = np.linspace(1,41,101)  
popt = np.array([ 0.00078803, -0.07833426])
plt.plot(Bs, q2T2_v3(Bs, *popt),  c='r')




# x = T2_driven_q1[:,0]
# y = T2_driven_q1[:,2]
# plt.plot(x, y, 'bo-')
# x = T2_driven_q2[:,0]
# y = T2_driven_q2[:,2]
# plt.plot(x, y, 'ro-')

x = T2_driven_q1[1:,0]
y = T2_driven_q1[1:,4]
plt.plot(x, y, 'bo--')
x = T2_driven_q2[:,0]
y = T2_driven_q2[:,4]
plt.plot(x, y, 'ro--')

plt.xlim(0.83,49)
# plt.ylim(0,) 
plt.xscale('log')
plt.yscale('log')
plt.xlabel('field (mT)')
plt.ylabel('T2 (us)')
plt.grid()




plt.subplot(235)
x = T2star_persist_q1[:,0]
y = T2star_persist_q1[:,2]
plt.plot(x, y, 'C0o-')
x = T2star_persist_q2[:,0]
y = T2star_persist_q2[:,2]
plt.plot(x, y, 'C1o-')

x = T2Hahn_persist_q1[:,0]
y = T2Hahn_persist_q1[:,2]
plt.plot(x, y, 'C0o--')
x = T2Hahn_persist_q2[:,0]
y = T2Hahn_persist_q2[:,2]
plt.plot(x, y, 'C1o--')

i=7
x = T2cpmg_persist_q1[:,0]
y = T2cpmg_persist_q1[:,i]
plt.plot(x, y, f'C0o-.')
x = T2cpmg_persist_q2[:,0]
y = T2cpmg_persist_q2[:,i]
plt.plot(x, y, f'C1o-.')

i=11
x = T2cpmg_persist_q1[:,0]
y = T2cpmg_persist_q1[:,i]
plt.plot(x, y, f'C0o:')
x = T2cpmg_persist_q2[:,0]
y = T2cpmg_persist_q2[:,i]
plt.plot(x, y, f'C1o:')
    

plt.xscale('log')
plt.yscale('log')
plt.xlabel('field (mT)')
plt.ylabel('T2 (us)')
plt.grid()


plt.rcParams.update({'lines.markersize':5})
plt.subplot(233)
i = 0
x = 2**np.arange(10)
y = T2cpmg_persist_q1[i,2:]
plt.plot(x, y, f'C0^-')
i = 3
x = 2**np.arange(10)
y = T2cpmg_persist_q1[i,2:]
plt.plot(x, y, f'C0v-')

# plt.plot( (16,512), ( 25* 16**0.5, 25* 512**0.5), 'k--')
# plt.plot( (32,512), ( 10* 32**(3/4), 10* 512**(3/4)  ), 'k--')
plt.plot( (32,512), ( 10* 32**(2/3), 10* 512**(2/3)  ), 'k--')
# plt.plot( (16,512), ( 25* 16**(2/3), 25* 512**(2/3) ), 'k--')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('N_pi')
plt.ylabel('T2 (us)')
plt.grid()


plt.subplot(236)
i = 0
x = 2**np.arange(10)
y = T2cpmg_persist_q2[i,2:]
plt.plot(x, y, f'C1^-')
i = 3
x = 2**np.arange(10)
y = T2cpmg_persist_q2[i,2:]
plt.plot(x, y, f'C1v-')

# plt.plot( (16,512), ( 25* 16**0.5, 25* 512**0.5), 'k--')
# plt.plot( (32,512), ( 10* 32**(3/4), 10* 512**(3/4)  ), 'k--')
plt.plot( (32,512), ( 10* 32**(2/3), 10* 512**(2/3)  ), 'k--')
# plt.plot( (16,512), ( 25* 16**(2/3), 25* 512**(2/3) ), 'k--')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('N_pi')
plt.ylabel('T2 (us)')
plt.grid()

# filename = 'FigureSupp10'  #'T2_all'
# plt.savefig(filename+'.pdf') 















