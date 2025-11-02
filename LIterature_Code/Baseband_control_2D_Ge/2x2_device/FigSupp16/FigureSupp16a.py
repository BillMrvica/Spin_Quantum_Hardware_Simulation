# this file is modified from Supp\CZ\CZ_leakage\Fig_2c_pulse_v0_simple.py
import numpy as np
import matplotlib.pyplot as plt

def voltage_to_J(voltage, j2v_params=dict(j_max=0.24e6, alpha=0.059 )  ):
    j_max, alpha = j2v_params['j_max'], j2v_params['alpha']
    '''
    voltage (float): fraction of max voltage. unit mV
    J = j_max when voltage == 1.0
    '''
    return j_max*np.exp(-alpha*(voltage))  # unit Hz


def J_to_voltage(J, j2v_params=dict(j_max=0.24e6, alpha=0.059 )  ):
    j_max, alpha = j2v_params['j_max'], j2v_params['alpha']
    J = np.asarray(J)  # unit Hz
    voltages = -np.log(J/j_max)/alpha

    return voltages  # unit mV

def cphase_function(duration, sample_rate, amplitude,
                    t_ramp=None, t_pulse=None, j_pulse_top=None, j_off=None,
                    j2v_params=dict(j_max=0.24e6, alpha=0.059 ) , shapeflag='hamming' ):
    sr = sample_rate/1e9 # convert to samples/ns
    # Use np.ceil, because np.arange count till, excluding last value.
    n_points = int(np.ceil(duration*sr))
    dc_start = int(np.ceil(t_ramp*sr))
    dc_stop  = int(np.ceil((t_ramp+t_pulse)*sr))
    offset = dc_stop - (t_ramp+t_pulse)*sr # offset in sample fraction due to rounding
    j_on_off = j_off/j_pulse_top
    # Tukey shaped pulse
    pulse_shape = np.zeros([n_points])
    if shapeflag == 'hann':
        #hann
        pulse_shape[0:dc_start] = (0.5+j_on_off/2) - (0.5-j_on_off/2)*np.cos(np.arange(dc_start)*np.pi/t_ramp)
        #hann
        pulse_shape[dc_start:dc_stop] = 1.0
        pulse_shape[dc_stop:n_points] = (0.5+j_on_off/2) - (0.5-j_on_off/2)*np.cos((np.arange(n_points-dc_stop)+offset)*np.pi/t_ramp+np.pi)

        v_fraction = J_to_voltage(pulse_shape*j_pulse_top, j2v_params  )
        
        #hann
        return v_fraction - v_fraction[0] , pulse_shape*j_pulse_top
    
    if shapeflag == 'hamming':
        #hamming
        pulse_shape[0:dc_start] = (0.54) - (0.46)*np.cos(np.arange(dc_start)*np.pi/t_ramp)
        pulse_shape[dc_start:dc_stop] = 1.0   
        #hamming
        pulse_shape[dc_stop:n_points] = (0.54) - (0.46)*np.cos((np.arange(n_points-dc_stop)+offset)*np.pi/t_ramp+np.pi)    
        #hamming
        v_fraction = J_to_voltage(pulse_shape*j_pulse_top, j2v_params  )
        return v_fraction -  J_to_voltage(j_off, j2v_params  ) , pulse_shape*j_pulse_top


# 2023-08-23 parameters for 2qRB


cz={
    't_padding': 16, 
    't_ramp':  23, #21.9, 
    't_pulse': 0,
    'gate_Joff': 14,  # gate_Joff and gate_Jon  should have same set of keys
    'gate_Jon': -76,
    't_ramp_in': 15, 
    't_ramp_out': 15, 
    'gate_voltages': -94, # relative values to _1100_s
}
_1100_s_vpB12 = 80
      


gate_Joff = cz['gate_Joff'] 
gate_Jon = cz['gate_Jon'] 

t_ramp = cz['t_ramp']
t_pulse = cz['t_pulse']

j2v_params=dict(j_max=0.24e6, alpha=0.059 )
t_gate = t_ramp*2 + t_pulse

j_off = voltage_to_J(gate_Joff, j2v_params )
j_on = voltage_to_J(gate_Jon, j2v_params )  

vgate, jgate = cphase_function(t_gate, 1e9, 1,
                    t_ramp=t_ramp, t_pulse=t_pulse, j_pulse_top=j_on, j_off=j_off,
                    j2v_params=dict(j_max=0.24e6, alpha=0.059 )  )


t_ramp_in = cz['t_ramp_in']
t_ramp_out = cz['t_ramp_out']



v_ramp_in = np.linspace(_1100_s_vpB12, _1100_s_vpB12+cz['gate_voltages'],  t_ramp_in)

v_ramp_out = np.linspace(_1100_s_vpB12+cz['gate_voltages'], _1100_s_vpB12,  t_ramp_out )

vgate = np.hstack(  [v_ramp_in, vgate+gate_Joff, v_ramp_out]   )





# plt.figure(figsize=(5,5))
# plt.plot(vgate, 'C0')
# ax = plt.gca().twinx()
# ax.plot(voltage_to_J(vgate), 'C1')

#%%
cz={
  't_padding': 16,
  't_ramp':  25, #24.84,
  't_pulse': 0,
  'gate_Joff': 40 , # gate_Joff and gate_Jon  should have same set of keys
  'gate_Jon':  -76,
  't_ramp_in': 15, 
  't_ramp_out': 15, 
  'gate_voltages': -40,
}
_1100_s_vpB12 = 80

gate_Joff = cz['gate_Joff'] 
gate_Jon = cz['gate_Jon'] 

t_ramp = cz['t_ramp']
t_pulse = cz['t_pulse']

j2v_params=dict(j_max=0.24e6, alpha=0.059 )
t_gate_hann = t_ramp*2 + t_pulse

j_off = voltage_to_J(gate_Joff, j2v_params )
j_on = voltage_to_J(gate_Jon, j2v_params )  

vgate_hann, jgate_hann = cphase_function(t_gate_hann, 1e9, 1,
                    t_ramp=t_ramp, t_pulse=t_pulse, j_pulse_top=j_on, j_off=j_off,
                    j2v_params=dict(j_max=0.24e6, alpha=0.059 ), shapeflag='hann'  )


t_ramp_in = cz['t_ramp_in']
t_ramp_out = cz['t_ramp_out']



v_ramp_in = np.linspace(_1100_s_vpB12, _1100_s_vpB12+cz['gate_voltages'],  t_ramp_in)[:-1]

v_ramp_out = np.linspace(_1100_s_vpB12+cz['gate_voltages'], _1100_s_vpB12,  t_ramp_out )

vgate_hann = np.hstack(  [v_ramp_in, vgate_hann+gate_Joff, v_ramp_out]   )



#%%
from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.2), Size.Fixed(0.7), Size.Fixed(1.2)]
v = [Size.Fixed(0.5), Size.Fixed(1.2), Size.Fixed(0.7), Size.Fixed(1.2)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)

ax11 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=1))
# ax13 = fig.add_axes(divider.get_position(),
#                     axes_locator=divider.new_locator(nx=1, ny=3))


# ax11.plot(vgate, 'C0')
# ax11.plot( np.arange(0, 15), vgate[0:15], 'k', linewidth=0.5   )
from scipy.signal.windows import general_hamming
import scipy.signal as signal
from scipy.fft import fft, fftshift, fftfreq
def mywin(alpha, t_gate, func):
    # alpha = 0.5
    if alpha is None:
        window = func( t_gate)
    else:
        window = func( t_gate, alpha)
    A = fft(window, 3000) / (len(window)/2.0)
    t = np.linspace(0, len(A)-1, len(A))
    freq = fftfreq(len(t), d=abs(t[1]-t[0]))
    return freq, abs(A)**2

def rec(t_gate):
    window = np.ones((t_gate,))
    A = fft(window, 3000) / (len(window)/2.0)
    t = np.linspace(0, len(A)-1, len(A))
    freq = fftfreq(len(t), d=abs(t[1]-t[0]))
    return freq, abs(A)**2

CQ1 = '#0D77BC'
CQ2 = '#DD5E27'

def lighter(color, percent):
    '''assumes color is rgb between (0, 0, 0) and (255, 255, 255)'''
    color = np.array(color)
    white = np.array([255, 255, 255])
    vector = white-color
    return (color + vector * percent)/255


# plt.figure(figsize=(12,6))
# plt.fill_between(  (0.044, 0.064), 1e-8, 2, color=lighter([13, 119, 188 ], 0.5)  )
# plt.fill_between(  (0.090, 0.112), 1e-8, 2, color=lighter([221, 94, 39 ], 0.5)  )
# plt.subplot(121)
ax = ax11

idx = 1+np.argmin(abs(mywin(0.54, t_gate, general_hamming)[0]*1e3 - 80))

t_gate = 46
t_gate_hann = 46


# plt.plot( *mywin(0.5, t_gate, general_hamming), 'k-', label=f'gate time {t_gate:.3g} ns ')
# plt.plot( *rec( t_gate_hann))   
ax.plot( rec( t_gate_hann)[0][0:idx]*1e3, 
          rec( t_gate_hann)[1][0:idx]/rec( t_gate_hann)[1][0]  , 
          'C0', linewidth=1 ) 
ax.plot( mywin(0.5, t_gate_hann, general_hamming)[0][0:idx]*1e3,
          mywin(0.5, t_gate, general_hamming)[1][0:idx]/mywin(0.5, t_gate, general_hamming)[1][0],
          'C1',  linewidth=1)
ax.plot( mywin(0.54, t_gate, general_hamming)[0][0:idx]*1e3, 
         mywin(0.54, t_gate, general_hamming)[1][0:idx]/mywin(0.54, t_gate, general_hamming)[1][0],
         'C2',  linewidth=1)


ax.plot( (43.7, 43.7), (1e-7, 2), 'r--', linewidth=1)

# plt.plot( (0.0426, 0.0426), (1e-5, 1), '-', c=CQ1)
# plt.plot( (0.0895, 0.0895), (1e-5, 1), '-', c=CQ2)


ax.set_yscale('log')
# ax.xlabel(' frequency (MHz)')
# plt.xlim(-0.12, 0.12)
ax.set_xlim(0, 80)
ax.set_ylim(1e-7, 2)
ax.set_xticks([0, 40, 80])
ax.set_yticks([1e-6, 1e-3, 1])

# filename = 'FigureSupp16a' # 'windowfunc_spectrum'
# plt.savefig(filename+'.pdf')
