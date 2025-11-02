# this file is modified from Fig2/Fig_2c_pulse.py
import numpy as np
import matplotlib.pyplot as plt

def voltage_to_J(voltage, j2v_params=dict(j_max=0.24e6, alpha=0.059 )  ):
    j_max, alpha = j2v_params['j_max'], j2v_params['alpha']
    '''
    voltage (float): fraction of max voltage. unit mV
    J = j_max when voltage == 1.0
    '''
    # print(j_max, alpha)
    return j_max*np.exp(-alpha*(voltage))  # unit Hz


def J_to_voltage(J, j2v_params=dict(j_max=0.24e6, alpha=0.059 )  ):
    j_max, alpha = j2v_params['j_max'], j2v_params['alpha']
    J = np.asarray(J)  # unit Hz
    voltages = -np.log(J/j_max)/alpha

    # voltages[voltages < 0] = 0
    return voltages  # unit mV

def cphase_function(duration, sample_rate, amplitude,
                    t_ramp=None, t_pulse=None, j_pulse_top=None, j_off=None,
                    j2v_params=dict(j_max=0.24e6, alpha=0.059 )  ):
    sr = sample_rate/1e9 # convert to samples/ns
    # Use np.ceil, because np.arange count till, excluding last value.
    n_points = int(np.ceil(duration*sr))
    dc_start = int(np.ceil(t_ramp*sr))
    dc_stop  = int(np.ceil((t_ramp+t_pulse)*sr))
    offset = dc_stop - (t_ramp+t_pulse)*sr # offset in sample fraction due to rounding
    j_on_off = j_off/j_pulse_top
    # Tukey shaped pulse
    pulse_shape = np.zeros([n_points])
    #hann
    # pulse_shape[0:dc_start] = (0.5+j_on_off/2) - (0.5-j_on_off/2)*np.cos(np.arange(dc_start)*np.pi/t_ramp)
    #hamming
    pulse_shape[0:dc_start] = (0.54) - (0.46)*np.cos(np.arange(dc_start)*np.pi/t_ramp)
    pulse_shape[dc_start:dc_stop] = 1.0
    #hann
    # pulse_shape[dc_stop:n_points] = (0.5+j_on_off/2) - (0.5-j_on_off/2)*np.cos((np.arange(n_points-dc_stop)+offset)*np.pi/t_ramp+np.pi)
    #hamming
    pulse_shape[dc_stop:n_points] = (0.54) - (0.46)*np.cos((np.arange(n_points-dc_stop)+offset)*np.pi/t_ramp+np.pi)
    # print(pulse_shape,j2v_params)
    v_fraction = J_to_voltage(pulse_shape*j_pulse_top, j2v_params  )
    # print(j_pulse_top)
    # print(v_fraction)
    
    
    #hann
    # return v_fraction - v_fraction[0]
    #hamming
    return v_fraction -  J_to_voltage(j_off, j2v_params  ) , pulse_shape*j_pulse_top


# 2023-08-23 parameters for 2qRB
cz={
    't_padding': 16, 
    't_ramp':  23, #21.9, 
    't_pulse': 0,
    'gate_Joff': -14,  # gate_Joff and gate_Jon  should have same set of keys
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


# plt.figure(figsize=(10,5))
# plt.subplot(121)
# plt.plot(vgate)

# plt.subplot(122)
# plt.plot(voltage_to_J(vgate))


plt.figure(figsize=(5,5))
plt.plot(vgate, 'C0')
ax = plt.gca().twinx()
ax.plot(voltage_to_J(vgate), 'C1')

#%%
from mpl_toolkits.axes_grid1 import Divider, Size

fig = plt.figure( figsize=(10,4))
h = [Size.Fixed(0.5), Size.Fixed(1.0), Size.Fixed(0.7), Size.Fixed(1.0)]
v = [Size.Fixed(0.5), Size.Fixed(1.0), Size.Fixed(0.7), Size.Fixed(1.0)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)

ax11 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=1))
# ax13 = fig.add_axes(divider.get_position(),
#                     axes_locator=divider.new_locator(nx=1, ny=3))


# ax11.plot(vgate, 'C0')
ax11.plot( np.arange(0, 15), vgate[0:15], 'k', linewidth=0.5   )
ax11.plot( np.arange(76-15, 76), vgate[76-15:76], 'k', linewidth=0.5)
ax11.plot( np.arange(14, 76-14), vgate[14:76-14], 'C0', linewidth=0.5)
ax11_twinx = ax11.twinx()
ax11_twinx.plot(voltage_to_J(vgate)/1e6, 'C2', linewidth=0.5)

# ax13.plot(voltage_to_J(vgate))
ax11.set_xticks([0, 30, 60])
# ax11.set_xlim([0, 73])
ax11.set_xlim([0, 75])

ax11.set_yticks([-80, 0, 80])
ax11.set_ylim([-85, 85])

ax11_twinx.set_ylim(0,)

ax = ax11
# ax.set_frame_on(True)
# ax.patch.set_visible(False)
# plt.setp(ax.spines["left"], visible=False)
# ax.spines["left"].set_visible(True)
# ax.spines["left"].set_linewidth(2)

# ax.spines["left"].set_color('C2')
# ax.spines["left"].set_edgecolor('C2')
# plt.show()

# filename = 'Figure2c_pulse' #'Fig_2c_pulse'
# plt.savefig(filename+'.pdf')








