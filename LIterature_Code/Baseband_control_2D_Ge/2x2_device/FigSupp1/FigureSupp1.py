# this file is modified from Supp/subns_ramp/subns_ramps.py
import os

import matplotlib.pyplot as pt
import numpy as np
import xarray as xr
from mpl_toolkits.axes_grid1 import Divider, Size


pt.rcParams.update({'axes.labelsize': 8,'xtick.labelsize': 8,'ytick.labelsize': 8,
                    'legend.fontsize': 8,'lines.markersize': 3,'font.size': 8,
                    'font.family': u'sans-serif','font.sans-serif': ['Arial'],'text.usetex': False})

format_string = 'pdf'
save = False
#%%

def get_awg_output(t, samples, analogue_shift, digital_filter_mode = 1):
    # pulse_response = xr.open_dataset(os.getcwd()+'/keysight_pulse_response_1.hdf5', engine='h5netcdf')
    fname = os.path.dirname(__file__) + f'/keysight_pulse_response_{digital_filter_mode}.hdf5'
    # print('File name: ', fname)
    pulse_response = xr.open_dataset(fname, engine='h5netcdf')
    t_response = pulse_response.coords['t'].data
    response = pulse_response['y'].data / 0.77
    sr = round(1/(t_response[1]-t_response[0]))

    t = np.linspace(t[0], t[-1]+1, len(t)*sr, endpoint=False)
    d = np.zeros(len(samples)*sr)
    d[::sr] = samples

    d = np.convolve(d, response)
    n_before = round(-t_response[0]*sr)
    n_after = round(t_response[-1]*sr)
    return t+analogue_shift, d[n_before: -n_after]


def render_ramp(t_start, t_stop, v_start, v_stop, samples):

    dv_dt = (v_stop - v_start) / (t_stop - t_start)
    i_start, dt_start = divmod(t_start, 1.0)
    i_stop, dt_stop = divmod(t_stop, 1.0)
    i_start = int(i_start)
    i_stop = int(i_stop)

    samples[i_start] += v_start*(1-dt_start)
    samples[i_stop] += v_stop*dt_stop - dv_dt*dt_stop
    samples[i_start+1:i_stop] += np.linspace(v_start + dv_dt*(1-dt_start),
                                             v_stop - dv_dt*dt_stop,
                                             i_stop-i_start-1, endpoint=False)

def quantize_amplitude(wave):
    min_vstep = 0.00038 # 0.38 mV
    wave = np.round(wave/min_vstep) * min_vstep
    return wave


def plot_samples(t, samples, color = 'k', linestyle = 'solid'):
    t2 = np.repeat(t, 2)[1:-1]
    samples2 = np.repeat(samples[:-1], 2)
    return pt.plot(t2, samples2, color)


t_max = 40
n_pulses = 2
t_ramp = 2.0
t_pulse = 10 #4.0
t_wait = 8.0
amplitude = 0.2
# amplitude = 0.02 / 0.126
# Shift for analogue signal to align with original ramp. Shift depends on ramp
# analogue_shift = -0.42 -0.16461214792859824 
analogue_shift =  -0.553 +0.00083289 #+0.03234041 
# 1: default digital filter; 0: no filter; 3: anti-ringing filter
digital_filter_mode = 1

# pt.close('all')
t = np.arange(t_max)

Pulses_ideal_sampled=pt.figure( figsize=(16,8))


HSPACE, VSPACE = 0.4, 0.4
h = [ s for _ in range(30) for s in [Size.Fixed(HSPACE), Size.Fixed(HSPACE)]  ]
v = [ s for _ in range(8) for s in [Size.Fixed(VSPACE), Size.Fixed(VSPACE)]  ]

divider = Divider(Pulses_ideal_sampled, (0, 0, 1, 1), h, v, aspect=False)
ax11 = Pulses_ideal_sampled.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=1, ny=1, nx1=6, ny1=6))

ax31_up = Pulses_ideal_sampled.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=8, ny=1, nx1=10, ny1=3))

ax31_dn = Pulses_ideal_sampled.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=8, ny=4, nx1=10, ny1=6))

ax51 = Pulses_ideal_sampled.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=12, ny=1, nx1=17, ny1=6))

ax71 = Pulses_ideal_sampled.add_axes(divider.get_position(),
                   axes_locator=divider.new_locator(nx=18, ny=1, nx1=23, ny1=6))



colors_1 = ['black', 'mediumblue']
colors_2 = ['gray', 'deepskyblue']
markers = ['^', 'v']
markerfmts = ['k^', 'bv']
for n, s in  enumerate([0,0.6]):
    wave = np.zeros(t_max)

    points = [(0.0, 0.0)]
    t0 = 1.0 + s
    # pt.title(f"first ramp starts at {t0:.2f} ns")

    for i in range(n_pulses):
        # add 1 ramped pulse
        points += [(t0, 0.0)]
        t0 += t_ramp
        points += [(t0, amplitude)]
        t0 += t_pulse
        points += [(t0, amplitude)]
        t0 += t_ramp
        points += [(t0, 0.0)]
        t0 += t_wait
    points += [(t_max-1, 0.0)]

    xy = np.array(points).T

    for i in range(len(points)-1):
        t_start, v_start = points[i]
        t_stop, v_stop = points[i+1]
        if t_stop > t_start:
            render_ramp(t_start, t_stop, v_start, v_stop, wave)

    wave = quantize_amplitude(wave)
    ta, out = get_awg_output(t, wave, analogue_shift, digital_filter_mode)


    pt.sca(ax11)
    pt.plot(xy[0], xy[1], colors_1[n], linestyle = (0, (1,1)),linewidth=1)
    l = plot_samples(t, wave, colors_1[n] )
    l[0].set_linewidth(1)
    pt.plot(ta, out, colors_2[n],linewidth=1)
    
    if n==0:
        pt.sca(ax31_dn)
        pt.plot(xy[0], xy[1], colors_1[n], linestyle = (0, (1,1)))
        # l = plot_samples(t, wave, colors_1[n])
        # l[0].set_linewidth(2)
        pt.plot(ta, out, colors_2[n])
        pt.tick_params(axis="y",direction="out", length=1.5, pad=3)
        pt.tick_params(axis="x",direction="out", length=1.5, pad=3)
        pt.xlim(1.95, 2.05)
        pt.ylim( 0.095, 0.105 )  
        pt.grid(True, which='both', axis='y', linestyle=':',linewidth=0.5)        
    else:
        pt.sca(ax31_up)
        pt.plot(xy[0], xy[1], colors_1[n], linestyle = (0, (1,1)))
        # l = plot_samples(t, wave, colors_1[n])
        # l[0].set_linewidth(2)
        pt.plot(ta, out, colors_2[n])
        pt.tick_params(axis="y",direction="out", length=1.5, pad=3)
        pt.tick_params(axis="x",direction="out", length=1.5, pad=3)
        pt.xlim(2.595, 2.605)
        pt.ylim(0.0995, 0.1005)        
        pt.grid(True, which='both', axis='y', linestyle=':',linewidth=0.5)        

    
    
    mx = np.logical_and(ta>xy[0][1], ta<xy[0][2]   )
    ta_mid = np.interp(amplitude/2, out[mx], ta[mx])
    print(ta_mid - (xy[0][1]+xy[0][2])/2)


for ax in [ax11, ]:
    pt.sca(ax)
    pt.grid(True, which='major', linestyle=':',linewidth=0.5)
    # pt.minorticks_on()
    # pt.locator_params(axis='x', nbins=6)
    pt.tick_params(axis="y",direction="out", length=1.5, pad=3)
    pt.tick_params(axis="x",direction="out", length=1.5, pad=3)
    pt.xticks([0,1,2,3,4,5])
    pt.xlim(0, 5)
    pt.ylim(-0.15 * amplitude, 1.15 * amplitude)
    pt.xlabel('Time (ns)')
    pt.ylabel('Voltage (V)')

# sarr = np.linspace(0, 1, 1001)[:-1]
sarr = np.linspace(0, 1, 1001)[:-1]
deviation = np.empty((len(sarr),2))
for n, s in  enumerate(sarr):
    wave = np.zeros(t_max)

    points = [(0.0, 0.0)]
    t0 = 1.0 + s
    # pt.title(f"first ramp starts at {t0:.2f} ns")

    for i in range(n_pulses):
        # add 1 ramped pulse
        points += [(t0, 0.0)]
        t0 += t_ramp
        points += [(t0, amplitude)]
        t0 += t_pulse
        points += [(t0, amplitude)]
        t0 += t_ramp
        points += [(t0, 0.0)]
        t0 += t_wait
    points += [(t_max-1, 0.0)]

    xy = np.array(points).T

    for i in range(len(points)-1):
        t_start, v_start = points[i]
        t_stop, v_stop = points[i+1]
        if t_stop > t_start:
            render_ramp(t_start, t_stop, v_start, v_stop, wave)

    wave = quantize_amplitude(wave)
    ta, out = get_awg_output(t, wave, analogue_shift, digital_filter_mode)

    mx = np.logical_and(ta>xy[0][1], ta<xy[0][2]   )
    ta_mid = np.interp(amplitude/2, out[mx], ta[mx])
    deviation[n,0] =  ta_mid - (xy[0][1]+xy[0][2])/2  

    mx = np.logical_and(ta>xy[0][3], ta<xy[0][4]   )
    ta_mid = np.interp(amplitude/2, out[mx][::-1], ta[mx][::-1])
    deviation[n,1] =  ta_mid - (xy[0][3]+xy[0][4])/2  

Nskip = 1
pt.sca(ax51)
# pt.plot(sarr, deviation, 'o')
pt.plot(sarr[::Nskip], deviation[::Nskip,0], '.', c= [ 255/256, 102/256, 0])

pt.ylim(-0.033, 0.033)
pt.xlim(0,1)
# pt.yticks([0,50,100])



pt.sca(ax71)
_, _, patches = ax71.hist(deviation[:,0], bins=20)
pt.setp(patches, facecolor='C0', edgecolor='k', lw=0.75)
pt.xticks(np.linspace(-0.03,0.03,7))
pt.xlim(-0.033, 0.033)
pt.yticks([0,50,100])
print(np.std(deviation,axis=0))
print( np.mean(deviation,axis=0) )

filename = 'FigureSupp1' # 'subns_ramps_v2'
pt.savefig(filename+'.pdf')
