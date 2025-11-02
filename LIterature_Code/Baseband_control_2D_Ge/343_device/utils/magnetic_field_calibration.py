# -*- coding: utf-8 -*-
"""
Created on Thu May 30 10:36:13 2024

@author: Francesco
"""


# this script illustrates the calibration of the magnetic field using CPMG measurements
# by tracking the Ge-73 spectral noise peak, we estimate the effective magnetic field perceived at the device using the precession frequency of Ge-73 well know in the literature

import array_343.figures.hopping_spins_paper.utils.analysis_tools as tools

#%%

#Q1 datasets at set field of 0.17 T
# we derive the field at the sample by correcting for a calibration factor that matches the revival effects
# 0.17 T * calibration_fraction = 0.117 T, as shown below.

uuid_CPMG = dict()
uuid_CPMG['1'] = 1716540442970283691
uuid_CPMG['2'] = 1716541446593283691
uuid_CPMG['4'] = 1716542770826283691


calibration_fraction = 0.69

f_Ge73 = 1.487659 * 0.17 * calibration_fraction
n = np.arange(1, 4) 

revival_times = 2* n / f_Ge73 # us


fig, axs = plt.subplots(1, 1, figsize=tools.cm2inch(9, 9))

for CPMG in uuid_CPMG:

   
    uuid = uuid_CPMG[CPMG]    
    dat = load_by_uuid(uuid)
    
    x, y, z = dat.m1_2.y(), dat.m1_2.x(), dat.m1_2() #2
    x_label, y_label, z_label = dat.m1_2.y.label, dat.m1_2.x.label, dat.m1_2.label
    x_unit, y_unit, z_unit = dat.m1_2.y.unit, dat.m1_2.x.unit, dat.m1_2.unit
    z = np.average(z, axis = 0)   
    
    z_max = np.average(z[:5])
    z_min = np.average(z[-5:])
    z_norm = (z - z_min) / (z_max - z_min)
    
    # correcting the time axis for the Y gate duration
    x90_duration = dat.snapshot['measurement']['sequence']['settings']['q1']['x90']['t_pulse']
    Y_duration = 2*x90_duration

    axs.plot((x + Y_duration)*1e-3, z_norm, label = f'CPMG-{CPMG}')
    
axs.axvline(revival_times[0], lw = 0.5, linestyle = '--', label = f'exp. rev. times at B = {0.17 * calibration_fraction:.3f} T', c = 'grey')
for revival in revival_times[1:]: 
    axs.axvline(revival, lw = 0.5, linestyle = '--')
    
axs.legend()
axs.set_xlabel('wait time $(\mu s)$')
axs.set_ylabel('$P_{\\uparrow, norm.}$')
axs.set_xlim(xmin=0)    

plt.tight_layout()
plt.show()

