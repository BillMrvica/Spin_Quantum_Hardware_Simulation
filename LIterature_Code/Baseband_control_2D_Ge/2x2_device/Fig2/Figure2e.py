# this file is modified from Fig2/2q_RB/2qRB_reorder_2023-08-22_21-41-44.py amd 2qRB_reorder_2023-08-22_21-41-44_CI.py
# for standard deviation of the data points, boostrap resampling, error bar of infidelity, see the data and scripts in TableSupp3

import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\Fig2','')
import sys
sys.path.insert(1, path_notebook)
import helper_functions
from helper_functions import hist_auto_su
#%% to plot histograms for readout check as function of sequence

import numpy as np
import pickle 
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size':12})


fdir = ''
all_sequences_data = pickle.load(  open(fdir + '2qRB_reorder_2023-08-22_21-41-44.pickle',  'rb')      )

sequence_length = all_sequences_data[0]['sequence_length']
hist_dict = all_sequences_data[0]['hist_dict']


autothreshold_rb_irb = {}
for i_interleaved in ['RB', 'IRB',  ]:
    num_sequences = all_sequences_data.size
    
    all_psb_q1q4_hists = []
    all_psb_q2q3_hists = []
    
    hist_range_q1q4 = (60,230)
    hist_range_q2q3 = (60,220)
    N_bins = 100
    
    psb14_vth = []
    psb23_vth = []
    p_autothreshold = np.zeros((num_sequences,4))
    for sequence_idx in range(1,num_sequences):
        L = len(all_sequences_data[sequence_idx][i_interleaved]['shots'])
        psb_q1q4 = np.hstack([all_sequences_data[sequence_idx][i_interleaved]['shots'][i][0] for i in range(L)])
        psb_q2q3 = np.hstack([all_sequences_data[sequence_idx][i_interleaved]['shots'][i][1] for i in range(L)])
        
        psb_q1q4_hists, psb14_binedges = np.histogram(psb_q1q4, bins=N_bins ,range=hist_range_q1q4)
        all_psb_q1q4_hists.append(psb_q1q4_hists)
    
        
        psb_q2q3_hists, psb23_binedges = np.histogram(psb_q2q3, bins=N_bins ,range=hist_range_q2q3)
        all_psb_q2q3_hists.append(psb_q2q3_hists)
    
    
    
        psb14_bincenters = (psb14_binedges[1:] + psb14_binedges[:-1])/2
        vth0 = hist_auto_su( psb14_bincenters, psb_q1q4_hists   )[3]
        psb14_vth += [ vth0 ]   
        
        psb23_bincenters = (psb23_binedges[1:] + psb23_binedges[:-1])/2
        vth1 = hist_auto_su( psb23_bincenters, psb_q2q3_hists   )[3]
        psb23_vth += [ vth1 ]    
        
        d_ss_trig = [psb_q1q4, psb_q2q3   ]
        inv0 = -1 if hist_dict['invert_result'][0] else 1
        inv1 = -1 if hist_dict['invert_result'][1] else 1    
        blist11 = np.logical_and( inv0*(d_ss_trig[0]-vth0)>0, inv1*(d_ss_trig[1]-vth1)>0  )
        blist10 = np.logical_and( inv0*(d_ss_trig[0]-vth0)>0, inv1*(d_ss_trig[1]-vth1)<=0  )
        blist01 = np.logical_and( inv0*(d_ss_trig[0]-vth0)<=0, inv1*(d_ss_trig[1]-vth1)>0  )
        blist00 = np.logical_and( inv0*(d_ss_trig[0]-vth0)<=0, inv1*(d_ss_trig[1]-vth1)<=0  )
        
        ss11 = np.count_nonzero(blist11)
        ss10 = np.count_nonzero(blist10)
        ss01 = np.count_nonzero(blist01)
        ss00 = np.count_nonzero(blist00)
        p_autothreshold[sequence_idx,:] = np.array([ss00, ss01, ss10, ss11 ] )
        
        if sequence_idx == 1:
            print(np.sum(psb_q1q4>144))
    
    all_psb_q1q4_hists = np.array(all_psb_q1q4_hists)
    all_psb_q2q3_hists = np.array(all_psb_q2q3_hists)
    
    hist_range14_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins), range(1,num_sequences)[::-1])
    hist_range23_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins), range(1,num_sequences)[::-1])
    
    VMAX=None
    plt.figure(figsize=(10,5))
    plt.subplot(121)
    plt.pcolormesh(hist_range14_mesh,sequence_idx_mesh,all_psb_q1q4_hists,vmax=VMAX)
    plt.plot( psb14_vth, sequence_idx_mesh[:,0]  )
    plt.plot( (151,)*2, (sequence_idx_mesh[0,0], sequence_idx_mesh[-1,0])  )
    
    plt.subplot(122)
    plt.pcolormesh(hist_range23_mesh,sequence_idx_mesh,all_psb_q2q3_hists,vmax=VMAX)
    plt.plot( psb23_vth, sequence_idx_mesh[:,0]  )
    plt.plot( (119,)*2, (sequence_idx_mesh[0,0], sequence_idx_mesh[-1,0])  )
    plt.show()

        
    autothreshold_rb_irb[i_interleaved] = np.array([p_autothreshold[ii,:]/np.sum(p_autothreshold[ii,:]) for ii in range(1, num_sequences)])

###


p00_rb = []
p00_irb = []
for i_gatestatement in range(1,len(all_sequences_data)):
    temp = all_sequences_data[i_gatestatement]['RB']['counts_all'] / np.sum(all_sequences_data[i_gatestatement]['RB']['counts_all'])
    p00_rb += [temp]
    
    temp = all_sequences_data[i_gatestatement]['IRB']['counts_all'] / np.sum(all_sequences_data[i_gatestatement]['IRB']['counts_all'])
    p00_irb += [temp]    
    
    
p00_rb = np.array(p00_rb)
p00_irb = np.array(p00_irb)


#%%

from projects.notebook_tools.notebook_tools import  fit_data
# def expdecay_alpha(x, A, x0, alpha, y0):
#     return A*np.exp(-(x/x0)**alpha) + y0

def expdecay_alpha(x, A, decayrate, alpha, y0):
    return A*(decayrate**(x**alpha)) + y0 

def expdecay(x, A, decayrate, y0):
    return A*(decayrate**x) + y0

xdata = np.array(sequence_length)


zdata_rb = p00_rb[:,0]   
# zdata_rb = autothreshold_rb_irb['RB'][:,0]
p0 = [0.65, 0.92,  0.25]
_, p1_rb = fit_data(xdata, zdata_rb, func=expdecay,p0=p0, plot=False, return_cov=True)



zdata_irb = p00_irb[:,0]   
# zdata_irb = autothreshold_rb_irb['IRB'][:,0]
p0 = [0.65, 0.85,  0.25]
_, p1_irb = fit_data(xdata, zdata_irb, func=expdecay,p0=p0, plot=False, return_cov=True)



decayrate_ref = p1_rb[1]
decayrate_interleaved = p1_irb[1]


d = 4 # 2 qubits
avg_err_ref = (1-decayrate_ref)*(d-1)/d
avg_interleaved_gate_infidelity = (1-decayrate_interleaved/decayrate_ref)*(d-1)/d

p1_rb_99 = np.array(p1_rb)
p1_rb_99[1] = (1-0.01*d/(d-1))
p1_rb_98 = np.array(p1_rb)
p1_rb_98[1] = (1-0.02*d/(d-1))

p1_irb_995 = np.array(p1_irb)
p1_irb_995[1] = (1-0.005*d/(d-1))*decayrate_ref
p1_irb_99 = np.array(p1_irb)
p1_irb_99[1] = (1-0.01*d/(d-1))*decayrate_ref
p1_irb_98 = np.array(p1_irb)
p1_irb_98[1] = (1-0.02*d/(d-1))*decayrate_ref


print('#######')
print('If enforce a simple exponential decay, alpha=1:')
print(f'reference sequence decay rate  = {p1_rb[1]:.4g}')
print(f'interleaved sequence decay rate = {p1_irb[1]:.4g}')
print(f'2Q Clifford infidelity = {avg_err_ref:.4g}')
print(f'interleaved CZ gate infidelity = {avg_interleaved_gate_infidelity:.4g}')
print(f'Infinite reference sequence will decay to {p1_rb[2]:.4g}')
print(f'Infinite interleaved sequence will decay to {p1_irb[2]:.4g}')

p0 = [0.65, 0.92, 1,  0.25] 
cov_rb_alpha, p1_rb_alpha = fit_data(xdata, zdata_rb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
p0 = [0.65, 0.92, 1,  0.25] 
cov_irb_alpha, p1_irb_alpha = fit_data(xdata, zdata_irb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
decayrate_ref_alpha = p1_rb_alpha[1] 
decayrate_interleaved_alpha = p1_irb_alpha[1] 
avg_err_ref_alpha = (1-decayrate_ref_alpha)*(d-1)/d
avg_interleaved_gate_infidelity_alpha = (1-decayrate_interleaved_alpha/decayrate_ref_alpha)*(d-1)/d

print('#######')
print('If not enforce alpha=1 but a free fitting parameter:')
print(f'reference sequence alpha  = {p1_rb_alpha[2]:.4g}')
print(f'interleaved sequence alpha = {p1_irb_alpha[2]:.4g}')
print(f'reference sequence decay rate  = {p1_rb_alpha[1]:.4g}')
print(f'interleaved sequence decay rate = {p1_irb_alpha[1]:.4g}')
print(f'2Q Clifford infidelity = {avg_err_ref_alpha:.4g}')
print(f'interleaved CZ gate infidelity = {avg_interleaved_gate_infidelity_alpha:.4g}')
print(f'Infinite reference sequence will decay to {p1_rb_alpha[3]:.4g}')
print(f'Infinite interleaved sequence will decay to {p1_irb_alpha[3]:.4g}')

# printing : 
# #######
# If enforce a simple exponential decay, alpha=1:
# reference sequence decay rate  = 0.9813
# interleaved sequence decay rate = 0.9725
# 2Q Clifford infidelity = 0.01403
# interleaved CZ gate infidelity = 0.006687
# Infinite reference sequence will decay to 0.2844
# Infinite interleaved sequence will decay to 0.283
# #######
# If not enforce alpha=1 but a free fitting parameter:
# reference sequence alpha  = 1.05
# interleaved sequence alpha = 0.9457
# reference sequence decay rate  = 0.9844
# interleaved sequence decay rate = 0.9667
# 2Q Clifford infidelity = 0.01167
# interleaved CZ gate infidelity = 0.01355
# Infinite reference sequence will decay to 0.2916
# Infinite interleaved sequence will decay to 0.278

#%%


xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))

idx=0
plt.figure()
plt.scatter(xdata, zdata_rb , c='k', s=25, label = 'RB')
plt.plot(xx, expdecay(xx, *p1_rb), 'k', label='RB ' )
# plt.plot(xx, expdecay_alpha(xx, *p1_rb_alpha) , label=f'RB alpha={p1_rb_alpha[2]:.3g}')
# plt.plot(xx, expdecay(xx, *p1_rb_99), 'k:'  )
# plt.plot(xx, expdecay(xx, *p1_rb_98), 'k:'  )



plt.scatter(xdata, zdata_irb , c='C2', s=25, label = 'IRB')
# plt.plot(xx, expdecay_alpha(xx, *p1_irb_alpha) , label=f'IRB alpha={p1_irb_alpha[2]:.3g}')
plt.plot(xx, expdecay(xx, *p1_irb), 'C2', label='IRB '  )
plt.plot(xx, expdecay(xx, *p1_irb_995), 'C2:'  )
plt.plot(xx, expdecay(xx, *p1_irb_99), 'C2:'  )
# plt.plot(xx, expdecay(xx, *p1_irb_98), 'C2:'  )



# plt.xscale('log')
plt.legend()




#%%
from mpl_toolkits.axes_grid1 import Divider, Size
fig = plt.figure( figsize=(6,4))
h = [Size.Fixed(0.5), Size.Fixed(3), Size.Fixed(0.7), Size.Fixed(1.0)]
v = [Size.Fixed(0.5), Size.Fixed(1.0), Size.Fixed(0.7), Size.Fixed(1.0)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)

ax11 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=1))

ax11.scatter(xdata, zdata_rb , c='k', s=5, label = 'RB')
ax11.plot(xx, expdecay(xx, *p1_rb), 'k', linewidth=0.5, label='RB ' )




ax11.scatter(xdata, zdata_irb , c='C2', s=5, label = 'IRB')

ax11.plot(xx, expdecay(xx, *p1_irb), 'C2', linewidth=0.5, label='IRB '  )
# ax11.plot(xx, expdecay(xx, *p1_irb_995), 'C2:'  )
# ax11.plot(xx, expdecay(xx, *p1_irb_99), 'C2:'  )

ax11.set_xlim(-5, 205 )
# ax11.set_xticks([0, 50, 100, 150, 200]) 
ax11.set_xticks( np.arange(0, 210, 25)) 
# ax11.set_ylim(0,1)
ax11.set_yticks([0.3, 0.5, 0.7, 0.9])   

# filename = 'Figure2e' # 'Fig_2qRB'
# plt.savefig(filename+'.pdf')


#%%  counting avg. premitive gate numbers in a Clifford

for idx in [1,2,20]:

    LL = sequence_length[idx-1] +1 
    L = len(all_sequences_data[idx]['RB']['gatelist'])
    countg = {'Gzpi2:0': 0, 'Gzpi2:1': 0,
                'Gxpi2:0': 0, 'Gxpi2:1': 0,
                'Gcphase:0:1': 0,
                'Gi:0:1': 0}
    for i in range(L):
        glist = all_sequences_data[idx]['RB']['gatelist'][i]
        for j in range(len(glist)):
            g = glist[j]
            countg[g] += 1
        
        
    countz = countg['Gzpi2:0'] + countg['Gzpi2:1']
    countx = countg['Gxpi2:0'] + countg['Gxpi2:1']
    countcz = countg  ['Gcphase:0:1']  
    counti = countg    ['Gi:0:1'] 
      
    print(f'number of randomization={L}')
    print(f'sequence length (including recovery Clifford) = {LL}')
    print(f'avg. number of Z_pi/2 per Clifford = {countz/(L*LL):.4g}' )
    print(f'avg. number of X_pi/2 per Clifford = {countx/(L*LL):.4g}' )    
    print(f'avg. number of CZ per Clifford = {countcz/(L*LL):.4g}' )     
    print(f'avg. number of idle per Clifford = {counti/(L*LL):.4g}' )     
    
#  1.627 CZ + 3.21 X90 + 5.36 Z90     
    




