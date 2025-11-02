# this file is modified from Supp\2q_RB_all\2qRB_reorder_2023-08-23_09-20-24_CI_supp.py
# confidence interval calculation of 2Q RB is at the bottom of the script
import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\TableSupp3','')
import sys
sys.path.insert(1, path_notebook)
import helper_functions
from helper_functions import hist_auto_su

from projects.notebook_tools.notebook_tools import  fit_data
#%% to plot histograms for readout check as function of sequence

import numpy as np
import pickle 
from matplotlib import pyplot as plt




plt.rcParams.update({'font.size':12})



fdir = ''
all_sequences_data = pickle.load(  open(fdir + '2023-08-23_09-20-24_2qRB_reorder.pickle',  'rb')      )

sequence_length = all_sequences_data[0]['sequence_length']
hist_dict = all_sequences_data[0]['hist_dict']
#all_sequences_data = gatestatement

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
        # data = all_sequences_data[sequence_idx][i_interleaved]['shots'][0]
        
        # psb_q1q4 = data[0]
        # psb_q2q3 = data[1]
        L = len(all_sequences_data[sequence_idx][i_interleaved]['shots'])
        psb_q1q4 = np.hstack([all_sequences_data[sequence_idx][i_interleaved]['shots'][i][0] for i in range(L)])
        psb_q2q3 = np.hstack([all_sequences_data[sequence_idx][i_interleaved]['shots'][i][1] for i in range(L)])
        
        psb_q1q4_hists, psb14_binedges = np.histogram(psb_q1q4, bins=N_bins ,range=hist_range_q1q4)
        all_psb_q1q4_hists.append(psb_q1q4_hists)
    
        
        psb_q2q3_hists, psb23_binedges = np.histogram(psb_q2q3, bins=N_bins ,range=hist_range_q2q3)
        all_psb_q2q3_hists.append(psb_q2q3_hists)
    
    
    
        psb14_bincenters = (psb14_binedges[1:] + psb14_binedges[:-1])/2
        # vth0 = hist_auto_su( psb14_bincenters, psb_q1q4_hists   )[3]
        vth0 = 151
        psb14_vth += [ vth0 ]   
        
        psb23_bincenters = (psb23_binedges[1:] + psb23_binedges[:-1])/2
        # vth1 = hist_auto_su( psb23_bincenters, psb_q2q3_hists   )[3]
        vth1 = 119
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
    # print(gatestatement[i_gatestatement]['counts_all'])
    # print(np.sum(gatestatement[i_gatestatement]['counts_all']))
    temp = all_sequences_data[i_gatestatement]['RB']['counts_all'] / np.sum(all_sequences_data[i_gatestatement]['RB']['counts_all'])
    p00_rb += [temp]
    
    temp = all_sequences_data[i_gatestatement]['IRB']['counts_all'] / np.sum(all_sequences_data[i_gatestatement]['IRB']['counts_all'])
    p00_irb += [temp]    
    
    
p00_rb = np.array(p00_rb)
p00_irb = np.array(p00_irb)

###
# idx=0
# plt.figure()
# plt.plot(sequence_length, p00_rb[:,idx], 'ob--')
# plt.plot(sequence_length, p00_irb[:,idx], '^b--')

# plt.plot(sequence_length, autothreshold_rb_irb['RB'][:,idx]  , 'ob-')
# plt.plot(sequence_length, autothreshold_rb_irb['IRB'][:,idx]  , 'ob-')

###
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
#######
# If enforce a simple exponential decay, alpha=1:
# reference sequence decay rate  = 0.9802
# interleaved sequence decay rate = 0.969
# 2Q Clifford infidelity = 0.01483
# interleaved CZ gate infidelity = 0.008598
# Infinite reference sequence will decay to 0.2868
# Infinite interleaved sequence will decay to 0.2835
# #######
# If not enforce alpha=1 but a free fitting parameter:
# reference sequence alpha  = 0.9881
# interleaved sequence alpha = 0.9536
# reference sequence decay rate  = 0.9794
# interleaved sequence decay rate = 0.9635
# 2Q Clifford infidelity = 0.01549
# interleaved CZ gate infidelity = 0.01212
# Infinite reference sequence will decay to 0.285
# Infinite interleaved sequence will decay to 0.2797

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




#%% collect and plot the data of 2Q RB data point distribution





import numpy as np
import pickle 
from matplotlib import pyplot as plt



plt.rcParams.update({'font.size':12})


fdir = ''
all_sequences_data = pickle.load(  open(fdir + '2023-08-23_09-20-24_2qRB_reorder.pickle',  'rb')      )

sequence_length = all_sequences_data[0]['sequence_length']
hist_dict = all_sequences_data[0]['hist_dict']
#all_sequences_data = gatestatement

autothreshold_rb_irb = {}
autothreshold_rb_irb_dist = {} 
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
    p_autothreshold_dist = np.zeros((num_sequences, ), dtype=object)
    for sequence_idx in range(1,num_sequences):
        # data = all_sequences_data[sequence_idx][i_interleaved]['shots'][0]
        
        # psb_q1q4 = data[0]
        # psb_q2q3 = data[1]
        L = len(all_sequences_data[sequence_idx][i_interleaved]['shots'])
        psb_q1q4 = np.hstack([all_sequences_data[sequence_idx][i_interleaved]['shots'][i][0] for i in range(L)])
        psb_q2q3 = np.hstack([all_sequences_data[sequence_idx][i_interleaved]['shots'][i][1] for i in range(L)])
        
        psb_q1q4_hists, psb14_binedges = np.histogram(psb_q1q4, bins=N_bins ,range=hist_range_q1q4)
        all_psb_q1q4_hists.append(psb_q1q4_hists)
    
        
        psb_q2q3_hists, psb23_binedges = np.histogram(psb_q2q3, bins=N_bins ,range=hist_range_q2q3)
        all_psb_q2q3_hists.append(psb_q2q3_hists)
    
    
    
        psb14_bincenters = (psb14_binedges[1:] + psb14_binedges[:-1])/2
        # vth0 = hist_auto_su( psb14_bincenters, psb_q1q4_hists   )[3]
        vth0 = 151
        psb14_vth += [ vth0 ]   
        
        psb23_bincenters = (psb23_binedges[1:] + psb23_binedges[:-1])/2
        # vth1 = hist_auto_su( psb23_bincenters, psb_q2q3_hists   )[3]
        vth1 = 119
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
        
        # if sequence_idx == 1:
        #     print(np.sum(psb_q1q4>144))
            ###
        temp_list = []
        for ii in range(L):
            temp_psb14= all_sequences_data[sequence_idx][i_interleaved]['shots'][ii][0]
            temp_psb23= all_sequences_data[sequence_idx][i_interleaved]['shots'][ii][1]
            temp_list += [ np.sum( np.logical_and( temp_psb14 < hist_dict['thresholds'][0], temp_psb23 > hist_dict['thresholds'][1] ) )/len(temp_psb14) ] 
        p_autothreshold_dist[sequence_idx ] = list(temp_list)
            
            
        ####
                
    all_psb_q1q4_hists = np.array(all_psb_q1q4_hists)
    all_psb_q2q3_hists = np.array(all_psb_q2q3_hists)
    
    hist_range14_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins), range(1,num_sequences)[::-1])
    hist_range23_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins), range(1,num_sequences)[::-1])
    
    VMAX=None
    # plt.figure(figsize=(10,5))
    # plt.subplot(121)
    # plt.pcolormesh(hist_range14_mesh,sequence_idx_mesh,all_psb_q1q4_hists,vmax=VMAX)
    # plt.plot( psb14_vth, sequence_idx_mesh[:,0]  )
    # plt.plot( (151,)*2, (sequence_idx_mesh[0,0], sequence_idx_mesh[-1,0])  )
    
    # plt.subplot(122)
    # plt.pcolormesh(hist_range23_mesh,sequence_idx_mesh,all_psb_q2q3_hists,vmax=VMAX)
    # plt.plot( psb23_vth, sequence_idx_mesh[:,0]  )
    # plt.plot( (119,)*2, (sequence_idx_mesh[0,0], sequence_idx_mesh[-1,0])  )
    # plt.show()
    

        
    autothreshold_rb_irb[i_interleaved] = np.array([p_autothreshold[ii,:]/np.sum(p_autothreshold[ii,:]) for ii in range(1, num_sequences)])
    autothreshold_rb_irb_dist[i_interleaved]  = p_autothreshold_dist[1:]
###


p00_rb = []
p00_irb = []
for i_gatestatement in range(1,len(all_sequences_data)):
    # print(gatestatement[i_gatestatement]['counts_all'])
    # print(np.sum(gatestatement[i_gatestatement]['counts_all']))
    temp = all_sequences_data[i_gatestatement]['RB']['counts_all'] / np.sum(all_sequences_data[i_gatestatement]['RB']['counts_all'])
    p00_rb += [temp]
    
    temp = all_sequences_data[i_gatestatement]['IRB']['counts_all'] / np.sum(all_sequences_data[i_gatestatement]['IRB']['counts_all'])
    p00_irb += [temp]    
    
    
p00_rb = np.array(p00_rb)
p00_irb = np.array(p00_irb)









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







p0 = [0.65, 0.92, 1,  0.25] 
cov_rb_alpha, p1_rb_alpha = fit_data(xdata, zdata_rb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
p0 = [0.65, 0.92, 1,  0.25] 
cov_irb_alpha, p1_irb_alpha = fit_data(xdata, zdata_irb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
decayrate_ref_alpha = p1_rb_alpha[1] 
decayrate_interleaved_alpha = p1_irb_alpha[1] 
avg_err_ref_alpha = (1-decayrate_ref_alpha)*(d-1)/d
avg_interleaved_gate_infidelity_alpha = (1-decayrate_interleaved_alpha/decayrate_ref_alpha)*(d-1)/d




# plot the data of 2Q RB data point distribution
xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))


for err_factor, titlestr in [ (1,'std of 128 randomizations'), (np.sqrt(128),'std of 128 randomizations/sqrt(128)') ]:
    idx=0
    plt.figure(figsize=(14,6))
    plt.suptitle(titlestr)
    plt.subplot(121)
    # plt.scatter(xdata, zdata_rb , c='k', s=25, label = 'RB')
    # for i_sequence, temp_list in enumerate(autothreshold_rb_irb_dist['RB']):
    #     plt.scatter((xdata[i_sequence],)*len(temp_list), temp_list  ,  s=25,  c='k' )
    
    plt.plot(xx, expdecay(xx, *p1_rb), 'k', label='RB ' )
    # plt.plot(xx, expdecay_alpha(xx, *p1_rb_alpha) , label=f'RB alpha={p1_rb_alpha[2]:.3g}')
    for i_sequence, temp_list in enumerate(autothreshold_rb_irb_dist['RB']):
        plt.errorbar(xdata[i_sequence]  , np.mean(np.array(temp_list)),
                     yerr = np.std(np.array(temp_list))/err_factor, c='k',
                     fmt ='o', markersize=3, linestyle='', linewidth=1)
    
    
    plt.subplot(122)
    # plt.scatter(xdata, zdata_irb , c='C2', s=25, label = 'IRB')
    # for i_sequence, temp_list in enumerate(autothreshold_rb_irb_dist['IRB']):
        # plt.scatter((xdata[i_sequence],)*len(temp_list), temp_list  ,  s=25,  c='C2' )
    plt.plot(xx, expdecay(xx, *p1_irb), 'C2', label='IRB '  )
    
    
    for i_sequence, temp_list in enumerate(autothreshold_rb_irb_dist['IRB']):
        plt.errorbar(xdata[i_sequence]  , np.mean(np.array(temp_list)),
                     yerr = np.std(np.array(temp_list))/err_factor, c='C2',
                     fmt ='o', markersize=3, linestyle='', linewidth=1)
    
    # plt.xscale('log')
    plt.legend()




idx=0
plt.figure(figsize=(14,6))
plt.title('scatter points randomized sequences')
plt.subplot(121)
for i_sequence, temp_list in enumerate(autothreshold_rb_irb_dist['RB']):
    plt.scatter((xdata[i_sequence],)*len(temp_list), temp_list  ,  s=5,  c='k' )
plt.plot(xx, expdecay(xx, *p1_rb), 'k', label='RB ' )



plt.subplot(122)
for i_sequence, temp_list in enumerate(autothreshold_rb_irb_dist['IRB']):
    plt.scatter((xdata[i_sequence],)*len(temp_list), temp_list  ,  s=5,  c='C2' )
plt.plot(xx, expdecay(xx, *p1_irb), 'C2', label='IRB '  )


plt.legend()



#%% load plot and print the 2Q RB bootstrap resampling  results
results = pickle.load( open('2023-08-23_09-20-24_2qRB_reorder_bootstrap_resample.pickle', 'rb') )['results']



plt.figure()
plt.title('histogram of avg. infidelity')
# plt.plot([results[i]['p1_rb'][1] for i in range(len(results))]   )
# plt.plot([results[i]['p1_irb'][1] for i in range(len(results))]   )
# plt.plot([results[i]['p1_irb_alpha'][3] for i in range(len(results))]   )


avg_err_ref_list = [ results[i]['avg_err_ref'] for i in range(len(results))]
u, v = np.histogram(avg_err_ref_list  , bins=50,  ) 
v = (v[1:]+v[:-1])/2 
plt.plot(v, u)


avg_interleaved_gate_infidelity_list = [ results[i]['avg_interleaved_gate_infidelity'] for i in range(len(results))]
u, v = np.histogram(avg_interleaved_gate_infidelity_list  , bins=50,  ) 
v = (v[1:]+v[:-1])/2 
plt.plot(v, u)
plt.xlabel('avg. gate infidelity')
plt.ylabel('counts')


avg_err_ref_sort = np.sort(np.array(avg_err_ref_list))
avg_interleaved_gate_infidelity_list = np.sort(np.array(avg_interleaved_gate_infidelity_list))


print('original data: 2Q Clifford avg infidelity = 0.01483')
print(f'bootstrap resampling average infidelity = {np.average(avg_err_ref_sort):.5g}', )
print('5% = ', avg_err_ref_sort[500] )
print('95% = ', avg_err_ref_sort[9500] )
print(f'-> 2Q Clifford avg infidelity = 0.01483 +- {abs(avg_err_ref_sort[9500]-avg_err_ref_sort[500])/2:.4g}, with 95% confidence interval')
# print(f'-> 2Q Clifford avg infidelity = 0.01403 +- {abs(avg_err_ref_sort[950]-avg_err_ref_sort[50])/2:.4g}, with 95% confidence interval')


print('original data: interleaved gate avg infidelity = 0.008598')
print(f'bootstrap resampling average infidelity = {np.average(avg_interleaved_gate_infidelity_list):.5g}', )
print('5% = ', avg_interleaved_gate_infidelity_list[500] )
print('95% = ', avg_interleaved_gate_infidelity_list[9500] )
print(f'-> interleaved gate avg infidelity = 0.008598 +- {abs(avg_interleaved_gate_infidelity_list[9500]-avg_interleaved_gate_infidelity_list[500])/2:.4g}, with 95% confidence interval')

# result of N_resample = 10000: 
# original data: 2Q Clifford avg infidelity = 0.01483
# bootstrap resampling average infidelity = 0.014822
# 5% =  0.014189333017929961
# 95% =  0.015467127585647805
# -> 2Q Clifford avg infidelity = 0.01483 +- 0.0006389, with 95% confidence interval
# original data: interleaved gate avg infidelity = 0.008598
# bootstrap resampling average infidelity = 0.0085726
# 5% =  0.007480128258760488
# 95% =  0.00966087401062296
# -> interleaved gate avg infidelity = 0.008598 +- 0.00109, with 95% confidence interval




plt.figure(figsize=(14,6))
plt.suptitle('histogram of fitting results (fit to super exponent)')
# plt.plot([results[i]['p1_rb'][1] for i in range(len(results))]   )
# plt.plot([results[i]['p1_irb'][1] for i in range(len(results))]   )
# plt.plot([results[i]['p1_irb_alpha'][3] for i in range(len(results))]   )

# avg_err_ref_list = [ results[i]['avg_err_ref'] for i in range(len(results))]
plt.subplot(121)
avg_err_ref_list = [ results[i]['avg_err_ref_alpha'] for i in range(len(results))]
u, v = np.histogram(avg_err_ref_list  , bins=50,  ) 
v = (v[1:]+v[:-1])/2 
plt.plot(v, u)

avg_err_ref_list = [ results[i]['avg_interleaved_gate_infidelity_alpha'] for i in range(len(results))]
u, v = np.histogram(avg_err_ref_list  , bins=50,  ) 
v = (v[1:]+v[:-1])/2 
plt.plot(v, u)
plt.xlabel('avg. gate infidelity')
plt.ylabel('counts')

plt.subplot(122)
idx_p = 2
alpha_rb_list = [ results[i]['p1_rb_alpha'][idx_p] for i in range(len(results))]
u, v = np.histogram(alpha_rb_list  , bins=50,  ) 
v = (v[1:]+v[:-1])/2 
plt.plot(v, u)
plt.xlabel('exponent  ' + r'$\alpha$')
plt.ylabel('counts')

tempvar95 = np.sort(np.array(alpha_rb_list))[int(0.95*len(alpha_rb_list))]
tempvar05 = np.sort(np.array(alpha_rb_list))[int(0.05*len(alpha_rb_list))] 
print( 'ref sequence alpha = 0.9881 +- ' + f'{(tempvar95-tempvar05)/2:.5g}'  )


alpha_irb_list = [ results[i]['p1_irb_alpha'][idx_p] for i in range(len(results))]
u, v = np.histogram(alpha_irb_list  , bins=50,  ) 
v = (v[1:]+v[:-1])/2 
plt.plot(v, u)
tempvar95 = np.sort(np.array(alpha_irb_list))[int(0.95*len(alpha_irb_list))]
tempvar05 = np.sort(np.array(alpha_irb_list))[int(0.05*len(alpha_irb_list))] 
print( 'interleaved sequence alpha = 0.9536 +- ' + f'{(tempvar95-tempvar05)/2:.5g}'  )



val_ref = [ results[i]['avg_err_ref_alpha'] for i in range(len(results))]
val_interleaved =  [ results[i]['avg_interleaved_gate_infidelity_alpha'] for i in range(len(results))]
val_ref_sort = np.sort(np.array(val_ref))
val_interleaved_sort = np.sort(np.array(val_interleaved))


print( 'reference sequence alpha  = 0.9881'  )
print('original data with free alpha: 2Q Clifford avg infidelity = 0.01549')
print(f'bootstrap resampling average infidelity = {np.average(val_ref_sort):.5g}', )
print('5% = ', val_ref_sort[int(0.05*len(avg_err_ref_list))] )
print('95% = ', val_ref_sort[int(0.95*len(avg_err_ref_list))] )
print(f'-> 2Q Clifford avg infidelity = 0.01549 +- {abs(val_ref_sort[int(0.95*len(avg_err_ref_list))]-val_ref_sort[int(0.05*len(avg_err_ref_list))])/2:.4g}, with 95% confidence interval')

print( 'interleaved sequence alpha = 0.9536' )
print('original data with free alpha: interleaved gate avg infidelity = 0.01212')
print(f'bootstrap resampling average infidelity = {np.average(val_interleaved_sort):.5g}', )
print('5% = ', val_interleaved_sort[int(0.05*len(avg_err_ref_list))] )
print('95% = ', val_interleaved_sort[int(0.95*len(avg_err_ref_list))] )
print(f'-> interleaved gate avg infidelity = 0.01212 +- {abs(val_interleaved_sort[int(0.95*len(avg_err_ref_list))]-val_interleaved_sort[int(0.05*len(avg_err_ref_list))])/2:.4g}, with 95% confidence interval')




#######
# ref sequence alpha = 0.9881 +- 0.056674
# interleaved sequence alpha = 0.9536 +- 0.047006
# reference sequence alpha  = 0.9881
# original data with free alpha: 2Q Clifford avg infidelity = 0.01549
# bootstrap resampling average infidelity = 0.015669
# 5% =  0.012650288939041299
# 95% =  0.018901698229179115
# -> 2Q Clifford avg infidelity = 0.01549 +- 0.003126, with 95% confidence interval
# interleaved sequence alpha = 0.9536
# original data with free alpha: interleaved gate avg infidelity = 0.01212
# bootstrap resampling average infidelity = 0.012059
# 5% =  0.006499520214681392
# 95% =  0.01771077485014816
# -> interleaved gate avg infidelity = 0.01212 +- 0.005606, with 95% confidence interval

