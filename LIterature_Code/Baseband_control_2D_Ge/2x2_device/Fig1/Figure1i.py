# this file is modified from q1andq2.py and q1andq2_CI.py

# -*- coding: utf-8 -*-

import os
path=str(os.getcwd())
path_notebook=str(os.getcwd())
path_notebook=path_notebook.replace('\\Fig1','')
import sys
sys.path.insert(1, path_notebook)
import helper_functions
from helper_functions import hist_auto_su

#%% to plot histograms for readout check as function of sequence

import numpy as np
import pickle 
from matplotlib import pyplot as plt


# from script_shuttling_15052023.helper_functions import  hist_auto_su

plt.rcParams.update({'font.size':12})

autothreshold_rb_irbs = []
fdir = ''
# fdir = r'C:\Users\TUD278427\env_qconstruct\spin-projects\stations\LD400Top_v3\script_shuttling_15052023\working\fromSander_RB\data_1q\\'

for fname in ['2023-08-30_13-37-48_1qRB_reorder.pickle',  '2023-08-30_05-20-12_1qRB_reorder.pickle' ]:
    all_sequences_data = pickle.load(  open(fdir + fname,  'rb')      )
    
    
    sequence_length = all_sequences_data[0]['sequence_length']
    hist_dict = all_sequences_data[0]['hist_dict']
    #all_sequences_data = gatestatement
    
    autothreshold_rb_irb = {}
    for i_interleaved in ['RB',   ]:
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
            vth0 = hist_auto_su( psb14_bincenters, psb_q1q4_hists   )[3]
            # vth0 = 159
            psb14_vth += [ vth0 ]   
            
            psb23_bincenters = (psb23_binedges[1:] + psb23_binedges[:-1])/2
            vth1 = hist_auto_su( psb23_bincenters, psb_q2q3_hists   )[3]
            # vth1 = 137
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
        
        # hist_range14_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins), range(1,num_sequences)[::-1])
        # hist_range23_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins), range(1,num_sequences)[::-1])
        
        sequence_idx_mesh = list(range(1,num_sequences))[::-1]
        hist_range14_mesh = np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins) 
        hist_range23_mesh = np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins)
        
        VMAX=None
        plt.figure(figsize=(10,5))
        plt.subplot(121)
        plt.pcolormesh(hist_range14_mesh,sequence_idx_mesh,all_psb_q1q4_hists,vmax=VMAX)
        plt.plot( psb14_vth, sequence_idx_mesh  )
        plt.gca().set_yticks(sequence_idx_mesh )
        plt.gca().set_yticklabels(  [str(sequence_length[i]) for i in range(len(sequence_length)) ]  )
        plt.plot( (166,)*2, (sequence_idx_mesh[0], sequence_idx_mesh[-1])  )
        
        plt.subplot(122)
        plt.pcolormesh(hist_range23_mesh,sequence_idx_mesh,all_psb_q2q3_hists,vmax=VMAX)
        plt.plot( psb23_vth, sequence_idx_mesh  )
        plt.gca().set_yticks(sequence_idx_mesh  )
        plt.gca().set_yticklabels(  [str(sequence_length[i]) for i in range(len(sequence_length)) ]  )
        plt.plot( (147,)*2, (sequence_idx_mesh[0], sequence_idx_mesh[-1])  )
        plt.show()
        
        # for i in range(1,num_sequences-1):
        #     plt.figure(figsize=(4,4))
        #     plt.title(i_interleaved + str(sequence_length[i]) )
        #     plt.plot( hist_range14_mesh[0,:],  all_psb_q1q4_hists[i,:]  )
        #     plt.plot( (psb14_vth[i],)*2,  (0, np.max(all_psb_q1q4_hists[i,:])) , 'r--' )
    
        # for i in range(1,num_sequences-1):
        #     plt.figure(figsize=(4,4))
        #     plt.title(i_interleaved + str(sequence_length[i]) )
        #     plt.plot( hist_range23_mesh[0,:],  all_psb_q2q3_hists[i,:]  )
        #     plt.plot( (psb23_vth[i],)*2,  (0, np.max(all_psb_q2q3_hists[i,:])) , 'r--' )
    
            
        autothreshold_rb_irb[i_interleaved] = np.array([p_autothreshold[ii,:]/np.sum(p_autothreshold[ii,:]) for ii in range(1, num_sequences)])
    autothreshold_rb_irbs +=  [autothreshold_rb_irb]
###

zdata_rbs = {}
zdata_rbs['q1'] =  autothreshold_rb_irbs[0]['RB'][:,0] + autothreshold_rb_irbs[0]['RB'][:,1]
zdata_rbs['q2'] =  autothreshold_rb_irbs[1]['RB'][:,0] + autothreshold_rb_irbs[1]['RB'][:,2]

###
#%% overlap q1 and q2

from projects.notebook_tools.notebook_tools import  fit_data

def expdecay_alpha(x, A, decayrate, alpha, y0):
    return A*(decayrate**(x**alpha)) + y0 

def expdecay(x, A, decayrate, y0):
    return A*(decayrate**x) + y0

xdata = np.array(sequence_length)

plt.figure()
for qubit_label in ['q1', 'q2']:
    zdata_rb = zdata_rbs[qubit_label]
    p0 = [0.4, 0.99,  0.5]
    _, p1_rb = fit_data(xdata, zdata_rb, func=expdecay,p0=p0, plot=False, return_cov=True)
    
    
    
    
    
    
    decayrate_ref = p1_rb[1]
    
    d = 2 # 1 qubits
    avg_err_ref = (1-decayrate_ref)*(d-1)/d
    
    
    p1_rb_99 = np.array(p1_rb)
    p1_rb_99[1] = (1-0.01*d/(d-1))
    
    print('#######')
    print('If enforce a simple exponential decay, alpha=1:')
    print(f'reference sequence decay rate  = {p1_rb[1]:.4g}')
    print(qubit_label + f' single-Q Clifford infidelity = {avg_err_ref:.4g}')
    print(f'Infinite sequence will decay to {p1_rb[2]:.4g}')
    
    
    
    p0 = [0.4, 0.99, 1,  0.5] 
    cov_rb_alpha, p1_rb_alpha = fit_data(xdata, zdata_rb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
    decayrate_ref_alpha = p1_rb_alpha[1] 
    avg_err_ref_alpha = (1-decayrate_ref_alpha)*(d-1)/d
    
    
    print('#######')
    print('If not enforce alpha=1 but a free fitting parameter:')
    print(f'reference sequence alpha  = {p1_rb_alpha[2]:.4g}')
    print(f'reference sequence decay rate  = {p1_rb_alpha[1]:.4g}')
    print(qubit_label + f' single-Q Clifford infidelity = {avg_err_ref_alpha:.4g}')
    print(f'Infinite sequence will decay to {p1_rb_alpha[3]:.4g}')
    
    

    
    xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
    
    idx=0
    
    plt.scatter(xdata, zdata_rb ,  s=25, label = 'RB')
    plt.plot(xx, expdecay(xx, *p1_rb),  label='RB ' )
    # plt.plot(xx, expdecay_alpha(xx, *p1_rb_alpha) , label=f'RB alpha={p1_rb_alpha[2]:.3g}')
    # plt.plot(xx, expdecay(xx, *p1_rb_99), 'k:'  )
    # plt.plot(xx, expdecay(xx, *p1_rb_98), 'k:'  )
    
    
    
    # plt.xscale('log')
    plt.legend()


#%%
from mpl_toolkits.axes_grid1 import Divider, Size


CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
colors = dict(q1=CQ1, q2=CQ2)

fig = plt.figure( figsize=(6,4))
h = [Size.Fixed(0.5), Size.Fixed(2.5), Size.Fixed(0.7), Size.Fixed(1.0)]
v = [Size.Fixed(0.5), Size.Fixed(1.0), Size.Fixed(0.7), Size.Fixed(1.0)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)

ax11 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=1))

ax13 = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=3))

axs = dict(q1=ax11, q2=ax13)

for qubit_label in ['q1', 'q2']:
    zdata_rb = zdata_rbs[qubit_label]
    p0 = [0.4, 0.99,  0.5]
    _, p1_rb = fit_data(xdata, zdata_rb, func=expdecay,p0=p0, plot=False, return_cov=True)
    
    
    
    
    
    
    decayrate_ref = p1_rb[1]
    
    d = 2 # 1 qubits
    avg_err_ref = (1-decayrate_ref)*(d-1)/d
    
    
    p1_rb_99 = np.array(p1_rb)
    p1_rb_99[1] = (1-0.01*d/(d-1))
    
    print('#######')
    print('If enforce a simple exponential decay, alpha=1:')
    print(f'reference sequence decay rate  = {p1_rb[1]:.4g}')
    print(qubit_label + f' single-Q Clifford infidelity = {avg_err_ref:.4g}')
    print(f'Infinite sequence will decay to {p1_rb[2]:.4g}')
    
    
    
    p0 = [0.4, 0.99, 1,  0.5] 
    cov_rb_alpha, p1_rb_alpha = fit_data(xdata, zdata_rb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
    decayrate_ref_alpha = p1_rb_alpha[1] 
    avg_err_ref_alpha = (1-decayrate_ref_alpha)*(d-1)/d
    
    
    print('#######')
    print('If not enforce alpha=1 but a free fitting parameter:')
    print(f'reference sequence alpha  = {p1_rb_alpha[2]:.4g}')
    print(f'reference sequence decay rate  = {p1_rb_alpha[1]:.4g}')
    print(qubit_label + f' single-Q Clifford infidelity = {avg_err_ref_alpha:.4g}')
    print(f'Infinite sequence will decay to {p1_rb_alpha[3]:.4g}')
    
    

    
    xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
    

    ax = axs[qubit_label]
    ax.scatter(xdata, zdata_rb , c=colors[qubit_label], s=8, label = 'RB')
    ax.plot(xx, expdecay(xx, *p1_rb),  'k', linewidth=0.5, label='RB ' )
    # ax.plot(xx, expdecay_alpha(xx, *p1_rb_alpha) , label=f'RB alpha={p1_rb_alpha[2]:.3g}')
    # ax.plot(xx, expdecay(xx, *p1_rb_99), 'k:'  )
    # ax.plot(xx, expdecay(xx, *p1_rb_98), 'k:'  )
    
    
    # ax.set_xlim(0, 80 )
    ax.set_xticks([0, 1000, 2000, 3000, 4000, 5000, 6000]) 
    # ax.set_xticks( np.arange(0, 210, 25)) 
    # ax.set_ylim(0,1)
    ax.set_yticks([0.5, 0.9])   
    
    # plt.xscale('log')
    # plt.legend()

# filename = 'Figure1i' #'Fig_1qRB'
# plt.savefig(filename+'.pdf')





# printing :
    
# #######
# If enforce a simple exponential decay, alpha=1:
# reference sequence decay rate  = 0.9993
# q1 single-Q Clifford infidelity = 0.0003284
# Infinite sequence will decay to 0.5456
# #######
# If not enforce alpha=1 but a free fitting parameter:
# reference sequence alpha  = 0.8332
# reference sequence decay rate  = 0.998
# q1 single-Q Clifford infidelity = 0.001012
# Infinite sequence will decay to 0.5197
# #######
# If enforce a simple exponential decay, alpha=1:
# reference sequence decay rate  = 0.9992
# q2 single-Q Clifford infidelity = 0.0003982
# Infinite sequence will decay to 0.5469
# #######
# If not enforce alpha=1 but a free fitting parameter:
# reference sequence alpha  = 0.8066
# reference sequence decay rate  = 0.9971
# q2 single-Q Clifford infidelity = 0.001469
# Infinite sequence will decay to 0.5253
# #######
# If enforce a simple exponential decay, alpha=1:
# reference sequence decay rate  = 0.9993
# q1 single-Q Clifford infidelity = 0.0003284
# Infinite sequence will decay to 0.5456
# #######
# If not enforce alpha=1 but a free fitting parameter:
# reference sequence alpha  = 0.8332
# reference sequence decay rate  = 0.998
# q1 single-Q Clifford infidelity = 0.001012
# Infinite sequence will decay to 0.5197
# #######
# If enforce a simple exponential decay, alpha=1:
# reference sequence decay rate  = 0.9992
# q2 single-Q Clifford infidelity = 0.0003982
# Infinite sequence will decay to 0.5469
# #######
# If not enforce alpha=1 but a free fitting parameter:
# reference sequence alpha  = 0.8066
# reference sequence decay rate  = 0.9971
# q2 single-Q Clifford infidelity = 0.001469
# Infinite sequence will decay to 0.5253







#%% collect the data of 1Q RB data point distribution

import numpy as np
import pickle 
from matplotlib import pyplot as plt


# from script_shuttling_15052023.helper_functions import  hist_auto_su

plt.rcParams.update({'font.size':12})

autothreshold_rb_irbs = []
num_sequences = 26
p_autothreshold_dist = np.zeros((num_sequences,2), dtype=object)
fdir = ''
# fdir = r'C:\Users\TUD278427\env_qconstruct\spin-projects\stations\LD400Top_v3\script_shuttling_15052023\working\fromSander_RB\data_1q\\'

for fname in ['2023-08-30_13-37-48_1qRB_reorder.pickle',  '2023-08-30_05-20-12_1qRB_reorder.pickle' ]:
# fname = '2023-08-30_13-37-48_1qRB_reorder.pickle'
    all_sequences_data = pickle.load(  open(fdir + fname,  'rb')      )
    
    
    sequence_length = all_sequences_data[0]['sequence_length']
    hist_dict = all_sequences_data[0]['hist_dict']
    #all_sequences_data = gatestatement
    
    autothreshold_rb_irb = {}
    
    for i_interleaved in ['RB',   ]:
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
            vth0 = hist_auto_su( psb14_bincenters, psb_q1q4_hists   )[3]
            # vth0 = 159
            psb14_vth += [ vth0 ]   
            
            psb23_bincenters = (psb23_binedges[1:] + psb23_binedges[:-1])/2
            vth1 = hist_auto_su( psb23_bincenters, psb_q2q3_hists   )[3]
            # vth1 = 137
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
                temp_data = all_sequences_data[sequence_idx][i_interleaved]['shots'][ii][0]
                temp_list += [ np.sum(temp_data < vth0)/len(temp_data) ]  # for q1 RB
            if fname == '2023-08-30_13-37-48_1qRB_reorder.pickle':
                p_autothreshold_dist[sequence_idx,0] = list(temp_list)    
            
            temp_list = []
            for ii in range(L):
                temp_data = all_sequences_data[sequence_idx][i_interleaved]['shots'][ii][1]
                temp_list += [ np.sum(temp_data > vth1)/len(temp_data)  ]  # for q2 RB            
            if fname == '2023-08-30_05-20-12_1qRB_reorder.pickle':
                p_autothreshold_dist[sequence_idx,1] = list(temp_list)
            ####
            
        all_psb_q1q4_hists = np.array(all_psb_q1q4_hists)
        all_psb_q2q3_hists = np.array(all_psb_q2q3_hists)
        
        # hist_range14_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins), range(1,num_sequences)[::-1])
        # hist_range23_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins), range(1,num_sequences)[::-1])
        
        sequence_idx_mesh = list(range(1,num_sequences))[::-1]
        hist_range14_mesh = np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins) 
        hist_range23_mesh = np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins)
        
        VMAX=None
        # plt.figure(figsize=(10,5))
        # plt.subplot(121)
        # plt.pcolormesh(hist_range14_mesh,sequence_idx_mesh,all_psb_q1q4_hists,vmax=VMAX)
        # plt.plot( psb14_vth, sequence_idx_mesh  )
        # plt.gca().set_yticks(sequence_idx_mesh )
        # plt.gca().set_yticklabels(  [str(sequence_length[i]) for i in range(len(sequence_length)) ]  )
        # plt.plot( (166,)*2, (sequence_idx_mesh[0], sequence_idx_mesh[-1])  )
        
        # plt.subplot(122)
        # plt.pcolormesh(hist_range23_mesh,sequence_idx_mesh,all_psb_q2q3_hists,vmax=VMAX)
        # plt.plot( psb23_vth, sequence_idx_mesh  )
        # plt.gca().set_yticks(sequence_idx_mesh  )
        # plt.gca().set_yticklabels(  [str(sequence_length[i]) for i in range(len(sequence_length)) ]  )
        # plt.plot( (147,)*2, (sequence_idx_mesh[0], sequence_idx_mesh[-1])  )
        # plt.show()
        
        # for i in range(1,num_sequences-1):
        #     plt.figure(figsize=(4,4))
        #     plt.title(i_interleaved + str(sequence_length[i]) )
        #     plt.plot( hist_range14_mesh[0,:],  all_psb_q1q4_hists[i,:]  )
        #     plt.plot( (psb14_vth[i],)*2,  (0, np.max(all_psb_q1q4_hists[i,:])) , 'r--' )
    
        # for i in range(1,num_sequences-1):
        #     plt.figure(figsize=(4,4))
        #     plt.title(i_interleaved + str(sequence_length[i]) )
        #     plt.plot( hist_range23_mesh[0,:],  all_psb_q2q3_hists[i,:]  )
        #     plt.plot( (psb23_vth[i],)*2,  (0, np.max(all_psb_q2q3_hists[i,:])) , 'r--' )
    
            
        autothreshold_rb_irb[i_interleaved] = np.array([p_autothreshold[ii,:]/np.sum(p_autothreshold[ii,:]) for ii in range(1, num_sequences)])
        # autothreshold_rb_irb_dist[i_interleaved] = np.array([p_autothreshold_dist[ii,:]/np.sum(p_autothreshold_dist[ii,:]) for ii in range(1, num_sequences)])
    
    autothreshold_rb_irbs +=  [autothreshold_rb_irb]


zdata_rbs = {}
zdata_rbs['q1'] =  autothreshold_rb_irbs[0]['RB'][:,0] + autothreshold_rb_irbs[0]['RB'][:,1]
zdata_rbs['q2'] =  autothreshold_rb_irbs[1]['RB'][:,0] + autothreshold_rb_irbs[1]['RB'][:,2]

zdata_rbs['q1_dist'] =  p_autothreshold_dist[1:,0]
zdata_rbs['q2_dist'] =  p_autothreshold_dist[1:,1]
###
#%% plot the data of 1Q RB data point distribution, for both qubits
CQ1 = '#0D77BC'
CQ2 = '#DD5E27'
colors = dict(q1=CQ1, q2=CQ2)
# 

def expdecay_alpha(x, A, decayrate, alpha, y0):
    return A*(decayrate**(x**alpha)) + y0 

def expdecay(x, A, decayrate, y0):
    return A*(decayrate**x) + y0

xdata = np.array(sequence_length)

plt.figure(figsize=(14,6))
plt.suptitle('scatter points randomized sequences')
for qubit_label in ['q1', 'q2']:
    zdata_rb = zdata_rbs[qubit_label]
    p0 = [0.4, 0.99,  0.5]
    _, p1_rb = fit_data(xdata, zdata_rb, func=expdecay,p0=p0, plot=False, return_cov=True)
    
    
    
    
    
    
    decayrate_ref = p1_rb[1]
    
    d = 2 # 1 qubits
    avg_err_ref = (1-decayrate_ref)*(d-1)/d
    
    
    p1_rb_99 = np.array(p1_rb)
    p1_rb_99[1] = (1-0.01*d/(d-1))
    
    print('#######')
    print('If enforce a simple exponential decay, alpha=1:')
    print(f'reference sequence decay rate  = {p1_rb[1]:.4g}')
    print(qubit_label + f' single-Q Clifford infidelity = {avg_err_ref:.4g}')
    print(f'Infinite sequence will decay to {p1_rb[2]:.4g}')
    
    
    
    p0 = [0.4, 0.99, 1,  0.5] 
    cov_rb_alpha, p1_rb_alpha = fit_data(xdata, zdata_rb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
    decayrate_ref_alpha = p1_rb_alpha[1] 
    avg_err_ref_alpha = (1-decayrate_ref_alpha)*(d-1)/d
    
    
    print('#######')
    print('If not enforce alpha=1 but a free fitting parameter:')
    print(f'reference sequence alpha  = {p1_rb_alpha[2]:.4g}')
    print(f'reference sequence decay rate  = {p1_rb_alpha[1]:.4g}')
    print(qubit_label + f' single-Q Clifford infidelity = {avg_err_ref_alpha:.4g}')
    print(f'Infinite sequence will decay to {p1_rb_alpha[3]:.4g}')
    
    

    
    xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
    
    idx=0
    plt.subplot(120+int(qubit_label[1]))
    
    # plt.scatter(xdata, zdata_rb ,  s=25,  c=colors[qubit_label], label = 'RB')
    for i_sequence, temp_list in enumerate(zdata_rbs[qubit_label+'_dist']):
        
        
        plt.scatter((xdata[i_sequence],)*len(temp_list), temp_list  ,  s=25,  c=colors[qubit_label], )
    plt.plot(xx, expdecay(xx, *p1_rb), c='k', )
    # plt.plot(xx, expdecay_alpha(xx, *p1_rb_alpha) , label=f'RB alpha={p1_rb_alpha[2]:.3g}')
    # plt.plot(xx, expdecay(xx, *p1_rb_99), 'k:'  )
    # plt.plot(xx, expdecay(xx, *p1_rb_98), 'k:'  )
    
    
    
    # plt.xscale('log')
    plt.legend()

for err_factor, titlestr in [ (1,'std of 128 randomizations'), (np.sqrt(128),'std of 128 randomizations/sqrt(128)') ]:
    plt.figure(figsize=(14,6))
    plt.suptitle(titlestr)
    
    for qubit_label in ['q1', 'q2']:
        zdata_rb = zdata_rbs[qubit_label]
        p0 = [0.4, 0.99,  0.5]
        _, p1_rb = fit_data(xdata, zdata_rb, func=expdecay,p0=p0, plot=False, return_cov=True)
        
    
        decayrate_ref = p1_rb[1]
        
        d = 2 # 1 qubits
        avg_err_ref = (1-decayrate_ref)*(d-1)/d
        
        
        p1_rb_99 = np.array(p1_rb)
        p1_rb_99[1] = (1-0.01*d/(d-1))
        
     
    
        
        xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
        
        idx=0
        plt.subplot(120+int(qubit_label[1]))
        
        # plt.scatter(xdata, zdata_rb ,  s=25,  c=colors[qubit_label], label = 'RB')
        for i_sequence, temp_list in enumerate(zdata_rbs[qubit_label+'_dist']):
            plt.errorbar(xdata[i_sequence]  , np.mean(np.array(temp_list)),
                         yerr = np.std(np.array(temp_list))/err_factor, c='k',
                         fmt ='o', markersize=3, linestyle='', linewidth=1)
        plt.plot(xx, expdecay(xx, *p1_rb), c='k', )



#%% q1 bootstrap resampling to get the confidence interval of the average clifford infidelity
from datetime import datetime

N_resample = 10000
results = []

# for fname in ['2023-08-30_13-37-48_1qRB_reorder.pickle',  '2023-08-30_05-20-12_1qRB_reorder.pickle' ]:
fdir = ''
fname = '2023-08-30_13-37-48_1qRB_reorder.pickle'
all_sequences_data = pickle.load(  open(fdir + fname,  'rb')      )

for i_resample in range(N_resample):
    # plt.close('all')
    # plt.figure(figsize=(15,5))
    

    
    sequence_length = all_sequences_data[0]['sequence_length']
    hist_dict = all_sequences_data[0]['hist_dict']

    
    autothreshold_rb_irb = {}
    
    for i_interleaved in ['RB',   ]:
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
            
    
            L = len(all_sequences_data[sequence_idx][i_interleaved]['shots'])
            Nshots = len(all_sequences_data[sequence_idx][i_interleaved]['shots'][1][0])
            re_arr = np.random.randint(0, L, L)
            # 
            re_arr_psbq1q4 = []
            re_arr_psbq2q3 = []
            for i_re_arr in re_arr:
                shots_arr = np.random.randint(0, Nshots, Nshots)
                re_arr_psbq1q4 += [ all_sequences_data[sequence_idx][i_interleaved]['shots'][i_re_arr][0][shots_arr]  ]
                re_arr_psbq2q3 += [ all_sequences_data[sequence_idx][i_interleaved]['shots'][i_re_arr][1][shots_arr]  ]
                
            psb_q1q4 = np.hstack([re_arr_psbq1q4 ])
            psb_q2q3 = np.hstack([re_arr_psbq2q3 ])
            
            psb_q1q4_hists, psb14_binedges = np.histogram(psb_q1q4, bins=N_bins ,range=hist_range_q1q4)
            all_psb_q1q4_hists.append(psb_q1q4_hists)
        
            
            psb_q2q3_hists, psb23_binedges = np.histogram(psb_q2q3, bins=N_bins ,range=hist_range_q2q3)
            all_psb_q2q3_hists.append(psb_q2q3_hists)
        
        
        
            psb14_bincenters = (psb14_binedges[1:] + psb14_binedges[:-1])/2
            vth0 = hist_auto_su( psb14_bincenters, psb_q1q4_hists   )[3]
            # vth0 = 159
            psb14_vth += [ vth0 ]   
            
            psb23_bincenters = (psb23_binedges[1:] + psb23_binedges[:-1])/2
            vth1 = hist_auto_su( psb23_bincenters, psb_q2q3_hists   )[3]
            # vth1 = 137
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
            # temp_list = []
            # for ii in range(L):
            #     temp_data = all_sequences_data[sequence_idx][i_interleaved]['shots'][ii][0]
            #     temp_list += [ np.sum(temp_data < vth0)/len(temp_data) ]  # for q1 RB
            # if fname == '2023-08-30_13-37-48_1qRB_reorder.pickle':
            #     p_autothreshold_dist[sequence_idx,0] = list(temp_list)    
            
            # temp_list = []
            # for ii in range(L):
            #     temp_data = all_sequences_data[sequence_idx][i_interleaved]['shots'][ii][1]
            #     temp_list += [ np.sum(temp_data > vth1)/len(temp_data)  ]  # for q2 RB            
            # if fname == '2023-08-30_05-20-12_1qRB_reorder.pickle':
            #     p_autothreshold_dist[sequence_idx,1] = list(temp_list)
            ####
            
        all_psb_q1q4_hists = np.array(all_psb_q1q4_hists)
        all_psb_q2q3_hists = np.array(all_psb_q2q3_hists)
        
        # hist_range14_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins), range(1,num_sequences)[::-1])
        # hist_range23_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins), range(1,num_sequences)[::-1])
        
        sequence_idx_mesh = list(range(1,num_sequences))[::-1]
        hist_range14_mesh = np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins) 
        hist_range23_mesh = np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins)
        
        VMAX=None

        # plt.subplot(131)
        # plt.pcolormesh(hist_range14_mesh,sequence_idx_mesh,all_psb_q1q4_hists,vmax=VMAX)
        # plt.plot( psb14_vth, sequence_idx_mesh  )
        # plt.gca().set_yticks(sequence_idx_mesh )
        # plt.gca().set_yticklabels(  [str(sequence_length[i]) for i in range(len(sequence_length)) ]  )

        
        # plt.subplot(132)
        # plt.pcolormesh(hist_range23_mesh,sequence_idx_mesh,all_psb_q2q3_hists,vmax=VMAX)
        # plt.plot( psb23_vth, sequence_idx_mesh  )
        # plt.gca().set_yticks(sequence_idx_mesh  )
        # plt.gca().set_yticklabels(  [str(sequence_length[i]) for i in range(len(sequence_length)) ]  )
        # plt.show()
        
    
            
        autothreshold_rb_irb[i_interleaved] = np.array([p_autothreshold[ii,:]/np.sum(p_autothreshold[ii,:]) for ii in range(1, num_sequences)])
        # autothreshold_rb_irb_dist[i_interleaved] = np.array([p_autothreshold_dist[ii,:]/np.sum(p_autothreshold_dist[ii,:]) for ii in range(1, num_sequences)])
    
    # autothreshold_rb_irbs +=  [autothreshold_rb_irb]
    
    
    zdata_rbs = {}
    zdata_rbs['q1'] =  autothreshold_rb_irb['RB'][:,0] + autothreshold_rb_irb['RB'][:,1]
    qubit_label = 'q1'
    
    
    ##
    xdata = np.array(sequence_length)
    zdata_rb = zdata_rbs[qubit_label]
    p0 = [0.4, 0.99,  0.5]
    _, p1_rb = fit_data(xdata, zdata_rb, func=expdecay,p0=p0, plot=False, return_cov=True)
    
    
    
    
    
    
    decayrate_ref = p1_rb[1]
    
    d = 2 # 1 qubits
    avg_err_ref = (1-decayrate_ref)*(d-1)/d
    
    
    p1_rb_99 = np.array(p1_rb)
    p1_rb_99[1] = (1-0.01*d/(d-1))
    
    # print('#######')
    # print('If enforce a simple exponential decay, alpha=1:')
    # print(f'reference sequence decay rate  = {p1_rb[1]:.4g}')
    # print(qubit_label + f' single-Q Clifford infidelity = {avg_err_ref:.4g}')
    # print(f'Infinite sequence will decay to {p1_rb[2]:.4g}')
    
    
    
    p0 = [0.4, 0.99, 1,  0.5] 
    cov_rb_alpha, p1_rb_alpha = fit_data(xdata, zdata_rb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
    decayrate_ref_alpha = p1_rb_alpha[1] 
    avg_err_ref_alpha = (1-decayrate_ref_alpha)*(d-1)/d
    
    
    # print('#######')
    # print('If not enforce alpha=1 but a free fitting parameter:')
    # print(f'reference sequence alpha  = {p1_rb_alpha[2]:.4g}')
    # print(f'reference sequence decay rate  = {p1_rb_alpha[1]:.4g}')
    # print(qubit_label + f' single-Q Clifford infidelity = {avg_err_ref_alpha:.4g}')
    # print(f'Infinite sequence will decay to {p1_rb_alpha[3]:.4g}')
    
    
    results += [ dict(xdata=xdata, zdata_rb=zdata_rb, 
                      p1_rb=p1_rb, decayrate_ref=decayrate_ref, avg_err_ref=avg_err_ref, 
                      p1_rb_alpha=p1_rb_alpha, decayrate_ref_alpha=decayrate_ref_alpha, avg_err_ref_alpha=avg_err_ref_alpha) ]
    
    xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
    
    idx=0
    CQ1 = '#0D77BC'
    CQ2 = '#DD5E27'
    colors = dict(q1=CQ1, q2=CQ2)
    
    # plt.subplot(133)
    
    # plt.scatter(xdata, zdata_rb ,  s=25,  c=colors[qubit_label]   )
    # plt.plot(xx, expdecay(xx, *p1_rb), c='k', )

    
    
    # plt.legend()
    if i_resample%1000 == 0:
        if i_resample==0:
            print('will take few minutes to do resampling')
        print( str(i_resample)+ '/' + str(N_resample) + '   ' + datetime.now().strftime('%Y-%m-%d_%H-%M-%S') )
    

#%%


# u, v = np.histogram(readdd , bins=hist_dict['Nbins'][i], range= (hist_dict['min_vals'][i], hist_dict['max_vals'][i]) ) 
# v = (v[1:]+v[:-1])/2 
plt.figure()
plt.title('histogram of qubit A avg gate infidelity')
# plt.plot([1-results[i]['p1_rb'][1] for i in range(len(results))]   )
# plt.plot([results[i]['p1_rb'][2] for i in range(len(results))]   )


# import pickle
# f = open('q1_RB_bootstrap_resample.pickle', 'wb')
# pickle.dump( dict(v=v, u=u, avg_err_ref_list=avg_err_ref_list, results=results ) , f  )
# f.close()
# results = pickle.load(  open('q1_RB_bootstrap_resample.pickle', 'rb' ))['results']

avg_err_ref_list = [ results[i]['avg_err_ref'] for i in range(len(results))]
u, v = np.histogram(avg_err_ref_list  , bins=50,  ) 
v = (v[1:]+v[:-1])/2 
plt.plot(v, u)



avg_err_ref_sort = np.sort(np.array(avg_err_ref_list))
print('\n')
print('qubit 1 original data: Clifford infidelity = 0.0003284.')
print(f'bootstrap resampling average infidelity = {np.average(avg_err_ref_sort):.5g}', )
print('5% = ', avg_err_ref_sort[500] )
print('95% = ', avg_err_ref_sort[9500] )
print(f'-> qubit 1 Clifford infidelity = 0.0003284 +- {abs(avg_err_ref_sort[9500]-avg_err_ref_sort[500])/2:.4g}, with 95% confidence interval')

# qubit 1 original data: Clifford infidelity = 0.0003284.
# bootstrap resampling average infidelity = 0.00032834
# 5% =  0.00029393543282679735
# 95% =  0.0003645989540071448
# -> qubit 1 Clifford infidelity = 0.0003284 +- 3.533e-05, with 95% confidence interval


#%% q2 bootstrap resampling to get the confidence interval of the average clifford infidelity
from datetime import datetime

N_resample = 10000
results = []

# for fname in ['2023-08-30_13-37-48_1qRB_reorder.pickle',  '2023-08-30_05-20-12_1qRB_reorder.pickle' ]:
fdir = ''
fname = '2023-08-30_05-20-12_1qRB_reorder.pickle'
all_sequences_data = pickle.load(  open(fdir + fname,  'rb')      )

for i_resample in range(N_resample):
    # plt.close('all')
    # plt.figure(figsize=(15,5))
    

    
    sequence_length = all_sequences_data[0]['sequence_length']
    hist_dict = all_sequences_data[0]['hist_dict']

    
    autothreshold_rb_irb = {}
    
    for i_interleaved in ['RB',   ]:
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
            
    
            L = len(all_sequences_data[sequence_idx][i_interleaved]['shots'])
            Nshots = len(all_sequences_data[sequence_idx][i_interleaved]['shots'][1][0])
            re_arr = np.random.randint(0, L, L)
            # 
            re_arr_psbq1q4 = []
            re_arr_psbq2q3 = []
            for i_re_arr in re_arr:
                shots_arr = np.random.randint(0, Nshots, Nshots)
                re_arr_psbq1q4 += [ all_sequences_data[sequence_idx][i_interleaved]['shots'][i_re_arr][0][shots_arr]  ]
                re_arr_psbq2q3 += [ all_sequences_data[sequence_idx][i_interleaved]['shots'][i_re_arr][1][shots_arr]  ]
                
            psb_q1q4 = np.hstack([re_arr_psbq1q4 ])
            psb_q2q3 = np.hstack([re_arr_psbq2q3 ])
            
            psb_q1q4_hists, psb14_binedges = np.histogram(psb_q1q4, bins=N_bins ,range=hist_range_q1q4)
            all_psb_q1q4_hists.append(psb_q1q4_hists)
        
            
            psb_q2q3_hists, psb23_binedges = np.histogram(psb_q2q3, bins=N_bins ,range=hist_range_q2q3)
            all_psb_q2q3_hists.append(psb_q2q3_hists)
        
        
        
            psb14_bincenters = (psb14_binedges[1:] + psb14_binedges[:-1])/2
            vth0 = hist_auto_su( psb14_bincenters, psb_q1q4_hists   )[3]
            # vth0 = 159
            psb14_vth += [ vth0 ]   
            
            psb23_bincenters = (psb23_binedges[1:] + psb23_binedges[:-1])/2
            vth1 = hist_auto_su( psb23_bincenters, psb_q2q3_hists   )[3]
            # vth1 = 137
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
            # temp_list = []
            # for ii in range(L):
            #     temp_data = all_sequences_data[sequence_idx][i_interleaved]['shots'][ii][0]
            #     temp_list += [ np.sum(temp_data < vth0)/len(temp_data) ]  # for q1 RB
            # if fname == '2023-08-30_13-37-48_1qRB_reorder.pickle':
            #     p_autothreshold_dist[sequence_idx,0] = list(temp_list)    
            
            # temp_list = []
            # for ii in range(L):
            #     temp_data = all_sequences_data[sequence_idx][i_interleaved]['shots'][ii][1]
            #     temp_list += [ np.sum(temp_data > vth1)/len(temp_data)  ]  # for q2 RB            
            # if fname == '2023-08-30_05-20-12_1qRB_reorder.pickle':
            #     p_autothreshold_dist[sequence_idx,1] = list(temp_list)
            ####
            
        all_psb_q1q4_hists = np.array(all_psb_q1q4_hists)
        all_psb_q2q3_hists = np.array(all_psb_q2q3_hists)
        
        # hist_range14_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins), range(1,num_sequences)[::-1])
        # hist_range23_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins), range(1,num_sequences)[::-1])
        
        sequence_idx_mesh = list(range(1,num_sequences))[::-1]
        hist_range14_mesh = np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins) 
        hist_range23_mesh = np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins)
        
        VMAX=None

        # plt.subplot(131)
        # plt.pcolormesh(hist_range14_mesh,sequence_idx_mesh,all_psb_q1q4_hists,vmax=VMAX)
        # plt.plot( psb14_vth, sequence_idx_mesh  )
        # plt.gca().set_yticks(sequence_idx_mesh )
        # plt.gca().set_yticklabels(  [str(sequence_length[i]) for i in range(len(sequence_length)) ]  )

        
        # plt.subplot(132)
        # plt.pcolormesh(hist_range23_mesh,sequence_idx_mesh,all_psb_q2q3_hists,vmax=VMAX)
        # plt.plot( psb23_vth, sequence_idx_mesh  )
        # plt.gca().set_yticks(sequence_idx_mesh  )
        # plt.gca().set_yticklabels(  [str(sequence_length[i]) for i in range(len(sequence_length)) ]  )
        # plt.show()
        
    
            
        autothreshold_rb_irb[i_interleaved] = np.array([p_autothreshold[ii,:]/np.sum(p_autothreshold[ii,:]) for ii in range(1, num_sequences)])
        # autothreshold_rb_irb_dist[i_interleaved] = np.array([p_autothreshold_dist[ii,:]/np.sum(p_autothreshold_dist[ii,:]) for ii in range(1, num_sequences)])
    
    # autothreshold_rb_irbs +=  [autothreshold_rb_irb]
    
    
    zdata_rbs = {}
    zdata_rbs['q2'] =  autothreshold_rb_irb['RB'][:,0] + autothreshold_rb_irb['RB'][:,2]
    qubit_label = 'q2'
    
    
    ##
    xdata = np.array(sequence_length)
    zdata_rb = zdata_rbs[qubit_label]
    p0 = [0.4, 0.99,  0.5]
    _, p1_rb = fit_data(xdata, zdata_rb, func=expdecay,p0=p0, plot=False, return_cov=True)
    
    
    
    
    
    
    decayrate_ref = p1_rb[1]
    
    d = 2 # 1 qubits
    avg_err_ref = (1-decayrate_ref)*(d-1)/d
    
    
    p1_rb_99 = np.array(p1_rb)
    p1_rb_99[1] = (1-0.01*d/(d-1))
    
    # print('#######')
    # print('If enforce a simple exponential decay, alpha=1:')
    # print(f'reference sequence decay rate  = {p1_rb[1]:.4g}')
    # print(qubit_label + f' single-Q Clifford infidelity = {avg_err_ref:.4g}')
    # print(f'Infinite sequence will decay to {p1_rb[2]:.4g}')
    
    
    
    p0 = [0.4, 0.99, 1,  0.5] 
    cov_rb_alpha, p1_rb_alpha = fit_data(xdata, zdata_rb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
    decayrate_ref_alpha = p1_rb_alpha[1] 
    avg_err_ref_alpha = (1-decayrate_ref_alpha)*(d-1)/d
    
    
    # print('#######')
    # print('If not enforce alpha=1 but a free fitting parameter:')
    # print(f'reference sequence alpha  = {p1_rb_alpha[2]:.4g}')
    # print(f'reference sequence decay rate  = {p1_rb_alpha[1]:.4g}')
    # print(qubit_label + f' single-Q Clifford infidelity = {avg_err_ref_alpha:.4g}')
    # print(f'Infinite sequence will decay to {p1_rb_alpha[3]:.4g}')
    
    
    results += [ dict(xdata=xdata, zdata_rb=zdata_rb, 
                      p1_rb=p1_rb, decayrate_ref=decayrate_ref, avg_err_ref=avg_err_ref, 
                      p1_rb_alpha=p1_rb_alpha, decayrate_ref_alpha=decayrate_ref_alpha, avg_err_ref_alpha=avg_err_ref_alpha) ]
    
    xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
    
    idx=0
    CQ1 = '#0D77BC'
    CQ2 = '#DD5E27'
    colors = dict(q1=CQ1, q2=CQ2)
    
    # plt.subplot(133)
    
    # plt.scatter(xdata, zdata_rb ,  s=25,  c=colors[qubit_label]   )
    # plt.plot(xx, expdecay(xx, *p1_rb), c='k', )

    
    
    # plt.legend()
    if i_resample%1000 == 0:
        if i_resample==0:
            print('will take few minutes to do resampling')        
        print( str(i_resample)+ '/' + str(N_resample) + '   ' + datetime.now().strftime('%Y-%m-%d_%H-%M-%S') )
    
#%%



plt.figure()
plt.title('histogram of qubit B avg gate infidelity')
# plt.plot([1-results[i]['p1_rb'][1] for i in range(len(results))]   )
# plt.plot([results[i]['p1_rb'][2] for i in range(len(results))]   )

# import pickle
# f = open('q2_RB_bootstrap_resample.pickle', 'wb')
# pickle.dump( dict(v=v, u=u, avg_err_ref_list=avg_err_ref_list, results=results ) , f  )
# f.close()

avg_err_ref_list = [ results[i]['avg_err_ref'] for i in range(len(results))]
u, v = np.histogram(avg_err_ref_list  , bins=50,  ) 
v = (v[1:]+v[:-1])/2 
plt.plot(v, u)

avg_err_ref_sort = np.sort(np.array(avg_err_ref_list))
print('\n')
print('qubit 2 original data: Clifford infidelity = 0.0003982.')
print(f'bootstrap resampling average infidelity = {np.average(avg_err_ref_sort):.5g}', )
print('5% = ', avg_err_ref_sort[500] )
print('95% = ', avg_err_ref_sort[9500] )
print(f'-> qubit 2 Clifford infidelity = 0.0003982 +- {abs(avg_err_ref_sort[9500]-avg_err_ref_sort[500])/2:.4g}, with 95% confidence interval')

# qubit 2 original data: Clifford infidelity = 0.0003982.
# bootstrap resampling average infidelity = 0.00039864
# 5% =  0.0003360649323924103
# 95% =  0.00046516669905061203
# ->  qubit 2 Clifford infidelity = 0.0003982 +- 6.455e-05, with 95% confidence interval