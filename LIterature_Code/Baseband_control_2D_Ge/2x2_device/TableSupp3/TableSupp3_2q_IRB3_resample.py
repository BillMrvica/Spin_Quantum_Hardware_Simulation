# this file is modified from Supp\2q_RB_all\2qRB_reorder_2023-08-23_09-20-24_CI_supp.py
# 2Q RB bootstrap resampling to get the confidence interval of the average clifford infidelity
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
from datetime import datetime

N_resample = 10000
results = []



fdir = ''
all_sequences_data = pickle.load(  open(fdir + '2023-08-23_09-20-24_2qRB_reorder.pickle',  'rb')      )

for i_resample in range(N_resample):
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
            Nshots = len(all_sequences_data[sequence_idx][i_interleaved]['shots'][1][0])
            # re_arr = np.arange(L, dtype=int) 
            re_arr = np.random.randint(0, L, L)

            re_arr_psbq1q4 = []
            re_arr_psbq2q3 = []
            for i_re_arr in re_arr:
                # shots_arr = np.arange(Nshots, dtype=int)
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
        
        all_psb_q1q4_hists = np.array(all_psb_q1q4_hists)
        all_psb_q2q3_hists = np.array(all_psb_q2q3_hists)
        
        hist_range14_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q1q4[0],hist_range_q1q4[1],N_bins), range(1,num_sequences)[::-1])
        hist_range23_mesh, sequence_idx_mesh = np.meshgrid(np.linspace(hist_range_q2q3[0],hist_range_q2q3[1],N_bins), range(1,num_sequences)[::-1])
        
        VMAX=None

    
        autothreshold_rb_irb[i_interleaved] = np.array([p_autothreshold[ii,:]/np.sum(p_autothreshold[ii,:]) for ii in range(1, num_sequences)])
    
    ###
    
    

    
    def expdecay_alpha(x, A, decayrate, alpha, y0):
        return A*(decayrate**(x**alpha)) + y0 
    
    def expdecay(x, A, decayrate, y0):
        return A*(decayrate**x) + y0
    
    xdata = np.array(sequence_length)
    
      
    zdata_rb = autothreshold_rb_irb['RB'][:,0]
    p0 = [0.65, 0.92,  0.25]
    _, p1_rb = fit_data(xdata, zdata_rb, func=expdecay,p0=p0, plot=False, return_cov=True)
    
    
     
    zdata_irb = autothreshold_rb_irb['IRB'][:,0]
    p0 = [0.65, 0.85,  0.25]
    _, p1_irb = fit_data(xdata, zdata_irb, func=expdecay,p0=p0, plot=False, return_cov=True)
    
    
    
    decayrate_ref = p1_rb[1]
    decayrate_interleaved = p1_irb[1]
    
    
    d = 4 # 2 qubits
    avg_err_ref = (1-decayrate_ref)*(d-1)/d
    avg_interleaved_gate_infidelity = (1-decayrate_interleaved/decayrate_ref)*(d-1)/d

    
    # print('#######')
    # print('If enforce a simple exponential decay, alpha=1:')
    # print(f'reference sequence decay rate  = {p1_rb[1]:.4g}')
    # print(f'interleaved sequence decay rate = {p1_irb[1]:.4g}')
    # print(f'2Q Clifford infidelity = {avg_err_ref:.4g}')
    # print(f'interleaved CZ gate infidelity = {avg_interleaved_gate_infidelity:.4g}')
    # print(f'Infinite reference sequence will decay to {p1_rb[2]:.4g}')
    # print(f'Infinite interleaved sequence will decay to {p1_irb[2]:.4g}')
    
    p0 = [0.65, 0.92, 1,  0.25] 
    cov_rb_alpha, p1_rb_alpha = fit_data(xdata, zdata_rb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
    p0 = [0.65, 0.92, 1,  0.25] 
    cov_irb_alpha, p1_irb_alpha = fit_data(xdata, zdata_irb, func=expdecay_alpha,p0=p0, plot=False, return_cov=True)
    decayrate_ref_alpha = p1_rb_alpha[1] 
    decayrate_interleaved_alpha = p1_irb_alpha[1] 
    avg_err_ref_alpha = (1-decayrate_ref_alpha)*(d-1)/d
    avg_interleaved_gate_infidelity_alpha = (1-decayrate_interleaved_alpha/decayrate_ref_alpha)*(d-1)/d
    
    # print('#######')
    # print('If not enforce alpha=1 but a free fitting parameter:')
    # print(f'reference sequence alpha  = {p1_rb_alpha[2]:.4g}')
    # print(f'interleaved sequence alpha = {p1_irb_alpha[2]:.4g}')
    # print(f'reference sequence decay rate  = {p1_rb_alpha[1]:.4g}')
    # print(f'interleaved sequence decay rate = {p1_irb_alpha[1]:.4g}')
    # print(f'2Q Clifford infidelity = {avg_err_ref_alpha:.4g}')
    # print(f'interleaved CZ gate infidelity = {avg_interleaved_gate_infidelity_alpha:.4g}')
    # print(f'Infinite reference sequence will decay to {p1_rb_alpha[3]:.4g}')
    # print(f'Infinite interleaved sequence will decay to {p1_irb_alpha[3]:.4g}')
    
    
    
    xx = np.linspace(min(xdata), max(xdata), 10*len(xdata))
    
    idx=0

    # plt.subplot(236)
    # plt.scatter(xdata, zdata_rb , c='k', s=25, label = 'RB')
    # plt.plot(xx, expdecay(xx, *p1_rb), 'k', label='RB ' )


    # plt.scatter(xdata, zdata_irb , c='C2', s=25, label = 'IRB')
    # plt.plot(xx, expdecay(xx, *p1_irb), 'C2', label='IRB '  )

    # plt.legend()
        
        
    results += [ dict(xdata=xdata, zdata_rb=zdata_rb, zdata_irb=zdata_irb,
                      p1_rb=p1_rb, p1_irb=p1_irb, 
                      decayrate_ref=decayrate_ref, decayrate_interleaved=decayrate_interleaved, 
                      avg_err_ref=avg_err_ref, avg_interleaved_gate_infidelity=avg_interleaved_gate_infidelity,
                      p1_rb_alpha=p1_rb_alpha, p1_irb_alpha=p1_irb_alpha, 
                      decayrate_ref_alpha=decayrate_ref_alpha, decayrate_interleaved_alpha=decayrate_interleaved_alpha, 
                      avg_err_ref_alpha=avg_err_ref_alpha, avg_interleaved_gate_infidelity_alpha=avg_interleaved_gate_infidelity_alpha, 
                      ) ]
    
    
    if i_resample%500 == 0:
        if i_resample==0:
            print('will take tens of minutes to do resampling')        
        print( str(i_resample)+ '/' + str(N_resample) + '   ' + datetime.now().strftime('%Y-%m-%d_%H-%M-%S') )

#%%  
# import pickle
# f = open('new_2023-08-23_09-20-24_2qRB_reorder_bootstrap_resample.pickle', 'wb')
# pickle.dump( dict(v=v, u=u, avg_err_ref_list=avg_err_ref_list, 
#                   avg_interleaved_gate_infidelity=avg_interleaved_gate_infidelity,
#                   results=results ) , f  )
# f.close()

#%%

