# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 23:19:31 2023

@author: TUD278427
"""
def data_setvaribles(d):
    setpoint_variables  = [key for key in d.arrays.keys() if getattr(d, key).is_setpoint] 
    return setpoint_variables, len(setpoint_variables)

def data2d_set_names(d):
    for key in data_setvaribles(d)[0]:
        setpoint1 = getattr(d, key).ndarray
        if len(setpoint1.shape) ==2:
            xname = key
        elif len(setpoint1.shape) ==1:
            yname = key
    return xname, yname    
    # for key in d.arrays.keys():
    #     if getattr(d, key).is_setpoint:
    #         setpoint1 = getattr(d, key).ndarray
    #         if len(setpoint1.shape) ==2:
    #             xname = key
    #         elif len(setpoint1.shape) ==1:
    #             yname = key
    # return xname, yname


def data2d_attr(d,  xyattr_list = ['ndarray', 'name', 'unit'], zname=None, zattr_list = ['ndarray', 'name', 'unit']):
    xname, yname = data2d_set_names(d)
    xattr = []
    yattr = []
    for att in xyattr_list:
        if att == 'ndarray':
            xattr += [getattr(getattr(d, xname), att)[0,:]]
        else:
            xattr += [getattr(getattr(d, xname), att),]
        yattr += [getattr(getattr(d, yname), att),] 
    zattr = []
    if zname is not None:
        for att in zattr_list:
            zattr += [getattr(getattr(d, zname), att) ]
    xattr = xattr[0] if len(xattr)==1 else xattr
    yattr = yattr[0] if len(yattr)==1 else yattr
    zattr = zattr[0] if len(zattr)==1 else zattr
    return xattr, yattr, zattr


def data1d_attr(d,  xattr_list = ['ndarray', 'name', 'unit'], zname=None, zattr_list = ['ndarray', 'name', 'unit']):
    xname = data_setvaribles(d)[0][0]
    xattr = []
    for att in xattr_list:
        xattr += [getattr(getattr(d, xname), att),]
    zattr = []
    if zname is not None:
        for att in zattr_list:
            zattr += [getattr(getattr(d, zname), att) ]
    xattr = xattr[0] if len(xattr)==1 else xattr
    zattr = zattr[0] if len(zattr)==1 else zattr
    return xattr, zattr


# def fullpath_datetime():
# try:
#     datetime.datetime.strptime('2023-07-1\\' +    '10-21-00', '%Y-%m-%d\\%H-%M-%S')
# except ValueError as Err:    
#     pass

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks

def hist_auto_su(y0s, y1s, expect_pos=None, kernal=np.array([1,7,21,35,35,21,7,1])/128, conv_mode='valid'):
    # suitable for the case of poor readout fidelity: find the highest two peaks
    # can also switch by giving the expected position
    # y0s = sensorval_ch2_set[0,0,:] 
    # y1s = hist_ch2[0,0,:]    
    # expect_pos: length must be 2

    #if y1s[14] ==  0.015 and y1s[15] ==  0.0825:
    #   print('!!!')


    ker = kernal #np.array([1,5,10,10,5,1])/32
    y1s_1 = np.convolve(y1s, ker, conv_mode)
    y0s_1 = y0s[len(ker)//2:-len(ker)//2+1]
    
    
    
    # np.argsort(y0s_1)
    peaks, properties = find_peaks(y1s_1)
    
    #if y1s[14] ==  0.015 and y1s[15] ==  0.0825:
    #    print(kernal)
    #    plt.figure()
    #    plt.plot(  y1s_1  )       
       
    peaks = np.sort(peaks)
    # if 0.0024 <= y1s[13] <=  0.0026 and 0.009 <= y1s[14] <=  0.011:
    #     print(peaks)
    #     # plt.figure()
    #     # plt.plot(  y1s_1  )     
    #     # plt.plot(peaks[0],  y1s_1[peaks[0]], 'ro'  )
    #     # plt.plot(peaks[1],  y1s_1[peaks[1]], 'ro'  )
    #     print(expect_pos)
        
        
        
    if len(peaks)<2:
        print('scipy.signal.find_peaks() fail. It found one or no peak.')
        p0_1, p1_1 = 0, 1
    elif len(peaks)==2:
        p0_1, p1_1 = peaks
    else:
        if expect_pos == None:
            heights = y1s_1[peaks]
            p0_1 = peaks[np.argsort(heights)[-1]]
            p1_1 = peaks[np.argsort(heights)[-2]]
         
        else:
            expect_pos_idx = [np.argmin(abs(y0s-expect_pos[i])) for i in range(2) ]

            # if len(expect_pos_idx) == 0:
            #     return None
            pmid = np.average(expect_pos_idx)-len(ker)//2
            pL = peaks[peaks<pmid]
            pR = peaks[peaks>=pmid]
            
            # print(len(expect_pos_idx))
            # print(expect_pos_idx)
            # print(peaks[peaks<pmid])            
            if len(pL) > 0:
                p0_1 = pL[np.argmin( abs(expect_pos_idx[0]-pL)  )]
            else:
                p0_1 = np.argmin(np.sum(y0s_1*y1s_1)/np.sum(y1s_1) - y0s_1)
                
                
            if len(pR) > 0:
                p1_1 = pR[np.argmin( abs(expect_pos_idx[1]-pR)  )]
            else:
                p1_1 = np.argmin(np.sum(y0s_1*y1s_1)/np.sum(y1s_1) - y0s_1)
            
    # print(len(peaks))    
    p0 = p0_1+ len(ker)//2
    p1 = p1_1+ len(ker)//2 
    vth = (y0s[p0] + y0s[p1])/2
    su = np.sum(y1s[y0s<vth])
    
        
    # plt.figure()
    # plt.plot(y0s, y1s)
    # # plt.plot(y0s, np.sum(hist_ch2[0,:,:],axis=0)/hist_ch2.shape[1])
    # plt.plot(y0s_1, y1s_1)
    # plt.scatter(x=y0s_1[p0_1], y=y1s_1[p0_1] , c='r',s=50, marker='o')
    # plt.scatter(x=y0s_1[p1_1], y=y1s_1[p1_1] , c='r',s=50, marker='o')
    # plt.scatter(x=y0s[p0], y=y1s[p0] , c='g',s=30, marker='o')
    # plt.scatter(x=y0s[p1], y=y1s[p1] , c='g',s=30, marker='o')
    
    # if 0.0012 <= y1s[12] <=  0.0013 and 0.0062 <= y1s[13] <=  0.0063:
        # print(p0, p1)
        # print(expect_pos)
        # print(y0s[p0], y0s[p1])
        # plt.figure()
        # plt.plot(  y1s_1  )     
        # plt.plot(peaks[0],  y1s_1[peaks[0]], 'ro'  )
        # plt.plot(peaks[1],  y1s_1[peaks[1]], 'ro'  )

    # if 0.0024 <= y1s[13] <=  0.0026 and 0.009 <= y1s[14] <=  0.011:
        # print(p0, p1)
        # print(expect_pos)
        # print(pL,pR)
        # print(y0s[p0], y0s[p1])
        # plt.figure()
        # plt.plot(  y1s_1  )     
        # plt.plot(peaks[0],  y1s_1[peaks[0]], 'ro'  )
        # plt.plot(peaks[1],  y1s_1[peaks[1]], 'ro'  )
        # print(expect_pos)
        
    return su, y0s[p0], y0s[p1], vth, p0, p1



def hist_auto_1d(sensorval_ch2_set, hist_ch2, xdata, zdata=None, plot_thresholding=False, kernal = np.array([1,7,21,35,35,21,7,1])/128, expect_pos_global=True, conv_mode='valid'):
    # sensorval_ch2_set = histlist[line_id][0]
    # hist_ch2 = histlist[line_id][1]
    
    su_autos = np.zeros(hist_ch2.shape[0:1])
    vths = np.zeros(hist_ch2.shape[0:1])
    peak01s = np.zeros((hist_ch2.shape[0],2)  )
    
    if expect_pos_global:
        expect_pos = hist_auto_su(sensorval_ch2_set[0,:], np.average(hist_ch2[:,:],axis=0), kernal=kernal, conv_mode=conv_mode)[1:3]
    else:
        expect_pos=None
    
    for j in range(hist_ch2.shape[0]):
        # su_auto, peak0, peak1, vth = hist_auto_su(sensorval_ch2_set[j,:], hist_ch2[j,:])[0:4]  
        # print(j)
        # if hist_ch2[j,41] ==  0.015 and hist_ch2[j,15] ==  0.0825:
        #     print(j, hist_ch2.shape[0])
        #     print(kernal)
        # if j==40:
        #     print(hist_ch2[j,:])
            
        su_auto, peak0, peak1, vth = hist_auto_su(sensorval_ch2_set[j,:], hist_ch2[j,:], 
                                                  expect_pos=expect_pos, kernal=kernal , conv_mode=conv_mode)[0:4] 
        su_autos[j] = su_auto
        vths[j] = vth
        peak01s[j,0] = peak0
        peak01s[j,1] = peak1
    
    if plot_thresholding:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,5))
        x0s = sensorval_ch2_set[0,:]
        x1s = xdata
        x2s = hist_ch2[:,:]
        
        # x0s_shifted1 = np.append(x0s, 2*x0s[-1] - x0s[-2])
        # x1s_shifted1 = np.append(x1s, 2*x1s[-1] - x1s[-2])
        # x0s_shifted2 = np.append( np.array([2*x0s[0] - x0s[1]]), x0s)
        # x1s_shifted2 = np.append( np.array([2*x1s[0] - x1s[1]]), x1s)
        # x0s_display = (x0s_shifted2 + x0s_shifted1)/2
        # x1s_display = (x1s_shifted2 + x1s_shifted1)/2
        # X0_display, X1_display = np.meshgrid(x0s_display, x1s_display, indexing='ij')
        # img = ax1.pcolormesh( X0_display, X1_display, x2s.transpose()  , shading='auto')
        
        img = ax1.pcolormesh( x0s, x1s, x2s , shading='auto')
        ax1.set_xlabel('sensorval_ch2_set (ns)')
        ax1.set_ylabel(' ()')
        # ax1.set_title( start_time+ '\n' + 'waittime vs ramptime')
        ax1.plot(vths[:] , x1s)
        ax1.plot(peak01s[:,0] , x1s)
        ax1.plot(peak01s[:,1] , x1s)

        ax2.plot(1-su_autos[:] , x1s)
        if zdata is not None:
            ax2.plot(zdata , x1s)
    return su_autos, peak01s[:,0], peak01s[:,1], vths





