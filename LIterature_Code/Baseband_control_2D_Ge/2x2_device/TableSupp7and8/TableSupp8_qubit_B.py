# this file is modified from Supp\SPAM\rho_povm_q2GST.py
import numpy as np

# main_dataset_XYI_length128_2023-08-23_01-09-15.html
# CPTP
rho= [[-0.0018049466031429837,0.04571634367422284],[0.9542836563257772,0.0001047001260929219]]
rho = np.array(rho)[::-1,:]
np.set_printoptions(precision=3, suppress=True)
# print(rho)


povm = []

tempvar = [[-0.0010972679908136768,0.03659566359206817],[0.9932464988936647,-0.0003102446134280004]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]



tempvar = [[0.0010972679908136768,0.9634043364079318],[0.0067535011063353245,0.0003102446134280004]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]    
    



sx = np.array([[0,1],[1,0]])
s0 = np.eye(2)



SPAM_estimate = np.zeros((2,2))
for i in range(2):
    tempvar = np.trace(povm[i] @ rho)
    # print(tempvar)
    SPAM_estimate[i,0] = tempvar


for i in range(2):
    tempvar = np.trace(povm[i] @ sx @ rho @ sx.conj().T)
    # print(tempvar)
    SPAM_estimate[i,1] = tempvar


np.set_printoptions(precision=5, suppress=True)
print('CPTP mode in GST estimate')
print('qubit B')
print(SPAM_estimate)

# [[0.94951 0.08033]
#  [0.05049 0.91967]]


np.average(np.diag(SPAM_estimate))


#%%
# HS_main_dataset_XYI_length128_2023-08-23_01-09-15.html
# H+S
rho= [[-0.0017326961767063977,0.046473845379061685],[0.9535261546209384,0.0005473832275447892]]
rho = np.array(rho)[::-1,:]
np.set_printoptions(precision=3, suppress=True)
# print(rho)


povm = []

tempvar = [[-0.0008880774152869858,0.022132024915769044],[0.977867975084231,-0.0005269123566124811]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]



tempvar = [[0.0008880774152869858,0.977867975084231],[0.022132024915769044,0.0005269123566124811]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]    
    



sx = np.array([[0,1],[1,0]])
s0 = np.eye(2)



SPAM_estimate = np.zeros((2,2))
for i in range(2):
    tempvar = np.trace(povm[i] @ rho)
    # print(tempvar)
    SPAM_estimate[i,0] = tempvar


for i in range(2):
    tempvar = np.trace(povm[i] @ sx @ rho @ sx.conj().T)
    # print(tempvar)
    SPAM_estimate[i,1] = tempvar


np.set_printoptions(precision=5, suppress=True)
print('H+S mode in GST estimate')
print('qubit B')
print(SPAM_estimate)
# [[0.93345 0.06655]
#  [0.06655 0.93345]]



np.average(np.diag(SPAM_estimate))


