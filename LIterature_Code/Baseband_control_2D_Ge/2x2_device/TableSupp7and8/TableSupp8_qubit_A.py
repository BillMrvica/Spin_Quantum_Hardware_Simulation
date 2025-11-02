# this file is modified from Supp\SPAM\rho_povm_q1GST.py
import numpy as np

#  main_XYI_length128_2023-08-23_01-33-58_re_thresholding.html
# CPTP
rho=  [[-0.013422417582771802,0.032463534134113126],[0.9675364658658869,0.004478388591321703]]
rho = np.array(rho)[::-1,:]
np.set_printoptions(precision=3, suppress=True)
# print(rho)


povm = []

tempvar = [[-0.01004862509276642,0.05506240599443024],[0.9998513333841017,0.006190106829317138]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]



tempvar = [[0.01004862509276642,0.9449375940055698],[0.00014866661589824082,-0.006190106829317138]]
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
print('qubit A')
print(SPAM_estimate)

# [[0.96905 0.0859 ]
#  [0.03095 0.9141 ]]


np.average(np.diag(SPAM_estimate))


#%%
#  HS_main_XYI_length128_2023-08-23_01-33-58_re_thresholding.html
# H+S
rho= [[-0.013542693473274667,0.03911583012540881],[0.9608841698745911,0.00648786991575501]]
rho = np.array(rho)[::-1,:]
np.set_printoptions(precision=3, suppress=True)
# print(rho)


povm = []

tempvar = [[-0.010565211952917831,0.018881441515646746],[0.9811185584843533,0.0070086090490420515]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]



tempvar = [[0.010565211952917831,0.9811185584843533],[0.018881441515646746,-0.0070086090490420515]]
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
print('qubit A')
print(SPAM_estimate)
# [[0.94332 0.05671]
#  [0.05668 0.94329]]



np.average(np.diag(SPAM_estimate))




















