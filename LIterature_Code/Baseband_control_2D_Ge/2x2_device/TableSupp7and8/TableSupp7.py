# this file is modified from Supp\SPAM\rho_povm_2023-08-20_15-35-39XYICphase_L8.py
import numpy as np

# CPTP
# line 18064
rho=[[-0.004205777489954665,-0.0013288483516887147,0.0006270737712282155,0.0013736575809230045],[-0.008937768364703745,-0.0006924772890706174,0.014143927616267177,0.0012435068768019445],[0.005096594694631841,0.02940941290988383,0.00038926781054240663,0.0007353621821927321],[0.9550730018929261,0.001274947987102306,0.0038296847149287636,-0.004449778517555351]]
rho = np.array(rho)[::-1,:]
np.set_printoptions(precision=3, suppress=True)
# print(rho)
# with np.set_printoptions(precision=3, suppress=True):
    # print(rho)

povm = []
# line 18192
tempvar = [[-0.0050569460408001685,-2.6964198311799255e-05,0.00015105317522456522,0.011466503335232192],[-0.008434877431959354,-0.0011582195838872032,0.07548983121397895,0.0014546501271592255],[-0.0002188023523507319,0.03431266261499205,0.0025162186943317177,0.00024698163556923164],[0.9817339549457175,-0.007750484880902476,0.007686513881790303,-0.0004402607506672732]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]


# line 18320
tempvar = [[0.00400209241488377,-0.003878740876454634,0.0013356915305164098,0.07493375101941058],[-4.845259092907308e-05,0.0004374478822894816,0.009442980900105674,-0.00036243119024548416],[-0.0002117339909825253,0.9482891073394426,-0.001243366119724633,0.0036626724994094535],[0.008964126316631876,0.00922205137097909,0.0009506621120936806,0.0025770382956718293]]   
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]    
    
# line 18448
tempvar =  [[-0.003580569308292279,-0.0009243478767941996,-0.002012070288439018,0.034972196348716206],[0.006914707775492587,0.001058756582163613,0.8931335105990561,-0.003135339647729018],[0.0004894386732586232,0.004926573623598635,-0.0016708842672949497,0.0015388806624159492],[0.00888908898097307,-0.00012564462409392525,-0.006955255206980893,-0.005409101612714552]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]   

# line 18576
tempvar = [[0.004635422934208678,0.004830052951560633,0.000525325582698043,0.8786275492966411],[0.0015686222473958395,-0.000337984880565891,0.021933677286859132,0.0020431207108152754],[-5.8902329925366106e-05,0.012471656421966792,0.0003980316926878646,-0.005448534797394634],[0.0004128297566775574,-0.0013459218659826891,-0.0016819207869030908,0.0032723240677099965]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]   



sx = np.array([[0,1],[1,0]])
s0 = np.eye(2)

UxQ1 = np.kron(sx,s0)
UxQ2 = np.kron(s0,sx)

SPAM_estimate = np.zeros((4,4))
for i in range(4):
    tempvar = np.trace(povm[i] @ rho)
    # print(tempvar)
    SPAM_estimate[i,0] = tempvar


for i in range(4):

    tempvar = np.trace(povm[i] @ UxQ2 @ rho @ UxQ2.conj().T)
    # print(tempvar)
    SPAM_estimate[i,1] = tempvar


for i in range(4):
    tempvar = np.trace(povm[i] @ UxQ1 @ rho @ UxQ1.conj().T)
    # print(tempvar)
    SPAM_estimate[i,2] = tempvar


for i in range(4):
    U = UxQ2 @ UxQ1
    tempvar = np.trace(povm[i] @ U @ rho @ U.conj().T)
    # print(tempvar)
    SPAM_estimate[i,3] = tempvar

#%%
print('CPTP mode in GST estimate')
np.set_printoptions(precision=5, suppress=True)
print('two-qubit SPAM')
print(SPAM_estimate)

# [[0.9396  0.06187 0.08647 0.01504]
#   [0.03669 0.90699 0.01267 0.08531]
#   [0.02144 0.00669 0.85407 0.05977]
#   [0.00227 0.02445 0.04679 0.83988]]


#q1 SPAM 
SPAM_q1 = np.array([
    [np.sum(SPAM_estimate[0:2,0:2])/2, np.sum(SPAM_estimate[0:2,2:4])/2],
 [np.sum(SPAM_estimate[2:4,0:2])/2, np.sum(SPAM_estimate[2:4,2:4])/2],])
print('')
print('two-qubit SPAM in qubit A subspace')
print(SPAM_q1)
# [[0.97258 0.09975]
#  [0.02742 0.90025]]

#q2 SPAM 
SPAM_q2 = np.array([
    [np.sum(SPAM_estimate[0::2,0::2])/2, np.sum(SPAM_estimate[0::2,1::2])/2],
 [np.sum(SPAM_estimate[1::2,0::2])/2, np.sum(SPAM_estimate[1::2,1::2])/2],])
print('')
print('two-qubit SPAM in qubit B subspace')
print(SPAM_q2)
# [[0.95079 0.07169]
#  [0.04921 0.92831]]
np.average(np.diag(SPAM_estimate))


#%% H+S


# line 18064
rho=[[0.0030434537224454986,-0.002768643830497375,0.001277117183097675,0.008308553049277856],[-0.030581150001569184,-0.003472171353389139,0.02132055364633989,6.304601322489837e-05],[0.0026012407759544837,0.03415499535990893,0.0048665129126561225,-0.0007205339079974782],[0.9362158979444734,-0.0020900507972824076,-0.01744150900323884,0.0002525803295357759]]
rho = np.array(rho)[::-1,:]
np.set_printoptions(precision=3, suppress=True)
# print(rho)
# with np.set_printoptions(precision=3, suppress=True):
    # print(rho)

povm = []
# line 18192
tempvar = [[0.0020349059957123303,-0.0006131823476615145,0.002684929202658397,0.005371491644768134],[0.00625140507289795,-0.00026328854025355064,0.04137817048005568,0.00857173239327421],[0.00989666165682122,0.026011446693786577,0.00040882777015622044,-0.000132483432924246],[0.9272388911813897,0.0035526979167460403,0.01720124266333021,-0.00016287065213742467]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]


# line 18320
tempvar = [[0.0010048613229086184,0.0017464988031179133,0.003496465378913742,0.04220960087898196],[0.0001228362714061871,0.003082056592135657,0.0012099143244186927,0.005938040062092479],[0.001982699025476018,0.9308644235389373,0.0018164665923344895,0.012996574589281116],[0.025716061257662015,0.00739440924273679,0.001258264795883481,-0.0009750548792228691]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]    
    
# line 18448
tempvar = [[-0.0006021729723049254,0.0002612536623860669,-0.0055868105271827005,0.02686000419465573],[-0.0063858592590425485,-0.0028931307155880973,0.9307880285161039,-0.00877317197407835],[-0.0018042684726760126,0.0010087714001012493,-0.001871625937059701,-0.00017997727258024093],[0.04134319588913915,-0.003195639343879781,-0.017142263870285883,0.001352707046743073]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]   

# line 18576
tempvar = [[-0.002437594346316027,-0.0013945701178424677,-0.0005945840543894248,0.9255589032815941],[1.161791473841111e-05,7.436266370599394e-05,0.0266238866794217,-0.005736600481288336],[-0.010075092209621241,0.04211535836717481,-0.0003536684254310088,-0.012684113883776622],[0.005701851671809172,-0.007751467815603046,-0.0013172435889277955,-0.00021478151538277936]]
tempvar = np.array(tempvar)[::-1,:]
# print(tempvar)
povm += [tempvar]   


sx = np.array([[0,1],[1,0]])
s0 = np.eye(2)

UxQ1 = np.kron(sx,s0)
UxQ2 = np.kron(s0,sx)

SPAM_estimate = np.zeros((4,4))
for i in range(4):
    tempvar = np.trace(povm[i] @ rho)
    # print(tempvar)
    SPAM_estimate[i,0] = tempvar

for i in range(4):
    tempvar = np.trace(povm[i] @ UxQ2 @ rho @ UxQ2.conj().T)
    # print(tempvar)
    SPAM_estimate[i,1] = tempvar

for i in range(4):
    tempvar = np.trace(povm[i] @ UxQ1 @ rho @ UxQ1.conj().T)
    # print(tempvar)
    SPAM_estimate[i,2] = tempvar
    
for i in range(4):
    U = UxQ2 @ UxQ1
    tempvar = np.trace(povm[i] @ U @ rho @ U.conj().T)
    # print(tempvar)
    SPAM_estimate[i,3] = tempvar
print('')
print('H+S mode in GST estimate')
np.set_printoptions(precision=5, suppress=True)
print('two-qubit SPAM')
print(SPAM_estimate)

# [[0.86927 0.05648 0.05843 0.0147 ]
#  [0.0562  0.87286 0.01084 0.05934]
#  [0.05942 0.0107  0.8737  0.05733]
#  [0.0151  0.05997 0.05703 0.86863]]


#q1 SPAM 
SPAM_q1 = np.array([
    [np.sum(SPAM_estimate[0:2,0:2])/2, np.sum(SPAM_estimate[0:2,2:4])/2],
 [np.sum(SPAM_estimate[2:4,0:2])/2, np.sum(SPAM_estimate[2:4,2:4])/2],])
print('')
print('two-qubit SPAM in qubit A subspace')
print(SPAM_q1)
# [[0.9274  0.07165]
#  [0.0726  0.92835]]

#q2 SPAM 
SPAM_q2 = np.array([
    [np.sum(SPAM_estimate[0::2,0::2])/2, np.sum(SPAM_estimate[0::2,1::2])/2],
 [np.sum(SPAM_estimate[1::2,0::2])/2, np.sum(SPAM_estimate[1::2,1::2])/2],])
print('')
print('two-qubit SPAM in qubit B subspace')
print(SPAM_q2)
# [[0.93042 0.0696 ]
#  [0.06958 0.9304 ]]


np.average(np.diag(SPAM_estimate))



