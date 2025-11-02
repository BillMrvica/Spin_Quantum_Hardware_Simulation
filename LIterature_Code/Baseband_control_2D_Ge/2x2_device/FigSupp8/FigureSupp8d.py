# this file is modified from Supp/unequal_time/bloch_slider_without_qutip_v3_q2x90.py
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np
import scipy.linalg
sx = np.array([[0,1],[1,0]],dtype=complex)
sy = np.array([[0,-1j],[1j,0]],dtype=complex)
sz = np.array([[1,0],[0,-1]],dtype=complex)
sm = np.array([[0,0],[1,0]],dtype=complex)




#%%
DEG = np.pi/180
fq2 = 0.070926*5/4   
fq3     = 0.06203767*5/4      # qubit frequency


theta = 44.7 * DEG       # q3 angle from sigma_z axis (toward sigma_x axis)
phi = 0 * DEG 

# initial state
psi00 = np.array([[1],[0]], dtype=complex) 
# psi00 = np.array([[0],[1]], dtype=complex)  


#####
outputs = [psi00]
tlists = []
hlists = []


dt = 0.01
H2 = (fq2/2) *  sz
dU2 = scipy.linalg.expm(  -2*np.pi*1j* H2 * dt  )
H3 = (fq3/2) * (np.cos(theta)*sz + np.sin(theta)*np.cos(phi)*sx  + np.sin(theta)*np.sin(phi)* sy)
dU3 = scipy.linalg.expm(  -2*np.pi*1j* H3 * dt  )

output = np.array(psi00)
  
wait_q2 = 3.42 + 2.288
wait_q3 = 4.86 + 1.542 


tlist = np.arange(0, wait_q3  , dt)
for j in range(len(tlist)):
    output = dU3 @ output
    outputs += [output    ]
tlists = np.hstack((tlists,tlist))  
hlists = hlists + [[np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]]*len(tlist)


tlist = np.arange( wait_q3 , wait_q3+wait_q2  , dt)
for j in range(len(tlist)):
    output = dU2 @ output
    outputs += [output]
tlists = np.hstack((tlists,tlist)) 
hlists = hlists + [[0,0,1]]*len(tlist)



hlists_plot = hlists[0::10]
tlists_plot = tlists[0::10]
outputs_plot = outputs[0::10]
#####

from my_bloch import Bloch
def expect(M, ss):
    if isinstance(ss, list):
        slist = []
        for ssx in ss:
            s = ssx.reshape(-1)
            slist += [  np.real(np.dot(s.conj(), M @ s))  ]
        return slist
    else:
        s = ss.reshape(-1)
        return np.real(np.dot(s.conj(), M @ s))
def state2vec(s):
    return expect(sx, s), expect(sy, s), expect(sz, s)

vx, vy, vz = expect(sx, psi00), expect(sy, psi00), expect(sz, psi00)
sphere_interactive = None

if sphere_interactive==None:
  
    # fig, axs = plt.subplots( figsize=(7,9))
    pltname = 'Interactive plot with slider'
    fig = plt.figure( pltname, figsize = (7,9))
    sphere_interactive=Bloch(fig)
else:
    sphere_interactive.clear()    
colorlist = list(plt.cm.Oranges(np.linspace(0.9,0,301)))
sphere_interactive.vector_color = ['g', colorlist[0]] 
sphere_interactive.point_color = colorlist
sphere_interactive.add_vectors(hlists_plot[0])
sphere_interactive.add_states( state2vec(psi00)  )
sphere_interactive.render()

#####


axcolor = 'lightgoldenrodyellow'
ax_t = plt.axes([0.1, 0.05, 0.8, 0.03], facecolor=axcolor)
s_t = Slider(ax_t, 't (ns)', tlists_plot[0], tlists_plot[-1], valinit=tlists_plot[0], valstep=0.02, color='C1')
s_t.label.set_size(15)


def update(val):
     
    tau = s_t.val
    idx = np.argmin(abs(tlists_plot-tau))    
    sphere_interactive.clear()

    sphere_interactive.add_vectors(hlists_plot[idx])
    sphere_interactive.add_states(  state2vec(outputs_plot[idx])  )
    
    if idx > 0:
        k = max(idx-300, 0)
        vx, vy, vz = expect(sx, outputs_plot[idx:k:-1]), expect(sy, outputs_plot[idx:k:-1]), expect(sz, outputs_plot[idx:k:-1]) 
        sphere_interactive.add_points([vx,vy,vz], 'm')
        # print(vx)
     
    sphere_interactive.render()

s_t.on_changed(update)
plt.show()


#%%
N_space = 30
hlists_plot = hlists[0::N_space]
tlists_plot = tlists[0::N_space]
outputs_plot = outputs[0::N_space]
vx, vy, vz = expect(sx, outputs_plot[:]), expect(sy, outputs_plot[:]), expect(sz, outputs_plot[:]) 

sphere_static = None
if sphere_static==None:
    pltname = f'Supp_q2x90_unequal'
    fig = plt.figure( pltname, figsize = (7,9))
    sphere_static=Bloch(fig)
else:
    sphere_static.clear()


sphere_static.vector_color = ['C0', 'C1']
sphere_static.add_vectors(hlists_plot[0])
sphere_static.add_vectors([0, 0, 1])

temp_var = np.array(hlists_plot)[:,0]<0.01
switch_pts = np.where(np.diff(  temp_var  ))[0]

colorlist = ['C0',]*(switch_pts[0]+1) \
            + ['C1',]*(len(temp_var)-switch_pts[0]-1) 
sphere_static.point_color = colorlist
sphere_static.add_points([vx,vy,vz], 'm')
sphere_static.render()

#%%
CQ1 = '#0D77BC'
CQ2 = '#DD5E27'


    
t_x = np.cumsum(np.array([0, 11.279361588134112,
2,
4.86,
2,
3.42,
11.279361588134112   ] ))   
    
    
y_x = np.array([ 0, 0, 
              1,
              1, 
              0, 
              0, 
              0  ])     


t = 0

plt.figure(figsize=(5,3))
plt.plot(  t + t_x , y_x  , 'k' )
# plt.plot(  t + t_x , y_x  , CQ2 )
t = ( t + t_x)[-1]     

plt.title('Qubit B unequal time gate voltage pulses')
plt.xlabel('time (ns)')
plt.ylabel('voltage pulse (a.u.)')





