# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:53:19 2021

@author: davem
"""


from scipy import integrate
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
#These import modules that I use a lot. You will pretty much ALWAYS use numpy and matploylib.pyplot
#The other ones are because I'm going to be doing a lot of differential equation solving, so I
#imported integrate. I probably don't need interpolate.


#%% This is a function based on experimental data that models the release of ATRA from PLGA-ATRA into synovium
def MP(t):
    y = 5.061*t/(2.003+t) + 0.004093*t - 0.03464
    
    return y

#Function to define the rate of ATRA release from particles. Fit is from in vitro release data, divisions are to get into moles/hr
# The final division in the function is to account for the fact that the release assay was performed with 10 mg of particles in 1 mL, while we
# Inject 2ug of particles
def k3dMP(t):
    dy = ((0.004093*(t/24)**2 + 0.0163966*(t/24) + 10.1536)/(((t/24) + 2.003)**2))/1000000/300/24*2/1000
    
    return dy

def k1dMP(t):
    dy = ((3.334*0.9215)/(((t/24) + 0.9215)**2) + 0.003773)/1000000/300/24*2/1000
    
    return dy

def k6dMP(t):
    dy = ((6.544*2.076)/(((t/24) + 2.076)**2) - 0.0212)/1000000/300/24*2/1000
    
    return dy
    

release = {'Rls1': k1dMP, 'Rls3': k3dMP, 'Rls6': k6dMP}
#%% This is the function that defines the mathematical model we are using to model the two compartments (peripheral blood and ankle synovium)
def TwoComp(Cs, t0, params, release):
    #define the function that is our model to produce derivative for integration
    [ksp, kps, ke, z, v2] = params
#    import the paramteres
    #Diff eq for peripheral blood concentration
    dC_PB = ksp*Cs[1]*z - ke*Cs[0] - kps*Cs[0]
   
    #Diff eq for synovium concentration
    dC_S = kps*Cs[0]/z-ksp*Cs[1] + release['Rls1'](t0)/v2
    
    #Diff eq for Random peripheral Organ Concentration
    dC_O = 0*Cs[0]/(z*5) - 0*Cs[2]
     
#    import the starting values and assign them to variable that make more sense
    dCs = [dC_PB, dC_S, dC_O]
    return dCs

#define the parameter values

ksp = 0.2
kps = 0.2
ke = 1.5
z = 0.001
v2 = 0.00002

params = [ksp, kps, ke, z, v2]
#Initial Values
C_PB0 = 0
C_S0 = 0
C_O0 = 0

Cs0 = [C_PB0, C_S0, C_O0]

ts = np.arange(0,504,0.1)
LL = np.zeros(len(ts))+100*10**-12

solve = integrate.odeint(TwoComp, Cs0, ts, args=(params,release,))

for i in range(len(solve[:,0])):
    if solve[i,0] < 0:
        solve[i,0] = 0
    if solve[i,1] < 0:
        solve[i,1] = 0

rls = k3dMP(ts)

font = {'family' : 'sans-serif',
        'sans-serif' : 'Arial',
        'size' : 26}
plt.rc('font', **font)

kernel_size = 20
kernel = np.ones(kernel_size)/kernel_size


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(ts/24,solve[:,0]*10**9,'b-',ts/24,LL*10**9,'b-.')
ax2.plot(ts/24,solve[:,1]*10**9,'r-')


ax1.set_xlabel('Time (Days)')
ax1.set_ylabel('Peripheral Blood (nM)',color = 'b')
ax2.set_ylabel('Synovium (nM)', color ='r')

fig, ax3 = plt.subplots()
ax3.plot(ts/24, solve[:,0]*10**9,'b-',ts/24, LL*10**9, 'k-.')
ax3.set_xlabel('Time (Days)')
ax3.set_ylabel('Peripheral Blood (nM)',color = 'b')
ax3.set_yscale('log')

fig, ax4 = plt.subplots()
ax4.plot(ts/24,solve[:,1]*10**9,'r-', ts/24, LL*10**9, 'k-.')
ax4.set_xlabel('Time (Days)')
ax4.set_ylabel('Synovium (nM)', color ='r')
ax4.set_yscale('log')

#plt.plot(ts/24,rls)

plt.show()


