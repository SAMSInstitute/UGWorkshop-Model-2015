#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Experimental system, based on Becks4Eq

Created on Thu Apr 30 11:34:27 2015

:author: Christopher Strickland
"""

from __future__ import division
import numpy as np
from scipy.integrate import ode
import scipy as sp
import matplotlib.pyplot as plt

##### Variables go here #####
N0 = 4.5e-5

D = 1.10
STOC_D = False #turn on or off noisy D (Gamma noise)
var_D = 0.5

OBS_NOISE = False #turn on or off observation noise


#time points to solve at
tpts = np.arange(0,1001)
#initial values
x0 = np.array([2.1e-7,4.56e-7,4.3e-5,N0])

#growth rates (1/sec)
#values below are in 1/day; divide by 24*3600 to get seconds
munr = 12.0
munc = 6.0
mupr =2.2
mupc = 2.2/8.0

#mass (g)
mr = 1.6e-12
mc = 8.2e-12
mp = 9.8e-10

#half-saturation constants (g/cc)
knr = 8.0e-6
knc = 8.0e-6
kpr = 1.0e-6
kpc = 1.0e-6

#death rates (1/sec)
#these values are in 1/day - divide by 24*3600 to get seconds
dr = 0.5/10.0
dc = 0.25/100.0
dp = 0.08/100.0

#yield coefficients
ypr = 0.12
ypc = 0.12
ynr = 0.1
ync = 0.1

#ODE function
def becks4Eq(t,x,N0,D):
    dx = np.zeros(4)
    
    #ODEs
    dx[0] = x[0]*(munr*x[3]/(x[3]+knr)-dr) - \
    mupr/ypr*mp/mr*x[0]/(kpr/mr+x[0])*x[2] - D*x[0]
    
    dx[1] = x[1]*(munc*x[3]/(x[3]+knc)-dc) - \
    mupc/ypc*mp/mc*x[1]/(kpc/mc+x[1])*x[2] - D*x[1]
    
    dx[2] = x[2]*(mupr*x[0]/(kpr/mr+x[0]) + \
    mupc*x[1]/(kpc/mc+x[1])-dp) - D*x[2]
    
    dx[3] = D*N0 - x[0]*munr*mr/ynr*x[3]/(knr+x[3]) - \
    x[1]*munc*mc/ync*x[3]/(knc+x[3]) - D*x[3]
    
    return dx
    
##### Solve procedure goes here #####
Rsol = []; Csol = []; Psol = []; Nsol = []
Rsol.append(x0[0])
Csol.append(x0[1])
Psol.append(x0[2])
Nsol.append(x0[3])
r = ode(becks4Eq).set_integrator('dopri5',nsteps=100000,verbosity=1)
r.set_initial_value(x0,0).set_f_params(N0,D)
for t in tpts[1:]:
    if STOC_D:
        D_t = sp.random.gamma(D**2/var_D,var_D/D)
    else:
        D_t = D
    r.set_f_params(N0,D_t)
    r.integrate(t)
    assert(r.successful())
    Rsol.append(r.y[0])
    Csol.append(r.y[1])
    Psol.append(r.y[2])
    Nsol.append(r.y[3])
    
##### Plot solution #####
plt.figure(1)
plt.subplot(221)
plt.plot(tpts,Rsol)
plt.title(r"Plot of $R$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$R$")
plt.subplot(222)
plt.plot(tpts,Csol)
plt.title(r"Plot of $C$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$C$")
plt.subplot(223)
plt.plot(tpts,Psol)
plt.title(r"Plot of $P$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$P$")
plt.subplot(224)
plt.plot(tpts,Nsol)
plt.title(r"Plot of $N$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$N$")
plt.show()