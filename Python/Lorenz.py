#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lorenz system

Created on Thu May 21 13:23:43 2015

:author: Christopher Strickland
"""

from __future__ import division
import numpy as np
import scipy as sp
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

##### Parameters go here #####

#time points to solve at
tpts = np.linspace(0,50,5001) #100,110,... for no transients
#initial values
x0 = np.array([10.,10.,10.])

####################
sigma = 10. #0.5 #10.
beta = 8./3. #1 #8./3.
rho = 0.5 #1. #14. #28.

####################

#Observation noise
OBS_NOISE = False #turn on or off observation noise
#Assume noise is gaussian
obs_mu = np.array([0,0,0]) #x,y,and z variables
obs_sig2 = np.array([.05,.05,.05])

##### ODE function #####
def LorenzODEs(t,x):
    dx = np.zeros(3)
    
    dx[0] = sigma*(x[1]-x[0])
    dx[1] = x[0]*(rho-x[2])-x[1]
    dx[2] = x[0]*x[1] - beta*x[2]
    
    return dx
    
##### Solve procedure goes here #####
Xsol = []; Ysol = []; Zsol = []
r = ode(LorenzODEs).set_integrator('dopri5',nsteps=100000,verbosity=1)
r.set_initial_value(x0,0)
for t in tpts:
    if t == 0:
        Xsol.append(x0[0]);Ysol.append(x0[1]);Zsol.append(x0[2])
        continue
    r.integrate(t)
    assert(r.successful())
    if OBS_NOISE:
        Xsol.append(max(r.y[0] + sp.random.normal(obs_mu[0],obs_sig2[0]),0))
        Ysol.append(max(r.y[1] + sp.random.normal(obs_mu[1],obs_sig2[1]),0))
        Zsol.append(max(r.y[2] + sp.random.normal(obs_mu[2],obs_sig2[2]),0))
    else:
        Xsol.append(r.y[0])
        Ysol.append(r.y[1])
        Zsol.append(r.y[2])
        
##### Plot solution #####
fig = plt.figure()
plt.subplot(221)
plt.plot(tpts,Xsol)
plt.title(r"Plot of $X$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$X$")
plt.ylim(-20,20)
plt.subplot(222)
plt.plot(tpts,Ysol)
plt.title(r"Plot of $Y$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$Y$")
plt.ylim(-30,30)
plt.subplot(223)
plt.plot(tpts,Zsol)
plt.title(r"Plot of $Z$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$Z$")
plt.ylim(0,50)
#attractor
ax = fig.add_subplot(2,2,4,projection='3d')
ax.plot(Xsol, Ysol, Zsol)
ax.set_title("Attractor")
ax.set_xlabel(r"X")
ax.set_ylabel(r"Y")
ax.set_zlabel(r"Z")
plt.show()