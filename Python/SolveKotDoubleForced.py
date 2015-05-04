#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Solve the Kot Forced Double-Monod Equations
-- based on the dimensionless equations from the Kot paper

Created on Mon May 04 10:01:39 2015

:author: Christopher Strickland
"""

from __future__ import division
import numpy as np
import scipy as sp
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

##### Variables go here #####

#time points to solve at
tpts = np.linspace(0,500,4001)
#initial values
x0 = np.array([1,1,1])

#Dilution rate
D = 0.1

#Si, mg/l
si = 115.0

#Maximum specific growth rate of prey and predator (1/hr)
mu1 = 0.5 #prey
mu2 = 0.2 #predator

#yield of prey per unit mass of substrate (dimensionless)
y1 = 0.4

#biomass yield of predator per unit mass of prey (dimensionsless)
y2 = 0.6

#half-saturation (Michaelis-Menten) constants for prey and predator
k1 = 8.0
k2 = 9.0

#Constants in dimensionless equations
A = mu1/D
a = k1/si
B = mu2/D
b = k2/y1/si

#Epsilon = 0.6 gives chaos
epsilon = 0.6

T = 100.0

#chaotic dynamics
omega = 5.0*np.pi/6.0
#other choices?
#omega = 2.*np.pi/D/T
#omega = 4.0*np.pi

##### ODE function #####
def KotODEs(t,x):
    dx = np.zeros(3)
    
    dx[0] = 1 + epsilon*np.sin(omega*t) - x[0] - A*x[0]*x[1]/(a+x[0])
    dx[1] = A*x[0]*x[1]/(a+x[0]) - x[1] - B*x[1]*x[2]/(b+x[1])
    dx[2] = B*x[1]*x[2]/(b+x[1]) - x[2]
    
    return dx
    
##### Solve procedure goes here #####
Xsol = []; Ysol = []; Zsol = []
Xsol.append(x0[0])
Ysol.append(x0[1])
Zsol.append(x0[2])
r = ode(KotODEs).set_integrator('dopri5',nsteps=100000,verbosity=1)
r.set_initial_value(x0,0)
for t in tpts[1:]:
#    if STOC_D:
#        D_t = sp.random.gamma(D**2/var_D,var_D/D)
#    else:
#        D_t = D
    r.integrate(t)
    assert(r.successful())
    Xsol.append(r.y[0])
    Ysol.append(r.y[1])
    Zsol.append(r.y[2])
    
##### Plot solution #####
fig = plt.figure()
plt.subplot(221)
plt.plot(tpts,Xsol)
plt.title(r"Plot of $X$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$R$")
plt.subplot(222)
plt.plot(tpts,Ysol)
plt.title(r"Plot of $Y$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$C$")
plt.subplot(223)
plt.plot(tpts,Zsol)
plt.title(r"Plot of $Z$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$P$")
#attractor
ax = fig.add_subplot(2,2,4,projection='3d')
ax.plot(Xsol, Ysol, Zsol)
ax.set_title("Attractor")
ax.set_xlabel(r"X")
ax.set_ylabel(r"Y")
ax.set_zlabel(r"Z")
plt.show()