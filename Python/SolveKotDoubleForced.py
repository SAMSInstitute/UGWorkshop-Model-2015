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
from scipy import stats
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

##### Parameters go here #####

#time points to solve at
tpts = np.linspace(0,500,4001)
#initial values
x0 = np.array([1,1,1])

####################

#Observation noise
OBS_NOISE = False #turn on or off observation noise
#Assume noise is gaussian
obs_mu = np.array([0,0,0]) #x,y,and z variables
obs_sig2 = np.array([.01,.01,.01])

####################

#Dilution rate
D = 0.15
#Most of Kot's stuff is based on D=0.1
#If you set epsilon=0, there is no forcing. In this case,
#   D=0.15 gives nice non-extinction steady states

####################

#For D=0.1:
#Epsilon = 0.6 gives chaos. >0.6 is really nice
#0.4 is cool. Several interacting oscillations.
#0.3 is two oscillations
#0.1 and 0.0 - limit cycle
epsilon = 0.0
STOC_EPS = False #turn on and off stochastic epsilon
if STOC_EPS:
    #create the random variable you want for epsilon here. E.g.
    #stats.truncnorm(min,max,mu,sig2)
    #stats.gamma(k,loc=0,scale=theta)
    eps_rv = stats.truncnorm(0,3,epsilon,0.15)

####################
#Note that this section does nothing if epsilon=0
#T = 100.0
T = 24

#chaotic dynamics, Fig. 6 when epsilon = 0.6
#omega = 5.0*np.pi/6.0 #T is not used here. This is equivalent to T=24.
omega = 2.*np.pi/D/T #nifty limit cycle w/ epsilon = 0.6, T=100!
                     #epsilon = 0.1 gives wave envelopes, T=100!
#omega = 4.0*np.pi

STOC_T = False #turn on and off stochastic forcing period
if STOC_T:
    #this function will replace the omega definition above
    def omega_f(T_arg):
        return 2.*np.pi/D/T_arg
    #create the random variable you want for T here
    #remember to require T>0
    T_rv = stats.truncnorm(0,150,50,10)

####################
#Other parameters

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

##### ODE function #####
def KotODEs(t,x,epsilon,omgea):
    dx = np.zeros(3)
    
    dx[0] = 1 + epsilon*np.sin(omega*t) - x[0] - A*x[0]*x[1]/(a+x[0])
    dx[1] = A*x[0]*x[1]/(a+x[0]) - x[1] - B*x[1]*x[2]/(b+x[1])
    dx[2] = B*x[1]*x[2]/(b+x[1]) - x[2]
    
    return dx
    
##### Solve procedure goes here #####
Xsol = []; Ysol = []; Zsol = []
r = ode(KotODEs).set_integrator('dopri5',nsteps=100000,verbosity=1)
r.set_initial_value(x0,0)
for t in tpts:
    if t == 0:
        Xsol.append(x0[0]);Ysol.append(x0[1]);Zsol.append(x0[2])
        continue
    if STOC_EPS:
        epsilon = eps_rv.rvs()
    if STOC_T:
        omega = omega_f(T_rv.rvs())
    r.set_f_params(epsilon,omega)
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

##### Output data #####
#uncomment for data generation
#np.savetxt('data.csv', np.asarray([Xsol, Ysol, Zsol]), delimiter=",")
        
##### Plot solution #####
fig = plt.figure()
plt.subplot(221)
plt.plot(tpts,Xsol)
plt.title(r"Plot of $X$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$R$")
plt.ylim(0,1.5)
plt.subplot(222)
plt.plot(tpts,Ysol)
plt.title(r"Plot of $Y$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$C$")
plt.ylim(0,1.5)
plt.subplot(223)
plt.plot(tpts,Zsol)
plt.title(r"Plot of $Z$ vs. time")
plt.xlabel(r"$t$")
plt.ylabel(r"$P$")
plt.ylim(0,1.5)
#attractor
ax = fig.add_subplot(2,2,4,projection='3d')
ax.plot(Xsol, Ysol, Zsol)
ax.set_title("Attractor")
ax.set_xlabel(r"X")
ax.set_ylabel(r"Y")
ax.set_zlabel(r"Z")
plt.show()