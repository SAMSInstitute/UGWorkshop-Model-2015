# -*- coding: utf-8 -*-
"""
Solve Beck's system of ODEs
Created on Wed Apr 29 15:12:24 2015

@author: Christopher Strickland
"""
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time
from __future__ import division

##### Variables go here #####
N0 = 4.5e-5
tpts = np.concatenate((np.array([0]),np.linspace(7000,7500,501)))
x0 = np.array([2.1e-7,4.56e-7,4.3e-5,4.5e-5])

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
def becks4Eq(x,t,N0,D):
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
for D in xrange(0.05,2.01,0.01):
    Usol = sp.integrate.odeint(becks4Eq,x0,tpts,(N0,D))
    plt.figure(1)
    plt.rc('text', usetex=True)
    plt.hold(True)
    plt.plot(D*np.ones(len(tpts)-1),Usol[1:,0],marker='.',markersize=2)
    plt.title(r"Bifurcation Diagram for $R$ vs $D$")
    plt.xlabel(r"Values for $D$")
    plt.ylabel(r"Rod Species $R$")
    plt.draw()
    
    plt.figure(2)
    plt.rc('text', usetex=True)
    plt.hold(True)
    plt.plot(D*np.ones(len(tpts)-1),Usol[1:,1],marker='.',markersize=2)
    plt.title(r"Bifurcation Diagram for $C$ vs $D$")
    plt.xlabel(r"Values for $D$")
    plt.ylabel(r"Cocci Species $C$")
    plt.draw()
    
    print(D)
    plt.sleep(0.05)
    
    