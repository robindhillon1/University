# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 17:37:47 2020

@author: robin
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl

#plt.rc('text',fontsize = 11 ,usetex=True) #In order to use LaTeX
#plt.rc('font', family='serif') #In order to use Serif (mathced font with LaTeX)
 #mpl.style.use('seaborn')

# Some unit conversion factors:
H0 = 70 # Hubble constant H = 70 km/s/Mpc
Mpc = 3.085677581e19 # km
km = 1.0
Gyr = 3.1536e16 # second

H_0 = (H0 * km * Gyr) / Mpc #Here H_0 is in the units of 1/Gyr

eps = 1e-4
def Friedmann(a,m,r,l):
    Omega_0 = m + r + l
    t = (1/ np.sqrt( (r/(a*a)) + (m/a) + (l * (a*a)) + (1-Omega_0) ) ) /H_0
#t_0 = (1/( (m/a) + (lam * (a*a)))**(0.5))/gyr # Matter Lambda Flat Universe Approximation
    return t



def Trapezoidal(a,b,m,r,l):
    n = 10000 # Step number
    deltaX = (b-a)/n # Step size

    AGE = 0

    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)

    for i in range(n):
        x[i] = a + i*deltaX
        y[i] = Friedmann(x[i],m,r,l)
        z[i] = (deltaX/2) * (2*np.sum(y) - y[0] - y[n-1])
        if (x[i] == 1 or 1-eps <= x[i] <= 1+eps):
            AGE = z[i]
            
    print('Age of the universe with m = %5.3f, r = %5.3f, lambda = %5.3f is %5.3f Gyr' %(m,r,l,AGE))
    
    return x,z,AGE

a,t,age = Trapezoidal(1e-10,10,0.27,8e-5,0.73)
a1,t1,age1 = Trapezoidal(1e-10,10,0.5,0.5,0)
a2,t2,age2 = Trapezoidal(1e-10,10,0.7,0.2,0.1)
a3,t3,age3 = Trapezoidal(1e-10,10,0.3,0.7,-0.1)

t_i = 0
t_f = 90
steps = 1000
t_e = np.linspace(t_i,t_f,steps)

a_e,age_e = EmptyUniverse(t_e)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.plot(t-age, a, color='blue', label = '$\Lambda$ CDM')
ax1.plot(t1-age1, a1, color='crimson', label = '$\Omega_{m,0} = 0.5,\Omega_{r,0} = 0.5$')
ax1.plot(t2-age2, a2, color='forestgreen', label = '$\Omega_{m,0} = 0.7,\Omega_{r,0} = 0.2, \Omega_{\Lambda,0} = 0.1 $')
ax1.plot(t_e-age_e,a_e,color = 'black',linestyle = ':',linewidth = 1)

ax1.set_xlim(-20,50)
ax1.set_ylim(0,5)

ax1.set_xlabel('$t$ (Gyr)')
ax1.set_ylabel('$a(t)$')
ax1.legend()

ax1.scatter(0,1,color = 'black',zorder = 10,s = 13)
plt.hlines(1,-60,0,color = 'black',linestyle = ':',linewidth = 0.5)
plt.vlines(0,0,1,color = 'black',linestyle = ':',linewidth = 0.5)