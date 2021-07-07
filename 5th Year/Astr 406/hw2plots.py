# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 15:57:12 2020

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt


def Fv(h, v, k, T):
    F = np.exp(-(h*v)/(k*T))
    return F


v = np.linspace(4.25E9, 10**17, 5000)
h = 6.63E-27
k = 1.38E-16
T = 8000

call = Fv(h, v, k, T)
print(call)
#plt.loglog(v, call)


def Fx2(v):
    F2 = v**2
    return F2


v = np.linspace(0,4.25E9,5000)

call2 = Fx2(v)
#plt.loglog(v,call2)



v2 = np.linspace(4.25E9,1.7E15,10000)
Tb = (1.81E19)*(v2**-2)*8000

plt.loglog(v2,Tb)
