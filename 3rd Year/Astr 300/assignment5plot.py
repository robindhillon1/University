# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 10:30:48 2018

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

def y(x):
    y = (0.5)*(1-(1/x))*(1+(1/x**2)) #Let c=1. 
    return y

def y2(x):
    y = (0.5)*(1-(1/x))*(1+16*(1/(x**2)))
    return y

x = np.linspace(0,10,1000)
plt.plot(x,y(x),label='L=$cr_s$')
plt.ylim(-5,5)
plt.xlim(0,10)
plt.plot(x,y2(x),label='L=$4cr_s$')
plt.xlabel('$r/r_s$')
plt.ylabel('$\Phi_{eff}$')
plt.title('Effective Potential')
plt.grid()
plt.legend()
plt.show()