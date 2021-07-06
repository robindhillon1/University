
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 22:37:34 2018

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

G = 1
Ma = 1

def v(r):
    v = np.sqrt(r**2/(r**2+1)**(3/2))
    return v

r = np.linspace(0,10,1000000)
plt.plot(r,v(r),label='L=$cr_s$')
plt.xlabel('R')
plt.ylabel('V(R)')
plt.title('Rotation Curve')
plt.grid()
plt.legend()
plt.show()
