# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 09:57:42 2018

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

def y(x):
    y = np.sqrt(1-(1/x)*np.arctan(x))
    return y

x = np.linspace(0,10,2000)
plt.plot(x,y(x))
plt.ylim(0,1)
plt.xlim(0,9)
plt.xlabel('r/r_0')
plt.ylabel('V/V_H')
plt.title('Relationship between V/V_H and r/r_0')
plt.grid()
plt.legend()
plt.show()