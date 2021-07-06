# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 18:39:07 2019

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 1, 100)

y = 1/(np.exp(1/x)-1)

plt.figure(1)
plt.plot(x, y)
plt.xlabel('1/T')
plt.ylabel('y')
plt.title('U as function of T')
plt.grid()

plt.figure(2)
x2 = np.linspace(-3,3,300)
y2 = (np.exp(x2)-np.exp(-x2))/(1+np.exp(x2)+np.exp(-x2))
y3 = np.tanh(x2)
plt.plot(x2, y2, label='3b')
plt.plot(x2, y3, ':', label='tanh(x)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('M as function of Beta')
plt.grid()
plt.legend()
