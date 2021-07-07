# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 12:59:03 2020

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

lamb = 1*10**(-6)
rho = np.linspace(-1*2*10**-3,2*10**-3,7000)
Rp = np.sqrt(rho**2+0.1**2)
y = np.cos((2*np.pi*(Rp-0.1))/lamb)

plt.plot(rho*1000, y)

plt.xlabel('Distance from axis, rho (mm)')
plt.ylabel('y')
plt.title('cos(2*pi*(R(rhp)-R(0))/lambda)')