# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 18:39:07 2019

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(255, 273, 200)
xx = np.linspace(273, 294, 200)

y = (334-(2.1)*(x-273)-(4.2*(273-255)))/(334+(2.1)*(x-273))
yy = (334-(2.1)*(xx-273)-(4.2*(273-255)))/(334+(2.1)*(xx-273))


delS = 4.2*np.log(273/255)-(334*(1-y))/273 + y*(4.2)*np.log(x/273) + (1-y)*(2.1)*np.log(x/273)

delSS = 4.2*np.log(273/255)-(334*(1-yy))/273 + yy*(4.2)*np.log(xx/273) + (1-yy)*(2.1)*np.log(xx/273)


plt.plot(x, delS)
plt.plot(xx, delSS, ':')
plt.xlabel('T2 (K)')
plt.ylabel('$\Delta$S')
plt.title('Change in Entropy of the system (J/k)')
plt.grid()