# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 13:50:51 2020

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc

x = np.linspace(-1,5,1000)
a=0.15  #Scale 
v1 = x*np.sqrt(a)
v2 = np.inf

S = sc.fresnel(v2)[0] + sc.fresnel(v1)[0]
C = sc.fresnel(v2)[1] + sc.fresnel(v1)[1]

I1 = (S**2+C**2)

plt.figure(1)
plt.plot(x,I1/np.max(I1))
plt.xlabel('x')
plt.ylabel('Normalized Intensity')
plt.title('Far Field Fresnel Diffraction of Knife Edge')

D = 5e-4
x = np.linspace(-4e-4,4e-4,10000)
lamb = 5e-7
z = 0.025
v1 = -1*np.sqrt(2/(z*lamb))*(D/2 + x)
v2 = np.sqrt(2/(z*lamb))*(D/2 - x)

S = sc.fresnel(v2)[0] - sc.fresnel(v1)[0]
C = sc.fresnel(v2)[1] - sc.fresnel(v1)[1]

I = (S**2+C**2)

plt.figure(2)
plt.plot(x*(1e3),I/np.mean(I))
plt.xlabel('x(mm)')
plt.ylabel('Normalized intensity on axis')
plt.title('Fresnel Diffraction single slit, z=0.025m')