# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 22:53:44 2020

@author: robin
"""

import matplotlib.pyplot as plt
import numpy as np

t = np.linspace(0,50,1000)

R = 0.7

num = R**(3/2)*(1-R**((3/2)*t))
den = 1-R**(3/2)
a = (1 + num/den)**2
b = 1 - R
Pout = a*b

plt.ioff()

plt.plot(t, Pout, label = R)

R = 0.8

num = R**(3/2)*(1-R**((3/2)*t))
den = 1-R**(3/2)
a = (1 + num/den)**2
b = 1 - R
Pout = a*b
plt.plot(t, Pout, label = R)

R = 0.9

num = R**(3/2)*(1-R**((3/2)*t))
den = 1-R**(3/2)
a = (1 + num/den)**2
b = 1 - R
Pout = a*b
plt.plot(t, Pout, label = R)

R = 0.95

num = R**(3/2)*(1-R**((3/2)*t))
den = 1-R**(3/2)
a = (1 + num/den)**2
b = 1 - R
Pout = a*b
plt.plot(t, Pout, label = R)
plt.show()
plt.legend()

plt.xlabel('time(s)')
plt.ylabel('Power')
plt.title('Power vs Time')


