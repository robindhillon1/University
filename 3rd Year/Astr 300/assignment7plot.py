# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 10:53:39 2018

@author: robin
"""
import numpy as np
import matplotlib.pyplot as plt

G = 1
Ma = 1

def func(t):
    x = 2*np.cos(t)
    y = np.cos((np.pi/2)*t)
    return x,y

t = np.arange(-10,10.1,0.1)
plt.plot(func(t)[0],func(t)[1])#,'ko',label='Initial position')
plt.plot(func(t)[0][0],func(t)[1][0],'ko')
plt.plot(func(0)[0],func(0)[1],'*',label='At t=0')
plt.xlabel('xPosition')
plt.ylabel('yPosition')
plt.title('Path of star from -10s<=t<=10s with interval of 0.1')
plt.legend()
plt.grid()
plt.show()
