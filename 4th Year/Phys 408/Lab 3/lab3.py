# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 20:20:13 2020

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

fname = 'D:/Uni/4th Year/Phys 408/2020/Lab 3/1umSpeedDataCalibration.txt'
data = np.loadtxt(fname)

data = data - np.mean(data)
N = len(data)
spec = abs(np.fft.fft(data))

spec = spec[2:N/2+1]
spec = spec/max(spec)
                     
x = np.arange(0, len(data), 1)


plt.plot(x, spec)
plt.xlabel('Frames')
plt.ylabel('Intensity')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.title('Fourier Transform - HeNe')
plt.xlim(100,400)
