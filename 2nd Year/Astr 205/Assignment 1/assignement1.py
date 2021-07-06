# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:44:42 2018

@author: robin
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

fname = 'assign1_2.csv'

#load the data from the file
data = np.loadtxt(fname, delimiter=',', comments='#')

# access the data
x = data

n, bins, patches = plt.hist(x, 100, normed = 1, facecolor = 'g', alpha=0.75)

guesses = (0.6, 0.4, 2.8701, 2.8701, 0, 0)
x_plot_initial = bins
y_plot = n
x_plot_0 = np.roll(x_plot_initial, -1)
x1_0 = 0.5*(x_plot_initial+x_plot_0)
print(len(x_plot_initial))
print(len(x1_0))
x_plot = x1_0[0: 100]
print(x_plot)
print(y_plot)

def fit_function(x, C_1, C_2, sigma_1, sigma_2, miu_1, miu_2):
    return  np.abs(C_1)/np.sqrt(2*math.pi*sigma_1**2)*np.exp(-(x-miu_1)**2/(2*sigma_1**2))+np.abs(C_2)/np.sqrt(2*math.pi*sigma_2**2)*np.exp(-(x-miu_2)**2/(2*sigma_2**2))

popt, pcov = curve_fit(fit_function, x_plot, y_plot, p0=guesses)

plt.plot(x_plot, fit_function(x_plot, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]), 'r-', label='fit')

print(popt)

plt.show()