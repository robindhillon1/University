# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 16:31:13 2019

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy.integrate import odeint

data = np.loadtxt('test.csv', delimiter=',', comments='#')

time = data[:, 0]
polish = data[:, 2]
rough = data[:, 3]
lacquer = data[:, 4]
unc = data[:,5]


r = 1.285/100
hP = 30.4/100
hR = 27.2/100
hL = 30.5/100

sigma = 5.67e-8

AreaP = 2*np.pi*r*hP
AreaR = 2*np.pi*r*hR
AreaL = 2*np.pi*r*hL

volP = np.pi*r**2*hP
volR = np.pi*r**2*hR
volL = np.pi*r**2*hL

density = 2710
C = 0.9

T0 = 20.4  # from the average of the thermometers

mP = density*volP
mR = density*volR
mL = density*volL

D = 2.57/100

timeD = np.arange(time[0], time[-1], 0.1)

def Temp_ODE(time, Temp, T0, mass, A, D, e, C):
    def model_dTdt(Temp, t):
        h = 1.32*((Temp-T0)/D)**(1/4)
        result = (1/(mass*C)) * ((A*e*sigma*(T0**4 - Temp**4))-(h*A*(Temp - T0)))
        return result
    
    TempResult = odeint(model_dTdt, polish[0], timeD) # can add h0
    return TempResult

# =============================================================================
# 
# def fxn(time, Temps, mass, A, D, e, C):
# 	eFit, eCov = curve_fit(Temp_ODE(time, T0, mass, A, D, e, C), time, Temps[0], p0=[e], maxfev = 10000, bounds = (0,1))
# 	print(eFit)
# =============================================================================

e = 1
a = Temp_ODE(timeD, polish[0], T0, mP, AreaP, D, e, C)

plt.plot(timeD, a[:,0], label='Fit')
plt.plot(time, polish, 'ko', markersize = '1.5', label='Polish')
plt.legend()
plt.grid()

# =============================================================================
# pFit = fxn(time, polish, mP, AreaP, D, 0.2, C)
# rFit = fxn(time, rough, mR, AreaR, D, 0.3, C)
# lFit = fxn(time, lacquer, mL, AreaL, D, 0.4, C)
# =============================================================================

# =============================================================================
# plt.plot(time, polish, 'ko', markersize = '1.5', label='Polish')
# plt.plot(time, rough, 'ro', markersize = '1.5', label='Rough')
# plt.plot(time,lacquer, 'yo', markersize = '1.5', label='Lacquer')
# plt.errorbar(time, polish, yerr = unc, c = 'k', marker= '.', ms = '1', linestyle = '', barsabove = True)
# plt.errorbar(time, rough, yerr = unc, c = 'r', marker= '.', ms = '1', linestyle = '', barsabove = True)
# plt.errorbar(time, lacquer, yerr = unc, c = 'y', marker= '.', ms = '1', linestyle = '', barsabove = True)
# plt.xlabel('Time(s)')
# plt.ylabel('Temperature(C)')
# plt.title('Temperature as a function of time')
# plt.legend()
# plt.grid()
# plt.show()
# 
# 
# =============================================================================
