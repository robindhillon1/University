# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 15:37:07 2020

@author: robin
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

H0 = (68/3.086E19)*(31536000*1E9)  # in 1/Gyr


def Friedmann(a, omega_m, omega_rad, omega_lamb):
    den = np.sqrt(omega_m/a + a**2*omega_lamb + omega_rad/(a*a) + (1-omega_lamb-omega_m-omega_rad))
    t = 1/(H0*den)
    return t


def bigbounce(a):
    omega_lamb = 1.8
    omega_m = 0.31
    den = np.sqrt(omega_m/a + a**2*omega_lamb + (1-omega_lamb-omega_m))
    t = 1/(H0*den)
    return t


def bigcrunch(a):
    omega_lamb = -0.31
    omega_m = 0.31
    den = np.sqrt(omega_m/a + a**2*omega_lamb + (1-omega_lamb-omega_m))
    t = 1/(H0*den)
    return t


aLim = float(input('Choose a value for the scale factor, a: '))
aBounce = np.linspace(0.552, aLim, 20001)
timeBounce = []

for i in range(len(aBounce)-1):
    time = integrate.quad(Friedmann, aBounce[0], aBounce[i+1], args=(0.31, 0, 1.8))
    timeBounce.append(time[0])

tPlot = [i-timeBounce[-1] for i in timeBounce]
tPlotN = [-1*i-timeBounce[-1] for i in timeBounce]

plt.plot(tPlot, aBounce[1:], 'k')
plt.plot(tPlotN, aBounce[1:], 'k')

aCrunch = np.linspace(0, aLim, 10001)
timeCrunch = []

for i in range(len(aCrunch)-1):
    time = integrate.quad(Friedmann, aCrunch[0], aCrunch[i+1], args=(0.31, 0, -0.31))
    timeCrunch.append(time[0])

tCrunch = [i-max(timeCrunch) for i in timeCrunch]
tCrunchN = [-1*i+max(timeCrunch) for i in timeCrunch]

plt.plot(tCrunch, aCrunch[1:])
plt.plot(tCrunchN, aCrunch[1:])
