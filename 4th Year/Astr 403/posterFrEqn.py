# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 14:00:43 2020

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sys

H0 = (68/3.086E19)*(31536000*1E9)  # in 1/Gyr


def Friedmann(a, omega_m, omega_rad, omega_lamb):
    den = np.sqrt(omega_m/a + a**2*omega_lamb + omega_rad/(a*a) + (1-omega_lamb-omega_m-omega_rad))
    t = 1/(H0*den)
    return t


# =============================================================================
# def EmptyUniverse(t):
#     a_e = H0*t
#     age_e = 1/H0
#     return a_e,age_e
# =============================================================================


aLim = float(input('Choose a value for the scale factor, a: '))

aBounce = np.linspace(0.552, aLim, 20001)
timeBounce = []

for i in range(len(aBounce)-1):
    time = integrate.quad(Friedmann, aBounce[0], aBounce[i+1], args=(0.31, 0, 1.8))
    timeBounce.append(time[0])

tPlot = [i-timeBounce[-1]+16.389 for i in timeBounce]
tPlotN = [-1*i-timeBounce[-1]+16.389 for i in timeBounce]

aCrunch = np.linspace(0, aLim, 10001)
timeCrunch = []

for i in range(len(aCrunch)-1):
    time = integrate.quad(Friedmann, aCrunch[0], aCrunch[i+1], args=(0.31, 0, -0.31))
    timeCrunch.append(time[0])

tCrunch = [i-max(timeCrunch)+25 for i in timeCrunch]
tCrunchN = [-1*i+max(timeCrunch)+25 for i in timeCrunch]


a = np.linspace(0, aLim, 10001)
allTime = []

params = {'legend.fontsize': 15}
plt.rcParams.update(params)

while True:
    omega_lamb = float(input('Input a value for density parameter of Lambda: '))
    omega_m = float(input('Input a value for density parameter of Matter: '))
    omega_rad = float(input('Input a value for density parameter of Radiation: '))

    for i in range(len(a)-1):
        time = integrate.quad(Friedmann, a[0], a[i+1], args=(omega_m, omega_rad, omega_lamb))
        allTime.append(time[0])

    print("Age of the universe: ", round(time[0], 2), "Gyr")

    sub = allTime[np.where(a == 1)[0][0] - 1]  # To normalize time with respect to t0.

    xTime = [i - sub for i in allTime]  # t-t0

    plt.ioff()  # stops automatically displaying plots. When plt.show() is called, then it plots.
    plt.plot(xTime, a[1:], label='$\Omega_m = %g,\Omega_r = %G, \Omega_{\Lambda} = %g $' % (omega_m, omega_rad, omega_lamb))
    plt.ylim(0, a[-1])
    #plt.xlim(-20,50)
    plt.xlabel('Time (t-t0) in Gyr', fontsize=15)
    plt.ylabel('Scale Factor, a', fontsize=15)
    plt.legend()

    allTime = []

    decision = str(input('To continue with more models, type y. To end, type n: '))

    if decision == 'y':
        continue
    elif decision == 'n':
        print('Ending program!')
        break
    while decision != 'y' or decision != 'n':
        decision = str(input('Please choose y to continue, or n to end! : '))
        if decision == 'y':
            break
        elif decision == 'n':
            print('Ending program!')
            sys.exit(0)

plt.plot(tPlot, aBounce[1:], 'k--')
plt.plot(tPlotN, aBounce[1:], 'k--', label='Big Bounce')
plt.plot(tCrunch, aCrunch[1:], 'r')
plt.plot(tCrunchN, aCrunch[1:], 'r', label='Big Crunch')
plt.title('Evolution of the scale factor for various Universe models', fontsize=20)
plt.legend()
plt.show()
#plt.xlim(-20,50)
