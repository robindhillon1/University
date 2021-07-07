# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 18:48:07 2020

@author: zhang
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sys
"""
h = 6.62607004*10**(-34)/(1.6022*10**(-10)*6.5823*10**(-25))
c = 2.9979*10**8/(2.9979*10**8)
k = 1.38*10**(-23)/(1.6022*10**(-10))
g = 2
T_0n = 1.95#1.68*10**(-4) #*k/(1.6022*10**(-10))
m_nu = (3.2*10**(-21)/c**2)/(1.7827*10**(-27))#0.02
#rho_crit = 8.7*10**(-27)/(2.0852*10**37)
rho_crit = 4870e6*1.602*10**(-19)/(2.0852*10**(37))
"""
h = 6.62607004*10**(-34)
c = 2.9979*10**8
k = 1.38*10**(-23)
T_0n = 1.95  # 1.68*10**(-4) #*k/(1.6022*10**(-10))
m_nu = 1.6*10**(-20)/c**2  # 0.02
rho_crit = 4870*10**6*1.602*10**(-19)

H0 = (68/3.086E19)*(31536000*1E9)  # in 1/Gyr


def Friedmann(a, omega_m, omega_rad, omega_lamb):
    den = np.sqrt(omega_m/a + a**2*omega_lamb + omega_rad/(a*a) + (1-omega_lamb-omega_m-omega_rad))
    t = 1/(H0*den)
    return t


def Friedmann2(a, omega_m, omega_rad, omega_lamb, neut, nFinal):
    den = np.sqrt(omega_m/a + a**2*omega_lamb + omega_rad/(a*a) + (1-omega_lamb-omega_m-omega_rad-nFinal) + neut*a**2)
    t = 1/(H0*den)
    return t


def Neutrino(x, m, T0, a):
    const = (8*(np.pi)*c/(h**3))*(T_0n*k/c)**3
    num = x**2*(np.sqrt((T_0n*k*x/c)**2+m_nu**2*c**2))
    den = np.exp(x*a) + 1
    drho = const*num/(den*rho_crit)
    return drho


params = {'legend.fontsize': 15}
plt.rcParams.update(params)

aLim = float(input('Choose a value for the scale factor, a: '))

a = np.linspace(0, aLim, 10001)

neutrinoValues = []
for i in range(len(a)-1):
    neu = integrate.quad(Neutrino, 0, 10000, args=(m_nu, T_0n, a[i]))
    neutrinoValues.append(neu[0])

nFinal = neutrinoValues[-1]
print(neu)

allTime = []
allTimeN = []

while True:
    omega_lamb = float(input('Input a value for density parameter of Lambda: '))
    omega_m = float(input('Input a value for density parameter of Matter: '))
    omega_rad = float(input('Input a value for density parameter of Radiation: '))

    for i in range(len(a)-1):
        time = integrate.quad(Friedmann, a[0], a[i+1], args=(omega_m, omega_rad, omega_lamb))
        allTime.append(time[0])

    print("Age of the universe not including neutrinos: ", round(time[0], 2), "Gyr")

    for i in range(len(a)-1):
        timeN = integrate.quad(Friedmann2, a[0], a[i+1], args=(omega_m, omega_rad, omega_lamb, neutrinoValues[i],  nFinal))
        allTimeN.append(timeN[0])

    print("Age of the universe  including neutrinos: ", round(timeN[0], 2), "Gyr")

    sub = allTime[np.where(a == 1)[0][0] - 1]  # To normalize time with respect to t0.
    xTime = [i - sub for i in allTime]  # t-t0

    subN = allTimeN[np.where(a == 1)[0][0] - 1]
    xTimeN = [i - sub for i in allTimeN]  # t-t0 for neutrinos

    plt.ioff()  # stops automatically displaying plots. When plt.show() is called, then it plots.
    plt.plot(xTime, a[1:], label='Age without neutrinos: {:.4f}Gyr'.format(time[0]))
    plt.plot(xTimeN, a[1:], 'k--', label='Age with neutrinos: {:.4f}Gyr'.format(timeN[0]))
    plt.ylim(0, a[-1])
    #plt.xlim(-20,50)
    plt.xlabel('Time (t-t0) in Gyr', fontsize=15)
    plt.ylabel('Scale Factor, a', fontsize=15)

    allTime = []
    allTimeN = []


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

plt.title(r'Benchmark Model including $m_{\nu} = 0.1eV/c^2$', fontsize=20)
plt.legend()
plt.show()

