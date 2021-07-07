# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 22:59:17 2019

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare

data = np.loadtxt('Week4Data.csv', delimiter=',', comments='#')

mass = data[:, 0]
vol = data[:, 1] 
volUnc = data[:, 2]
rebounds = data[:, 5]

volList = vol.tolist()

P1 = 1.013e5

P = [P1]

T = 293

n = 1.44*10**(-4)

for i in range(len(vol)):
	if i == 6:
		break
	
	pres = P[i]*volList[i]/volList[i+1]
	P.append(pres)


# =============================================================================
# plt.plot(volList,P, 'ko')
# 
# vol = np.array(volList)
# 
# 
# plt.xlabel('Volume')
# plt.ylabel('Pressure')
# plt.title('Pressure as Volume is varied by adding Masses')
# plt.grid()
# =============================================================================


def func(m, P0, R):
	T = 293
	#P0 = 1.013e5
	A = np.pi*(7.15e-3)**2
	n = 0.00017857142857142857  # at STP, 4/22.4e3 = n. THis gives R = 8.27
	g = 9.81
	y = (1/(n*R*T))*(P0-(m*g)/A + 1/A)
	return y


vol2 = 1/(vol*10**(-6))
pres2 = P

Pres = np.array(P)
popt, pcov = curve_fit(func, mass, vol2, p0=[100000,8.31], maxfev=10000000)
print("R from Hanging Volumes: "+str(popt))

plt.plot(mass, vol2, 'ko', label = 'Data')
plt.plot(mass, func(mass, *popt), 'r:', linewidth = 2, label = 'Fit')
plt.xlabel('Mass (kg)')
plt.ylabel('1/Volume (m^(-3))')
plt.title('1/Volume vs Mass')
plt.grid()
plt.legend()


rebounds = 1/(rebounds*10**(-6))
popt2, pcov2 = curve_fit(func, mass[1:6], rebounds[1:6], p0=[100000,8.31], maxfev=10000000)
print("R from Rebound Volumes: "+str(popt2),'\n')

R = (popt[1]+popt2[1])/2
P = (popt[0]+popt2[0])/2

print("Avg Value of R: {}, and rounding it gives: {:3.3}".format(R,R))
print("\n")
print("Avg Value of P0: {}, and rounding it gives: {:3.5}".format(P,P))
print("\n")
# =============================================================================
# vRatio = (vol/vol[0])
# vDiff = (vol - vol[0])*10**(-6)
# 
# lnRatio = np.log(vRatio)
# 
# x = ((P[0]*np.pi*r**2 - mass*g)/(n*T))*(vDiff)
# 
# plt.figure(2)
# plt.plot(x, lnRatio, 'ko')
# =============================================================================

#xx = linear(vol2, *popt)

#R = xx[0]*vol[0]*(10**(-6))/(n*T)
#print('R = ', R)


# =============================================================================
# data = np.loadtxt('Week2Data.csv', delimiter=',', comments='#')
# 
# mass = data[:, 0]
# massUnc = data[:, 1]
# dV = data[:, 2]
# uDV = data[:, 3]
# vFinal = data[:, 4]
# vInit = data[:, 5]
# radius = data[:, 6]
# uRadius = data[:, 7]
# 
# n = 4.2446299064738245e-05
# 
# vRatio = vFinal/vInit
# vDiff = vFinal - vInit
# 
# lnRatio = np.log(vRatio)
# x = (P[0]-mass*g/np.pi*radius**2)*(vDiff/n*T)
# 
# plt.figure(2)
# plt.plot(x, lnRatio)
# =============================================================================


# =============================================================================
# chi2 = chisquare(P, linear(volList, *popt))
# print("\nChisquare: ",chi2)
# =============================================================================
