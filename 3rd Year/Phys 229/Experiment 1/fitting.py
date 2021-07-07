# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 16:02:28 2019

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare

data = np.loadtxt('800HzResonances.csv', delimiter=',', comments='#')

freq = data[:, 0]
freqUnc = data[:, 1]
amp = data[:, 2]
ampUnc = data[:, 3]


def combined(f, a1, fr1, w1, w2, a2, fr2, m, b):
	#gauss = a1*np.exp(-(f-fr1)**2/(2*w1**2)) + (m*f + b)
	lorentz = (w1*a1)/(w1**2+(f-fr1)**2)**0.5 + (m*f + b)
	lorentz2 = (w2*a2)/(w2**2+(f-fr2)**2)**0.5 + (m*f + b)
	return lorentz + lorentz2

plt.figure()

popt3, pcov3 = curve_fit(combined, freq, amp, p0=[2, 787, 5, 5, 3, 798, -1, 1], maxfev=10000000)
print("Lorentz+Lorentz Fit Parameters: "+str(popt3))

plt.subplot(212)

plt.plot(freq, amp, 'ko', markersize = 1.5, label = 'Data')
plt.plot(freq, combined(freq, *popt3), 'r:', linewidth = 2, label = 'Fit')
plt.errorbar(freq, amp, yerr= ampUnc, c = 'black', marker= '.', linestyle = '')

plt.axvline(popt3[1], linestyle = '--', color = 'green')
plt.text(786,1,'Fitted $f_r$: \n<--{:3.5} Hz'.format(popt3[1]))#,rotation=90)
plt.axvline(popt3[5], linestyle = '--', color = 'green')
plt.text(799,0.5,'Fitted $f_r$: \n<--{:3.5} Hz'.format(popt3[5]))#,rotation=90)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (V)')
plt.title('Resonant Frequencies characterized by Lorentz Functions')
plt.grid()
plt.legend()

chi1 = chisquare(amp, combined(freq, *popt3))
print("\nChiqsuare for Lorentz-Lorentz: ",chi1)


def combined2(f, a1, fr1, w1, w2, a2, fr2, m, b):
	gauss = a1*np.exp(-(f-fr1)**2/(2*w1**2)) + (m*f + b)
	lorentz = (w2*a2)/(w2**2+(f-fr2)**2)**0.5 + (m*f + b)
	return gauss + lorentz

popt3, pcov3 = curve_fit(combined2, freq, amp, p0=[2, 787, 5, 5, 3, 798, -1, 1], maxfev=10000000)
print("\nLorentz+Gaussian Fit Parameters: "+str(popt3))

plt.subplot(211)

plt.plot(freq, amp, 'ko', markersize = 1.5, label = 'Data')
plt.plot(freq, combined2(freq, *popt3), 'r:', linewidth = 2, label = 'Fit')
plt.errorbar(freq, amp, yerr= ampUnc, c = 'black', marker= '.', linestyle = '')
plt.axvline(popt3[1], linestyle = '--', color = 'green')
plt.text(786,1,'Fitted $f_r$ (Gaussian) \n <--{:3.5} Hz'.format(popt3[1]))
plt.axvline(popt3[5], linestyle = '--', color = 'green')
plt.text(799,0.5,'Fitted $f_r$ (Lorentz) \n <--{:3.5} Hz'.format(popt3[5]))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (V)')
plt.title('Resonant Frequencies characterized by Lorentz-Gaussian Functions')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('finalPlot.pdf')

chi2 = chisquare(amp, combined2(freq, *popt3))
print("\nChisquare for Lorentz-Gaussian: ",chi2)