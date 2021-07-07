# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import matplotlib.pyplot as plt
import numpy as np

fname = '29cmData.csv'
data = np.loadtxt(fname, delimiter=',',usecols=(3,4,9,10))
                     
timeR = 1000*data[:,0]
ampR = data[:,1]

timePeaks = 1000*data[:,2]
ampPeaks = data[:,3]

plt.figure(1)
plt.plot(timeR, ampR, linewidth = 0.3, label = 'Ramp')
plt.plot(timePeaks, ampPeaks, '.-', markersize = 0.5, linewidth = 1, label = 'Tranmission')
plt.xlabel('Time(ms)')
plt.ylabel('Voltage(V)')
plt.title('Transmission Spectrum at Cavity Length of 29cm')
#plt.legend(frameon=False)
plt.show()

cavLen = [8,11,13,16,20]
cavFin = [191,175,193,193,203]
cavWid = [151,143,132,125,117] #us

expFin = [61.2,61.2,61.2,61.2,61.2]
expWid = [32.64,44.88,53.04,65.28] #ns

plt.figure(2)
plt.plot(cavLen, cavFin, 'o-', label = 'Measured')
plt.plot(cavLen, expFin, label = 'Expected')
plt.xlabel('Length (cm)')
plt.ylabel('Finesse')
plt.title('Finesse as function of length')
plt.grid()
plt.legend()

calcWid = [9.82,7.79,5.98,4.86,3.69]#Mhz
cavWidExp = [30.6,22.3,18.8, 15.3,12.3]

plt.figure(3)
plt.plot(cavLen, calcWid, 'o-', label='Experimental')
plt.plot(cavLen, cavWidExp, '*-', label = 'Expected')
plt.xlabel('Length (cm)')
plt.ylabel('Linewidth (MHz)')
plt.title('Linewidth as function of Length')
plt.grid()
plt.legend()
