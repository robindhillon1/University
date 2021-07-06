# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 12:56:56 2018

@author: robin
"""

import math
import numpy as np
import matplotlib.pyplot as plt

fname = 'V_edit.txt'
data = np.genfromtxt(fname, comments='#')
                     
#f = open('V_edit.txt','r')                    
#file_contents = f.read()
#print(file_contents)
#f.close()

dateJD = data[:,0]
vmag = data[:,1]
photoerr = data[:,2]
phase = data[:,3] #Runs between 0 and 1.
#Represents the fraction of the period when the data were taken 

plt.figure(1)
plt.clf()
plott = plt.scatter(dateJD,vmag,s=5,marker=".",color="red")
plt.xlabel('JD')
plt.ylabel('V-Magnitude')
plt.title('Light Curve of Eclipsing Binary')
plt.gca().invert_yaxis()

period =  8.1113 #period (must be known already!) (2767.766470-2763.713520),1460.04032
phase = dateJD / period # divide by period to convert to phase
phase = phase % 1 # take fractional part of phase only (i.e. discard whole number part)

plt.figure(2)
plt.clf()
plott2 = plt.scatter(phase,vmag,s=10,marker=".",color="red")
plt.xlabel('Orbital Phase')
plt.ylabel('V-Magnitude')
plt.title('Light Curve of Eclipsing Binary')
plt.gca().invert_yaxis()