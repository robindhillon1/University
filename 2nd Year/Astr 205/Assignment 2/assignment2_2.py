# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 15:16:46 2018

@author: robin
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

fname = 'HW2_Q2_data.txt'
#load the data from the file
data = np.genfromtxt(fname, comments='#')

vmag = data[:,0]
imag = data[:,1]
dx = data[:,2]
dy = data[:,3]

plt.figure(1)
plt.clf()   
plt.scatter((vmag-imag),vmag,s=0.5,marker=".",color="red")
plt.xlabel('V-I')
plt.ylabel('V')
plt.title('Colour-Magnitude Diagram(CMD)')
plt.gca().invert_yaxis()
#This is an HR diagram

plt.figure(2)
plt.clf()
alls = plt.scatter(dx,dy,s=0.5,marker=".",color="blue")
rangtuc = np.where((dx**2+dy**2)<0.3**2)
rangcld = np.where(((dx+0.6)**2+(dy+0.2)**2)<0.15**2)
sepr = plt.scatter(dx[rangtuc],dy[rangtuc],s=0.5,marker=".",color="red")
plt.xlabel('dx(Milli-arcseconds/year)')
plt.ylabel('dy(Milli-arcseconds/year)')
plt.title('Proper Motion of the Stars')

lgnd = plt.legend((alls,sepr),('All Stars','47 Tuc Stars'),scatterpoints=1,fontsize=8)
lgnd.legendHandles[0]._sizes=[6]
lgnd.legendHandles[1]._sizes=[6]

plt.figure(3)
plt.clf()
alls2 = plt.scatter(vmag[rangtuc]-imag[rangtuc],vmag[rangtuc],s=0.5,marker=".",color="red")
sepr2 = plt.scatter(vmag[rangcld]-imag[rangcld],vmag[rangcld],s=0.5,marker=".",color="blue")
plt.xlabel('V-I')
plt.ylabel('V')
plt.title('CMD of 47 Tuc stars and SMC')
plt.gca().invert_yaxis()

lgnd = plt.legend((alls2,sepr2),('47 Tuc Stars','SMC Stars'),scatterpoints=1,fontsize=8)
lgnd.legendHandles[0]._sizes=[6]
lgnd.legendHandles[1]._sizes=[6]