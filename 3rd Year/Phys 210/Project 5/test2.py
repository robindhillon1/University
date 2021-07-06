# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 15:17:40 2018

@author: robin
"""

# WHAT THE PROJECT DOES

#BONUS: Solving the Heat Equation at different temperatures

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.animation as animation

# constants
L = 1 # size of "box"
D = 1e-3
N = 100 # number of grid points in one-d

dt = 1e-2
dx = L/N

tmax = 10
steps = int(tmax/dt)+1

# choose some times to make plots:
# given in step numbers
plotsteps = np.array([0.001, 0.01, 0.1, 1, 10, 100, 1000])
plotsteps /= dt
plotsteps = plotsteps.astype(int)
# create initial conditions:
C = np.zeros((N,N))

# Set particles in a blob in the center:
C[N//2, N//2] = 10
k = dt/dx/dx*D

Cp = np.zeros((N, N))
isum = np.sum(C)

pcol = []

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

sumsteps = np.array([0, 0.1, 1, 10, 100])
sumsteps /= dt
sumsteps = sumsteps.astype(int)

esumm = []

for i in range(steps):
# do all the points, except the boundaries
	#Cp[j+1, k+1] = C[j+1, k+1]+k*(C[j, k+1]+C[j+2, k+1]-2*C[j+1, k+1])
	Cp_xLe = np.roll(C, -1, axis=0)
	Cp_xRi = np.roll(C, 1, axis=0)
	Cp_yUp = np.roll(C, 1, axis=1)
	Cp_yDw = np.roll(C, -1, axis=1)
	Cp_xnew = C + k*(Cp_xLe+Cp_xRi-2*C)
	Cp_new = Cp_xnew + k*(Cp_yUp+Cp_yDw-2*C)

	Cp_new[0,:] = Cp_new[1,:]  # first row, all columns
	Cp_new[-1,:] = Cp_new[-2,:]  # last row, all columns.
	Cp_new[:,0] = Cp_new[:,1]  # first column, all rows
	Cp_new[:,-1] = Cp_new[:,-2]  # last column, all rows
		# swap C and Cp so that they don't end up as the same array:
	C, Cp_new = Cp_new, C
	
	pcol.append((ax2.pcolormesh(C.copy(), cmap='jet'),))

	#if i in plotsteps:
	#	ax1.plot(C[:, N//2], label="t = %g" % (i*dt))
	
	if i in sumsteps:
		esumm.append(np.sum(C))
		# write at i, this is the sum to the text file.
#if the boundary conditions are done well, the number of particles
# ie, the integral of C, should be constant. Check:
# =============================================================================
# esum = np.sum(C)
# print("initial and final integrals of concentration:", isum, esum, "\n")
# plt.xlabel('Position')
# plt.ylabel('Concentration')
# plt.grid()
# plt.legend()
# plt.show()
# =============================================================================

# =============================================================================

plt.figure()
C_slice = C[:,50]  # 49th row, all columns
time = np.linspace(0.01,10,100)
plt.plot(time[::3], C_slice[::3],'ko',markersize=2.5)

#def gauss(x, sigma, miu):
#    return  1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-miu)**2/(2*sigma**2))

def gauss(x, a, sigma, miu):
    return  a*np.exp(-(x-miu)**2/(2*sigma**2))

popt, pcov = curve_fit(gauss, time[::3], C_slice[::3], p0=[0.6, 2.8, 0])
print("Fit parameters: ", *popt)
plt.plot(time[::3], gauss(time[::3], *popt), 'r-', label='fit')
plt.grid()
plt.legend()
plt.show()

# =============================================================================

anim = animation.ArtistAnimation(fig2, pcol, interval=10, repeat=False)
anim.save('converges.mp4')

# =============================================================================
f = open('totals.txt','w')
for i in range(len(esumm)):
	f.write('At t={}s, number = {}\n.'.format(sumsteps[i], esumm[i]))
f.close()