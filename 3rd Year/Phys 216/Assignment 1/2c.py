# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 16:15:43 2019

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt


# ===============================3(b)==========================================

# Initials conditions and constants we'll be using to solve the ODE:
b = 0.5  # drag constant in kg/s
g = -9.8  # gravity in m/s^2 y direction. 0 in x direction.
m = 100  # mass in kg

v = t = 0

vlist = []
tlist = []
ylist = []

dt = 1.56  # time step
y = 30000

while y > 0:
	v += g*dt
	vlist.append(v)

	t += dt
	tlist.append(t)

	y += v*dt
	ylist.append(y)

print("\nImpact velocity(m/s) without drag: {:3.3f}".format(v))
print("Time taken(s) without drag: {:3.3f}".format(t))
print("==============================================")	

#Interpolation 

if ylist[-1] < 0:
	deltat = -ylist[-2]*(tlist[-1]-tlist[-2])/(ylist[-1]-ylist[-2])
	print("\nNew delta-t after interpolation: {:3.3f}".format(deltat))
	
	yold = ylist[-1]
	y = ylist[-2] + v*deltat
	ylist[-1] = y
	print('Final position before and after interpolation respectively: {:3.3f}, {:3.3f}'.format(yold, ylist[-1]))

	v = vlist[-2] + g*deltat
	vlist[-1] = v
	
	t = tlist[-2] + deltat
	tlist[-1] = t

print("\nImpact velocity(m/s) without drag: {:3.3f}".format(v))
print("Time taken(s) without drag: {:3.3f}".format(t))
print("==============================================")




# ===============================3(c)==========================================

v = t = 0
vlist = []
tlist = []
ylist = []
y = 30000

while y > 0:
	v += ((b*v**2/m) + g)*dt  # gravity has been initialized as = -9.8 before
	vlist.append(v)
	
	t += dt
	tlist.append(t)
	
	y += v*dt
	ylist.append(y)

#Interpolation 

if ylist[-1] < 0:
	deltat = -ylist[-2]*(tlist[-1]-tlist[-2])/(ylist[-1]-ylist[-2])
	print("\nNew delta-t after interpolation: {:3.3f}".format(deltat))
	
	yold = ylist[-1]
	y = ylist[-2] + v*deltat
	ylist[-1] = y
	print('Final position before and after interpolation respectively: {:3.3f}, {:3.3f}'.format(yold, ylist[-1]))

	v = vlist[-2] + ((b*vlist[-2]**2/m) + g)*deltat
	vlist[-1] = v
	
	t = tlist[-2] + deltat
	tlist[-1] = t

print("\nImpact velocity(m/s) with drag: {:3.3f}".format(v))
print("Time taken(s) with drag: {:3.3f}".format(t))
print("==============================================")


# ====================================3(d)=====================================

H = 8000
rEarth = 6370000
y = 30000
v = t = 0

vlist = []
tlist = []
ylist = []
alist = []

def drag(z,v):
	b = 0.5*np.exp(-z/H)
	g = 9.8*(1+(z/rEarth))**(-2)
	a = ((b*v**2/m) - g)
	return a

while y > 0:
	a = drag(y,v)
	alist.append(a)
	
	v += a*dt
	vlist.append(v)
	
	t += dt
	tlist.append(t)
	
	y += v*dt
	ylist.append(y)

# Interpolation below:
	
if ylist[-1] < 0:
	deltat = -ylist[-2]*(tlist[-1]-tlist[-2])/(ylist[-1]-ylist[-2])
	print("\nNew delta-t after interpolation: {:3.3f}".format(deltat))
	
	yold = ylist[-1]
	y = ylist[-2] + v*deltat
	ylist[-1] = y
	print('Final position before and after interpolation respectively: {:3.3f}, {:3.3f}'.format(yold, ylist[-1]))
	
	a = drag(y,vlist[-2])
	alist[-1] = a
	
	v = vlist[-2] + a*deltat
	vlist[-1] = v
	
	t = tlist[-2] + deltat
	tlist[-1] = t


print("\nImpact velocity(m/s) with drag and corrected parameters: {:3.3f}".format(v))
print("Time taken(s) with drag and corrected parameters: {:3.3f}".format(t))
print("==============================================")

plt.figure(1)
plt.plot(tlist, vlist, label= '(d) Scaled b,g')
plt.xlabel('Time(s)')
plt.ylabel('Velocity(m/s)')
plt.title('Velocity vs Time - With Drag and Corrected parameters')
plt.grid()
plt.legend()

plt.figure(2)
plt.plot(tlist, ylist)
plt.xlabel('Time(s)')
plt.ylabel('Altitude (m)')
plt.title('Altitude vs Time - With Drag and Corrected parameters')
plt.grid()

plt.figure(3)
plt.plot(tlist, alist)
plt.xlabel('Time(s)')
plt.ylabel('Acceleration (m/s^2)')
plt.title('Acceleration vs Time - With Drag and Corrected parameters')
plt.grid()

