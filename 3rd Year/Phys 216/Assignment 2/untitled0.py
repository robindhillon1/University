# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 20:11:04 2019

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt


# ===============================3(b)==========================================

# Initials conditions and constants we'll be using to solve the ODE:
b = 0.5  # drag constant in kg/s
g = -9.8  # gravity in m/s^2 y direction. 0 in x direction.
m = 100  # mass in kg

v = t = vPC = 0

vlist = []
tlist = []
ylist = []
ylistPC = []
vlistPC = []

vFinList = []
yFinList = []

dt = 1.58  # time step
y = yPC = 30000

while y > 0:
	vm = v
	v = vm + g*dt
	vlist.append(v)

	t += dt
	tlist.append(t)
	
	ym = y
	y = ym + v*dt
	ylist.append(y)
	
	vPC = vm + g*dt  # same as before since g is the same in both cases.
	vlistPC.append(vPC)
	
	yPC = ym + vPC*dt
	ylistPC.append(yPC)
	
	vFinal = (v+vPC)/2
	vFinList.append(vFinal)
	yFinal = (y+yPC)/2
	yFinList.append(yFinal)
	
print("=============3(a)================")
print("Impact velocity(m/s) without drag: {:3.3f}".format(v))
print("Time taken(s) without drag: {:3.3f}".format(t))

#Interpolation 

deltat = -ylist[-2]*(tlist[-1]-tlist[-2])/(ylist[-1]-ylist[-2])
print("\nNew delta-t after interpolation: {:3.3f}".format(deltat))

yold = ylist[-1]
y = ylist[-2] + v*deltat
ylist[-1] = y
print('Final position before and after interpolation respectively: {:3.3f}, {:3.3f}'.format(yold, ylist[-1]))

vold = vlist[-1]
v = vlist[-2] + g*deltat
vlist[-1] = v

told = tlist[-1]
t = tlist[-2] + deltat
tlist[-1] = t

print("\nImpact velocity(m/s) without drag: {:3.3f}".format(v))
print("Time taken(s) without drag: {:3.3f}".format(t))


print("\n--> Now, using the Predictor-Corrector or Modified Euler:")
print("\nImpact velocity(m/s) without drag: {:3.3f}".format(vPC))
print("Time taken(s) without drag: {:3.3f}".format(t))

deltatPC = -yFinList[-2]*(told-tlist[-2])/(yFinList[-1]-yFinList[-2])
print("\nNew delta-t for PC after interpolation: {:3.3f}".format(deltatPC))

yoldPC = yFinList[-1]
yPC = yFinList[-2] + vFinal*deltatPC
yFinList[-1] = yPC
print('Final position before and after interpolation respectively: {:3.3f}, {:3.3f}'.format(yoldPC, yFinList[-1]))

v = vFinList[-2] + g*deltatPC
vlist[-1] = v

t = tlist[-2] + deltatPC
tlist[-1] = t

print("\nImpact velocity(m/s) without drag: {:3.3f}".format(v))
print("Time taken(s) without drag: {:3.3f}".format(t))
print("==============================================")













# =============================================================================
# print("\n=============3(b)================")
# 
# H = 8000
# rEarth = 6370000
# y = yPC = 30000
# v = t = vPC = 0
# m = 100
# 
# vlist = []
# vlistPC = []
# tlist = []
# ylist = []
# ylistPC = []
# alist = []
# 
# vFinList = []
# yFinList = []
# alistPC = []
# 
# dt = 7.54
# 
# def drag(z,v):
# 	b = 0.5*np.exp(-z/H)
# 	g = 9.8*(1+(z/rEarth))**(-2)
# 	a = ((b*v**2/m) - g)
# 	return a
# 
# while (y and yPC) > 0:
# 	a = drag(y,v)
# 	alist.append(a)
# 	
# 	vm = v
# 	v = vm + a*dt  #predicted
# 	vlist.append(v)
# 	
# 	t += dt
# 	tlist.append(t)
# 	
# 	ym = y
# 	y = ym + v*dt  #predicted
# 	ylist.append(y)
# #=====================================================================
# 	aPC = drag(y,v)
# 	alistPC.append(aPC)
# 	
# 	vPC = vm + aPC*dt
# 	vlistPC.append(vPC)
# 	
# 	yPC = ym + v*dt
# 	ylistPC.append(yPC)
# 	
# 	vFinal = (v+vPC)/2
# 	vFinList.append(vFinal)
# 	yFinal= (y+yPC)/2
# 	yFinList.append(yFinal)
# 
# # Interpolation below:
# deltat = -ylist[-2]*(tlist[-1]-tlist[-2])/(ylist[-1]-ylist[-2])
# print("New delta-t after interpolation: {:3.3f}".format(deltat))
# 
# yold = ylist[-1]
# y = ylist[-2] + v*deltat
# ylist[-1] = y
# print('Final position before and after interpolation respectively: {:3.3f}, {:3.3f}'.format(yold, ylist[-1]))
# 
# a = drag(y,vlist[-2])
# alist[-1] = a
# 
# vold = v
# v = vlist[-2] + a*deltat
# vlist[-1] = v
# 
# told = tlist[-1]
# t = tlist[-2] + deltat
# tlist[-1] = t
# 
# print("\nImpact velocity(m/s) with drag and corrected parameters: {:3.3f}".format(v))
# print("Time taken(s) with drag and corrected parameters: {:3.3f}".format(t))
# 
# 
# print("\n--> Now, using the Predictor-Corrector or Modified Euler:")
# 
# deltatPC = -yFinList[-2]*(told-tlist[-2])/(yFinList[-1]-yFinList[-2])
# print("\nNew delta-t after interpolation: {:3.3f}".format(deltatPC))
# 
# yoldPC = yFinList[-1]
# yPC = yFinList[-2] + vFinal*deltatPC
# yFinList[-1] = yPC
# print('Final position before and after interpolation respectively: {:3.3f}, {:3.3f}'.format(yoldPC, yFinList[-1]))
# 
# a = drag(yPC,vFinList[-2])
# alistPC[-1] = a
# 	
# v = vFinList[-2] + a*deltatPC
# vlist[-1] = v
# 
# t = tlist[-2] + deltatPC
# tlist[-1] = t
# 
# print("\nImpact velocity(m/s) with drag and corrected parameters: {:3.3f}".format(v))
# print("Time taken(s) with drag and corrected parameters: {:3.3f}".format(t))
# print("==============================================")
# # =============================================================================
# # 
# # plt.figure(1)
# # plt.plot(tlist, vlist, 'r-', linewidth = 2, label= 'Euler')
# # plt.plot(tlist, vFinList, 'k:', linewidth = 2, label= 'Modified Euler')
# # plt.xlabel('Time(s)')
# # plt.ylabel('Velocity(m/s)')
# # plt.title('Velocity vs Time - With Drag and Corrected parameters')
# # plt.grid()
# # plt.legend()
# # 
# # plt.figure(2)
# # plt.plot(tlist, ylist, 'r-', linewidth = 2, label = 'Euler')
# # plt.plot(tlist, yFinList, 'k:', linewidth = 2, label = 'Modified Euler')
# # plt.xlabel('Time(s)')
# # plt.ylabel('Altitude (m)')
# # plt.title('Altitude vs Time - With Drag and Corrected parameters')
# # plt.grid()
# # plt.legend()
# # 
# # plt.figure(3)
# # plt.plot(tlist, alist, 'r-', linewidth = 2, label = 'Euler')
# # plt.plot(tlist, alistPC, 'k:', linewidth = 2, label = 'ModifiedEuler')
# # plt.xlabel('Time(s)')
# # plt.ylabel('Acceleration (m/s^2)')
# # plt.title('Acceleration vs Time - With Drag and Corrected parameters')
# # plt.grid()
# # plt.legend()
# # =============================================================================
# 
# =============================================================================
