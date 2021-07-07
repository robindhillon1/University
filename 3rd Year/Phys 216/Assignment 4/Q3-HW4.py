# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 16:00:22 2019

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

azi = (round(0.432*360))*(np.pi/180)
alt = (round(0.432*70))*(np.pi/180)
lat = 49*(np.pi/180)
w = 7.29e-5

v0 = 200
vx = v0*np.cos(alt)*np.sin(azi)
vy = v0*np.cos(alt)*np.cos(azi)
vz = v0*np.sin(alt)

m = 10
H = 8000
rEarth = 6370000
t = 0

dt = 0.03
x = y = 0
z = 4

tlist = []
xlist = []
ylist = []
zlist = []
alist = []

alistPC = []
alistPCx = []
alistPCy = []

xlistPC = []
ylistPC = []
zlistPC = []

vxFinList = []
vyFinList = []
vzFinList = []
vMagList = []

xFinList = []
yFinList = []
zFinList = []

ranges = []

def drag(z,b):
	drag = b*np.exp(-z/H)
	return drag

def grav(z):
	g = 9.8*(1+(z/rEarth))**(-2)
	return g

while z > 0:
	vMag = np.sqrt((vx**2+vy**2+vz**2))
	vMagList.append(vMag)
	
	ax = 2*w*(vy*np.sin(lat) - np.cos(lat)*(vz+grav(z)*t)) + x*w**2 - drag(z, 0.043)*vx*vMag/m
	ay = -2*w*vx*np.sin(lat) + (np.sin(lat)*w**2)*(-z*np.cos(lat) + y*np.sin(lat)) - drag(z, 0.043)*vy*vMag/m
	az = 2*w*vx*np.cos(lat) - grav(z) + (np.sin(lat)*w**2)*(-y*np.sin(lat) + z*np.cos(lat)) - drag(z, 0.043)*vMag*vz/m
	alist.append(az)
	
	zm = z*1
	z = zm + vz*dt
	zlist.append(z)
	
	xm = x*1
	x = xm + vx*dt
	xlist.append(x)
	
	ym = y*1
	y = ym + vy*dt
	ylist.append(y)
	
	vz0 = vz*1
	vz = vz0 + az*dt
	
	vx0 = vx*1
	vx = vx0 + ax*dt
	
	vy0 = vy*1
	vy = vy0 +ay*dt
	
	t += dt
	tlist.append(t)
	
#=====================================================================
	vMag = np.sqrt((vx**2+vy**2+vz**2))
	vMagList.append(vMag)
	
	aPCx = 2*w*(vy*np.sin(lat) + np.cos(lat)*(grav(z)*t-vz)) + x*w**2 - drag(z, 0)*vx*vMag/m
	aPCy = -2*w*vx*np.sin(lat) + (np.sin(lat)*w**2)*(-z*np.cos(lat) + y*np.sin(lat)) - drag(z, 0)*vy*vMag/m
	aPCz = 2*w*vx*np.cos(lat) - grav(z) + (np.cos(lat)*w**2)*(-y*np.sin(lat) + z*np.cos(lat)) - drag(z, 0)*vMag*vz/m
	alistPC.append(aPCz)
	alistPCx.append(aPCx)
	alistPCy.append(aPCy)
	
	vPCx = vx0 + aPCx*dt
	vPCy = vy0 + aPCy*dt
	vPCz = vz0 + aPCz*dt
	
	zPC = zm + vPCz*dt
	xPC = xm + vPCx*dt
	yPC = ym + vPCy*dt
	
	vzFinal = (vz+vPCz)/2
	vzFinList.append(vzFinal)
	
	vxFinal = (vx+vPCx)/2
	vxFinList.append(vxFinal)
	
	vyFinal = (vy+vPCy)/2
	vyFinList.append(vyFinal)
	
	vMag = np.sqrt((vx**2+vy**2+vz**2))
	vMagList.append(vMag)
	
	zFinal= (z+zPC)/2
	zFinList.append(zFinal)
	
	xFinal= (x+xPC)/2
	xFinList.append(xFinal)
	
	yFinal= (y+yPC)/2
	yFinList.append(yFinal)
	
	rngs = np.sqrt(xFinal**2+yFinal**2)
	ranges.append(rngs)

# INTERPOLATION
	
deltatPC= -zFinList[-2]/vzFinList[-2]
print("\nNew delta-t after interpolation: {:3.6f}\n".format(deltatPC))

zoldPC = zFinList[-1]
zPC = zFinList[-2] + vzFinList[-2]*deltatPC
zFinList[-1] = zPC

xPC = xFinList[-2] + vxFinList[-2]*deltatPC
xFinList[-1] = xPC

yPC = yFinList[-2] + vyFinList[-2]*deltatPC
yFinList[-1] = yPC

ranges[-1] = np.sqrt(xPC**2+yPC**2)

vz = vzFinList[-2] + alistPC[-2]*deltatPC
vzFinList[-1] = vz

vx = vxFinList[-2] + alistPCx[-2]*deltatPC
vxFinList[-1] = vx

vy = vyFinList[-2] + alistPCy[-2]*deltatPC
vyFinList[-1] = vy

aPCz = ((-drag(zFinList[-2], 0.043)*vMag*vzFinList[-2]/m) - grav(zFinList[-2]))
alistPC[-1] = aPCz

aPCx = -drag(zFinList[-2], 0.043)*vxFinList[-2]*vMag/m
aPCy = -drag(zFinList[-2], 0.043)*vyFinList[-2]*vMag/m

t = tlist[-2] + deltatPC
tlist[-1] = t


print('Final position before and after interpolation respectively: {:3.3f}, {:3.3f}\n'.format(zoldPC, zPC))

print('(a) Final range (m) = {:3.5}, at x = {:3.5f} and y = {:3.5f}.\n'.format(ranges[-1], xFinList[-1], yFinList[-1]))

print('(b) Altitude = {:3.3f}, and azimuth = {:3.3f}\n'.format(alt, azi))

print('(c) Flight Time of Pumpkin (s): {:3.5f}\n'.format(tlist[-1]))

print('(d) Velocity vector at time of landing (cm/s): <vx,vy,vz> = <{:3.5f},{:3.5f},{:3.5f}>'.format(100*vxFinList[-1], 100*vyFinList[-1], 100*vzFinList[-1]))

plt.figure(1)
plt.plot(ranges, zFinList)
plt.title('Range vs Height')
plt.xlabel('Range (m)')
plt.ylabel('Height (m)')
plt.grid()

plt.figure(2)
plt.plot(xFinList, yFinList)
plt.title('x-y Path of Pumpkin')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid()