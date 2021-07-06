# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 19:08:37 2018

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('vel.txt',skip_header=5)

HJD = data[:,0]
v_p = data[:,1]
v_s = data[:,2]
phase = data[:,3] #runs from -0.5 to 0.5
#v_p = v_p[np.argsort(phase)]
#phase=np.sort(phase)
plt.figure(1)
plt.clf()
plott = plt.scatter(phase,v_p,s=50,marker=".",color="red",label='Primary')
plott2 = plt.scatter(phase,v_s,s=50,marker=".",color="blue",label='Secondary')
plt.xlabel('Orbital Phase')
plt.ylabel('Radial Velocity(km/s)')
plt.title('Radial Velocity Curve')
plt.grid(True)

from scipy import optimize

def test_func(x, a, b,offset,c):
    return a*np.sin(b*x+offset)+c

testphase=np.linspace(-0.5,0.5,100)

params, params_covariance = optimize.curve_fit(test_func, phase, v_p, p0=[50.99,1,0,80])
print(params)
sin1 = plt.plot(testphase,test_func(testphase,*params),'r-',label='BestFit-P')

def test_func2(y, d, e,offset2,f):
    return d*np.sin(e*y+offset2)+f

testphase=np.linspace(-0.5,0.5,100)

params2, params_covariance2 = optimize.curve_fit(test_func, phase, v_s, p0=[50.99,1,0,80])
print(params2)
sin2 = plt.plot(testphase,test_func(testphase,*params2),'b-',label='BestFit-S')
plt.legend()

print("\n")

#Calculations
G = 6.67*10**-11
period = 8.1113*24*3600
vp = params[0]
print("Velocity of Primary Stars(km/s): "+str(vp))
ap = (vp*period)/(2*np.pi)
print("Semi-major axis of Primary Star(km): "+str(ap)+"\n")

vs = params2[0]
print("Velocity of Secondary Star(km/s): "+str(vs))
ass = (vs*period)/(2*np.pi)
print("Semi-major axis of Secondary Star(km): "+str(ass)+"\n")

comvel = (params[3]+params2[3])/2
print("The center of mass velocity of the system(km/s): "+str(comvel)+"\n")

massS = ((vp*10**3+vs*10**3)**3*period)/(2*np.pi*G*((vs/vp)+1))
print("Mass of Secondary Star(kg): "+str(massS))

massP = (vs/vp)*massS
print("Mass of Primary Star(kg): "+str(massP)+"\n")

phase1 = 0.710033
phase2 = 0.722614
phase3 = 0.724023
 
t1 = phase1*period
t2 = phase2*period
t3 = phase3*period

v = vp+vs

Rs = 0.5*v*(t2-t1)
print("Radius of Secondary Star(km): "+str(Rs))

Rp = 0.5*v*(t3-t1)
print("Radius of Primary Star(km): "+str(Rp))