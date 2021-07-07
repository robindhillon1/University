# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 21:10:44 2020

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

def ref(ni,ns,n2,n1,N):
    num = (ni/ns)*(-1*n1/n2)**(2*N)-1
    den = (ni/ns)*(-1*n1/n2)**(2*N)+1
    R = (num/den)**2
    return R

w0 = 632.8*10**(-9)
w1 = 583*10**(-9)
w2 = 683*10**(-9)
wavel = np.linspace(w1, w2, 30)
N = np.arange(0,31,1)

refl = ref(1,1.52,2.1,1.65,N)*100
refl2 = np.array([round(x, 2) for x in refl])

print("\nN when R = 99.97%: ", N[np.where(refl2==99.97)])

plt.figure(1)
plt.plot(N, refl)
plt.xlabel('N')
plt.ylabel('Reflectances (%)')
plt.title('Reflectivity as a function of N')
plt.grid()
plt.show()


##B

n2 = 2.1
n1 = 1.65
ns = 1.52
# =============================================================================
# coef = 337#3.33*10**(-9)
# delta = (wavel/w0)*(0.5*np.pi)
# 
# def func2(N):
#     A_n = ((np.cos(delta))**2-(n1/n2)*(np.sin(delta))**2)**N
#     D_n = ((np.cos(delta))**2-(n2/n1)*(np.sin(delta))**2)**N
#     gam_n =  ((np.cos(delta))*(np.sin(delta))*(1/coef)*(1/n2 + 1/n1))**N
#     beta = (coef*(np.cos(delta))*(np.sin(delta))*(n1+n2))**N
#     return [A_n, D_n, gam_n, beta]
# 
# [A_n, D_n, gam_n, beta] = func2(1)
# 
# E = (A_n - (ns*D_n))/coef
# F = ((coef**2)*(ns*gam_n) - beta)
# G = (A_n + (ns*D_n))/coef
# H = (coef**2)*(ns*gam_n) + beta
# 
# R = (((E**2 + F**2)/(G**2 + H**2))**19)*100
# =============================================================================

ra = (1-n1)/(1+n1)
rb = (n1-n2)/(n1+n2)
rbPr = -1*ra

Refl = np.zeros(30)
phi = np.pi*(wavel/w0)

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R = ((num/den))

ra = -1*rb
num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))

ra = (n2-ns)/(n2+ns)
rb = (ns-1)/(ns+1)
num = ra**2+rb**2+2*ra*rb*np.cos(phi)
den = (ra**2)*(rb**2)+1+2*ra*rb*np.cos(phi)
R += ((num/den))


plt.figure(2)
plt.plot(wavel, (R*100)-2.7)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflectances (%)')
plt.title('Reflectivity as a function of wavelength')
plt.grid()
plt.show()


W = 1
Pow = [1]
T = 0.0003
TPow = []

while W > 1*10**(-4):
    Tpow = T*W
    TPow.append(Tpow)
    W = W - Tpow

time = np.linspace(0,2*10**(-4),len(TPow))

plt.figure(3)
plt.plot(time*(1*10**3), TPow)
plt.xlabel('Time(ms)')
plt.ylabel('Intensity(W)')
plt.title('Intensity of light exiting a cavity')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid()
plt.show()


plt.figure(4)
x = np.linspace(0, 1, 100)
y = 1 - np.exp(-2*x/0.404)
yy = (1-np.cos(13.5*x))/(1-np.cos(1+13.5*x))
plt.plot(x,y)
plt.plot(x, yy)







    