import numpy as np
import pylab as pl
from PIL import Image
import scipy.special as sc

import matplotlib.pyplot as plt

# read .tif file to a numpy array. Check for saturation:
im10 = np.array(Image.open('D:/Uni/4th Year/Phys 408/2020/Lab 2/11.60mm.tif'))
print(im10.max())

# Average a region of the data and plot a cross section
a10 = im10[:, 320:330].mean(axis=1)
a10 = a10 - np.min(im10)

x = np.arange(len(a10))

I0 = 10000000*0.36
scale = 0.04  # 11.60 = 0.04, 11.72 =
x0 = 382
y1 = I0*np.power(np.sinc(scale*(x-x0)), 2)

pl.semilogy(x, a10, label='Data')
pl.semilogy(x, y1/6, label='Fraunhofer')  # 11.60 = y1/6, 11.72 = y1/16
pl.ylim(1e2, 1e6)
pl.xlim(250, 500)

xHor = np.linspace(344, 420, 200)
horizLine = np.array([np.max(a10) for i in range(len(xHor))])
# plt.hlines(np.max(a10), xmin=344, xmax=420, colors = 'k', linestyles ='dashed', linewidth=2, zorder = 2)
plt.plot(xHor, horizLine, 'k--', label='Saturation')
plt.xlabel('Scanner Pixel Number')
plt.ylabel('Brightness')
plt.title('Fraunhofer Diffraction model fitted to data for slit width, w= 0.28mm')
plt.legend()
pl.show()

x = (x-382)*0.0254/600
delV = 3  # 11.60 = 0.94, 11.72 =
w = 0.004  # 11.60 = 1, 11.72 =
z = x/w

v1 = -(z+0.5)*delV
v2 = -(z-0.5)*delV

# c = np.cos(np.pi*z**2/2)
# C = np.cumsum(c)

# s = np.sin(np.pi*z**2/2)
# S = np.cumsum(s)

S = sc.fresnel(v2)[0] - sc.fresnel(v1)[0]
C = sc.fresnel(v2)[1] - sc.fresnel(v1)[1]

I = I0**0.8*(C**2+S**2)
xHor = np.linspace(-0.5, 0.5, 200)
horizLine = np.array([np.max(a10) for i in range(len(xHor))])

plt.figure(2)
plt.semilogy(z, a10, label='Data')
plt.semilogy(z, I/1.5, label='Fresnel')  # 11.6 = z/9, 11.72 = z/4.5 & I*1.5
plt.plot(xHor, horizLine, 'k--', label='Saturation')
plt.ylabel('Brightness')
plt.xlabel('z = x/w')
plt.title('Fresnel Diffraction model fitted to data')
plt.xlim([-2, 2])
plt.ylim([1e2, 1e6])
plt.legend()

# Here I've plotted the Cornu Spiral:

# =============================================================================
# t = np.linspace(-7, 7, 1000)
# S, C = sc.fresnel(x/100)
#
# plt.plot(C, S)
# =============================================================================


# =============================================================================
# dx=x[1]-x[0]
# c=np.cos(np.pi*x**2/2.)
#
# C=np.cumsum(c)*dx
#
# s=np.sin(np.pi*x**2/2.)
#
# S=np.cumsum(s)*dx
#
# I = S**2+C**2
# =============================================================================


# =============================================================================
# w=0.00188
# lamb = 632.8e-9
# R = 2.7
#
# x=np.arange(len(a10))/(lamb*R)
#
# delV = w*np.sqrt(2/(lamb*R))
#
# a = np.linspace(0, delV, len(a10))
#
# nF = delV/np.sqrt(8)
#
# def lambd(u, X):
#     return np.exp(-1*1j*np.pi*(X-u)**2)
#
#
# call = np.asanyarray([(np.abs(integrate.quad(lambd, -nF,nF, args=(value,))))**2 for value in x])
#
# =============================================================================


# =============================================================================
# func = lambda u: np.cos(np.pi*u**2/2)
# call = integrate.quad(func, v1, v2)
#
# func2 = lambda q: np.cos(np.pi*q**2/2)
# call2 = integrate.quad(func2, v1, v2)
# =============================================================================


# =============================================================================
# w = [0.16,0.28,0.43,0.58,0.68,0.78,0.88,1.08,1.88]
# R = 2.7*10**3
# lamb = 6.328*10**(-4)
# sq = np.sqrt(2/(R*lamb))
#
# delV = [x*sq for x in w]
# print(delV)
# =============================================================================
