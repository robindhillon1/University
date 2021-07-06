# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 20:10:42 2018

@author: robin
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

fname = 'HW2_Q1_Full_Sky_Catalogue.csv'
#load the data from the file
data = np.genfromtxt(fname, comments='#')

h = data[:,0]
m = data[:,1]
s = data[:,2]

deg = data[:,3]
arcmin = data[:,4]
arcsec = data[:,5] 
mv = data[:,6]

RAs = (h*3600)+(m*60)+s
dec = (deg*3600)+(arcmin*60)+arcsec

RA = (RAs/(3600*24))*360 #
declin = dec/3600

visb = np.where(mv<=6)
# np.shape(visb) gives: Number of visible stars: 5074
visbh = h[visb]
visbm = m[visb]
visbs = s[visb]

visbdeg = deg[visb]
visbarcmin = arcmin[visb]
visbarcsec = arcsec[visb]

visbRAs = (visbh*3600)+(visbm*60)+visbs
visbdec = (visbdeg*3600)+(visbarcmin*60)+visbarcsec

vRA = (visbRAs/(3600*24))*360#degrees
vdeclin = visbdec/3600

#Circumpolar stars are stars that never set: 90-L. Latitude of Vancouver is 49,
#so the circumpolar stars are: 90-49=41.  

circd = vdeclin[vdeclin>=41]
circRA = vRA[vdeclin>=41]
#Number of circumpolar visible stars: 860

fig = plt.figure(num=1, figsize=(10, 5))
plt.clf()
ax = fig.add_subplot(111, projection="mollweide")

alls = ax.scatter(RA*np.pi/180 -np.pi, declin*np.pi/180, s=0.01, marker=".", color="red")
visbs = ax.scatter(vRA*np.pi/180 -np.pi, vdeclin*np.pi/180, s=1, marker="*", color="black")
circs = ax.scatter(circRA*np.pi/180 -np.pi, circd*np.pi/180, s=1, marker="o", color="blue")

lgnd = ax.legend((alls,visbs,circs),('All Stars','Visible Stars','Visible Circumpolar Stars'),scatterpoints=1,fontsize=8)
lgnd.legendHandles[0]._sizes=[6]
lgnd.legendHandles[1]._sizes=[6]
lgnd.legendHandles[2]._sizes=[6]

ax.grid(True)
plt.xlabel('Right-Ascension(360 degrees = 24hr)')
plt.ylabel('Declination(Degrees)')
plt.title('RA vs Dec; Positions of Stars')
plt.show()
