# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 21:05:20 2018

@author: robin
"""
import numpy as np
import matplotlib.pyplot as plt

f = open('clusterUBV.txt','r')                    
file_contents = f.read()
print(file_contents)

data = np.genfromtxt('clusterUBV.txt',skip_footer=6)
data2 = np.genfromtxt('UBV_intrinsic_colour.txt')
data3 = np.genfromtxt('UBV_intrinsic_ms.txt')
data4 = np.genfromtxt('isochrone_3.16e7yrs.txt')
data5 = np.genfromtxt('isochrone_1.00e8yrs.txt')
data6 = np.genfromtxt('isochrone_3.16e8yrs.txt')

ID = data[:,0]
source = data[:,1]
V = data[:,2]
BV = data[:,3]
UB = data[:,4]
N = data[:,5]

BVint = data2[:,0]
UBint = data2[:,1]

Mv = data3[:,0]
BVo = data3[:,1]

Biso1 = data4[:,0]
Viso1 = data4[:,1]

Biso2 = data5[:,0]
Viso2 = data5[:,1]

Biso3 = data6[:,0]
Viso3 = data6[:,1]

plt.figure(1)
plt.clf()
plt.scatter(BV-0.037,V-0.111,s=5,marker=".",color="red",label='Cluster CMD')
plt.plot(BVo,Mv+5.65,color="black",label='Fiducial Main Sequence')#distance using distance modulus
plt.xlabel('B-V')
plt.ylabel('V')
plt.title('Color Magnitude Diagram of Cluster')
plt.gca().invert_yaxis()
plt.grid(True)
plt.legend()

#rangBV = BV[(BV>=0) & (BV<=0.721828)]
#modelV = 5.65*rangBV+6.5
#linemodel = np.where(V[np.where(rangBV)]>modelV)
#plott9 = plt.plot(rangBV-0.037,modelV-0.111,color="blue",linewidth = 0.5)

plt.figure(2)
plt.clf()

def F(x):
    m = 5.65
    b = 6.5
    return m*x+b

x = BV
y = V

ind = np.where((x>=0)&(x<=0.75))
y2 = y[ind]
x2 = x[ind]

binrange = np.where(y2<F(x2))
bifq = len(binrange[0])/len(x2)
print('\nBinary Frequency: '+str(bifq)+'\n')

test = np.linspace(0,0.75,100)

plt.scatter(x,y,s=5,marker='.',color='red',label='Cluster CMD')
plt.plot(test,F(test),'k-',lw=2,label='Model Range')

#plt.plot(x[ind],y[ind],'y.')
plt.plot(x2[binrange],y2[binrange],'b.',label='Binaries')
plt.xlabel('B-V')
plt.ylabel('V')
plt.title('Color Magnitude Diagram of Cluster')
plt.gca().invert_yaxis()
plt.grid(True)
plt.legend()

plt.figure(3)
plt.clf()
plt.scatter(BV,UB,s=5,marker=".",color="red",label='Cluster CCD')
plt.xlabel('B-V')
plt.ylabel('U-B')
plt.title('Color-Color Diagram of Cluster')
plt.plot(BVint+0.037,UBint+(0.72*0.037),label='UBV Intrinsic Colour')
plt.gca().invert_yaxis()

d = 10*(10**(0.2*5.65))
print('\nDistance is(pc): '+str(d)+'\n')
plt.grid(True)
plt.legend()

plt.figure(4)
plt.clf()
plt.scatter(BV,V,s=5,marker=".",color="red",label='Cluster CMD')
plt.scatter(BV-0.037,V-0.111,s=5,marker=".",color="blue",label='De-reddened, extinction-corrected CMD')
plt.xlabel('B-V')
plt.ylabel('V')
plt.title('Color Magnitude Diagram of Cluster')
plt.gca().invert_yaxis()
plt.grid(True)
plt.legend()

plt.figure(5)
plt.clf()
plt.scatter(BV-0.037,V-0.111,s=5,marker=".",color="blue",label='Cluster CMD')
plt.plot((Biso1-Viso1),Viso1+5.65,color="green",lw = 0.5,label='Isochrone: 3.16e7yrs')
plt.plot((Biso2-Viso2),Viso2+5.65,color="red",lw = 0.5,label='Isochrone: 1.00e8yrs')
plt.plot((Biso3-Viso3),Viso3+5.65,color="purple",lw = 0.5,label='Isochrone: 3.16e8yrs')
plt.xlabel('B-V')
plt.ylabel('V')
plt.title('Color Magnitude Diagram of Cluster')
plt.gca().invert_yaxis()
plt.grid(True)
plt.legend()