# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 13:02:43 2018

@author: robin
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('standard_solar_model.txt',skip_header=20,skip_footer=2)

M_Ms = data[:,0]
R_Rs = data[:,1]
T_data = data[:,2]
rho_data = data[:,3]
P_data = data[:,4]

Rs = 6.96*10**(8)
Msun = 1.989*10**30
rho_c = (3*Msun)/(np.pi*(Rs**3))
rho_cCGS = ((3*Msun)/(np.pi*(Rs**3)))/1000
print('\nCentral Density(kg/m^3) using Linear Model: '+"{:.3e}".format(rho_c))

G = 6.67*10**(-11)

P_c = (5*G*Msun**2)/(4*np.pi*(Rs**4))
print('Central Pressure(Pa) using Linear Model: '+"{:.3e}".format(P_c))

u = 0.6
k = 1.3807*10**(-23)
m_h = 1.67*10**(-27)

T_c = (5*G*Msun*u*m_h)/(12*Rs*k)
print('Central Temperature(K) using Linear Model: '+"{:.3e}".format(T_c))

rho_cSolM = 1.505*10**5
P_cSolM = 2.338*10**16
T_cSolM = 1.548*10**7

print('\nCentral Density(kg/m^3) using Standard Solar Model: '+"{:.3e}".format(rho_cSolM))
print('Central Pressure(Pa) using Standard Solar Model: '+"{:.3e}".format(P_cSolM))
print('Central Temperature(K) using Standard Solar Model: '+"{:.3e}".format(T_cSolM)+'\n')

def rho(r):
    return rho_c*(1-(r/Rs))

def M(r):
    return 4*np.pi*rho_c*(((r**3)/3)-(r**4)/(4*Rs))

def P(r):
    return P_c - (((4*np.pi*G*(rho_c)**2)/Rs)*(((Rs*r**2)/6)-((7*r**3)/36)+(r**4/(16*Rs))))

def T(r):
    return (P(r)*u*m_h)/(rho(r)*k)

r = np.linspace(0,Rs,1284)

plt.figure(1)
plt.clf()
plt.plot(R_Rs,M_Ms,color="black",label='Standard Solar Model')
plt.plot(r/Rs,M(r)/Msun,color="red",label='Linear Solar Model')
plt.title('Mass')
plt.xlabel('R/Rsun')
plt.ylabel('M/Msun')
plt.grid(True)
plt.legend()

plt.figure(2)
plt.clf()
plt.plot(R_Rs,rho_data/(rho_cSolM/1000),color="black",label='Standard Solar Model')
plt.plot(r/Rs,rho(r)/rho_c,color="red",label='Linear Solar Model')
plt.title('Density')
plt.xlabel('R/Rsun')
plt.ylabel('rho/rho_c')
plt.grid(True)
plt.legend()

plt.figure(3)
plt.clf()
plt.plot(R_Rs,P_data/(10*P_cSolM),color="black",label='Standard Solar Model')
plt.plot(r/Rs,P(r)/P_c,color="red",label='Linear Solar Model')
plt.title('Pressure')
plt.xlabel('R/Rsun')
plt.ylabel('P/P_c') 
plt.grid(True)
plt.legend()

plt.figure(4)
plt.clf()
plt.plot(R_Rs,T_data/T_cSolM,color="black",label='Standard Solar Model')
plt.plot(r/Rs,T(r)/T_c,color="red",label='Linear Solar Model')
plt.title('Temperature')
plt.xlabel('R/Rsun')
plt.ylabel('T/T_c')
plt.grid(True)
plt.legend()
