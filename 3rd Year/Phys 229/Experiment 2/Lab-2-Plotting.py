#!/usr/bin/env python
# coding: utf-8

# In[32]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import scipy.constants


# ## Theory
# 
# Here we define the theoretical models to integrate and fit

# In[33]:


# Constants:

baseTemp = 20 + 273.15 # K
heatCapacity = 900 # J/(kgK)
sig = scipy.constants.sigma # Stefan Boltzmann constant  5.67e−8 W*m−2*K−4


# In[34]:


sig


# In[35]:


def Temp_ODE(time, Temp_0, mass, A, D, e, heatCapacity=heatCapacity):
    # mass is mass of rod in kg
    # A is surface area
    # D is diameter
    # e is emissivity of material, this is the quantity we are trying to fit
    # 
    def theory_dTdt(Temp, t):
        h = 1.32*((Temp-baseTemp)/D)**(1/4)
        result = (1/(mass*heatCapacity)) * ((-A*e*sig*(Temp**4 - baseTemp**4))-(h*A*(Temp - baseTemp)))
        return result
    
    TempResult = odeint(theory_dTdt, Temp_0, time, mxordn=20, mxords=20, h0=0.001)
    return TempResult[:,0]
    


# In[36]:


# Following is a wrapper for scipy curve_fit where we define a small function for the Temp_ODE above with 
# all relevant parameters inserted
def fit_Temp_ODE(time, Temp, dTemp, mass, A, D, heatCapacity=heatCapacity):
    def fitfn(t_, e):
        return Temp_ODE(t_, Temp[0], mass, A, D, e, heatCapacity)
    
    e_fit, e_cov = curve_fit(fitfn, time, Temp, p0=(0), sigma=dTemp, absolute_sigma=True, 
                             maxfev=1000, bounds=(0,1), method='dogbox') #bounds=np.array([(0, -np.inf), (1, np.inf)])) 
    return (e_fit, e_cov)


# ## Experiment
# Here we have all our data and experimental measurements

# In[37]:


def area(h, d):
    return h*np.pi*d+np.pi*(d/2)**2


# In[38]:


# Polished:

mP = 0.427 # kg +/- 0.001 kg  Mass
lP = 0.304 # m +/- 0.001 m  Length
dP = 0.027 # m +/- 0.001 m  Diameter
aP = area(lP, dP)

# Rough:

mR = 0.382 # kg +/- 0.001 kg  Mass
lR = 0.272 # m +/- 0.001 m  Length
dR = 0.027 # m +/- 0.001 m  Diameter
aR = area(lR, dR)

# Lacquered

mL = 0.429 # kg +/- 0.001 kg  Mass
lL = 0.305 # m +/- 0.001 m  Length
dL = 0.027 # m +/- 0.001 m  Diameter
aL = area(lL, dL)


# In[39]:


data = pd.read_csv('day1.csv').iloc[1:]


# In[40]:


data


# In[ ]:





# In[41]:


time = data['Time (min)'].values
uTime = data['ut'].values
TPolished = data['Polished'].values + 273.15
TPolished[0] = TPolished[0]+2
TRough = data['Rough'].values + 273.15
TLacquered = data['Lacquered'].values + 273.15
uTemp = data['dTemp (dec C)'].values*2
#uTemp = np.ones(len(TPolished))


# Because we collected more data earlier, if we were to use all our data to fit the curve it would skew the curve fit algorithm to fit the first section rather than the entire curve, thus we balance and only use the max time before fitting

# In[42]:


time2 = np.arange(0, 4081, 120)

ll = []
for i, val in enumerate(time):
    if val in time2:
        ll.append(i)
        
ll = np.array(ll)

time2 = data['Time (min)'].values[ll]
uTime2 = data['ut'].values[ll]
TPolished2 = data['Polished'].values[ll] + 273.15
TRough2 = data['Rough'].values[ll] + 273.15
TLacquered2 = data['Lacquered'].values[ll] + 273.15
uTemp2 = data['dTemp (dec C)'].values[ll]*2


# ## Data Analysis
# Here we use the data and theory to analyze and fit parameters

# In[43]:


lP*np.pi*(dP/2)**2/(lR*np.pi*(dR/2)**2)


# In[44]:


test_data1 = Temp_ODE(time, TPolished[0], mP, aP, dP, e=0, heatCapacity=heatCapacity)
test_data2 = Temp_ODE(time, TRough[0], mR, aR, dR, e=0, heatCapacity=heatCapacity)
test_data3 = Temp_ODE(time, TLacquered[0], mL, aL, dL, e=0, heatCapacity=heatCapacity)


# The below cell uses my fitting ODE function to fit ODES, both with the new sparse evenly weighted dataset as well as the old skewed dataset

# In[45]:


# new sparse dataset, evenly weighted on all ends
e_P, e_P_cov = fit_Temp_ODE(time2, TPolished2, uTemp2, mP, aP, dP)
e_R, e_R_cov = fit_Temp_ODE(time2, TRough2, uTemp2, mR, aR, dR)
e_L, e_L_cov = fit_Temp_ODE(time2, TLacquered2, uTemp2, mL, aL, dL)

# old skewed dataset, with a lot of initial points but not a lot at the end. currently commented out.
#e_P, e_P_cov = fit_Temp_ODE(time, TPolished, uTemp, mP, aP, dP)
#e_R, e_R_cov = fit_Temp_ODE(time, TRough, uTemp, mR, aR, dR)
#e_L, e_L_cov = fit_Temp_ODE(time, TLacquered, uTemp, mL, aL, dL)


# In[29]:


#print(e_P)
#print(e_R)
#print(e_L)


# In[30]:


# Covariances are 
#print(e_P_cov)
#print(e_R_cov)
#print(e_L_cov)


# In[31]:


fig = plt.figure(figsize=(15,10))

plt.errorbar(time, TPolished, xerr=uTime, yerr=uTemp, c='r', label='Polished', capsize=2, fmt='.',ms=3,)
plt.errorbar(time, TRough, xerr=uTime, yerr=uTemp, c='g', label='Rough', capsize=2, fmt='.', ms=3,)
plt.errorbar(time, TLacquered, xerr=uTime, yerr=uTemp, c='b', label='Lacquered', capsize=2, fmt='.', ms=3,)
#plt.plot(time, test_data1, label='P')
#plt.plot(time, test_data2, label='R')
#plt.plot(time, test_data3, label='L')

plt.plot(time, Temp_ODE(time, TPolished[0], mP, aP, dP, e=e_P[0], heatCapacity=heatCapacity,), c='r', ls=':',
         label=('Polished Fit: e={:.2f} +/- {:.2f}'.format(e_P[0], np.sqrt(e_P_cov[0, 0]))))
plt.plot(time, Temp_ODE(time, TRough[0], mR, aR, dR, e=e_R[0], heatCapacity=heatCapacity,), c='g', ls=':',
         label=('Rough Fit: e={:.2f} +/- {:.2f}'.format(e_R[0], np.sqrt(e_R_cov[0, 0]))))
plt.plot(time, Temp_ODE(time, TLacquered[0], mL, aL, dL, e=e_L[0], heatCapacity=heatCapacity,), c='b', ls=':',
         label=('Lacquered Fit: e={:.2f} +/- {:.2f}'.format(e_L[0], np.sqrt(e_L_cov[0, 0]))))

plt.grid()
plt.legend()

plt.title('Temperature and Emissivity of Rods')
plt.ylabel('Temperature (K)')
plt.xlabel('Time (s)')

#plt.show()
plt.savefig('Fit.pdf', dpi=800)
