# Note: I use Spyder on my own laptop. Did not use the lab computer.
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

L = float(input('Please enter a value for the inductance, L: '))
C = float(input('Please enter a value for the capacitance, C: '))

w0 = 1/np.sqrt(L*C)  # omega_0
freq = np.linspace(0.1*w0/(2*np.pi), w0/np.pi, 100)  # given frequency range
m = [5, 1, 0.2]  # given values of m
yAf = np.array([0, 0])  # initial charge and current.
tAf = np.arange(0, 60, np.pi/w0/20)  # time array
V0 = 1  # given value of v0

R1 = 2*m[0]*np.sqrt(L/C)  # value of resistance at m = 5
R2 = 2*m[1]*np.sqrt(L/C)  # value of resistance at m = 1
R3 = 2*m[2]*np.sqrt(L/C)  # value of resistance at m = 0.2
R_all = [R1, R2, R3]


def dIdtAF(yAf, tAf, w, R):
    # This function takes the array with intial charge and current, as well as
    # the time array, the omega (angular frequency), and resistance, as
    # parameters. We will use these to solve the differential equation below,
    # dI, by using the odeint function.
    q = yAf[0]  # charge
    I = yAf[1]  # current = dq/dt
    dI = V0*np.cos(w*tAf)/L - q/(L*C) - I*R/L  # the differential equation
    return [I, dI]  # return the current and the result(s) of the DE above


res1 = []  # empty list for each resistor/resistance which I'll be appending to
res2 = []  # could have preassigned them using np.zeros() or such, but the prof
res3 = []  # said the change is miniscule and unnecessary in this case.

for i in R_all:  # iterating over the 3 resistances for respective m
    for j in freq:  # iterating over all the frequencies
        w = 2*np.pi*j  # angular frequency, omega
        call = odeint(dIdtAF, yAf, tAf, args=(w, i))  # calling odeint to
        # solve differential equation above
        y = call[:, 1]  # for readablility, assigning y to the result returned.
        imax = np.max(y[200:])  # taking the max of the values after steady
        # state reached. Using the values after 200 because by doing so, we
        # neglect the first few values where the solution might not be steady.
        if i == R_all[0]:
            res1.append(imax)  # appending to corresponding resistance
        if i == R_all[1]:
            res2.append(imax)
        if i == R_all[2]:
            res3.append(imax)

plt.figure(1)  # creating a figure
plt.subplot(311)  # subplots for each resistance
plt.plot(freq, res1)  # plotting each resonance for corresponding R with freq
plt.xlabel("Frequency, Hz")
plt.ylabel("Current, I(t)")
plt.title("L={}H,C={}F,R={}Ohms".format(L, C, R_all[0]))
# Here I assumed the units are in C, F, and Ohms. We were not told if they were
# mH, nH etc. which I why I made this general assumption.
plt.grid()

plt.subplot(312)
plt.plot(freq, res2)
plt.xlabel("Frequency, Hz")
plt.ylabel("Current, I(t)")
plt.title("L={}H,C={}F,R={}Ohms".format(L, C, R_all[1]))
plt.grid()

plt.subplot(313)
plt.plot(freq, res3)
plt.xlabel("Frequency, Hz")
plt.ylabel("Current, I(t)")
plt.title("L={}H,C={}F,R={}Ohms".format(L, C, R_all[2]))
plt.tight_layout()  # allows us to see the titles clearly of each subplot
plt.grid()
plt.savefig('resonance.pdf')
plt.show()

# =============================================================================
# Here, we answer question 3. For readability, use another for loop. Could have
# done it all in one nested loop, but again, that looks messy and decreases
# readability, as the TAs said. Hence, I'll stick to this as it's cleaner.

I_1 = []  # empty lists for each resistance that I'll be appending to.
I_2 = []
I_3 = []
for i in R_all:  # iterating over the resistances
    call = odeint(dIdtAF, yAf, tAf, args=(w0, i))  # calling odeint to solve
    # DE. However, now we keep w0 fixed, instead of using w for each frequency
    # The steps below are the same as before.
    y = call[:, 1]
    if i == R_all[0]:
        I_1.append(y)
    if i == R_all[1]:
        I_2.append(y)
    if i == R_all[2]:
        I_3.append(y)

plt.figure(2)
plt.subplot(311)
plt.plot(tAf, I_1[0])  # Instead of frequency, now we plot time on the x-axis.
plt.xlabel("Time(s)")
plt.ylabel("Current, I(t)")
plt.title("L={}H,C={}F,R={}Ohms".format(L, C, R_all[0]))
plt.grid()

plt.subplot(312)
plt.plot(tAf, I_2[0])
plt.xlabel("Time(s)")
plt.ylabel("Current, I(t)")
plt.title("L={}H,C={}F,R={}Ohms".format(L, C, R_all[1]))
plt.grid()

plt.subplot(313)
plt.plot(tAf, I_3[0])
plt.xlabel("Time(s)")
plt.ylabel("Current, I(t)")
plt.title("L={}H,C={}F,R={}Ohms".format(L, C, R_all[2]))
plt.tight_layout()  # allows us to see the titles clearly of each subplot.
plt.grid()
plt.savefig('transients.pdf')
plt.show()
