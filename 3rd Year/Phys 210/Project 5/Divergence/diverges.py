# Project 5: In this project, we will explore the 2-Dimensional Diffusion
# equation; where concentration of dye molecules in a square dish obeys the
# said equation. We know that the concentration is a function of 3 variables:
# space (x and y), and time. We will experiment with this setup and solve the
# equation numerically. Then we'll plot the concentrations at different times,
# create an animation showing the evolution of C(x,y,t), which is the Diffusion
# equation, as time increases showing the convergence. Also, I will find the
# sums of the number of particles at different times.

# Lastly, I will also experiment with the size of dt and see at which value the
# solution diverges, and how the aforementioned phenomena are affected.

# BONUS: I will be solving the Heat Equation at different temperatures, as well
# as the Method of Relaxation. Please see the script 'bonus.py' for a detailed
# description of the setups.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.animation as animation

# constants
L = 1  # size of "box"
D = 1e-3
N = 100  # number of grid points in one-d

dt = 2.6e-2
dx = L/N

tmax = 10  # max time
steps = int(tmax/dt)+1

animTime = np.arange(0, 100.1, 0.1)  # time interval of 0 to 100, steps of 0.1
# For the animation. There will be 1000 frames, which is what we want.

# choose some times to make plots:
# given in step numbers
plotsteps = np.array([0.001, 0.01, 0.1, 1, 10, 100, 1000])
plotsteps /= dt
plotsteps = plotsteps.astype(int)
# create initial conditions:
C = np.zeros((N, N))  # initial array of zeros

# Set particles in a blob in the center:
C[N//2, N//2] = 10  # central point (origin)
k = dt/dx/dx*D

isum = np.sum(C)  # initial sum

pcol = []  # empty list that I'll be appending to for the animation

fig1, ax1 = plt.subplots()  # here are the different figures for the plots
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

sumsteps = np.array([0, 0.1, 1, 10, 100])  # Will calculate the sums of C at
# these times.
sumsteps /= dt
sumsteps = sumsteps.astype(int)

esumm = []  # empty list to append the sums to.
times = []  # for sqrt(time)
stds = []  # for standard deviation

xvalues = np.linspace(0, 1, 100)  # x values for the gaussian function


def gauss(x, a, sigma):
    """
    This function will be used to find the best fit parameters for sigma, which
    are the standard deviations that I'll be plotting against Square Root of
    Time later on. Takes the x values as parameters, as well as the parameters
    that will be fitted and returned a value of; in this case, the amplitude a
    and sigma. Mean is simply L/2. Returns fit parameter a and sigma.
    """
    return a*np.exp(-(x-L/2)**2/(2*sigma**2))  # here L/2 is the mean


for i in range(steps):
    # Here I do all the points, except the boundaries. np.roll is effecient and
    # reduces the time taken when comparing to nested for loops.
    Cp_xLe = np.roll(C, -1, axis=0)  # roll left
    Cp_xRi = np.roll(C, 1, axis=0)  # roll right
    Cp_yUp = np.roll(C, 1, axis=1)  # roll up
    Cp_yDw = np.roll(C, -1, axis=1)  # roll down
    Cp_xnew = C + k*(Cp_xLe+Cp_xRi-2*C)
    Cp_new = Cp_xnew + k*(Cp_yUp+Cp_yDw-2*C)  # the 2D Diffusion equation

    # Take care of the boundaries here
    Cp_new[0, :] = Cp_new[1, :]  # first row, all columns
    Cp_new[-1, :] = Cp_new[-2, :]  # last row, all columns.
    Cp_new[:, 0] = Cp_new[:, 1]  # first column, all rows
    Cp_new[:, -1] = Cp_new[:, -2]  # last column, all rows

    # swap C and Cp_new so that they don't end up as the same array:
    C, Cp_new = Cp_new, C
    C_slice = C[:, 50]  # slice middle column for standard deviation fits

    if i in plotsteps:  # plots the concetration vs position
        ax1.plot(C[:, N//2], label="t = %g" % (i*dt))

    if i in sumsteps:  # calculates the sum at specified times
        esumm.append(np.sum(C))

    if (i*dt) <= 10:  # extracts standard deviation at specific times, which
        # will be plotted later.
        popt, pcov = curve_fit(gauss, xvalues, C_slice, maxfev=10000)
        times.append(np.sqrt(i*dt))  # appending to preassigned lists
        stds.append(abs(popt[1]))

    if i < len(animTime):  # 100s / 0.1 = 1000 frames.
        pcol.append((ax2.pcolormesh(C.copy(), cmap='jet'),))

esum = np.sum(C)
print("initial and final integrals of concentration:", isum, esum, "\n")
ax1.set_xlabel('Position')
ax1.set_ylabel('Concentration')
ax1.set_title('Concentration vs Position')
ax1.grid()
ax1.legend()
fig1.savefig('concentrations.pdf')

# =============================================================================

theory = np.linspace(0, 10, 1001)
theoryLine = np.sqrt(2*D*theory)  # theory function which we plot below

plt.plot(times[::34], stds[::34], 'bo', markersize=2.5, label='data')
plt.plot(np.sqrt(theory), theoryLine, 'k-', label='Theory')
ax3.set_xlabel('Square Root of Time')
ax3.set_ylabel('Standard Deviation')
ax3.set_title('Evolution of Standard Deviation Vs. Time')
ax3.grid()
ax3.legend()
fig3.savefig('diffusion.pdf')

# =============================================================================

anim = animation.ArtistAnimation(fig2, pcol, interval=30, repeat=False)
ax2.set_xlabel('x(m)')
ax2.set_ylabel('y(m)')
ax2.set_title('Evolution of Diffusion Equation, C(x,y,t) from t=[0,100]s')
anim.save('diverges.mp4')

# =============================================================================

f = open('totals.txt', 'w')  # opening and writing the sums to text file.
for i in range(len(esumm)):
	f.write('At t={}s, sum = {}\n'.format(sumsteps[i]*dt, esumm[i]))
f.close()

# The textfile 'diverges.txt' is attached.
