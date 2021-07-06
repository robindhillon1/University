import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Initials conditions and constants we'll be using to solve the ODE:
c = 0.65  # drag constant in kg/s
g = 9.81  # gravity in m/s^2 y direction. 0 in x direction.
m = 0.1  # mass in kg
v0 = 10  # initial velocity in m/s
theta = 50*np.pi/180  # angle above horizontal in radians
vx0 = v0*np.cos(theta)  # Initial velocity in x direction in m/s
vy0 = v0*np.sin(theta)  # Initial velocity in y direction in m/s


def dvdtAF(yAF, tF):
    # We will be numerically integrating with odeint. The function returns a
    # list of the values of the positions and velocites.
    vx = yAF[2]  # velocity in the x direction (third element of yAF)
    vy = yAF[3]  # velocity in the y direction (fourth element of yAF)
    return np.array([vx, vy, -(c*vx)/m, -g-(c*vy)/m])  # g = 0 in x direction


t0 = 0  # intial time
tf = 1.2  # final time
tstep = 0.001  # This is the for the calulcations at the end. This gives more
# points and hence better accuracy.
tAF = np.arange(t0, tf, tstep)  # creating the 1D time array for calculations
tAFpl = tAF[::45]  # same array as above, but the one to plot the 20 points
y0AF = np.array([0, 0, vx0, vy0])  # initial conditions for position, velocity.
# intial x and y positions are 0 (origin).

yM = odeint(dvdtAF, y0AF, tAF)  # Intergate the function using odeint and plot
pos = np.where(yM[:, 1] >= 0)  # Position at or above ground.
plt.plot(yM[:, 0][pos][::45], yM[:, 1][pos][::45], 'ro', label='Numerical')

vT = m*g/c  # Terminal velocity


def x(t):
    # function for the position of projectile at some time t, in x direction
    xpos = (v0*vT/g)*np.cos(theta)*(1-np.exp(-g*t/vT))
    return xpos


def y(t):
    # function for the position of projectile at some time t, in y direction
    ypos = (vT/g)*(v0*np.sin(theta)+vT)*(1-np.exp(-g*t/vT)) - vT*t
    return ypos


select = np.where(y(tAF) >= 0)  # Selecting the positive y-values.
plt.plot(x(tAF)[select], y(tAF)[select], 'b-', label='Analytical')
plt.xlabel('Position(x)')
plt.ylabel('Position(y)')
plt.title('Trajectory of the projectile in the xy-plane')

# Answers to the 4 questions:

dimpact = max(yM[:, 0][pos])  # Distance to impact.
print("Distance to impact: {:4.3}m".format(dimpact))

ymax = max(yM[:, 1])  # max height
print("Max Height: {:4.4}m".format(ymax))

tflight = max(tAF[pos])  # Time of flight
print("Time of flight: {:4.3}s".format(tflight))

vxf = yM[:, 2][pos][-1]  # final x-velocity at impact
vyf = yM[:, 3][pos][-1]  # final y-velocity at impact
vFinal = np.sqrt(vxf**2+vyf**2)  # final impact velocity
print("Final impact velocity: {:4.4}m/s".format(vFinal))

# BONUS

plt.axhline(y=ymax, color='k', linestyle='--', label='Max Height')
plt.axvline(x=dimpact, color='m', linestyle='--', label='Distance to impact')
plt.grid()
plt.legend()
plt.savefig('ballistic.pdf')
plt.show()
