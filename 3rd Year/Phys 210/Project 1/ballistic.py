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
vT = m*g/c  # Terminal velocity


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


def x(t):
    # function for the position of projectile at some time t, in x direction
    xpos = (v0*vT/g)*np.cos(theta)*(1-np.exp(-g*t/vT))
    return xpos


def y(t):
    # function for the position of projectile at some time t, in y direction
    ypos = (vT/g)*(v0*np.sin(theta)+vT)*(1-np.exp(-g*t/vT)) - vT*t
    return ypos


yM = odeint(dvdtAF, y0AF, tAF)  # Intergate the function using odeint and plot
pos = np.where(yM[:, 1] >= 0)  # Position at or above ground for Numerical.

# Answers to the 4 questions: I will be answering them using both the numerical
# and analytical methods, and compare them for accuracy. Professor approved. 
xNume, yNume = yM[:, 0][pos][::45], yM[:, 1][pos][::45]
# I didn't have to do this; I could have simply
# plotted them without assigning, but the professor said he prefers this as it
# is more organized. Trade-off between more lines and organization, but he said
# he'd prefer this and more lines are okay here. Same for anaylitical method.

select = np.where(y(tAF) >= 0)  # Position at or above ground for Analytical.
xAnaly, yAnaly = x(tAF)[select], y(tAF)[select]
dimpact = max(yM[:, 0][pos])  # Distance to impact using Numerical method.
dimpact2 = max(xAnaly)  # Distance to impact using Analytical.
print("Distance to impact: Numerical = {:4.3}m".format(dimpact))
print("Comparing, distance to impact: Analytical = {:4.3}m".format(dimpact2))

ymax = max(yM[:, 1])  # max height using Numerical.
ymax2 = max(yAnaly)  # max heigh using Analytical
print("\nMax Height: Numerical = {:4.4}m".format(ymax))  # Using Numerical
print("Max Height: Analytical = {:4.4}m".format(ymax2))  # Using Analytical

tflight = max(tAF[pos])  # Time of flight using Numerical.
tflight2 = max(tAF[select])  # Time of flight using Analytical
print("\nTime of flight: Numerical = {:4.3}s".format(tflight))
print("Time of flight: Analytical = {:4.3}s".format(tflight2))

# As we can see on the output screen, using very small increments for the
# numerical model, our results are identical to that of the analytical model.
# If out increments were larger (tstep > 0.001), our numerical method would be
# less accurate. I just printed both models for comparison.
vxf = yM[:, 2][pos][-1]  # final x-velocity at impact
vyf = yM[:, 3][pos][-1]  # final y-velocity at impact
vFinal = np.sqrt(vxf**2+vyf**2)  # final impact velocity
print("\nFinal impact velocity: {:4.4}m/s".format(vFinal))

plt.plot(xNume, yNume, 'ro', label='Numerical')
plt.plot(xAnaly, yAnaly, 'b-', label='Analytical')
plt.xlabel('Position(x)')
plt.ylabel('Position(y)')
plt.title('Trajectory of the projectile in the xy-plane')
plt.grid()
plt.legend()
plt.savefig('ballistic.pdf')
plt.show()