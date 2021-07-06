# Project 4:
# This code will generate the animation of collisions of (400) gas particles.
# Then, the plots/histrograms of the probability distributions of speed and
# energy will be generated according to the Maxwell-Boltzmann Distribution of
# Speed and the Boltzmann Distribution provided, respectively. Once these are
# plotted, I will fit the data using curve_fit to get an approximate value for
# the temperature (at which the particles have been interacting at). It will be
# obtained from fitting the Maxwell-Boltzmann Distribution of Speed, and the
# resulting temperature will be used to plot the analytical solution of the
# Boltzmann Distribution, as required.

# Bonus: A seperate file, called Bonus, contains most of the similar code but
# also contains code that will calculate the temperature from fitting the
# Boltzmann Distribution function and compare it with the temperature obtained
# from fitting the Maxwell-Boltzmann Distribution of speed. Note: we are only
# supposed to plot the analytical solution to Boltzmann Distribution using the
# temperature obtained from the first fit. We are not asked to fit the
# Boltzmann Distribution. Nonetheless, I will fit them both and compare them in
# the Bonus.py file; I will compare the temperature obtained from fitting/based
# on the velocity and the tempeature based on the energy.

# Also, for pep8, I use the lab computers as well as the following website to
# check: http://pep8online.com/

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from itertools import combinations
from scipy.optimize import curve_fit

pRadius = 0.0015  # radius of particle, in meters, m
pMass = 2.672e-26  # mass of an O2 molecule in kg

k_B = 1.38e-23   # Boltzmann constant

npoint = 400  # number of particles
nframe = 1000  # number of steps
xmin, xmax, ymin, ymax = 0, 1, 0, 1  # max and min values of the axes
Dt = 0.00002  # time step
fig, ax = plt.subplots()
plt.xlim(xmin, xmax)  # aforementioned max and min values used here
plt.ylim(ymin, ymax)
plt.title('Collisions of particles')
plt.xlabel('x')
plt.ylabel('y')


def update_point(num):
    """
    This function takes in the nframe variable assinged above, which is the
    number of steps, as the parameter, and updates the preassigned velocities
    vx and vy. The collisions are elastic for the purpose of this project and
    every time two particles in a gas collide, the new velocities will be
    calculated according to the conservation of energy and momentum. With these
    new velocities we will create an animation that shows the movement and
    collisions of particles. Then, the velocities will be plotted below.
    """
    global x, y, vx, vy  # assings these variables to be global. So they can be
    # used outside the function, and are not local any more.
    indx = np.where((x < xmin) | (x > xmax))
    indy = np.where((y < ymin) | (y > ymax))
    vx[indx] = -vx[indx]
    vy[indy] = -vy[indy]
    xx = np.asarray(list(combinations(x, 2)))  # pair indices
    yy = np.asarray(list(combinations(y, 2)))
    dd = (xx[:, 0]-xx[:, 1])**2+(yy[:, 0]-yy[:, 1])**2  # distance squared bet-
    # ween the pairs.
    collision = np.where(dd <= ((2*pRadius)**2))  # condition that finds where
    # the distance between pairs is equal to or less than 2 times the radius
    # of each particle, allowing a collision.
    pInd2 = np.asarray(list(combinations(pInd, 2)))  # using combinations,
    # creates pairs of the colliding particles.
    collInd = pInd2[collision]  # the indices of the colliding particles

    i, j = [collInd[:, 0], collInd[:, 1]]  # particle(s) 1(i) and 2(j)
    # for readability, I will break up the equation for the new velocities.
    num = ((vx[i]-vx[j])*(x[i]-x[j])+(vy[i]-vy[j])*(y[i]-y[j]))  # numerator &
    deno = ((x[i]-x[j])**2+(y[i]-y[j])**2)  # denominator of the equation(s)

    dx = Dt*vx  # As mentioned by Colby, we call the dx and dy before we
    dy = Dt*vy  # calculate the new velocities; account for initial velocities
    # The 4 lines of code below calculate the new velocities after collision(s)
    vx[i] -= (num/deno)*(x[i]-x[j])
    vx[j] -= (num/deno)*(x[j]-x[i])
    vy[i] -= (num/deno)*(y[i]-y[j])
    vy[j] -= (num/deno)*(y[j]-y[i])

    x = x + dx
    y = y + dy
    data = np.stack((x, y), axis=-1)
    im.set_offsets(data)


pInd = np.arange(0, npoint)  # indices of the 400 particles

x = np.random.random(npoint)  # random 400 numbers
y = np.random.random(npoint)

c = np.empty_like(x, dtype='str')  # Returns new array that has the same shape
# and type as the given array; in this case, x.
c[np.where(x > 0.5)] = 'red'  # Assigning colors to the points, red for x>0.5
c[np.where(x <= 0.5)] = 'blue'  # Blue for x<=0.5; about half red, half blue

vx = -500.*np.ones(npoint)
vy = np.zeros(npoint)
vx[np.where(x <= 0.5)] = -vx[np.where(x <= 0.5)]
im = ax.scatter(x, y, c=c)
im.set_sizes([10])

animation = animation.FuncAnimation(fig, update_point, nframe, interval=10,
                                    repeat=False)
# calling and creating the animation
animation.save('collisions.mp4')  # saving the animation

vel = np.sqrt(vx**2+vy**2)  # velocity of particles in the end for use in fit
energy = 0.5*pMass*(vel**2)  # energy of particles in the end


def MBvelDistri(vel, T):
    """
    The purpose of this function is to ultimately fit the Maxwell-Boltzmann
    velocity distribution, hence the function name, and obtain an approximate
    value for the temperature. Takes the velocities as parameter as well as the
    temperature T, which it will fit and return a value for when we use
    curve_fit.
    """
    f_distri = ((pMass*vel)/(k_B*T))*np.exp(-0.5*pMass*vel**2/(k_B*T))
    return f_distri  # return data corresponding to equation


plt.clf()
fig2 = plt.figure()
plt.subplot(211)
nv, binsv, patchesv = plt.hist(vel, bins=50, density=True, facecolor='y',
                               alpha=0.75, label='Numerical Values (v)')
# Here I plotted the histogram with 50 bins of the velocity/speed values.
# The Density = True normalizes the data, facecolor assings to color to bins,
# and alpha decides the opcaity. 1= opaque, 0 = transparent.
midPointsV = np.array(binsv[:-1]) + 0.5*(binsv[1]-binsv[0])  # midpoints of the
# (velocity) bins, but excluding the last bin.
Temp1, pcov = curve_fit(MBvelDistri, midPointsV, nv, p0=[250])  # curve fitting
# the function above, and obtaining the temperature approximation. Plots below.

plt.plot(midPointsV, MBvelDistri(midPointsV, Temp1), 'c-', label='Fit of f(v)')
plt.title('Maxwell_Boltzmann Distribution of Speed')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Probability Distribution, f(v)')
plt.grid()
plt.legend()

plt.subplot(212)
nE, binsE, patchesE = plt.hist(energy, bins=50, density=True, facecolor='g',
                               alpha=0.75, label='Numerical values (E)')
midPointsE = np.array(binsE[:-1]) + 0.5*(binsE[1]-binsE[0])  # midpoints of the
# (energy) bins, but excluding the last bin.

g_distri = (1/(k_B*Temp1))*(np.exp(-1*midPointsE/(k_B*Temp1)))
# As mentioned in the instructions, we simply plot the analytical solution g(E)
# We don't fit it. I will nonetheless fit it in the bonus and compare the Temps
# Plots below.
plt.plot(midPointsE, g_distri, 'k-', label='Analytical g(E)')
plt.title('Boltzmann Distribution')
plt.xlabel('Energy (J)')
plt.ylabel('Probability Distribution, g(E)')
plt.grid()
plt.legend()
plt.tight_layout()  # allows to see titles and labels clearly of each subplot.
plt.savefig('distributions.pdf')  # saving as pdf file.
plt.show()

f = open('collisions.txt', 'w')  # opening and allowing to write to file.
f.write('Temperature obtained from fitting f(v): {} K'.format(Temp1[0]))
f.close()  # closing the file after writing
