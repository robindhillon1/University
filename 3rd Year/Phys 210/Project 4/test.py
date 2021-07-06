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
xmin, xmax, ymin, ymax = 0, 1, 0, 1
Dt = 0.00002  # time step
fig, ax = plt.subplots()
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)


def update_point(num):  # num is the nframe as parameter
	global x, y, vx, vy
	#print(num)
	indx = np.where((x < xmin) | (x > xmax))
	indy = np.where((y < ymin) | (y > ymax))
	vx[indx] = -vx[indx]
	vy[indy] = -vy[indy]
	xx = np.asarray(list(combinations(x, 2)))
	yy = np.asarray(list(combinations(y, 2)))
	dd = (xx[:, 0]-xx[:, 1])**2+(yy[:, 0]-yy[:, 1])**2  # distance sqayred between pairs

	collision = np.where(dd <= ((2*pRadius)**2))

	pInd2 = np.asarray(list(combinations(pInd, 2)))

	collInd = pInd2[collision]

	i,j = [collInd[:,0], collInd[:,1]]
	
	x1 = x[i]
	x2 = x[j]
	y1 = y[i]
	y2 = y[j]
	v1x = vx[i]
	v2x = vx[j]
	v1y = vy[i]
	v2y = vy[j]
	
	newv = ((v1x-v2x)*(x1-x2) + (v1y-v2y)*(y1-y2))/((x2-x1)**2 + (y2-y1)**2)
	newx = newv*(x1-x2)
	newy = newv*(y1-y2)
	
	dx = Dt*vx
	dy = Dt*vy
	
	vx[i] -= newx
	vx[j] += newx
	vy[i] -= newy
	vy[j] += newy
	
	
	x = x + dx
	y = y + dy
	data = np.stack((x, y), axis=-1)
	im.set_offsets(data)

pInd = np.arange(0,npoint)

x = np.random.random(npoint)
y = np.random.random(npoint)

c = np.empty_like(x,dtype = 'str')
c[np.where(x > 0.5)] = 'red'
c[np.where(x <= 0.5)] = 'blue'

vx = -500.*np.ones(npoint)
vy = np.zeros(npoint)
vx[np.where(x <= 0.5)] = -vx[np.where(x <= 0.5)]
im = ax.scatter(x, y ,c=c)
im.set_sizes([10])

animation = animation.FuncAnimation(fig,update_point, nframe, interval=10,
									repeat=False)


animation.save('collisions.mp4')

vel = np.sqrt(vx**2+vy**2)   # velocity of particles in the end
energy = 0.5*pMass*(vel**2)   # energy of particles in the end

# plot the histograms
def MBvelDistri(vel, T):
    f_distri = ((pMass*vel)/(k_B*T))*np.exp(-0.5*pMass*vel**2/(k_B*T))
    return f_distri


plt.clf()
fig2 = plt.figure(2)
plt.subplot(211)
nv, binsv, patchesv = plt.hist(vel, bins=50, density = True, facecolor='y',
							   alpha=0.75, label='Numerical Values (v)')

middle = np.array(binsv[:-1]) +0.5*(binsv[1]-binsv[0])
Temp1, pcov = curve_fit(MBvelDistri,middle,nv,p0=[250])

plt.plot(middle,MBvelDistri(middle,Temp1),'c-',label='Fit of f(v)')
plt.title('Probability Distribution for Velocity')
plt.xlabel('v (m/s)')
plt.ylabel('f(v)')
plt.grid()
plt.legend()

plt.subplot(212)
nE, binsE, patchesE = plt.hist(energy, bins=50, density = True, facecolor='g', alpha=0.75, label='Numerical Values (E)')

middle2 = np.array(binsE[:-1]) +0.5*(binsE[1]-binsE[0])

g_distri = (1/(k_B*Temp1))*(np.exp(-1*middle2/(k_B*Temp1)))

plt.plot(middle2,g_distri,'k-',label = 'Theoretical fit (E)')
plt.title('Probability Distribution for Energy')
plt.xlabel('E (J)')
plt.ylabel('g(E)')
plt.grid()
plt.legend()
plt.tight_layout()  # allows us to see the titles clearly of each subplot.
plt.savefig('distributions.pdf')
plt.show()

f = open('collisions.txt', 'w')
f.write('Temperature obtained from fitting: {} K'.format(Temp1[0]))
f.close()
