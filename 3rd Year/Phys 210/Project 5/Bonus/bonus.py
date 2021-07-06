# For the bonus, I will solve the Heat Equation which is very similar to the
# diffusion equation. The setup will be the following: there's a metallic room
# which is heated from the right, left, and the top (bird's eye view). These 3
# walls are held at a constant temperature of 400 Kelvin. The 'bottom' wall is
# non-existent. It's a room with 3 walls, and the 4th wall is not there. The
# temperature at this bottom non-existent wall is 200K. The room itself is at
# an initial temperature of 300K. Over time, I will investigate the convergence
# of the temperature in the room. So that the code runs fast and we see the
# trend, I will just run it up to 10 seconds.

# After this is done, I will solve the Laplace equation using the method of
# relaxation. The boundary conditions will be a bit similar; they will be
# described below.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# the constants below are the same as the real project
L = 1
D = 1e-3
N = 100

dt = 1e-2  # time interval for convergence
dx = L/N

k = dt/dx/dx*D

tmax = 10  # so that we get quick results and code doesn't take forever to run

T_sides, T_room, T_bot = 400, 300, 200.  # initial temperatures in K

T = np.zeros((101, 101))
T[0, :] = T_sides  # assigning the initial temperatures; boundary conditions
T[-1, :] = T_bot
T[:, 0] = T_sides
T[:, -1] = T_sides
T[1: -1, 1: -1] = T_room

nframe = 200  # number of frames for animation. Since I used ArtistAnimation, I
# will use FuncAnimation below.
fig, ax = plt.subplots()
im = plt.imshow(T)  # initial state
steps = int(tmax/dt)+1  # same as the project


def animations(steps):
    """
    This function calculates the distribution/convergence of the temperature in
    a room of 3 walls that are being heated at a constant temperature. As time
    increases, the state of the room is calculated (temperature distribution).
    In this case, 10 seconds.
    """
    global T
    for i in range(steps):
        Tx_left = np.roll(T, -1, axis=0)  # same procedure as in the project.
        Tx_right = np.roll(T, 1, axis=0)
        Ty_Up = np.roll(T, 1, axis=1)
        Ty_Down = np.roll(T, -1, axis=1)
        Tx_new = T + k*(Tx_left+Tx_right-2*T)
        T_new = Tx_new + k*(Ty_Up+Ty_Down-2*T)

        T_new[0, :] = T_sides
        T_new[-1, :] = T_bot
        T_new[:, 0] = T_sides
        T_new[:, -1] = T_sides
        T = T_new
    im.set_array(T)  # needed for FuncAnimation


anim = animation.FuncAnimation(fig, animations, nframe, interval=50,
                               repeat=False)
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.colorbar()  # displays the temperature bar
plt.title('Temperature Distribution in Room after 10s')
anim.save('bonusConvergence.mp4')
plt.close()
# =============================================================================

# A very similar branch of physics that we can explore is the Laplace Equation.
# Below, I will solve the Laplace equation with 4 boundary conditions, 2D, like
# we have done so in this project. Please see the boundary conditions and code
# below. This is the method of relaxation.

N = 35  # Number of points where we compute the solution
stp = 1  # Step size

x = np.linspace(0, 10, N)  # Locations where we'll solve the equation
y = np.linspace(0, 10, N)
V = np.zeros((N, N))  # initial array for potential

# The boundary values have been initialized in the for loop below:
plt.subplot(211)
for i in range(N):
    V[i, 0] = np.cos(np.pi*x[i]/x[-1])  # boundary conditions
    V[i, -1] = np.cos(np.pi*x[i]/x[-1])
    V[0, i] = np.sin(np.pi*y[i]/y[-1])
    V[-1, i] = np.sin(np.pi*y[i]/y[-1])
    V[i, i] = np.cos(np.pi*y[i]/y[-1])
    V[i, -i] = np.sin(np.pi*x[i]/x[-1])

plt.imshow(V)  # initial state of the potential
plt.xlabel('x')
plt.ylabel('y')
plt.title('Laplace Equation Before Relaxation')

plt.subplot(212)
# Iterating the soluton(s), 'relaxation', and plotting result after convergence
for i in range(N):
    for j in range(N-1):  # exluding endpoints
        for k in range(N-1):
            V[j, k] = 0.25*(V[j-stp, k]+V[j+stp, k]+V[j, k+stp]+V[j, k-stp])
            # Above is the equation for Method of Relaxation in 2D.

plt.imshow(V)  # after relaxation
plt.xlabel('x')
plt.ylabel('y')
plt.title('Laplace Equation After Relaxation')
plt.tight_layout()  # allows to see titles and labels clearly of each subplot.
plt.savefig('relaxation.pdf')
