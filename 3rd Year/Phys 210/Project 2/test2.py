import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

T = [0.01, 0.1, 1, 2, 3, 4, 5, 10, 100]  #given temperatures


def new_config(grInit):
    # This function takes the initial grid as a parameter, selects an index
    # (electron) randomly, flips its sign, and returns the new configuration.
    grNew = np.copy(grInit)
    i = np.random.randint(0, 50)
    j = np.random.randint(0, 50)
    grNew[i, j] = -1*grNew[i, j]
    return grNew, grInit, i, j  # returns new grid, old, row, column.


def local_energy(grInit, i, j):
    # Instead of calculating the whole energy of the initial grid and the new
    # grid, we can just find the local energy of the chosen electron and find
    # the local energy. This saves us time.
    e = (grInit[(i+1) % 50, j]+grInit[i, (j+1) % 50]+grInit[(i-1) % 50, j] +
         grInit[i, (j-1) % 50])
    # If the electron is on a boundary of the grid, %50 accounts for that.
    return grInit[i, j]*e  # returns the local energy


# =============================================================================
# def mainfunc(T):
#     # This is the function that does the majority of the work. We get a
#     # temperature as the paramter, generate a random grid with values (1,-1)
#     # and run it 600,000 times. See below for the specific tasks of the code.
#     grInit = np.random.choice((1, -1), (50, 50))  # aforementioned grid
#     for k in range(600000):
#         grNew, grInit, ii, jj = new_config(grInit)  # new grid, old grid, row,
#         # column. The row and column are for calling the local energy function
#         e_local = local_energy(grInit, ii, jj)  # local energy around electron
#         if e_local < 0:  # professor told me that if energy if less than 0,
#             # then we accept the new grid, which is what we do below. However,
#             # we have to use a negative 2 instead of positive 2 in the else
#             # statement below. If if use e_local > 0, then we would use +2 in
#             # the else statement, as told.
#             grInit = grNew
#         else:
#             probability = np.exp((-2*e_local)/T)  # probability for evolution
#             # of system.
#             a = np.random.uniform(0, 1)
#             if a < probability:  # if probability is higher, we accept new grid
#                 grInit = grNew
#     return grInit
# 
# 
# moments = np.zeros((9, 5))  # Pre-assinging moments which we'll be modifying as
# # we progress further.
# for i, t in enumerate(T):  # t in enumerate(T) takes out the temperatures one
#     # by one, which we use to call the mainfunc() above.
#     for j in range(5):  # to get 5 independent measurements of moments.
#         grInit = mainfunc(t)
#         moments[i, j] = abs(np.sum(grInit))  # moments which we'll take max of.
#     moments_final = np.max(moments, axis=1)  # these max moments we plot
# 
# avg_moment = np.mean(moments_final)  # the average moment. ~= 1135 in this run.
# plt.figure()
# plt.semilogx(T, moments_final, 'ro-')  # semilog since we want T in log scale.
# plt.xlabel('Temperature')
# plt.ylabel('Magnetic Moment')
# plt.title('Magnetization vs Temperature of Ising Model')
# plt.grid()
# plt.savefig('Tcurie.pdf')
# plt.show()
# 
# # The curie temperature is approximately where magnetization changes from
# # approximately 0 (high temperatures) to not zero (low temperatures). In the
# # figure, we can see that the phase/drop occurs in temperature range between
# # T = 0 and T = 10. Hence, the curie temperature is within this range. I would
# # estimate that the curie temperature is around T = 2 to 4, from the figure.
# =============================================================================


def animations(grInit,Temp):
    # In this function, we will get an initial grid and temperature as the
    # parameters. The new lines are identical from the mainfunc() function.
    grids = []  # empty list that we'll be appending the updated grids to.
    for i in range(600000):
        grNew, grInit, ii, jj = new_config(grInit)  # new configuration of grid
        e_local = local_energy(grInit, ii, jj)  # local energy at electron
        if i % 1000 == 0:  # accept network image every 1000 times.
            grids.append(grInit)  # appending it to the predefined grids above
        if e_local < 0:  # this if, else follows directly from mainfunc
            grInit = grNew
        else:
            probability = np.exp((-2*e_local)/Temp)
            a = np.random.uniform(0, 1)
            if a < probability:
                grInit = grNew
    fig = plt.figure()  # New figures for the animations
    ims = []  # emplty list, which I'll be appending the grids to. 
    for i in grids: # iterating through the list of grids
        ims.append((plt.pcolormesh(i, cmap='seismic'),))  # plotting the base
        # with chosen colormap. 
    imani = animation.ArtistAnimation(fig, ims, interval=50, repeat=False)
    return imani  # here we return the updated list of grids

Temps = [0.1,2.5,100]

grInit = np.random.choice((1,-1),(50,50))
for i in Temps:
    mov = animations(grInit,i)
    mov.save('temp_{}.mp4'.format(i))
    plt.show()   