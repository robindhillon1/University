import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 50
grInit = np.random.choice((1, -1), (N, N))

T = [0.01, 0.1, 1, 2, 3, 4, 5, 10, 100]  #given temperatures


def energy(grInit):
    # This function calculates the energy of the system once we accomodate for
    # the nearest neighbours. Using np.roll, we don't double count and we get
    # the desired results. 
    grR = np.roll(grInit, 1, axis = 1)  # Shifted 1 to the right
    grU = np.roll(grInit, -1, axis = 0)  # Shifted 1 up
    H = -1*np.sum(grInit*(grR+grU))
    return H


# =============================================================================
# def energy(grInit):
#     # This function calculates energy of the system once we accomodate for
#     # the nearest neighbours. Using np.roll, we don't double count and we get
#     # the desired results. Just accomodating for a horizontal and vertical
#     # shift accounts for the 4 neighbours.
#     grR = np.roll(grInit, 1, axis=1)  # Shifted 1 to the right
#     grU = np.roll(grInit, -1, axis=0)  # Shifted 1 up
#     H = -1*np.sum(grInit*(grR+grU))
#     return H
# =============================================================================


def new_config(grInit):
    grNew = np.copy(grInit)
    i = np.random.randint(0, 50)
    j = np.random.randint(0, 50)
    grNew[i,j] = -1*grNew[i, j]
    return grNew, grInit, i, j  # returns new grid, old, row, column. 


def local_energy(grInit,i,j):
    e = grInit[(i+1)%50,j]+grInit[(i-1)%50,j]+grInit[i,(j+1)%50]+grInit[i,(j-1)%50]
    # If the electron is on a boundary of the grid, %50 accounts for that. 
    return grInit[i,j]*e
    

moments=np.zeros((9,5))  # Final moments. Moments which we'll be taking the max of. 
for i, t in enumerate(T):
    for j in range(5):
        grInit = np.random.choice((1, -1), (N, N))
        for k in range(600000):
            grNew, grInit, ii, jj = new_config(grInit)  #new,old,row,col
            e_local = local_energy(grInit,ii,jj)  # local energy around electron
            if e_local<0:
                grInit = grNew
            else:
                probability = np.exp((-2*e_local)/t)
                a = np.random.uniform(0,1)
                if a < probability:  # condition if probability is higher
                    grInit = grNew
        moments[i,j] = abs(np.sum(grInit))  # moments which we'll take max of.   
    moments_final = np.max(moments,axis=1)  # these moments are the ones we plot

avg_moment = np.mean(moments_final)
print(moments_final)

plt.figure()
plt.plot(T,moments_final,'ro-')
plt.xlabel('Temperature')
plt.ylabel('Magnetic Moment')
plt.grid()
plt.savefig('Tcurie.pdf')
plt.show()

#################################################### Animations section. 

# do a for loop for 600 hunder iters, then iside 1000 times main function. 
# the ncheck if temp is the one for the movies, then add the grid. 
def anim_grids(grInit,energy1,T):   # numITer
    grids = []  # empty list that we'll be appending the updated grids to.
    grInit = np.random.choice((1, -1), (N, N))
    for i in range(600000):
        grNew, grInit, ii, jj = new_config(grInit)  # new configuration of grid
        e_local = local_energy(grInit,ii,jj)  # local energy at electron
        if e_local < 0:
            grInit = grNew
            if i%1000 == 0:  # accept every 1000 times. 
                grids.append(grInit)
        else:
            probability = np.exp((-2*e_local)/T)
            a = np.random.uniform(0,1)
            if a < probability:
                grInit = grNew
                if i%1000 == 0:
                    grids.append(grInit)
    return grids  # here we return the updated list of grids

Temps = [0.1,2.5,100]


def animations(grids):
    fig = plt.figure()  # New figures for the animations
    ims = []  # emplty list, which I'll be appending the grids to. 
    for i in grids: # iterating through the list of grids
        ims.append((plt.pcolormesh(i, cmap='seismic'),))  # plotting the base
        # with chosen colormap. 
    imani = animation.ArtistAnimation(fig, ims, interval=50, repeat=False)
    return imani


grInit = np.random.choice((1,-1),(N,N))
energy1 = energy(grInit)

for i in Temps:
    grids = anim_grids(grInit,energy1,i)
    mov = animations(grids)
    mov.save('temp{}.mp4'.format(i))
    plt.show()   