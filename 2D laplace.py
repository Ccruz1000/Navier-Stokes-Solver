
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


# 3D Plotting

def plot2D(x, y, p):
    fig = plt.figure(figsize=(11, 7), dpi=100)
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.viridis,
            linewidth=0, antialiased=False)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    plt.show()

#Laplace 2-D Function

def laplace2d(p, y, dx, dy, l1norm_target):
    l1norm = 1
    pn = np.empty_like(p)

    while l1norm > l1norm_target:
        pn = p.copy() #Create a new list and saves old list if while condition is not satisfied
        p[1:-1,1:-1] = ((dy ** 2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) + dx ** 2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) / (2 *
                                                        (dx ** 2 + dy ** 2))) # Discretized equation to solve for p

        #Boundary Conditions
        p[:, 0] = 0 #p = 0 @ x = 0
        p[:,-1] = y #p = y @ x = 2
        p[0, :] = p[1, :] # dp/dy = 0 @ y = 0
        p[-1,:] = p[-2, :] # dp/dy = 0 @ y = 1
        l1norm = (np.sum(np.abs(p[:]) - np.abs(pn[:]))) / (np.sum(np.abs(pn[:]))) # Compares copy of old p 'pn' to new p
    return p

#physical parameters
L = 1
c = 1e-4 # L1 Target
#space parameters
N = 200
M = 200
dx = L / N
dy = L / M

#Initial Conditions
p = np.zeros((N, M))

##plotting aids
x = np.linspace(0, 2, N)
y = np.linspace(0, 1, M)

# boundary conditions
p[:, 0] = 0  # p = 0 @ x = 0
p[:, -1] = y  # p = y @ x = 2
p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

# 2D Plot with 3D projection

p = laplace2d(p,y,dx,dy, c)

plot2D(x, y, p)