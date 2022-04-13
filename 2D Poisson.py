

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


# Define Parameters

#physical parameters
L = 1 #If both the length of the x and y are the same
#space parameters
N = 200
M = 200
nt = 100 # number of timesteps
dx = L / (N-1)
dy = L / (M-1)

# Initialization
p = np.zeros((M, N))
pd = np.zeros((M, N))
b = np.zeros((M, N))
x = np.linspace(L, N)
y = np.linspace(L, M)

# Source to be added to laplace equation
b[int(M / 4), int(N / 4)] = 100
b[int(3 * M / 4), int(3 * N / 4)] = -100

for it in range(nt):

    pd = p.copy()

    p[1:-1,1:-1] = (((pd[1:-1, 2:] + pd[1:-1, :-2]) * dy**2 +
                    (pd[2:, 1:-1] + pd[:-2, 1:-1]) * dx**2 -
                    b[1:-1, 1:-1] * dx**2 * dy**2) /
                    (2 * (dx**2 + dy**2)))

    p[0, :] = 0
    p[M-1, :] = 0
    p[:, 0] = 0
    p[:, N-1] = 0

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

plot2D(x, y, p)