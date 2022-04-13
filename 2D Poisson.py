
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def plot2D(x, y, p):
    fig = plt.figure(figsize=(11, 7), dpi=100)
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.viridis,
            linewidth=0, antialiased=False)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    plt.show()

# physical parameters

lx = 1
ly = 1

# space parameters
N = 200
M = 200
nt = 100 # Total time elapsed
dx = lx / (N - 1)
dy = ly / (M - 1)

# Initialization
p = np.zeros((M, N))
pd = np.zeros((M, N))
b = np.zeros((M, N))
x = np.linspace(0, lx, N)
y = np.linspace(0, ly, M)

# Source to be added to laplace equation
b[int(M / 4), int(N / 4)] = 100
b[int(3 * M / 4), int(3 * N / 4)] = -100

def poisson2D(p, pd, b, dx, dy):

    for it in range(nt):

        pd = p.copy()

        p[1:-1, 1:-1] = (((pd[1:-1, 2:] + pd[1:-1, :-2]) * dy ** 2 +
                      (pd[2:, 1:-1] + pd[:-2, 1:-1]) * dx ** 2 -
                      b[1:-1, 1:-1] * dx ** 2 * dy ** 2) /
                     (2 * (dx ** 2 + dy ** 2)))

        p[0, :] = 0
        p[M -1, :] = 0
        p[:, 0] = 0
        p[:, N - 1] = 0
    return p

p = poisson2D(p, pd, b, dx, dy)

plot2D(x, y, p)