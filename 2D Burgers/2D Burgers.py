#2D PDEs:
# du/dt + udu/dx + vdu/dy = nu(d2u/dx2 + d2u/dy2)
# dv/dt +udv/dx + vdv/dy = nu(d2v/dx2 + d2v/dy2)

import numpy as np
import matplotlib.pyplot as plt

def step(t, t_end, step_value):
    if t < t_end:
        value = step_value
    else:
        value = 0
    return value


#Define the function to calculate the first derivative
def conv(dx, u, N):
    N = u.shape[0]
    dudx1 = np.zeros((len(u),2))
    dudx1[0, :] = (u[1, :] - u[0, :]) / dx
    dudx1[N - 1, :] = (u[-2, :] - u[-1, :]) / dx


    for j in np.arange(0, 1):
        for i in np.arange(1, N - 1):
            dudx1[i, j] = (u[i, j] - u[i-1, j]) / dx
    return dudx1


def diff(dx, u, N):
    N = u.shape[0]
    dudx2 = np.zeros((len(u), 2))
    dudx2[0, :] = (u[2, :] + u[0, :] - 2 * u[1, :]) / (dx ** 2)
    dudx2[N - 1, :] = 0

    for j in np.arange(0, 1):
        for i in np.arange(1, N - 1):
            dudx2[i, j] = (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / (dx ** 2)
    return dudx2


#time variables
t = 0
t_final = 0.75
dt = 1e-4

#physical parameters
L = 1
b = 1e-3

#space parameters
N = 200
M = 200
dx = L / N
dy = L / M

#Initialize solution arrays
u1 = np.zeros((N + 1, 2))
u1_old = np.zeros((N + 1, 2))
v1 = np.zeros((N + 1, 2))
v1_old = np.zeros((N + 1, 2))

x = np.zeros(N + 1)
y = np.zeros(M + 1)


for i in np.arange(N):
    x[i + 1] += x[i] + L / N

for i in np.arange(M):
    y[i + 1] += y[i] + L / M

while t < t_final:
    t += dt


    u1[0, 0] = step(t, 0.2, 1)

    #Define discretization for x component
    dudx = conv(dx, u1_old, N)
    dudy = conv(dy, u1_old, N)
    du2dx2 = diff(dx, u1_old, N)
    du2dy2 = diff(dy, u1_old, N)

    #Define discretizatiom for y component
    dvdx = conv(dx, v1_old, N)
    dvdy = conv(dy, v1_old, N)
    dv2dx2 = diff(dx, v1_old, N)
    dv2dy2 = diff(dy, v1_old, N)

    u1[1: -1, 1: -1] = u1_old[1: -1, 1: -1] - dt * u1_old[1: -1, 1: -1] * dudx[1: -1, 1: -1] - dt * v1_old[1: -1, 1: -1] * dudy[1: -1, 1: -1] + b * dt * du2dx2[1: -1, 1: -1] + b * dt * du2dy2[1: -1, 1: -1]
    v1[1: -1, 1: -1] = v1_old[1: -1, 1: -1] - dt * u1_old[1: -1, 1: -1] * dvdx[1: -1, 1: -1] - dt * v1_old[1: -1, 1: -1] * dvdy[1: -1,1: -1] + b * dt * dv2dx2[1: -1, 1: -1] + b * dt * dv2dy2[1: -1, 1: -1]
    u1_old = u1
    v1_old = v1
    print(u1)

plt.plot(x, u1[:, 0], label= "1D Burger's")
plt.plot(x, v1[:, 0])
plt.legend()
plt.grid()
plt.title('dt= ' + str(dt) + ' / N= ' + str(N))
plt.show()



