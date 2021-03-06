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
def duj(dx, u, N):
    N = u.shape[0]
    dudx = np.zeros_like(u)
    dudx[0] = (u[1] - u[0]) / dx

    for i in np.arange(1, N-1):
        dudx[i] = (u[i] - u[i-1]) / dx
    return dudx

def dui(dy, u, N):
    N = u.shape[0]
    dudy = np.zeros_like(u)
    dudy[0] = (u[1] - u[0]) / dy

    for i in np.arange(1, N-1):
        dudy[i] = (u[i] - u[i-1]) / dy
    return dudy

def d1x_diff(dx, u, N):
    N = u.shape[0]
    dudx2 = np.zeros_like(u)
    dudx2[0] = (u[2] + u[0] - 2 * u[1]) / (dx ** 2)

    for i in np.arange(1, N-1):
        dudx2[i] = (u[i+1] - 2 * u[i] + u[i-1]) / (dx ** 2)

    return dudx2

def d1y_diff(dy, u, N):
    N = u.shape[0]
    dudy2 = np.zeros_like(u)
    dudy2[0] = (u[2] + u[0] - 2 * u[1]) / (dy ** 2)

    for i in np.arange(1, N-1):
        dudy2[i] = (u[i+1] - 2 * u[i] + u[i-1]) / (dy ** 2)

    return dudy2

#time variables
t = 0
t_final = 0.75
dt = 1e-4

#physical parameters
L = 1
b = 1e-3

#space parameters
N = 400
M = 400
dx = L / N
dy = L / M

#Initialize solution arrays
u1x = np.zeros(N + 1)
u1x_old = np.zeros(N + 1)
du1dx = np.zeros(N + 1)
du2dx = np.zeros(N + 1)
du1dy = np.zeros(M + 1)
du2dy = np.zeros(M + 1)

x = np.zeros(N + 1)
y = np.zeros(M + 1)

for i in np.arange(N):
    x[i + 1] += x[i] + L / N

for i in np.arange(M):
    y[i + 1] += y[i] + L / M

while t < t_final:
    t += dt


    q1[0, i] = step(t, 0.2, 1)

    dq1dx = d1x_conv(dx, q1_old, N)
    dq2dx = d1x_diff(dx, q1_old, N)

    q1[1: -2] = q1_old[1: -2] + b * dt * dq2dx[1: -2] - q1_old[1: -2] * dt * dq1dx[1: -2]
    q1[-1] = q1[-2]     #zero gradient boundary layer
    q1_old = q1

plt.plot(x, q1, label= "1D Burger's")
plt.legend()
plt.grid()
plt.title('dt= ' + str(dt) + ' / N= ' + str(N))
plt.show()



