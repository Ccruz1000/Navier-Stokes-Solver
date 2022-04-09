import numpy as np
import matplotlib.pyplot as plt

def step(t, t_end, step_value):
    if t < t_end:
        value = step_value
    else:
        value = 0
    return value


#Define the function to calculate the first derivative
def d1_conv(dx, q, N):
    N = q.shape[0]
    dqdx = np.zeros_like(q)
    dqdx[0] = (q[1] - q[0]) / dx

    for i in np.arange(1, N-1):
        dqdx[i] = (q[i] - q[i-1]) / dx
    return dqdx

def d1_diff(dx, q, N):
    N = q.shape[0]
    dqdx2 = np.zeros_like(q)
    dqdx2[0] = (q[2] + q[0] - 2 * q[1]) / (dx ** 2)

    for i in np.arange(1, N-1):
        dqdx2[i] = (q[i+1] - 2 * q[i] + q[i-1]) / (dx ** 2)

    return dqdx2


#time variables
t = 0
t_final = 0.75
dt = 1e-4

#physical parameters
L = 1
b = 0.07

#space parameters
N = 200
dx = L / N

#Initialize solution arrays
q1 = np.zeros(N + 1)
q1_old = np.zeros(N + 1)
dq1dx = np.zeros(N + 1)
dq2dx = np.zeros(N + 1)

x = np.zeros(N + 1)

for i in np.arange(N):
    x[i + 1] += x[i] + L / N

while t < t_final:
    t += dt

    q1[0] = step(t, 0.2, 1)

    dq1dx = d1_conv(dx, q1_old, N)
    dq2dx = d1_diff(dx, q1_old, N)

    q1[1: -2] = q1_old[1: -2] + b * dt * dq2dx[1: -2] - q1_old[1: -2] * dt * dq1dx[1: -2]
    q1[-1] = q1[-2]     #zero gradient boundary layer
    q1_old = q1

plt.plot(x, q1, label= "1D Burger's")
plt.legend()
plt.grid()
plt.show()



