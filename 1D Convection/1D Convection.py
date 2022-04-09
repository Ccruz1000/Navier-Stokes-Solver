import numpy as np
import matplotlib.pyplot as plt

def step(t, t_end, step_value):
    if t < t_end:
        value = step_value
    else:
        value = 0
    return value


#Define the function to calculate the first derivative
def d1_conv (dx, q, N):
    N = q.shape[0]
    dqdx = np.zeros_like(q)
    dqdx[0] = (q[1] - q[0]) / dx

    for i in np.arange(1, N-1):
        dqdx[i] = (q[i] - q[i-1]) / dx
    return dqdx

#time variables
t = 0
t_final = 0.75
dt = 1e-4

#physical parameters
L = 1
c = 1

#space parameters
N = 200
dx = L / N

#Initialize solution arrays
q1 = np.zeros(N + 1)
q1_old = np.zeros(N + 1)
dq1dx = np.zeros(N + 1)

x = np.zeros(N + 1)

for i in np.arange(N):
    x[i + 1] += x[i] + L / N

while t < t_final:
    t += dt

    q1[0] = step(t, 0.2, 1)

    dq1dx = d1_conv(dx, q1_old, N)

    q1[1: -1] = q1_old[1: -1] - c * dt * dq1dx[1: -1]

plt.plot(x, q1, label='1D convection')
plt.legend()
plt.grid()
plt.show()


