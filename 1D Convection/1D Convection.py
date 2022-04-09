import numpy as np
import matplotlib as plt

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
dx = 1 / N

#Initialize solution arrays
q1 = np.zeros(N + 1)
q1_old = np.zeros(N + 1)
dq1dx = np.zeros(N + 1)


