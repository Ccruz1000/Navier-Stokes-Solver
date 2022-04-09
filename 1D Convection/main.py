import numpy as np
import matplotlib as plt

def d1_conv (dx, q, N):
    N = q.shape[0]
    dqdx = np.zeros_like(q)
    dqdx[0] = (q[1] - q[0]) / dx

    for i in np.arange(1, N-1):
        dqdx[i] = (q[i] - q[i-1]) / dx
    return dqdx

t = 0
t_final = 0.75
dt = 1e-4

N = 200
dx = 1 / N

