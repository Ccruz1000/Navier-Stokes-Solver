# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 15:47:27 2022

@author: Rayan
"""
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 00:17:59 2022

@author: Rayan
"""

import numpy as np
import matplotlib.pyplot as plot

#Inlet Conditions
ui = 2.3
Tin = 25
visc = 1.6e-5
x = 0.3048
pr = 0.71

#Physical Properties
a = visc / pr
q = 2100.
re = (ui * x) / visc
k = 0.0255

#Boundary Layer
y = 0.03
si = 4.92 * (re ** -0.5)
so = (si / y) * 100

#Grid Generation
imax = 65
jmax = 65

#Variable Initiation
dx = x / imax
dy = y / jmax
N = -6 / imax
h = np.arange(28, 22, N)
u = np.ones((imax, jmax))
v = np.zeros((1, jmax))
x = np.zeros((1, jmax))
T = np.zeros((imax, jmax))

#Boundary Conditions
u[:, :] = 1
u[:, 0] = 0
u[0, :] = ui
u[:, jmax - 1] = ui
u[0, 0] = 0

v[0, :] = 0
v[:, 0] = 0
v[0, 0] = 0
T[0,:]=Tin
T[:,jmax-1]=Tin

for o in range(0, imax-1):
    for j in range(1, jmax - 2):
        for i in range(0, imax - 2):
            u[i, j-1] = (u[i-1, j-1]) + (visc * (u[i-1, j ] - 2 * u[i-1, j-1] + u[i-1, j - 2]) * (dx)) / ((u[i-1, j-1] * dy ** 2)) - ((u[i-1, j ] - u[i-1, j - 2]) * (v[i-1, j-1] * dx) / (2 * u[i-1, j-1] * dy))
            v[i+1, j] = v[i+1 , j - 1] - ((0.5 * dy * (u[i+1, j] - u[i, j] + u[i+1 , j - 1] - u[i, j - 1])) / dx)

for m in range(0, imax-1):
    for j in range(1, jmax - 2):
        for i in range(0, imax - 2):
            T[i, 0] = (T[i, 1] + q / h[0, i])
            T[i + 1, j] = (T[i, j]) + (a * (T[i, j + 1] - 2 * T[i, j] + T[i, j - 1]) * dx) / (u[i, j] * dy ** 2) + (
                        (T[i, j + 1] - T[i, j - 1]) * (v[i, j] * dx) / (u[i, j] * 2 * dy))
            T[imax, 0] = (T[imax - 1, 0] - T[imax, 1])

Nu_length = 50
re1 = (u[Nu_length, 2] * dx * Nu_length) / visc
cfx1 = 0.664 * np.sqrt(1 / re1)
cfx = (2 * visc * (u[Nu_length, 2] - u[Nu_length, 1])) / (dy * ui ** 2)

nu1 = 0.453 * pr ** (1 / 3) * np.sqrt(re1)
Nu1 = (q * (dx * Nu_length)) / (k * (T[Nu_length, 1] - T[Nu_length, 2]))

plot.plot(x, u)

plot.show()
