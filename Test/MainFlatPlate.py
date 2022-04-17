import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


# Define Domain
num_x = 20  # Number of X points
num_y = 20  # Number of Y points
x_end = 0.3048  # Final Point in x
x_start = 0.05  # First point in x (absolute value, shifted to left on x axis)
y_end = 0.05  # Final Point in y
y_start = 0
MAXIT = 1e4  # Number of time steps allowed before stopping program

# Inflow Conditions

M_inf = 4  # Free stream mach number
T_inf = 288.16  # Free stream temperature (K)
p_inf = 101325  # Free stream pressure (pa)
inflow = Primitives(0, 0, p_inf,T_inf)
inflow.u = M_inf*inflow.a # m/s
Re_inf = inflow.calculateReynoldsNumber(x_end)

