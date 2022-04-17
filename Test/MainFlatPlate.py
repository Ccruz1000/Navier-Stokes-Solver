import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Plate length
lhori = 0.00001  # m

# Courant number
K = 0.8

# grid size and max iterations
nx = 70
ny = 70
maxiter = 10000

# Inflow conditions
Minf = 4
pinf = 101325  # Pa
Tinf = 288.16  # Kelvin
inflow = Primitives(0, 0, pinf, Tinf)
inflow.u = Minf * inflow.a  # m/s
Reinf = inflow.calculateReynoldsNumber(lhori)

# boundary layer size & vertical domain size
delta = 5 * lhori / np.sqrt(Reinf)
lvert = 5 * delta

# grid
x, y = np.meshgrid(np.linspace(0, lhori, nx), np.linspace(0, lvert, ny))

# Set initial conditions
# Set all intial values to inflow values. Conditions on boundaries updated by solveMacCormack()
primitives = Primitives(inflow.u * np.ones(ny, nx), inflow.v * np.ones(ny, nx), inflow.p * np.ones(ny, nx), \
                        inflow.T * np.ones(ny, nx))

# Solve two wall temperature conditions
Tw_Tinf = 1.0  # constant wall temperature
constantTw = solveMacCormack(primitives, inflow, Tw_Tinf, K, x, y, maxiter)
