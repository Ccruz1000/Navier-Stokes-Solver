import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def ddxb(f, dx):
    [ny, nx, n3] = np.shape(f)
    inX = [1, nx - 1]
    dfdx = np.zeros(ny, nx, n3)
    dfdx[:, inX, :] = (f[:, inX , :] - f[:, inX-1, :]) / dx
    dfdx[:, 0, :] = dfdx[:, 1, :]

    return dfdx


def ddxc(f, dx):
    [ny, nx] = np.shape(f)
    inX = (1, nx - 2)
    dfdx = np.zeros(ny, nx)
    dfdx[:, inX ] = (f[:, inX+1] - f[:, inX - 1]) / (2 * dx)
    dfdx[:, 0] = (f[:, 1] - f[:, 0]) / (dx);
    dfdx[:, nx - 1] = (f[:, nx - 1] - f[:, nx - 2]) / (dx)

    return dfdx


def ddxf(f, dx):
    [ny, nx, n3] = np.shape(f)
    inX = [0, nx - 2]
    dfdx = np.zeros(ny, nx, n3)
    dfdx[:, inX , :] = (f[:, inX+1, :] - f[:, inX , :]) / dx
    dfdx[:, nx - 1, :] = dfdx[:, nx - 2, :]

    return dfdx


def ddyb(f, dy):
    [ny, nx, n3] = np.shape(f)
    inY = [1, ny - 1]
    dfdy = np.zeros(ny, nx, n3)
    dfdy[inY , :, :] = (f[inY , :, :] - f[inY - 1, :, :]) / dy
    dfdy[0, :, :] = dfdy[1, :, :]

    return dfdy


def ddyc(f, dy):
    [ny, nx] = np.shape(f)
    inY = [1, ny - 2]
    dfdy = np.zeros(ny, nx)
    dfdy[inY , :] = (f[inY+1, :] - f[inY - 1, :]) / (2 * dy)
    dfdy[0, :] = (f[1, :] - f[0, :]) / (dy)
    dfdy[ny - 1, :] = (f[ny - 1, :] - f[ny - 2, :]) / (dy)
    return dfdy


def ddyf(f, dy):
    [ny, nx, n3] = np.shape(f)
    inY = [0, ny - 2]
    dfdy = np.zeros(ny, nx, n3)
    dfdy[inY, :, :] = (f[inY+1, :, :] - f[inY , :, :]) / dy
    dfdy[ny - 1, :, :] = dfdy[ny - 2, :, :]
    return dfdy


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
