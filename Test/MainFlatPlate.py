import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def ddxb(f, dx):
    [ny, nx, n3] = np.shape(f)
    inX = [1, nx - 1]
    dfdx = np.zeros(ny, nx, n3)
    dfdx[:, inX, :] = (f[:, inX, :] - f[:, inX - 1, :]) / dx
    dfdx[:, 0, :] = dfdx[:, 1, :]

    return dfdx


def ddxc(f, dx):
    [ny, nx] = np.shape(f)
    inX = (1, nx - 2)
    dfdx = np.zeros(ny, nx)
    dfdx[:, inX] = (f[:, inX + 1] - f[:, inX - 1]) / (2 * dx)
    dfdx[:, 0] = (f[:, 1] - f[:, 0]) / (dx);
    dfdx[:, nx - 1] = (f[:, nx - 1] - f[:, nx - 2]) / (dx)

    return dfdx


def ddxf(f, dx):
    [ny, nx, n3] = np.shape(f)
    inX = [0, nx - 2]
    dfdx = np.zeros(ny, nx, n3)
    dfdx[:, inX, :] = (f[:, inX + 1, :] - f[:, inX, :]) / dx
    dfdx[:, nx - 1, :] = dfdx[:, nx - 2, :]

    return dfdx


def ddyb(f, dy):
    [ny, nx, n3] = np.shape(f)
    inY = [1, ny - 1]
    dfdy = np.zeros(ny, nx, n3)
    dfdy[inY, :, :] = (f[inY, :, :] - f[inY - 1, :, :]) / dy
    dfdy[0, :, :] = dfdy[1, :, :]

    return dfdy


def ddyc(f, dy):
    [ny, nx] = np.shape(f)
    inY = [1, ny - 2]
    dfdy = np.zeros(ny, nx)
    dfdy[inY, :] = (f[inY + 1, :] - f[inY - 1, :]) / (2 * dy)
    dfdy[0, :] = (f[1, :] - f[0, :]) / (dy)
    dfdy[ny - 1, :] = (f[ny - 1, :] - f[ny - 2, :]) / (dy)
    return dfdy


def ddyf(f, dy):
    [ny, nx, n3] = np.shape(f)
    inY = [0, ny - 2]
    dfdy = np.zeros(ny, nx, n3)
    dfdy[inY, :, :] = (f[inY + 1, :, :] - f[inY, :, :]) / dy
    dfdy[ny - 1, :, :] = dfdy[ny - 2, :, :]
    return dfdy


# Calculate E flux array from primitive variables
# To maintain second-order accuracy, the x-derivative terms appearing in E
# are differenced in the opposite direction to that used for dE/dx, while
# the y-derivative terms are approximated using central differences.
# Likewise, the y-derivative terms appearing in F are differenced in the
# opposite direction to that used for dF/dy, while the x-derivative terms
# in F are central differenced

def calculateE(primitives, dx, dy, direction):
    # primitives
    [u, v, p, T] = primitives.deal()
    [mu, lam, k] = primitives.getMuLambdaK()
    r = primitives.r
    Et = primitives.Et

    # gradients
    if direction == 'forward':
        dudx = ddxf(u, dx)
        dvdx = ddxf(v, dx)
        dTdx = ddxf(T, dx)
    elif direction == 'backward':
        dudx = ddxb(u, dx)
        dvdx = ddxb(v, dx)
        dTdx = ddxb(T, dx)
    else:
        ValueError('direction not allowed')
    dudy = ddyc(u, dy)
    dvdy = ddyc(v, dy)

    # stresses and heat fluxes
    txx = lam * (dudx + dvdy) + 2 * mu * dudx
    # tyy = lam * (dudx + dvdy)+ 2 * mu * dvdy
    txy = mu * (dudy + dvdx)
    qx = -k * dTdx

    E = np.zeros(np.shape(r, 1), np.shape(r, 2), 4)
    E[:, :, 0] = r * u
    E[:, :, 1] = r * u ** 2 + p - txx
    E[:, :, 2] = r * u * v - txy
    E[:, :, 3] = (Et + p) * u - u * txx - v * txy + qx
    return E


def calculateF(primitives, dx, dy, direction):
    [u, v, p, T] = primitives.deal()
    [mu, lam, k] = primitives.getMuLambdaK()
    r = primitives.r
    Et = primitives.Et
    if direction == "forward":

        dudy = ddyf(u, dy)
        dvdy = ddyf(v, dy)
        dTdy = ddyf(T, dy)

    elif direction == "backward":
        dudy = ddyf(u, dy)
        dvdy = ddyf(v, dy)
        dTdy = ddyf(T, dy)

    else:
        ValueError("direction not allowed")
    dudx = ddxc(u, dx)
    dvdx = ddxc(v, dx)
    tyy = lam * (dudx + dvdy) + 2 * mu * (dvdy)
    txy = mu * (dudy + dvdx)
    qy = -k * dTdy
    F = np.zeros(np.shape(r, 1), np.shape(r, 2), 4)
    F[:, :, 0] = r * v
    F[:, :, 1] = r * u * v - txy
    F[:, :, 2] = r * v ** 2 + p - tyy
    F[:, :, 3] = (Et + p) * v - u * txy - v * tyy + qy

    return F


def calculateU(primitives):
    U = np.zeros(np.shape(primitives.u, 1), np.shape(primitives.u, 2), 4)
    U[:, :, 0] = primitives.r
    U[:, :, 1] = U[:, :, 0] * primitives.u
    U[:, :, 2] = U[:, :, 0] * primitives.v
    U[:, :, 3] = primitives.Et
    return U


# Calculate primitive variables from solution vector
def decodeSolutionVector(U):
    r = U[:, :, 0]
    u = U[:, :, 1] / r
    v = U[:, :, 2] / r
    Et = U[:, :, 3]
    e = Et / r - 0.5 * (u ** 2 + v ** 2)
    cv = Primitives.R / (Primitives.gm - 1)
    T = e / cv
    p = r * Primitives.R * T
    primitives = Primitives(u, v, p, T)
    return primitives


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
