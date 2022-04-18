import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import time
from Primitives import *


def ddxb(f, dx):
    ny, nx, n3 = np.shape(f)
    inX = np.linspace(1, nx, num=nx, endpoint=False)
    dfdx = np.zeros(ny, nx, n3)
    dfdx[:, inX, :] = (f[:, inX, :] - f[:, inX - 1, :]) / dx
    dfdx[:, 0, :] = dfdx[:, 1, :]

    return dfdx


def ddxc(f, dx):
    ny, nx = np.shape(f)
    inX = np.linspace(1, nx-1, num=nx-1, endpoint=False)
    dfdx = np.zeros(ny, nx)
    dfdx[:, inX] = (f[:, inX + 1] - f[:, inX - 1]) / (2 * dx)
    dfdx[:, 0] = (f[:, 1] - f[:, 0]) / dx
    dfdx[:, nx - 1] = (f[:, nx - 1] - f[:, nx - 2]) / dx

    return dfdx


def ddxf(f, dx):
    ny, nx, n3 = np.shape(f)
    inX = np.linspace(0, nx-1, num=nx-1, endpoint=False)
    dfdx = np.zeros(ny, nx, n3)
    dfdx[:, inX, :] = (f[:, inX + 1, :] - f[:, inX, :]) / dx
    dfdx[:, nx - 1, :] = dfdx[:, nx - 2, :]

    return dfdx


def ddyb(f, dy):
    ny, nx, n3 = np.shape(f)
    inY = np.linspace(1, ny, num=ny, endpoint=False)
    dfdy = np.zeros(ny, nx, n3)
    dfdy[inY, :, :] = (f[inY, :, :] - f[inY - 1, :, :]) / dy
    dfdy[0, :, :] = dfdy[1, :, :]

    return dfdy


def ddyc(f, dy):
    ny, nx = np.shape(f)
    inY = np.linspace(1, ny-1, num=ny-1, endpoint=False)
    dfdy = np.zeros(ny, nx)
    dfdy[inY, :] = (f[inY + 1, :] - f[inY - 1, :]) / (2 * dy)
    dfdy[0, :] = (f[1, :] - f[0, :]) / dy
    dfdy[ny - 1, :] = (f[ny - 1, :] - f[ny - 2, :]) / dy
    return dfdy


def ddyf(f, dy):
    ny, nx, n3 = np.shape(f)
    inY = np.linspace(0, ny-1, num=ny-1, endpoint=False)
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
    tyy = lam * (dudx + dvdy) + 2 * mu * dvdy
    txy = mu * (dudy + dvdx)
    qy = -k * dTdy
    F = np.zeros(np.shape(r, 1), np.shape(r, 2), 4)
    F[:, :, 0] = r * v
    F[:, :, 1] = r * u * v - txy
    F[:, :, 2] = r * v ** 2 + p - tyy
    F[:, :, 3] = (Et + p) * v - u * txy - v * tyy + qy

    return F


def calculateU(primitives):
    numx, numy = np.shape(primitives.u)
    U = np.zeros((numx, numy, 4))
    U[:, :, 0] = primitives.r
    U[:, :, 1] = U[:, :, 0] * primitives.u
    U[:, :, 2] = U[:, :, 0] * primitives.v
    U[:, :, 3] = primitives.Et
    return U


def postStressAndHeatFlux(primitives, x, y):
    dx = x[0, 1] - x[0, 0]
    dy = y[1, 0] - y[0, 0]
    [u, v, [], T] = primitives.deal()
    [mu, lam, k] = primitives.getMuLambdaK()
    dudx = ddxc(u, dx)
    dudy = ddyc(u, dy)
    dvdx = ddxc(v, dx)
    dvdy = ddyc(v, dy)
    txx = lam * (dudx + dvdy) + 2 * mu * dudx
    tyy = lam * (dudx + dvdy) + 2 * mu * dvdy
    txy = mu * (dudy + dvdx)
    dTdx = ddxc(T, dx)
    dTdy = ddyc(T, dy)
    qx = -k * dTdx
    qy = -k * dTdy
    return txx, tyy, txy, qx, qy


# Solve Navier-Stokes equations for supersonic flow over a flat plate using MacCormack method.
def solveMacCormack(primitives, inflow, Tw_Tinf, K, x, y, maxiter):
    # mesh data
    ny, nx = np.shape(x)
    dx = x[0, 1] - x[0, 0]
    dy = y[1, 0] - y[0, 0]
    inX = 1, nx - 2  # interior x
    inY = 1, ny - 2  # interior y

    # set boundary conditions
    primitives = updateBoundaryConditions(primitives, inflow, Tw_Tinf)

    # time march
    i = 0
    converged = False

    while not converged and i < maxiter:
        start_time = time.time()
        print(i)
        i += 1

        # time step size
        dt = primitives.calculateTimeStep(dx, dy, K)

        # solution vector
        U = calculateU(primitives)

        # flux vectors
        E = calculateE(primitives, dx, dy, 'backward')
        F = calculateF(primitives, dx, dy, 'backward')

        # forward finite difference predictor
        dUdt_predictor = -ddxf(E, dx) - ddyf(F, dy)

        # Predictor step and corrector vector calculations
        U2 = U
        U2[inY, inX, :] = U[inY, inX, :] + dt * dUdt_predictor[inY, inX, :]  # interior points only
        primitives2 = decodeSolutionVector(U2)
        E2 = calculateE(primitives2, dx, dy, 'forward')
        F2 = calculateF(primitives2, dx, dy, 'forward')

        # backward finite difference corrector
        dUdt_corrector = -ddxb(E2, dx) - ddyb(F2, dy)

        # MacCormack solution step
        dUdt = 0.5 * (dUdt_predictor + dUdt_corrector)
        U[inY, inX, :] = U[inY, inX, :] + dt * dUdt[inY, inX, :]  # interior points only
        primitives = decodeSolutionVector(U)
        primitives = updateBoundaryConditions(primitives, inflow, Tw_Tinf)

        # check density convergence
        rCurrent = primitives.r
        if i > 1:
            deltaR = np.max(np.max(np.abs(rCurrent - rLast)))
            if deltaR < 1e-8:
                converged = True
            print('Iteration: ' + i + '| delta rho:' + deltaR)
            # fprintf(1, 'Iteration:%5d | delta rho: %8e\n', i, deltaR)
        rLast = rCurrent
    #  runtime = time.time() - start_time

    # Mass Flow Check
    massIn = np.trapz(y[:, 0], primitives.u[:, 0] * primitives.r[:, 0])
    massOut = np.trapz(y[:, -1], primitives.u[:, -1] * primitives.r[:, -1])
    massDiffCheck = 100 * abs(massIn - massOut) / massIn
    print('Mass inflow matches mass outflow within ')
    # fprintf(1, 'Mass inflow matches mass outflow within %.3f%%.\n', massDiffCheck)
    # fprintf(1, 'Runtime: %.2f seconds.\n', runTime)

    return primitives, massDiffCheck, converged, i


def decodeSolutionVector(U):
    r = U[:, :, 0]
    u = U[:, :, 1] / r
    v = U[:, :, 2] / r
    Et = U[:, :, 3]
    e = Et / r - .5 * (u ** 2 + v ** 2)
    cv = Primitives.R / (Primitives.gm - 1)
    T = e / cv
    p = r * Primitives.R * T
    primitives = Primitives(u, v, p, T)
    return primitives


def updateBoundaryConditions(primitivesIn, inflow, Tw_Tinf):
    [u, v, p, T] = primitivesIn.deal()
    [Vinf, _, pinf, Tinf] = inflow.deal()
    u[:, 0] = Vinf
    v[:, 0] = 0
    p[:, 0] = pinf
    T[:, 0] = Tinf

    u[-1, :] = Vinf
    v[-1, :] = 0
    p[-1, :] = pinf
    T[-1, :] = Tinf

    u[:, -1] = 2 * u[:, -2] - u[:, -3]
    v[:, -1] = 2 * v[:, -2] - v[:, -3]
    p[:, -1] = 2 * p[:, -2] - p[:, -3]
    T[:, -1] = 2 * T[:, -2] - T[:, -3]
    u[0, :] = 0
    v[0, :] = 0
    p[0, :] = 2 * p[1, :] - p[2, :]

    if Tw_Tinf > 0:
        T[0, :] = Tinf * Tw_Tinf
    else:
        T[0, :] = T[1, :]

    u[0, 0] = 0
    v[0, 0] = 0
    p[0, 0] = pinf
    T[0, 0] = Tinf
    primitivesOut = Primitives(u, v, p, T)
    return primitivesOut


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
delta = 5 * lhori / np.sqrt(Reinf[0])
lvert = 5 * delta

# grid
x, y = np.meshgrid(np.linspace(0, lhori, nx), np.linspace(0, lvert, ny))

# Set initial conditions

# Set all intial values to inflow values. Conditions on boundaries updated by solveMacCormack()
primitives = Primitives(inflow.u * np.ones((ny, nx)), inflow.v * np.ones((ny, nx)), inflow.p * np.ones((ny, nx)),
                        inflow.T * np.ones((ny, nx)))

# Solve two wall temperature conditions
Tw_Tinf = 1.0 # constant wall temperature
constantTw = solveMacCormack(primitives, inflow, Tw_Tinf, K, x, y, maxiter)
