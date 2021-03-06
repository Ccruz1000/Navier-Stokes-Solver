import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import time
from Primitives import *


def ddxb(f, dx):
    ny, nx, n3 = np.shape(f)
    inX = np.linspace(1, nx, num=nx, endpoint=False, dtype=int)
    dfdx = np.zeros((ny, nx, n3))
    dfdx[:, inX, :] = (f[:, inX, :] - f[:, inX - 1, :]) / dx
    dfdx[:, 0, :] = dfdx[:, 1, :]

    return dfdx


def ddxc(f, dx):
    ny, nx, n3 = np.shape(f)
    inX = np.linspace(1, nx-1, num=nx-1, endpoint=False, dtype=int)
    dfdx = np.zeros((ny, nx, n3))
    dfdx[:, inX, :] = (f[:, inX + 1, :] - f[:, inX - 1, :]) / (2 * dx)
    dfdx[:, 0, :] = (f[:, 1, :] - f[:, 0, :]) / dx
    dfdx[:, nx - 1, :] = (f[:, nx - 1, :] - f[:, nx - 2, :]) / dx

    return dfdx


def ddxf(f, dx):
    ny, nx, n3 = np.shape(f)
    inX = np.linspace(0, nx-1, num=nx-1, endpoint=False, dtype=int)
    dfdx = np.zeros((ny, nx, n3))
    dfdx[:, inX, :] = (f[:, inX + 1, :] - f[:, inX, :]) / dx
    dfdx[:, nx - 1, :] = dfdx[:, nx - 2, :]

    return dfdx


def ddyb(f, dy):
    ny, nx, n3 = np.shape(f)
    inY = np.linspace(1, ny, num=ny, endpoint=False, dtype=int)
    dfdy = np.zeros((ny, nx, n3))
    dfdy[inY, :, :] = (f[inY, :, :] - f[inY - 1, :, :]) / dy
    dfdy[0, :, :] = dfdy[1, :, :]

    return dfdy


def ddyc(f, dy):
    ny, nx, n3 = np.shape(f)
    inY = np.linspace(1, ny-1, num=ny-1, endpoint=False, dtype=int)
    dfdy = np.zeros((ny, nx, n3))
    dfdy[inY, :, :] = (f[inY + 1, :, :] - f[inY - 1, :, :]) / (2 * dy)
    dfdy[0, :, :] = (f[1, :, :] - f[0, :, :]) / dy
    dfdy[ny - 1, :, :] = (f[ny - 1, :, :] - f[ny - 2, :, :]) / dy
    return dfdy


def ddyf(f, dy):
    ny, nx, n3 = np.shape(f)
    inY = np.linspace(0, ny-1, num=ny-1, endpoint=False, dtype=int)
    dfdy = np.zeros((ny, nx, n3))
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
    ny, nx, n3 = np.shape(r)
    E = np.zeros((ny, nx, 4))
    E[:, :, 0] = r[:, :, 0] * u[:, :, 0]
    E[:, :, 1] = r[:, :, 0] * u[:, :, 0] ** 2 + p[:, :, 0] - txx[:, :, 0]
    E[:, :, 2] = r[:, :, 0] * u[:, :, 0] * v[:, :, 0] - txy[:, :, 0]
    E[:, :, 3] = (Et[:, :, 0] + p[:, :, 0]) * u[:, :, 0] - u[:, :, 0] * txx[:, :, 0] - v[:, :, 0] * txy[:, :, 0] + qx[:, :, 0]
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
    ny, nx, n3 = np.shape(r)
    F = np.zeros((ny, nx, 4))
    F[:, :, 0] = r[:, :, 0] * v[:, :, 0]
    F[:, :, 1] = r[:, :, 0] * u[:, :, 0] * v[:, :, 0] - txy[:, :, 0]
    F[:, :, 2] = r[:, :, 0] * v[:, :, 0] ** 2 + p[:, :, 0] - tyy[:, :, 0]
    F[:, :, 3] = (Et[:, :, 0] + p[:, :, 0]) * v[:, :, 0] - u[:, :, 0] * txy[:, :, 0] - v[:, :, 0] * tyy[:, :, 0] + qy[:, :, 0]

    return F


def calculateU(primitives):
    numy, numx, n3 = np.shape(primitives.u)
    U = np.zeros((numx, numy, 4))
    U[:, :, 0] = primitives.r[:, :, 0]
    U[:, :, 1] = U[:, :, 0] * primitives.u[:, :, 0]
    U[:, :, 2] = U[:, :, 0] * primitives.v[:, :, 0]
    U[:, :, 3] = primitives.Et[:, :, 0]
    return U


def postStressAndHeatFlux(primitives, x, y):
    dx = x[0, 1] - x[0, 0]
    dy = y[1, 0] - y[0, 0]
    [u, v, _, T] = primitives.deal()
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

    start_time = time.time()
    while not converged and i < maxiter:
        i += 1

        # time step size
        dt = primitives.calculateTimeStep(dx, dy, K)
        # dt = 9e-15

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
            if deltaR < 1e-16:
                converged = True
            # elif np.isnan(deltaR):
            #     raise ValueError('Delta R failed to converge')
            print('Iteration: ' + str(i) + ' | delta rho:' + str(deltaR))
            # fprintf(1, 'Iteration:%5d | delta rho: %8e\n', i, deltaR)
        rLast = rCurrent
    runtime = time.time() - start_time
    # print('Runtime: ' + str(runtime) ' Seconds')
    # Mass Flow Check
    # massIn = np.trapz(y[:, 0], primitives.u[:, 0] * primitives.r[:, 0])
    # massOut = np.trapz(y[:, -1], primitives.u[:, -1] * primitives.r[:, -1])
    # massDiffCheck = 100 * abs(massIn - massOut) / massIn
    # print('Mass inflow matches mass outflow within ')
    # fprintf(1, 'Mass inflow matches mass outflow within %.3f%%.\n', massDiffCheck)
    # fprintf(1, 'Runtime: %.2f seconds.\n', runTime)

    return primitives,  converged, deltaR, runtime #, i  massDiffCheck,


def decodeSolutionVector(U):
    r = U[:, :, 0]
    u = U[:, :, 1] / r
    v = U[:, :, 2] / r
    Et = U[:, :, 3]
    numx = np.shape(r[0])
    numy = np.shape(r[1])
    r = r.reshape((r.shape[0], r.shape[1], 1))
    u = u.reshape((u.shape[0], u.shape[1], 1))
    v  = v.reshape((v.shape[0], v.shape[1], 1))
    Et  = Et.reshape((Et.shape[0], Et.shape[1], 1))
    e = Et / r - .5 * (u ** 2 + v ** 2)
    prim_val = Primitives(1, 1, 1, 1)
    cv = prim_val.R / (prim_val.gm - 1)
    T = e / cv
    p = r * prim_val.R * T
    primitives = Primitives(u, v, p, T)
    number = 1
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
# K = np.linspace(0.6, 0.8, 200)
# K = 0.8
# K = [0.6271356783919598, 0.6301507537688442, 0.6613065326633166, 0.6623115577889447, 0.6653266331658292, 0.699497487437186, 0.700502512562814, 0.7427135678391961, 0.7437185929648241, 0.7447236180904523, 0.7929648241206031, 0.7939698492462313, 0.7979899497487437]
# K = [0.6271356783919598, 0.6301507537688442, 0.6613065326633166, 0.6623115577889447, 0.6653266331658292, 0.699497487437186, 0.700502512562814, 0.7437185929648241, 0.7447236180904523, 0.7929648241206031, 0.7939698492462313, 0.7979899497487437]
# K = [0.6271356783919598, 0.6301507537688442, 0.6613065326633166, 0.6623115577889447, 0.6653266331658292, 0.699497487437186, 0.700502512562814, 0.7437185929648241, 0.7447236180904523, 0.7929648241206031, 0.7939698492462313, 0.7979899497487437]
# K = [0.6271356783919598, 0.6301507537688442, 0.6613065326633166, 0.6623115577889447, 0.6653266331658292, 0.699497487437186, 0.700502512562814, 0.7437185929648241, 0.7447236180904523, 0.7929648241206031, 0.7939698492462313, 0.7979899497487437]
K = [0.7437185929648241]
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
primitives = Primitives(inflow.u * np.ones((ny, nx, 1)), inflow.v * np.ones((ny, nx, 1)), inflow.p * np.ones((ny, nx, 1)),
                        inflow.T * np.ones((ny, nx, 1)))

# Solve two wall temperature conditions
Tw_Tinf = 1.0 # constant wall temperature
courant = []
for k in K:
    constantTw = solveMacCormack(primitives, inflow, Tw_Tinf, k, x, y, maxiter)
    if not np.isnan(constantTw[-2]):
        courant.append(k)
        print('Runtime: ' + str(constantTw[-1]) + ' seconds')
        plt.figure(1)
        plt.contourf(x, y, constantTw[0].u[:, :, 0])
        plt.colorbar()
        plt.show()
print(courant)
# Succesful Courants first attempt iter = 50
# [0.6271356783919598, 0.6301507537688442, 0.6613065326633166, 0.6623115577889447, 0.6653266331658292, 0.699497487437186, 0.700502512562814, 0.7427135678391961, 0.7437185929648241, 0.7447236180904523, 0.7929648241206031, 0.7939698492462313, 0.7979899497487437]
# Succesful Courants second attempt iter = 500
# [0.6271356783919598, 0.6301507537688442, 0.6613065326633166, 0.6623115577889447, 0.6653266331658292, 0.699497487437186, 0.700502512562814, 0.7437185929648241, 0.7447236180904523, 0.7929648241206031, 0.7939698492462313, 0.7979899497487437]
# Successful courants third attempt iter = 1000
# [0.6271356783919598, 0.6301507537688442, 0.6613065326633166, 0.6623115577889447, 0.6653266331658292, 0.699497487437186, 0.700502512562814, 0.7437185929648241, 0.7447236180904523, 0.7929648241206031, 0.7939698492462313, 0.7979899497487437]
# Successful courants fourth attempt iter = 5000
# [0.6271356783919598, 0.6301507537688442, 0.6613065326633166, 0.6623115577889447, 0.6653266331658292, 0.699497487437186, 0.700502512562814, 0.7437185929648241, 0.7447236180904523, 0.7929648241206031, 0.7939698492462313, 0.7979899497487437]
# Successful courants fifth attempt iter = 10000
# [0.6271356783919598, 0.6301507537688442, 0.6613065326633166, 0.6623115577889447, 0.6653266331658292, 0.699497487437186, 0.700502512562814, 0.7437185929648241, 0.7447236180904523, 0.7929648241206031, 0.7939698492462313, 0.7979899497487437]
