import matplotlib.pyplot as plt
import numpy as np

# Import Packages
import numpy as np
import matplotlib.pyplot as plt
import scipy


# Calculate viscosity based on Sutherlands Law
def DYNVIS(T, mu_0, T0):
    mu = mu_0 * (T / T0) ** (3 / 2) * (T0 + 110) / (T + 110)
    return mu


# Calculate thermal conductivity from Prandtl number
def THERMC(Pr, c_p, mu):
    k = mu * c_p / Pr
    return k


# Calculate x-direction of normal stress
def TAUXX(u, v, lam, mu, DX, DY, call_case):
    num_y, num_x = np.shape(u)
    # Initialize arrays
    du_dx = np.zeros_like(u)
    dv_dy = np.zeros_like(u)
    # Calculate derivative of u wrt x and v wrt to y
    # Calculate du_dx
    if call_case == 'Predict_E':
        for i in range(1, num_x-1):
            for j in range(0, num_y-1):
                du_dx[j, i] = (u[j, i] - u[j, i - 1]) / DX  # Backward differencing
        du_dx[:, 0] = (u[:, 1] - u[:, 0]) / DX  # Forward difference at i = 0
    elif call_case == 'Correct_E':
        for i in range(0, num_x - 2):
            for j in range(0, num_y-1):
                du_dx[j, i] = (u[j, i + 1] - u[j, i]) / DX  # Forward differencing
        du_dx[:, num_x - 1] = (u[:, num_x - 1] - u[:, num_x - 2]) / DX  # Backward difference at i = numx
    # Compute dv_dy
    for i in range(0, num_x-1):
        for j in range(1, num_y - 2):
            dv_dy[j, i] = (v[j + 1, i] - v[j - 1, i]) / (2 * DY)  # Central difference
        dv_dy[0, :] = (v[1, :] - v[0, :]) / DY  # Forward at i = 0
        dv_dy[num_y - 1, :] = (v[num_y - 1, :] - v[num_y - 2, :]) / DY  # Backward at j = numy
    tau_xx = lam * (du_dx + dv_dy) * 2 * mu * du_dx
    return tau_xx


# Calculate normal stress in X-Y direction
def TAUXY(u, v, mu, DX, DY, call_case):
    num_y, num_x  = np.shape(u)
    # Initialize arrays
    du_dy = np.zeros_like(u)
    dv_dx = np.zeros_like(u)
    # Calculate derivative of u wrt x and v wrt to y
    # Calculate du_dy
    if call_case == 'Predict_E' or 'Correct_E':
        for i in range(0, num_x-1):
            for j in range(1, num_y - 2):
                du_dy[j, i] = (u[j + 1, i] - u[j - 1, i]) / (2 * DY)  # Backward differencing
        du_dy[0, :] = (u[1, :] - u[0, :]) / DY  # Forward difference at j = 0
        du_dy[num_y - 1, :] = (u[num_y - 1, :] - u[num_y - 2, :]) / DY  # Backward at j = numy
        # The calculation of dv_dx is different for each case
        if call_case == 'Predict_E':
            for i in range(1, num_x-1):
                for j in range(0, num_y-1):
                    dv_dx[j, i] = (v[j, i] - v[j, i - 1]) / DX  # Backward differencing
            dv_dx[:, 0] = (v[:, 1] - v[:, 0]) / DX  # Forward difference at i = 0
        else:
            for i in range(0, num_x - 2):
                for j in range(0, num_y-1):
                    dv_dx[j, i] = (v[j, i + 1] - v[j, i]) / DX  # Forward difference
                dv_dx[:, num_x - 1] = (v[:, num_x - 1] - v[:, num_x - 2]) / DX  # Backward difference at i = numx
    elif call_case == 'Predict_F' or 'Correct_F':
        # The calculation of dv_dx is the same for both cases
        for i in range(1, num_x - 2):
            for j in range(0, num_y-1):
                dv_dx[j, i] = (v[j, i + 1] - v[j, i - 1]) / (2 * DX)  # Central Difference
        dv_dx[:, 0] = (v[:, 1] - v[:, 0]) / DX  # Forward difference at i = 1
        dv_dx[:, num_x - 1] = (v[:, num_x - 1] - v[:, num_x - 2]) / DX  # Backward differance at i = numx

        # The calculation of du_dy is different for each case
        if call_case == 'Predict_F':
            for i in range(0, num_x-1):
                for j in range(1, num_y-1):
                    du_dy[j, i] = (u[j, i] - u[j - 1, i]) / DY  # Backward difference
                du_dy[0, :] = (u[1, :] - u[0, :]) / DY  # Forward at j = 1
        else:
            for i in range(0, num_x-1):
                for j in range(0, num_y - 2):
                    du_dy[j, i] = (u[j + 1, i] - u[j, i]) / DY  # Forward difference
                du_dy[num_y - 1, :] = (u[num_y - 1, :] - u[num_y - 2, :]) / DY  # Backward difference at j = numy
    tau_xy = mu * (du_dy + dv_dx)
    return tau_xy


# Calculate normal stress in Y direction
def TAUYY(u, v, lam, mu, DX, DY, call_case):
    # Initialize arrays
    numy, numx = np.shape(v)
    du_dx = np.zeros_like(v)
    dv_dy = np.zeros_like(v)
    # Calculate dv_dy
    if call_case == 'Predict_F':
        for i in range(0, numx-1):
            for j in range(1, numy-1):
                dv_dy[j, i] = (v[j, i] - v[j - 1, i]) / DY  # Backward difference
        dv_dy[0, :] = (v[1, :] - v[0, :]) / DY  # Forward difference at j = 0
    elif call_case == 'Correct_F':
        for i in range(0, numx-1):
            for j in range(0, numy - 2):
                dv_dy[j, i] = (v[j + 1, i] - v[j, i]) / DY  # Forward difference
        dv_dy[numy - 1, :] = (v[numy - 1, :] - v[numy - 2, :]) / DY  # Backward difference at j = numy
    # Calculate du_dx
    for i in range(1, numx - 2):
        for j in range(0, numy-1):
            du_dx[j, i] = (u[j, i + 1] - u[j, i - 1]) / (2 * DX)  # Central difference
    du_dx[:, 0] = (u[:, 1] - u[:, 0]) / DX  # Forward difference at i = 1
    du_dx[:, numx - 1] = (u[:, numx - 1] - u[:, numx - 2]) / DX  # Backward difference at i = numx
    # Calculate tau_yy
    tau_yy = lam * (du_dx + dv_dy) + 2 * mu * dv_dy
    return tau_yy


# Calculate x component of heat flux vector
def QX(T, k, DX, call_case):
    # Initialize arrays
    numy, numx = np.shape(T)
    dT_dx = np.zeros_like(T)
    if call_case == 'Predict_E':
        for i in range(1, numx-1):
            for j in range(0, numy-1):
                dT_dx[j, i] = (T[j, i] - T[j, i - 1]) / DX  # Backward difference
        dT_dx[:, 0] = (T[:, 1] - T[:, 0]) / DX  # Forward difference at i = 1
    elif call_case == 'Correct_E':
        for i in range(0, numx - 2):
            for j in range(0, numy-1):
                dT_dx[j, i] = (T[j, i + 1] - T[j, i]) / DX  # Forward difference
        dT_dx[:, numx - 1] = (T[:, numx - 1] - T[:, numx - 2]) / DX  # Backward difference at i = numx
    q_x = -1 * k * dT_dx
    return q_x


# Calculate y component of heat flux vector
def QY(T, k, DY, call_case):
    # Initialize arrays
    numy, numx = np.shape(T)
    dT_dy = np.zeros_like(T)
    if call_case == 'Predict_F':
        for i in range(0, numx-1):
            for j in range(1, numy-1):
                dT_dy[j, i] = (T[j, i] - T[j - 1, i]) / DY  # Backward difference
        dT_dy[0, :] = (T[1, :] - T[0, :]) / DY  # Forward difference at j = 1
    elif call_case == 'Correct_F':
        for i in range(0, numx-1):
            for j in range(0, numy - 2):
                dT_dy[j, i] = (T[j + 1, i] - T[j, i]) / DY  # Forward difference
        dT_dy[numy - 1, :] = (T[numy - 1, :] - T[numy - 2, :]) / DY  # Backward dfference at j = numy
    q_y = -1 * k * dT_dy
    return q_y


# Calculate flux vector E from primitive flow field variables
def Primitive2E(rho, u, p, v, T, mu, lam, k, c_v, dx, dy, call_case):
    E1 = rho * u  # Continuity
    tau_xx = TAUXX(u, v, lam, mu, dx, dy, call_case)
    E2 = rho * u ** 2 + p - tau_xx  # X-momentum
    tau_xy = TAUXY(u, v, mu, dx, dy, call_case)
    E3 = rho * u * v - tau_xy  # Y-momentum
    q_x = QX(T, k, dx, call_case)
    E5 = (rho * (c_v * T + (u ** 2 + v ** 2) / 2) + p) * u - u * tau_xx - v * tau_xy + q_x  # From energy
    return E1, E2, E3, E5


# Calculate flux vector F from primitive flow field variables
def Primitive2F(rho, u, p, v, T, mu, lam, k, c_v, DX, DY, call_case):
    F1 = rho * v  # Continuity
    tau_yx = TAUXY(u, v, mu, DX, DY, call_case)
    F2 = rho * u * v - tau_yx  # From X-momentum
    tau_yy = TAUYY(u, v, lam, mu, DX, DY, call_case)
    F3 = rho * v ** 2 + p - tau_yy  # Y-momentum
    q_y = QY(T, k, DY, call_case)
    F5 = (rho * (c_v * T + (u ** 2 + v ** 2) / 2) + p) * v - u * tau_yx - v * tau_yy + q_y
    return F1, F2, F3, F5


# Decode primitive flow field variables
def U2Primitive(U1, U2, U3, U5, c_v):
    rho = U1
    u = U2 / U1
    v = U3 / U1
    T = (U5 / U1 - ((U2 / U1) ** 2 + (U3 / U1) ** 2) / 2) / c_v
    return rho, u, v, T


# Apply Boundary conditions
def BC(rho, u, v, p, T, rho_inf, u_inf, p_inf, T_inf, T_w_T_inf, R, x):
    numy, numx = np.shape(rho)

    # Case 1:
    T[0, 0] = T_inf
    p[0, 0] = p_inf
    rho[0, 0] = rho_inf

    # Case 2:

    u[1: numy-1, 0] = u_inf
    # v(1: numy-1, 0) = 0
    p[1: numy-1, 0] = p_inf
    T[1: numy-1, 0] = T_inf
    rho[1: numy-1, 0] = rho_inf

    # For Upper boundary

    u[numy-1, 1: numx-1] = u_inf
    # v[numy-1, 1: numx-1] = 0
    p[numy-1, 1: numx-1] = p_inf
    T[numy-1, 1: numx-1] = T_inf
    rho[numy-1, 1: numx-1] = rho_inf

    # There are four cases -> Before plate, on plate, inlet, outlet/outlet
    # Case 4 -> Outlet: Temperature is calculated based on extrapolation from adjacent interior points,
    # Density is computed from equation of state
    u[1:numy - 2, numx - 1] = 2 * u[1:numy - 2, numx - 2] - u[1:numy - 2, numx - 3]
    v[1:numy - 2, numx - 1] = 2 * v[1:numy - 2, numx - 2] - v[1:numy - 2, numx - 3]
    p[1:numy - 2, numx - 1] = 2 * p[1:numy - 2, numx - 2] - p[1:numy - 2, numx - 3]
    T[1:numy - 2, numx - 1] = 2 * T[1:numy - 2, numx - 2] - T[1:numy - 2, numx - 3]
    rho[1:numy - 2, numx - 1] = p[1:numy - 2, numx - 1] / (R * T[1:numy - 2, numx - 1])

    # Case 3 - > Free stream for inlet

    if np.all(x):
        T[0, 1: numx-1] = T[1, 1: numx-1]
    else:
        T[0, 1: numx-1] = T_w_T_inf * T_inf

    p[0, 1: numx-1] = 2 * p[1, 1: numx-1] - p[2, 1: numx-1]
    rho[0, 1: numx-1] = p[0, 1: numx-1] / (R * T[0, 1:numx-1])

    return rho, u, v, p, T


# Check convergence
def CONVER(rho_old, rho):
    # Check for a converged solution. The converged criterion is that the change in density between
    # time steps is lower than 1e-14
    if not np.isreal(np.all(rho)):
        raise ValueError('The calculation has failed. A complex number has been detected')
    elif np.max(np.max(np.abs(rho_old - rho))) < 1e-8:
        converged = True
    else:
        converged = False
    return converged


# Check validity of numerical solution by confirming conservation of mass
# def MDOT(rho, u, rho_inf, u_inf, y_end, DY):
#     m_dot_in = -rho_inf * u_inf * (y_end - DY) - DY * rho_inf * u_inf / 2
#     numx = np.size([rho, 2])
#     y = 0, dy, y_end
#     m_dot_out = np.trapz(y, rho[:, numx-1] * u[:, numx-1])
#
#     deviation = np.abs((m_dot_in + m_dot_out) / m_dot_in) * 100
#     if deviation < 1:
#         continuity = True
#     else:
#         continuity = False
#     return continuity



# Problem Setup
M_inf = 4  # Free stream mach number
T_inf = 288.16  # Free stream temperature (K)
p_inf = 101325  # Free stream pressure (pa)
MAXIT = 1e4  # Number of time steps allowed before stopping program
t = 1  # Time step counter
time = 0  # Initial time (s)
dt = 1e-2  # time step

# Define Domain
num_x = 20  # Number of X points
num_y = 20  # Number of Y points
x_end = 0.3048  # Final Point in x
x_start = 0.05  # First point in x (absolute value, shifted to left on x axis)
y_end = 0.05  # Final Point in y
y_start = 0

y_refine = 2  # Value to refine mesh in y direction
x_refine = 2  # Value to refine mesh in x direction

y_coarse = np.linspace(y_start ** (1 / y_refine), y_end ** (1 / y_refine), num_y)  # Mesh divisions in Y
x_coarse = np.linspace(-x_start ** (1 / x_refine), x_end ** (1 / x_refine), num_x)  # Mesh divisions in x

# Refine Mesh
x = np.sign(x_coarse) * x_coarse ** x_refine
y = np.sign(y_coarse) * y_coarse ** y_refine

# Calculate dx and dy
dx = np.zeros_like(x)
dy = np.zeros_like(y)
for i in range(len(x)):
    if i == 0:
        dx[i] = x[i + 1] - x[i]
    else:
        dx[i] = x[i] - x[i - 1]
for i in range(len(y)):
    if i == 0:
        dy[i] = y[i + 1] - y[i]
    else:
        dy[i] = y[i] - y[i - 1]

X, Y = np.meshgrid(x, y)  # Mesh

# Constant Parameters
T_w_T_inf = 1  # Ratio of wall temperature to free stream
gamma = 1.4  # Air specific heat ratio
R = 287  # Air gas constant (J/kgK)
a_inf = np.sqrt(gamma * R * T_inf)  # Free stream speed of sound (m/s)
mu_0 = 1.7894e-4  # Reference viscosity (kg/ms)
T_0 = 288.16  # Reference temperature (K)
Pr = 0.71  # Prandtl Number
c_v = R / (gamma - 1)  # Specific heat at constant volume
c_p = gamma * c_v  # Specific heat at constant pressure
rho_inf = p_inf / (R * T_inf)  # Free stream density (kg/m^3)
Re = rho_inf * M_inf * a_inf * x_end / DYNVIS(T_inf, mu_0, T_0)
CFL = 0.6 # Courant-Friedrichs-Lewy number (should be between 0.5 and 0.8)
# Define initial conditions
p = np.ones((num_y, num_x)) * p_inf  # Initialize pressure to be ambient pressure
rho = np.ones((num_y, num_x)) * rho_inf  # Initialize density to be ambient density
T = np.ones((num_y, num_x)) * T_inf  # Initialize temperature to be ambient temperature
u = np.ones((num_y, num_x)) * M_inf * a_inf  # Initialize u velocity array
u[0, :] = 0  # No slip boundary condition
v = np.zeros((num_y, num_x))  # Initialize v velocity array
T[0, :] = T_w_T_inf * T_inf  # Constant wall temperature boundary condition
mu = DYNVIS(T, mu_0, T_0)  # Initialize dynamic viscosity
lam = - 2 / 3 * mu  # Second viscosity (based on Stokes' Hypothesis [kg/m*S]
k = THERMC(Pr, c_p, mu)

# Other required variables
U1_p = np.zeros((num_y, num_x))
U2_p = np.zeros((num_y, num_x))
U3_p = np.zeros((num_y, num_x))
U5_p = np.zeros((num_y, num_x))
rho_p = np.zeros((num_y, num_x))
u_p = np.zeros((num_y, num_x))
v_p = np.zeros((num_y, num_x))
T_p = np.zeros((num_y, num_x))
p_p = np.zeros((num_y, num_x))

# Begin loop
converged = False
rho_old = rho

while not converged and t <= MAXIT:
    print(t)
    # Time step needed to satisfy CFL stability criterion
    v_prime = np.max(np.max(4 / 3 * mu[1:num_y - 2, 1: num_x - 2] ** 2 * gamma / (Pr * rho[1:num_y - 2, 1:num_x-2])))
    delta_t_CFL = 1/(np.abs(u[1:num_y-2,1:num_x-2])/dx[1:num_y-2] + np.abs(v[1:num_y-2,1:num_x-2])/dy[1:num_x-2] + np.sqrt(gamma*R*T[1:num_y-2,1:num_x-2])*np.sqrt(1/dx[1:num_y-2]**2 + 1/dy[1:num_x-2]**2) + 2*v_prime*(1/dx[1:num_y-2]**2 + 1/dy[1:num_x-2]**2))
    dt = CFL * np.min(np.min(delta_t_CFL))

    # Apply MacCormack's Technique to compute solution vector U
    U1 = rho  # Continuity
    U2 = rho * u  # X-momentum
    U3 = rho * v  # Y-momentum
    U5 = rho * (c_v * T + (u ** 2 + v ** 2) / 2)  # Energy

    # Predict Step
    for i in range(len(dx)):
        for j in range(len(dy)):
            E1, E2, E3, E5 = Primitive2E(rho, u, p, v, T, mu, lam, k, c_v, dx[i], dy[j], 'Predict_E')
            F1, F2, F3, F5 = Primitive2F(rho, u, p, v, T, mu, lam, k, c_v, dx[i], dy[j], 'Predict_F')

    # Predict flow field properties at next time step for interior points using forward difference
    for i in range(1, num_x - 2):
        for j in range(1, num_y - 2):
            U1_p[j, i] = U1[j, i] - dt * ((E1[j, i + 1] - E1[j, i]) / dx[i] + (F1[j + 1, i] - F1[j, i]) / dy[j])
            U2_p[j, i] = U2[j, i] - dt * ((E2[j, i + 1] - E2[j, i]) / dx[i] + (F2[j + 1, i] - F2[j, i]) / dy[j])
            U3_p[j, i] = U3[j, i] - dt * ((E3[j, i + 1] - E3[j, i]) / dx[i] + (F3[j + 1, i] - F3[j, i]) / dy[j])
            U5_p[j, i] = U5[j, i] - dt * ((E5[j, i + 1] - E5[j, i]) / dx[i] + (F5[j + 1, i] - F5[j, i]) / dy[j])

    # Predict flow field variables needed for calculations
    # Calculate pressure with ideal gas law
    for i in range(1, num_x - 2):
        for j in range(1, num_y - 2):
            rho_p[j, i], u_p[j, i], v_p[j, i], T_p[j, i] = U2Primitive(U1_p[j, i], U2_p[j, i], U3_p[j, i], U5_p[j, i], c_v)
            p_p[j, i] = rho_p[j, i] * R * T_p[j, i]


    # Apply Boundary Conditions
    rho_p, u_p, v_p, p_p, T_p = BC(rho_p, u_p, v_p, p_p, T_p, rho_inf, M_inf * a_inf, p_inf, T_inf, T_w_T_inf, R, False)
    # for remaining properties
    mu_p = DYNVIS(T_p, mu_0, T_0)
    lam_p = -2 / 3 * mu_p
    k_p = THERMC(Pr, c_p, mu_p)

    # End of predictor step, now move to corrector step
    # Calculate flux vectors for E and F
    for i in range(len(dx)):
        for j in range(len(dy)):
            E1_p, E2_p, E3_p, E5_p = Primitive2E(rho_p, u_p, p_p, v_p, T_p, mu_p, lam_p, k_p, c_v, dx[i], dy[j],
                                                 'Correct_E')
            F1_p, F2_p, F3_p, F5_p = Primitive2F(rho_p, u_p, p_p, v_p, T_p, mu_p, lam_p, k_p, c_v, dx[i], dy[j],
                                                 'Correct_F')
    # Now we can predict flow field properties at next step for interior points using backwards differencing
    for i in range(1, num_x - 2):
        for j in range(1, num_y - 2):
            U1[j, i] = 1/2 * (U1[j, i] + U1_p[j, i] - dt * ((E1_p[j, i] - E1_p[j, i - 1]) / dx[i] + (
                    F1_p[j, i] - F1_p[j - 1, i]) / dy[j]))
            U2[j, i] = 1 / 2 * (U2[j, i] + U2_p[j, i] - dt * ((E2_p[j, i] - E2_p[j, i - 1]) / dx[i] + (
                    F2_p[j, i] - F2_p[j - 1, i]) / dy[j]))
            U3[j, i] = 1 / 2 * (U3[j, i] + U3_p[j, i] - dt * ((E3_p[j, i] - E3_p[j, i - 1]) / dx[i] + (
                    F3_p[j, i] - F3_p[j - 1, i]) / dy[j]))
            U5[j, i] = 1 / 2 * (U5[j, i] + U5_p[j, i] - dt * ((E5_p[j, i] - E5_p[j, i - 1]) / dx[i] + (
                    F5_p[j, i] - F5_p[j - 1, i]) / dy[j]))

    # Finally decode flow variables
    for i in range(1, num_x - 2):
        for j in range(1, num_y - 2):
            rho[j, i], u[j, i], v[j, i], T[j, i] = U2Primitive(U1[j, i], U2[j, i], U3[j, i], U5[j, i], c_v)
            p[j, i] = rho[j, i] * R * T[j, i]


    # rho[1:num_y - 2, 1:num_x - 2], u[1:num_y - 2, 1:num_x - 2], v[1:num_y - 2, 1:num_x - 2], \
    # T[1:num_y - 2, 1:num_x - 2] = U2Primitive(U1[1:num_y - 2, 1:num_x - 2], U2[1:num_y - 2, 1:num_x - 2],
    #                                           U3[1:num_y - 2, 1:num_x - 2], U5[1:num_y - 2, 1:num_x - 2], c_v)
    # # Use ideal gas law to calculate pressure
    # p[1:num_y - 2, 1:num_x - 2] = rho[1:num_y - 2, 1:num_x - 2] * R * T[1:num_y - 2, 1:num_x - 2]

    # Apply Boundary Conditions
    rho, u, v, p, T = BC(rho, u, v, p, T, rho_inf, M_inf*a_inf, p_inf, T_inf, T_w_T_inf, R, False)
    # Finally for remaining properties
    mu = DYNVIS(T, mu_0, T_0)
    lam = -2 / 3 * mu
    k = THERMC(Pr, c_p, mu)

    # This is the end of the corrector step. The flow field properties are known at each grid point for t + 1

    # Check convergence
    converged = CONVER(rho_old, rho)
    rho_old = rho
    # Change time
    time = time + dt
    t = t + 1
# if converged:
#     for j in range(len(dy)):
#         continuity = MDOT(rho, u, rho_inf, M_inf * a_inf, y_end, dy[j])
# else:
#     raise ValueError('Calculation Failed: The max number of iterations has been reached before '
#                      'achieving continuity')
# if continuity:
#     print('The calculation has finished succesfully')
# else:
#     raise ValueError('The solution has converged to an invalid result')

M = np.sqrt(u ** 2 + v ** 2) / np.sqrt(gamma * R * T)
plt.figure(1)
plt.contour(X, Y, M)
plt.title('Mach Contour')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.show()

fig = plt.figure(figsize=(11,7), dpi=100)
# plotting the pressure field as a contour
plt.contourf(X, Y, M, alpha=0.5, cmap=cm.viridis)
plt.colorbar()
# plotting the pressure field outlines
plt.contour(X, Y, M, cmap=cm.viridis)
# plotting velocity field
plt.quiver(X, Y, u, v)
plt.xlabel('X')
plt.ylabel('Y');
plt.show()
