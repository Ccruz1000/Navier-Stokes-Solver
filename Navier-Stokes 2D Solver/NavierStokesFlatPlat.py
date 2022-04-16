import matplotlib.pyplot as plt
import numpy as np


# Import Packages
import numpy as np
import matplotlib.pyplot as plt
import scipy


# Calculate viscosity based on Sutherlands Law
def DYNVIS(T, mu_0, T0):
    mu = mu_0 * (T / T0) ** (3 / 2)*(T0 + 110) / (T + 110)
    return mu


# Calculate thermal conductivity from Prandtl number
def THERMC(Pr, c_p, mu):
    k = mu*c_p/Pr
    return k


# Calculate x-direction of normal stress
def TAUXX(u, v, lam, mu, DX, DY, call_case):
    num_x, num_y = np.shape(u)
    # Initialize arrays
    du_dx = np.zeros_like(u)
    dv_dy = np.zeros_like(u)
    # Calculate derivative of u wrt x and v wrt to y
    # Calculate du_dx
    if call_case == 'Predict_E':
        for i in range(1, num_x):
            for j in range(0, num_y):
                du_dx[j, i] = (u[j, i] - u[j, i - 1]) / DX  # Backward differencing
        du_dx[:, 0] = (u[:, 1] - u[:, 0]) / DX  # Forward difference at i = 0
    elif call_case == 'Correct_E':
        for i in range(1, num_x):
            for j in range(0, num_y):
                du_dx[j, i] = (u[j, i + 1] - u[j, i]) / DX  # Forward differencing
        du_dx[:, num_x] = (u[:, num_x] - u[:, num_x - 1]) / DX  # Backward difference at i = numx
    # Compute dv_dy
    for i in range(0, num_x):
        for j in range(1, num_y):
            dv_dy[j, i] = (v[j + 1, i] - v[j - 1, i]) / (2 * DY)    # Central difference
        dv_dy[0, :] = (v[1, :] - v[0, :]) / DY  # Forward at i = 0
        dv_dy[num_y, :] = (v[num_y, :] - v[num_y - 1, :]) / DY  # Backward at j = numy
    tau_xx = lam * (du_dx + dv_dy) * 2 * mu * du_dx
    return tau_xx


def TAUXY(u, v, lam, mu, DX, DY, call_case):
    num_x, num_y = np.shape(u)
    # Initialize arrays
    du_dy = np.zeros_like(u)
    dv_dx = np.zeros_like(u)
    # Calculate derivative of u wrt x and v wrt to y
    # Calculate du_dy
    if call_case == 'Predict_E' or 'Correct_E':
        for i in range(0, num_x):
            for j in range(1, num_y - 1):
                du_dy[j, i] = (u[j + 1, i] - u[j - 1, i]) / (2 * DY)  # Backward differencing
        du_dy[0, :] = (u[1, :] - u[0, :]) / DY  # Forward difference at j = 0
        du_dy[num_y, :] = (u[num_y, :] - u[num_y - 1, :]) / DY  # Backward at j = numy
        # The calculation of dv_dx is different for each case
        if call_case == 'Predict_E':
            for i in range(1, num_x):
                for j in range(0, num_y):
                    dv_dx[j, i] = (v[j, i] - u[j, i - 1]) / DX  # Backward differencing
            dv_dx[:, 0] = (v[:, 1] - u[:, 0]) / DX  # Forward difference at i = 0
        else:
            for i in range(0, num_x - 1):
                for j in range(0, num_y):
                    dv_dx[i, j] = (v[j, i + 1] - v[j, i]) / DX  # Forward difference
                dv_dx[:, num_x] = (v[:, num_x] - v[:, num_x - 1]) / DX  # Backward difference at i = numx
    if call_case == 'Predict_F' or 'Correct_F':
        # The calculation of dv_dx is the same for both cases
        for i in range(0, num_x):
            for j in range(1, num_y):
                dv_dx[j, i] = (v[j, i + 1] - v[j, i - 1]) / (2 * DX)    # Central Difference
        dv_dx[:, 0] = (v[:, 1] - v[:, 1]) / DX  # Forward difference at i = 1
        dv_dx[:, num_x] = (v[:, num_x] - v[:, num_x - 1]) / DX  # Backward differance at i = numx

        # The calculation of du_dy is different for each case
        if call_case == 'Predict_F':
            for i in range(0, num_x):
                for j in range(1, num_y):
                    du_dy[j, i] = (u[j, i] - u[j - 1, i]) / DY  # Backward difference
                du_dy[0, :] = (u[1, :] - u[0, :]) / DY  # Forward at j = 1
        else:
            for i in range(0, num_x):
                for j in range(0, num_y - 1):
                    du_dy[j, i] = (u[j + 1, i] - u[j, i]) / DY  # Forward difference
                du_dy[num_y, :] = (u[num_y, :] - u[num_y, :]) / DY  # Backward difference at j = numy
    tau_xy = mu * (du_dy + dv_dx)
    return tau_xy



# Calculate flux vector E from primitive flow field variables
def Primitive2E(rho, u, p, v, T, mu, lam, k, c_v, dx, dy, call_case):
    E1 = rho * u    # Continuity
    tau_xx = TAUXX(u, v, lam, mu, dx, dy, call_case)
    E2 = rho * u ** 2 + p - tau_xx  # X-momentum
    tau_xy = TAUXY(u, v, lam, mu, dx, dy, call_case)
    E3 = rho * u * v - tau_xy # Y-momentum
    # Ended on line 27 of primitive2e


def main():
    # Problem Setup
    M_inf = 4   # Free stream mach number
    T_inf = 288.16  # Free stream temperature (K)
    p_inf = 101325  # Free stream pressure (pa)
    MAXIT = 1e4  # Number of time steps allowed before stopping program
    t = 1   # Time step counter
    time = 0    # Initial time (s)


    # Define Domain
    num_x = 65  # Number of X points
    num_y = 65  # Number of Y points
    x_end = 0.3048    # Final Point in x
    x_start = 0.05   # First point in x (absolute value, shifted to left on x axis)
    y_end = 0.05    # Final Point in y
    y_start = 0

    y_refine = 2    # Value to refine mesh in y direction
    x_refine = 2    # Value to refine mesh in x direction

    y_coarse = np.linspace(y_start**(1/y_refine), y_end**(1/y_refine), num_y)    # Mesh divisions in Y
    x_coarse = np.linspace(-x_start**(1/x_refine), x_end**(1/x_refine), num_x)    # Mesh divisions in x

    # Refine Mesh
    x = np.sign(x_coarse)*x_coarse**x_refine
    y = np.sign(y_coarse)*y_coarse**y_refine

    # Calculate dx and dy
    dx = np.zeros_like(x)
    dy = np.zeros_like(y)
    for i in range(len(x)):
        if i == 0:
            dx[i] = x[i + 1] - x[i]
        else:
            dx[i] = x[i] - x[i-1]
    for i in range(len(y)):
        if i == 0:
            dy[i] = y[i + 1] - y[i]
        else:
            dy[i] = y[i] - y[i-1]

    X, Y = np.meshgrid(x, y)    # Mesh

    # Constant Parameters
    T_w_T_inf = 1   # Ratio of wall temperature to free stream
    gamma = 1.4  # Air specific heat ratio
    R = 287  # Air gas constant (J/kgK)
    a_inf = np.sqrt(gamma*R*T_inf)  # Free stream speed of sound (m/s)
    mu_0 = 1.7894e-4    # Reference viscosity (kg/ms)
    T_0 = 288.16    # Reference temperature (K)
    Pr = 0.71   # Prandtl Number
    c_v = R / (gamma - 1)   # Specific heat at constant volume
    c_p = gamma * c_v   # Specific heat at constant pressure
    rho_inf = p_inf / (R * T_inf)   # Free stream density (kg/m^3)
    Re = rho_inf * M_inf * a_inf * x_end / DYNVIS(T_inf, mu_0, T_0)

    # Define initial conditions
    p = np.ones((num_x, num_y)) * p_inf  # Initialize pressure to be ambient pressure
    rho = np.ones((num_x, num_y)) * rho_inf  # Initialize density to be ambient density
    T = np.ones((num_x, num_y)) * T_inf  # Initialize temperature to be ambient temperature
    u = np.ones((num_x, num_y)) * M_inf * a_inf  # Initialize u velocity array
    v = np.zeros((num_x, num_y))    # Initialize v velocity array
    for i in range(num_x):
        if x[i] > 0:
            T[0, i] = T_w_T_inf * T_inf  # Constant wall temperature boundary condition
            u[0, i] = 0  # No slip boundary condition
    mu = DYNVIS(T, mu_0, T_0)  # Initialize dynamic viscosity
    lam = - 2 / 3 * mu  # Second viscosity (based on Stokes' Hypothesis [kg/m*S]
    k = THERMC(Pr, c_p, mu)

    # Other required variables
    U1_p = np.zeros((num_x, num_y))
    U2_p = np.zeros((num_x, num_y))
    U3_p = np.zeros((num_x, num_y))
    u5_p = np.zeros((num_x, num_y))
    rho_p = np.zeros((num_x, num_y))
    u_p = np.zeros((num_x, num_y))
    v_p = np.zeros((num_x, num_y))
    T_p = np.zeros((num_x, num_y))
    p_p = np.zeros((num_x, num_y))

    # Begin loop
    converged = False
    rho_old = rho

    while not converged and t <= MAXIT:
        # Apply MacCormack's Technique to compute solution vector U
        U1 = rho    # Continuity
        U2 = rho * u  # X-momentum
        U3 = rho * v    # Y-momentum
        U5 = rho * (c_v * T + (u ** 2 + v ** 2) / 2)    # Energy

        # Predict Step



    plt.scatter(y, dy)
    plt.show()


if __name__== "__main__":
  main()