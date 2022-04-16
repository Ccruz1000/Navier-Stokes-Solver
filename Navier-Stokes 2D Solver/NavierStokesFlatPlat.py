# Import Packages
import numpy as np
import matplotlib.pyplot as plt
import scipy

# Calculate viscosity based on Sutherlands Law
def DYNVIS(T, mu_0, T0):
    mu = mu_0 * (T / T0) ** (3 / 2)*(T0 + 110) / (T + 110)
    return mu

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
    rho_inf = p_inf / (R * T_inf)   # Freestream density (kg/m^3)
    Re = rho_inf * M_inf * a_inf * x_end / DYNVIS(T_inf, mu_0, T_0)

    # Define initial conditions
    p = np.ones((x, y)) * p_inf # Initialize pressure to be ambient pressure
    rho = np.ones((x, y)) * rho_inf  # Initialize density to be ambient density
    T = np.ones((x, y)) * T_inf # Initialize temperature to be ambient temperature




    plt.scatter(y, dy)
    plt.show()


if __name__== "__main__":
  main()