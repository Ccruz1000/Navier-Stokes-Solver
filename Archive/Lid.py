import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
nx = 41  # number of x points
ny = 41  # number of y points
x_length = 1.0  # domain size in x
y_length = 1.0  # domain size in y
nt = 500  # number of iterations
dt = 0.001
mu = 0.1  # kinematic viscosity
rho = 1.0  # density
top_u  = 1  # x velocity at top boundary

n_poisson = 50  # define number of poission iterations

def main():
    element_length = x_length
    x = np.linspace(0.0, x_length, nx)  # Generate x points array
    y = np.linspace(0.0, y_length, ny)  # Generate y points array

    X, Y = np.meshgrid(x, y)  # Generate computational domain

    # Initialize arrays to store previous values
    u_prev = np.zeros_like(X)
    v_prev = np.zeros_like(Y)
    p_prev = np.zeros_like(X)

    def central_difference_x(f):
        diff = np.zeros_like(f)
        diff[1:-1, 1:-1] = (f[1:-1, 2:] - f[1:-1, 0:-2]) / (2 * dx)
        return diff

    def central_difference_y(f):
        diff = np.zeros_like(f)
        diff[1:-1, 1:-1] = (f[2:, 1:-1] - f[0:-2, 1:-1]) / (2 * dy)
        return diff

    def laplace(f):
        diff = np.zeros_like(f)
        diff[1:-1, 1:-1]  = (f[1:-1, 0:-2] + f[0:-2, 1:-1] - 4 * f[1:-1, 1:-1] + f[1:-1, 2:] + f[2:, 1:-1]) / (dx ** 2)
        return diff

    for _ in tqdm(range(nt)):
        d_u_prev_dx = central_difference_x(u_prev)
        d_u_prev_dy = central_difference_y(u_prev)
        d_v_prev_dx = central_difference_x(v_prev)
        d_v_prev_dy = central_difference_y(v_prev)
        laplace_u_prev = laplace(u_prev)
        laplace_v_prev = laplace(v_prev)


        u_tent = (u_prev + dt * (u_prev * d_u_prev_dx + v_prev * d_u_prev_dy) + mu * laplace_u_prev)
        v_tent = (v_prev + dt * (u_prev * d_v_prev_dx + v_prev * d_v_prev_dy) + mu * laplace_v_prev)

        # Velocity Boundary Conditions
        u_tent[0, :] = 0.0
        u_tent[:, 0] = 0.0
        u_tent[: -1] = 0.0
        u_tent[-1, :] = top_u
        v_tent[0, :] = 0.0
        v_tent[:, 0] = 0.0
        v_tent[: -1] = 0.0
        v_tent[-1, :] = 0.0

        d_u_tent_dx = central_difference_x(u_tent)
        d_v_tent_dy = central_difference_y(v_tent)

        # Compute a pressure correction by solving pressure-poisson equation
        rhs = (rho / dt) * (d_u_tent_dx + d_v_tent_dy)

        for _ in range(n_poisson):
            p_next = np.zeros_like(p_prev)
            p_next[1:-1, 1:-1] = 1/4 * (p_prev[1:-1, 0:-2] + p_prev[0:-2, 1:-1] + p_prev[1:-1, 2:] + p_prev[2:, 1:-1] - dx**2 * rhs[1:-1, 1:-1])
            p_next[:, -1] = p_next[:, -2]
            p_next[0, :] = p_next[1, :]
            p_next[:, 0] = p_next [:, 1]
            p_next[-1, :] = 0.0

            p_prev = p_next

        d_p_next_dx = central_difference_x(p_next)
        d_p_next_dy = central_difference_y(p_next)

        # Correct velocity so fluid stays incompressible
        u_next = (u_tent - dt/rho * d_p_next_dx)
        v_next = (v_tent - dt/rho * d_p_next_dy)
        # Velocity Boundary Conditions
        u_next[0, :] = 0.0
        u_next[:, 0] = 0.0
        u_next[: -1] = 0.0
        u_next[-1, :] = top_u
        v_next[0, :] = 0.0
        v_next[:, 0] = 0.0
        v_next[: -1] = 0.0
        v_next[-1, :] = 0.0

        # Advance in time
        u_prev = u_next
        v_prev = v_next
        p_prev = p_next

        plt.style.use("dark_background")
        plt.figure()
        plt.contourf(X[::2, ::2], Y[::2, ::2], p_next[::2, ::2], cmap="coolwarm")
        plt.colorbar()

        plt.quiver(X[::2, ::2], Y[::2, ::2], u_next[::2, ::2], v_next[::2, ::2], color="black")
        # plt.streamplot(X[::2, ::2], Y[::2, ::2], u_next[::2, ::2], v_next[::2, ::2], color="black")
        plt.xlim((0, 1))
        plt.ylim((0, 1))
        plt.show()

if __name__ == '__main__':
    main()
