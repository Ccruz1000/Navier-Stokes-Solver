# Imported Packages
import numpy as np
import matplotlib.pyplot as plt
import sys
np.set_printoptions(threshold=sys.maxsize)

# User Defined Function


def mesh_generate(geom, coord):
    if geom == 'airfoil':
        # Import airfoil
        airfoil_points = np.loadtxt(coord, skiprows=1, dtype='float')  # Import airfoil coordinates
        chord = airfoil_points[0, 0]
        # Close trailing edge
        if airfoil_points[0].all != airfoil_points[-1].all:
            closed_point = np.array([airfoil_points[0]])  # Create numpy array with first point to close
            airfoil_points = np.append(airfoil_points, closed_point, axis=0)  # Append original array with new array

        plt.plot(airfoil_points[:, 0], airfoil_points[:, 1], marker='x')
        plt.ylim(-0.5, 0.5)
        plt.show()
    elif geom == 'circle':
        # Define circle
        x_0, y_0, r = 0.0, 0.0, 0.5  # Define circle x-center, y-center and radius
        # Define X values to divide mesh domain
        coarse_div_x = 20
        fine_div_x = 75
        coarse_start_x = -3
        fine_start_x = -1.1 * r + x_0
        fine_end_x = 1.1 * r + x_0
        coarse_end_x = 3
        # Define Y values to divide mesh domain
        coarse_div_y = 10
        fine_div_y = 36
        coarse_start_y = -1.5
        fine_start_y = -1.1 * r + y_0
        fine_end_y = 1.1 * r + y_0
        coarse_end_y = 1.5
        # Concatenate domain to include coarse and fine meshing
        x = np.concatenate(
            [np.linspace(coarse_start_x, coarse_end_x, coarse_div_x), np.linspace(fine_start_x, fine_end_x,
                                                                                  fine_div_x),
             np.linspace(fine_end_x, coarse_end_x, coarse_div_x)])
        y = np.concatenate(
            [np.linspace(coarse_start_y, coarse_end_y, coarse_div_y), np.linspace(fine_start_y, fine_end_y,
                                                                                  fine_div_y),
             np.linspace(fine_end_y, coarse_end_y, coarse_div_y)])

        X, Y = np.meshgrid(x, y)  # Define complete domain mesh

        # Define points within circle
        pts_in = (X - x_0) ** 2 + (Y - y_0) ** 2 < r ** 2
        Xin = X[pts_in]
        Yin = Y[pts_in]
        # Define points outside of circle
        pts_out = (X - x_0) ** 2 + (Y - y_0) ** 2 > r ** 2
        Xout = X[pts_out]
        Yout = Y[pts_out]
        plt.scatter(Xout, Yout, s=0.02)
        # plt.scatter(Xin, Yin, c='r', s=0.02)
        plt.show()
    else:
        print('Check your input geometry')

    # TODO Determine how to return boundaries (ie object domain, fluid domain, walls, inlet, and outlet
    # TODO add airfoil functionality
