# Imported Packages
import numpy as np
import matplotlib.pyplot as plt
import sys
np.set_printoptions(threshold=sys.maxsize)

# User Defined Function


# Define Shape
x_0, y_0, r = 0.0, 0.0, 0.5  # Define circle x-center, y-center and radius
# Define X values to divide mesh domain
coarse_div_x = 40
fine_div_x = 150
coarse_start_x = -3
fine_start_x = -1.1*r + x_0
fine_end_x = 1.1*r + x_0
# Define Y values to divide mesh domain
coarse_div_y = 20
fine_div_y = 75
coarse_end_x = 3
coarse_start_y = -1.5
fine_start_y = -1.1*r + y_0
fine_end_y = 1.1*r + y_0
coarse_end_y = 1.5
# Concatenate domain to include coarse and fine meshing
x = np.concatenate([np.linspace(coarse_start_x, coarse_end_x, coarse_div_x), np.linspace(fine_start_x, fine_end_x,
                    fine_div_x), np.linspace(fine_end_x, coarse_end_x, coarse_div_x)])
y = np.concatenate([np.linspace(coarse_start_y, coarse_end_y, coarse_div_y), np.linspace(fine_start_y, fine_end_y,
                    fine_div_y), np.linspace(fine_end_y, coarse_end_y, coarse_div_y)])

X, Y = np.meshgrid(x, y)  # Define complete domain mesh

# Define points within circle

pts_out = (X - x_0)**2 + (Y - y_0)**2 > r**2
pts_in = (X - x_0)**2 + (Y - y_0)**2 < r**2
pts_on = (X - x_0)**2 + (Y - y_0)**2 == r**2
Xin = X[pts_in]
Yin = Y[pts_in]
Xout = X[pts_out]
Yout = Y[pts_out]
Xon = X[pts_on]
Yon = Y[pts_on]
# plt.scatter(Xout, Yout, s=0.01, c='r')
# plt.scatter(Xin , Yin, s=0.01, c='b')
plt.scatter(X, Y, s=0.01)
plt.show()
