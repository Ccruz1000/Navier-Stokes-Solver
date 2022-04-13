# Imported Packages
import numpy as np
import matplotlib.pyplot as plt

# User Defined Function

# Define Boundary
nx = 50  # Number of  points in x
ny = 50  # Number of points in y
x_start = -3  # Channel start in X
x_end = 3  # Channel end in X
y_start = -3  # Channel start in y
y_end = 3  # Channel end in y
x = np.linspace(x_start, x_end, nx)  # Discretize in x
y = np.linspace(y_start, y_end, ny)  # Discretize in y

# Define Shape
x_0, y_0, r = 0.0, 0.0, 0.5  # Define circle x-center, y-center and radius
X, Y = np.meshgrid(x, y)  # Define cavity mesh

# Define points within circle
pts = (X - x_0)**2 + (Y - y_0)**2 > r**2
X = X[pts]
Y = Y[pts]
plt.scatter(X, Y, s=0.1, c='r')
plt.show()
