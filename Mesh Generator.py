# Imported Functions
import numpy as np
import pygmsh
import gmsh
import meshio
import optimesh
import matplotlib.pyplot as plt

# User Defined Functions


# Chanel Parameters
Start = -1  # Start of rectangular channel
End = 4  # End of rectangular channel
Top = 1  # Top of rectangular channel
Bottom = -1  # Bottom of rectangular channel
resolution = 1e-2  # Mesh resolution
# Load airfoil from txt selig file
airfoil = np.loadtxt("naca0012h_selig.txt", skiprows=1, dtype='float')
# plt.plot(airfoil[:, 0], airfoil[:, 1])
# plt.ylim((-0.5, 0.5))
# plt.xlim((-0.1, 1.1))
# plt.show()
# # Initialize empty geometry using the built in kernel in gmsh
geometry = pygmsh.geo.Geometry()
# Fetch model we would like to add data to
model = geometry.__enter__()
points = []  # Initialize empty list to store airfoil points in
# Loop through airfoil points to generate points
for i in range(len(airfoil)):
    points.append(model.add_point((airfoil[i, 0], airfoil[i, 1], airfoil[i, 2]), mesh_size=resolution))
boundary_points = [model.add_point((Start, Bottom, 0), mesh_size=2*resolution),
                   model.add_point((End, Bottom, 0), mesh_size=6*resolution),
                   model.add_point((End, Top, 0), mesh_size=6*resolution),
                   model.add_point((Start, Top, 0), mesh_size=2*resolution)]
# Add lines between all points creating the rectangle
channel_lines = [model.add_line(boundary_points[i], boundary_points[i+1])
                 for i in range(-1, len(boundary_points)-1)]
airfoil_lines = [model.add_line(points[i], points[i+1])
                 for i in range(-1, len(points)-1)]
# Create a line loop and plane surface for meshing
channel_loop = model.add_curve_loop(channel_lines)
airfoil_loop = model.add_curve_loop(airfoil_lines)
plane_surface = model.add_plane_surface(channel_loop, holes=[airfoil_loop])
# Call gmsh kernal before adding physical entities
model.synchronize()
volume_marker = 6
model.add_physical([plane_surface], "Volume")
model.add_physical([channel_lines[0]], "Inflow")
model.add_physical([channel_lines[2]], "Outflow")
model.add_physical([channel_lines[1], channel_lines[3]], "Walls")
model.add_physical(airfoil_loop, "Obstacle")
# Generate geometry, and save as gmsh file
geometry.generate_mesh(dim=2)

gmsh.write("mesh.msh")
gmsh.clear()
geometry.__exit__()

# HELLO RAYAN

# resolution = 0.01
# # Chanel Parameters
# L = 2.2
# H = 0.41
# C = [0.2, 0.2, 0]
# r = 0.05
#
# # Initialize empty geometry using the built in kernel in GMSH
# geometry = pygmsh.geo.Geometry()
# # Fetch model we would like to add data to
# model = geometry.__enter__()
# # Add Circle
# circle = model.add_circle(C, r, mesh_size=resolution)
# # Add points with finer resolution on the left side
# points = [model.add_point((0, 0, 0), mesh_size=resolution),
#           model.add_point((L, 0, 0), mesh_size=5*resolution),
#           model.add_point((L, H, 0), mesh_size=5*resolution),
#           model.add_point((0, H, 0), mesh_size=resolution)]
# # Add lines between all points creating the rectangle
# channel_lines = [model.add_line(points[i], points[i+1])
#                  for i in range(-1, len(points)-1)]
# # Create a line loop and plane surface for meshing
# channel_loop = model.add_curve_loop(channel_lines)
# plane_surface = model.add_plane_surface(channel_loop, holes=[circle.curve_loop])
# # Call gmsh kernal before adding physical entities
# model.synchronize()
# volume_marker = 6
# model.add_physical([plane_surface], "Volume")
# model.add_physical([channel_lines[0]], "Inflow")
# model.add_physical([channel_lines[2]], "Outflow")
# model.add_physical([channel_lines[1], channel_lines[3]], "Walls")
# model.add_physical(circle.curve_loop.curves, "Obstacle")
# # Generate geometry, and save as gmsh file
# geometry.generate_mesh(dim=2)
#
# gmsh.write("mesh.msh")
# gmsh.clear()
# geometry.__exit__()
