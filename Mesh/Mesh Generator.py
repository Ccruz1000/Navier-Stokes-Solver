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

# Create boundary layer by extrusion or as a mesh constraint?
by_extrusion = False

# Generate curved mesh?
order2 = True

# Load airfoil from txt selig file
airfoil = np.loadtxt("naca0012h_selig.txt", skiprows=1, dtype='float')

# Initialize empty geometry using the built in kernel in gmsh
geometry = pygmsh.geo.Geometry()
# Fetch model we would like to add data to
model = geometry.__enter__()
points = []  # Initialize empty list to store airfoil points in
curve = []  # Initialize empty list to store airfoil curve in

# Loop through airfoil points to generate points
for i in range(len(airfoil)):
    points.append(model.add_point((airfoil[i, 0], airfoil[i, 1], airfoil[i, 2]), mesh_size=resolution))

# Create curve from airfoil points
# curve.append(model.addSpline(points))

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
