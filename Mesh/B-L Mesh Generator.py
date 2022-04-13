# Imported functions
import gmsh
import numpy as np
import pygmsh

# User defined functions


# Import airfoil coordinates
coordinate_file = 'naca0012h_selig.txt'
airfoil_points = np.loadtxt(coordinate_file, skiprows=1, dtype='float')

# Define channel (will change to C mesh later
Start = -1  # Start of rectangular channel
End = 4  # End of rectangular channel
Top = 1  # Top of rectangular channel
Bottom = -1  # Bottom of rectangular channel

# Define mesh resolution
chord = 1
fine_res = chord*0.01
coarse_res = chord*0.15

# Define distance from trailing edge that closed point occurs
dist_TE = 1e-2

# Create boundary layer by extrusion or as a mesh constraint?
by_extrusion = False

# Generate curved mesh?
order2 = True

# Initialize empty geometry using the built in kernel in gmsh
geometry = pygmsh.geo.Geometry()

# Fetch model we would like to add data to
model = geometry.__enter__()
points = []  # Initialize empty list to store airfoil points in
curve = []  # Initialize empty list to store airfoil curve in

# Add airfoil points to model
for i in range(len(airfoil_points)):
    points.append(model.add_point((airfoil_points[i, 0], airfoil_points[i, 1], airfoil_points[i, 2]), fine_res))

# Create spline from airfoil points
curve.append(model.add_spline(points))
# Close trailing edge (sharp trailing edge, can add rounded later)
pt = model.add_point(((dist_TE + airfoil_points[0, 0]), 0.0, 0), fine_res)
curve.append(model.add_line(points[-1], pt))
curve.append(model.add_line(pt, points[0]))

# Add points for channel
boundary_points = [model.add_point((Start, Bottom, 0), mesh_size=2*fine_res),
                   model.add_point((End, Bottom, 0), mesh_size=coarse_res),
                   model.add_point((End, Top, 0), mesh_size=coarse_res),
                   model.add_point((Start, Top, 0), mesh_size=2*fine_res)]

# Add lines for channel
channel_lines = [model.add_line(boundary_points[i], boundary_points[i+1])
                 for i in range(-1, len(boundary_points)-1)]

# Add line loop for channel
channel_loop = model.add_curve_loop(channel_lines)

# Add line loop for airfoil
airfoil_loop = model.add_curve_loop(curve)

# Add plane surface for meshing
plane_surface = model.add_plane_surface(channel_loop, holes=[airfoil_loop])
recombined_surface = model.set_recombined_surfaces(plane_surface)

# Call gmsh kernal before adding physical entities
model.synchronize()
volume_marker = 6
model.add_physical([plane_surface], "Volume")
model.add_physical([channel_lines[0]], "Inflow")
model.add_physical([channel_lines[2]], "Outflow")
model.add_physical([channel_lines[1], channel_lines[3]], "Wall")
model.add_physical(airfoil_loop, "Airfoil")
# Generate geometry, and save as gmsh file
geometry.generate_mesh(dim=2)

gmsh.write("mesh.msh")
gmsh.clear()
geometry.__exit__()
