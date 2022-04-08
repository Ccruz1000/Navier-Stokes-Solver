import gmsh
import sys
import numpy as np

gmsh.initialize(sys.argv)
gmsh.model.add("Test")


# Define coordinate file
coordinate_file = 'naca0012h_selig.txt'

# Import airfoil coordinates
airfoil_points = np.loadtxt(coordinate_file, skiprows=1, dtype='float')

chord = airfoil_points[0, 0]  # Airfoil Chord Length
width = 5 * chord
height = 2.5 * chord

# Define resolution
fine_res = 0.01 * chord
coarse_res = 0.3 * chord

# Define distance from trailing edge that closed point occurs
dist_TE = 1e-2

# Add origin
origin = gmsh.model.occ.add_point(0, 0, 0, fine_res)

points = []  # Initialize empty list to store airfoil points in
curve = []  # Initialize empty list to store airfoil curve in

# Import airfoil points
for i in range(len(airfoil_points)):
    points.append(gmsh.model.occ.addPoint(airfoil_points[i, 0], airfoil_points[i, 1], airfoil_points[i, 2], fine_res))

# Add points
boundary_points = [gmsh.model.occ.add_point(0, -height, 0, coarse_res),
                   gmsh.model.occ.add_point(0, height, 0, coarse_res),
                   gmsh.model.occ.add_point(width, height, 0, coarse_res),
                   gmsh.model.occ.add_point(width, -height, 0, coarse_res)]

# Create lines for channel
channel_lines = [gmsh.model.occ.add_line(boundary_points[i], boundary_points[i+1])
                 for i in range(-1, len(boundary_points)-1)]

# Add loop for channel
channel_loop = gmsh.model.occ.addCurveLoop(channel_lines)

# Add plane surface for channel_loop
channel_surface = gmsh.model.occ.addPlaneSurface([channel_loop])

# Add airfoil curve
curve.append(gmsh.model.occ.add_spline(points))
pt = gmsh.model.occ.add_point((dist_TE + airfoil_points[0, 0]), 0.0, 0.0, fine_res)
curve.append(gmsh.model.occ.add_line(points[-1], pt))
curve.append(gmsh.model.occ.add_line(pt, points[0]))

# Add line loop for airfoil
airfoil_loop = gmsh.model.occ.add_curve_loop(curve)

# Create plane surface for airfoil
s = gmsh.model.occ.addPlaneSurface([airfoil_loop])

# Remove airfoil surface from main surface
gmsh.model.occ.cut([(2, channel_surface)], [(2, s)])
gmsh.model.occ.synchronize()

