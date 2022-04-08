# Imported packages
import warnings
import gmsh
import numpy as np
# User defined functions

# Ignore warnings
warnings.filterwarnings("ignore")

# Initialize gmsh
gmsh.initialize()

# Create rectangular channel for airfoil
gmsh.model.add('NACA Channel')
L, W, r = 2.5, 0.41, 0.25  # Add rectangular channel size
channel = gmsh.model.occ.addRectangle(-1, -0.25, 0, L, W)  # Add rectangular channel
circle = gmsh.model.occ.add_circle(0, 0, 0, r)  # Add circle
fluid = gmsh.model.occ.cut([(2, channel)], [(2, circle)])  # Subtract circle from channel
gmsh.model.occ.synchronize()  # Synchronize CAD model
# Tag entities
volumes = gmsh.model.getEntities(dim=2)
assert(volumes == fluid[0])
fluid_marker = 11
gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], fluid_marker)
gmsh.model.setPhysicalName(volumes[0][0], fluid_marker, "Fluid Volume")
# Find all entities to add physical groups
surfaces = gmsh.model.occ.getEntities(dim=2)
inlet_marker, outlet_marker, wall_marker, obstacle_marker = 1, 3, 5, 7
walls = []
obstacles = []

