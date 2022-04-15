from Numpy_Meshgrid import *
from Find_Neighbor import find_neighbor

# Define Initial Conditions
coordinate_file = 'naca0012h_selig.txt'
geometry = 'circle'
mesh_file = 'Structured Rectangular Mesh.msh'
# mesh_file = 'Structured Tri_Quad.msh'
# mesh_generate(geometry, coordinate_file)
find_neighbor(mesh_file)
