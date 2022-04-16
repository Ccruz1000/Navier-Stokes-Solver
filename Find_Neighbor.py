# Imported Packages
import gmsh
import numpy as np
import sys

# User Defined Functions

# Define function to find neighbors
def find_neighbor(file):
    gmsh.initialize(sys.argv)

    gmsh.open(file)
    gmsh.model.getCurrent()

    # Obtain list of nodes
    nodes = gmsh.model.mesh.getNodes()
    node_index = nodes[0]   # Get list of node index
    node_coord = np.array_split(nodes[1], np.arange(3, len(nodes[1]), 3))   # Get nodes coordinates

    # Obtain list of elements
    element_list = []
    elements = gmsh.model.mesh.getElements()
    # Obtain lists of nodes

    print(elements)





    gmsh.finalize()
