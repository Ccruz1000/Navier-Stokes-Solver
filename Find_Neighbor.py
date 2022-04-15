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
    # Find tetrahedrons, and face nodes
    print('--- Getting elements and edge nodes')
    quad, _ = gmsh.model.mesh.getElementsByType(3)  # Find quadrangles
    tri, _ = gmsh.model.mesh.getElementsByType(2)   # Find triangles
    line, _ = gmsh.model.mesh.getElementsByType(1)  # Find lines
    quad_nodes = gmsh.model.mesh.getElementEdgeNodes(3)  # Find quad edges
    tri_nodes = gmsh.model.mesh.getElementEdgeNodes(2)  # Find tri edges
    line_nodes = gmsh.model.mesh.getElementEdgeNodes(1)  # Find line nodes

    print('--- Compuing edge x element incidence')
    edges = []
    fxt = {}
    # for i in range(0, len(quad_nodes), 4):
    #     f = tuple(sorted(quad_nodes[i: i + 4]))
    #     edges.append(f)
    #     t = quad_nodes[i // 4]
    #     if not f in fxt:
    #         fxt[f].append(t)
    #     else:
    #         fxt[f].append(t)

    print(quad[2], 'Quad')
    # print(tri, 'Tri')
    # print(line, 'Line')
    print(quad_nodes[2], 'Quad Edge')
    # print(tri_nodes, 'Tri Edge')
    # print(line_nodes, 'Line Edge')




    gmsh.finalize()
