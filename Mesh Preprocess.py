# Open mesh file
lines_list = open('Structured Tri_Quad.msh').read().splitlines()
mesh_format = []  # Initialize array to store mesh format type
physical_name = []  # Initialize array to store physical names
nodes = []  # Initialize array to store nodes
element = []    # Initialize array to store elements
cntr = 0
for line in lines_list:
    # Extract Mesh Format
    if line == '$MeshFormat':
        mesh_format.append(lines_list[cntr + 1])
    # Extract Physical Names
    if line == '$PhysicalNames':
        name_number = int(lines_list[cntr + 1])
        for i in range(name_number + 1):
            physical_name.append(lines_list[i + cntr + 1])
    # Extract Nodes
    if line == '$Nodes':
        node_number = int(lines_list[cntr + 1])
        for i in range(node_number):
            nodes.append(lines_list[i + cntr + 2])
    # Extract Elements
    if line == '$Elements':
        element_number = int(lines_list[cntr + 1])
        for i in range(element_number):
            element.append(lines_list[i + cntr + 2])
    cntr += 1

print(lines_list)
print(mesh_format, 'Mesh')
print(name_number, 'Name number')
print(physical_name, 'Name')
print(node_number, 'Node Number')
print(nodes, 'Node')
print(element_number, 'Element Number')
print(element, 'Element')
