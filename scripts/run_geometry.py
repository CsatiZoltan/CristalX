#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import json
import os
from os.path import join
from math import pi

import numpy as np
import matplotlib.pyplot as plt

from grains.geometry import TriMesh


# Directory where we load from and write to
script_dir = os.path.dirname(__file__)  # absolute path to the directory the script is in
data_dir = join(script_dir, 'data')

# Load the previously saved mesh data
mesh_file = join(data_dir, '1_mesh_extended.npz')
with np.load(mesh_file, allow_pickle=True) as mesh:
    if np.alltrue([variable_name in {'nodes', 'elements', 'element_groups', 'node_groups'} for
                   variable_name in mesh]):
        nodes = mesh['nodes']
        elements = mesh['elements']
        element_groups = mesh['element_groups']
        node_groups = mesh['node_groups']
    else:
        raise Exception('Mesh must contain variables "nodes", "elements", "element_groups" and '
                        '"node_groups".')
# Retrieve groups, which were stored in dictionaries (https://stackoverflow.com/a/40220343/4892892)
element_groups = element_groups.item()
node_groups = node_groups.item()

# Create a Mesh object so that we perform manipulations on it
mesh = TriMesh(nodes, elements)
# Associate the element and node groups to it
for group_name, elements in element_groups.items():
    mesh.create_cell_set(group_name, elements)
for group_name, nodes in node_groups.items():
    mesh.create_vertex_set(group_name, nodes)

# In some cases, the node numbering is not counter-clockwise, as required by many finite element
# software.
mesh.change_vertex_numbering('ccw', inplace=True)

# For the simulation, we need to work on the actual domain. In this example, 1 mm corresponds to
# 29.55 pixels.
mesh.scale(1/29.55, inplace=True)

# Rotate the mesh so that it will be easier to deal with it later
mesh.rotate(-pi/2, inplace=True)

# Plot the mesh to see if that is what we wanted
mesh.plot(vertex_legends=True)
plt.show()

# Save the modified mesh to a different file
np.savez_compressed(join(data_dir, '1_mesh_extended_scaled.npz'), nodes=mesh.vertices,
                    elements=mesh.cells, element_groups=mesh.cell_sets,
                    node_groups=mesh.vertex_sets)
