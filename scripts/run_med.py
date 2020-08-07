#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import json
import os
from os.path import join
from math import sqrt, pi

from grains import HAS_MED
if HAS_MED:
    from grains.med import *
else:
    raise ImportError('MEDCoupling is required for this script to run.')


# Directory where we load from and write to
script_dir = os.path.dirname(__file__)  # absolute path to the directory the script is in
data_dir = join(script_dir, 'data')

# Use a sample mesh shipped with the code
mesh_file = join(data_dir, '1_mesh_extended.med')
msh = read_mesh(mesh_file)

# Obtain the mesh topology
nodes, node_groups = get_nodes(msh)
elems, elem_groups = get_elements(msh, 'group')

# In some cases, the node numbering is not counter-clockwise, as required by many finite element
# software.
elems = change_node_numbering(elems, nodes, orientation='ccw')

# Rotate the mesh so that it will be easier to deal with it later
nodes = rotate_mesh(nodes, -pi/2)

# Save the .inp file
write_inp(join(data_dir, '1_mesh_extended_med.inp'), nodes, elems, elem_groups, node_groups)

# Compute the diameter of each element group. There are multiple definitions for the diameter.
# Here, we use the equivalent diameter, which requires the area of the element group. The area
# is computed based on its mesh. Note that the diameters are given in the same unit as what is
# used in the mesh. For the example mesh above, the unit is pixel.
d = {}
for group_name, group_elems in elem_groups.items():
    group_area = element_group_area(group_name, elem_groups, elems, nodes)
    group_diameter = sqrt(4*group_area/pi)
    d[group_name] = group_diameter
# Write it to a text file so that we can use this information by another module that does not
# have access to MEDCoupling.
diameters_file = os.path.join(data_dir, '1_diameters_extended.json')
with open(diameters_file, 'w') as f:
    json.dump(d, f, indent=4)
