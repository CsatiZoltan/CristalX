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

# Save the mesh both as a numpy array and as an Abaqus .inp file
np.savez_compressed(join(data_dir, '1_mesh_extended.npz'), nodes=nodes, elements=elems,
                    element_groups=elem_groups, node_groups=node_groups)
