# -*- coding: utf-8 -*-
"""
Demonstrates how to read material data from and write material into an Abaqus input file.
It also shows an example to scale an existing geometry by modifying the coordinates of the
nodes in the mesh.
"""

import os.path as path

from grains.abaqus import Geometry, Material

# # Create material
# Material.remove('data/1_mesh.inp', 'data/material_removed.inp')
# # Read material from input file
# mat = Material()
# mat.read('data/1_mesh.inp')
# # Add new material to the existing database
# mat.add_material('MAT-1')
# mat.add_linearelastic('MAT-1', 1, 0)
# # Save the modifications to a new file
# mat.write()
# mat.show()
# print(mat)

# Fetch the geometry from an input file
geom = Geometry()
geom.read('data/1_elastoplastic.inp')
# Shrink the original geometry and save the modified file
geom.scale(1/2)
geom.write()
