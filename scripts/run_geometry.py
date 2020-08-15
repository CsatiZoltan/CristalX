#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from os import path
from os.path import join

import numpy as np
import matplotlib.pyplot as plt

from grains import HAS_OCCT
from grains.geometry import polygonize, overlay_regions, splinegonize, regions2step, plot_splinegons
from grains.simulation import change_domain
if HAS_OCCT:
    from OCC.Display.SimpleGui import init_display
else:
    raise ImportError('PythonOCC is required for this script to run.')

# Directory where we read from and write to
script_dir = path.dirname(__file__)  # absolute path to the directory the script is in
data_dir = join(script_dir, 'data')

# Choose whether you want to perform the simulation on the original segmented image or on the
# image obtained by extending the original segmented image
use_extended = True  # change to `False` if desired

# Use a sample image shipped with the code (https://stackoverflow.com/a/36476869/4892892)
if use_extended:
    image_matrix = join(data_dir, '1_labelimage_extended.npy')
    if not path.isfile(image_matrix):  # extended image not yet created
        segmented_image_matrix = join(data_dir, '1_labelimage.npy')
        image_matrix = np.load(segmented_image_matrix)
        image_matrix = change_domain(image_matrix, 0.4395, 0, 0, 0, 500)
        image_matrix = change_domain(image_matrix, 0, 0.5568, 0, 0, 501)
        np.save(join(data_dir, '1_labelimage_extended.npy'), image_matrix)
        image_matrix = join(data_dir, '1_labelimage_extended.npy')
else:
    image_matrix = join(data_dir, '1_labelimage.npy')
test_image = np.load(image_matrix)

# Polygon representation of the label image
# polygons = polygonize(test_image, connectivity=1, look_around=2, close=True)
# # Plot the polygonized regions superposed on the original image
# ax = overlay_regions(test_image, {15: polygons[15]})
# plt.show()

# Splinegon representation of the label image
splinegons, _ = splinegonize(test_image, connectivity=1, look_around=2, degree_min=3,
                             degree_max=3, continuity='C0', tol=1)
# Plot the splinegonized regions
plot_splinegons(list(splinegons.values()), color=(0, 0, 1))
# Write the geometry to a STEP file
regions2step(list(splinegons.values()), path.join(data_dir, 'microstructure.stp'))
