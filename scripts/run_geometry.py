#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os.path as path

import numpy as np
import matplotlib.pyplot as plt

from grains import HAS_OCCT
from grains.geometry import polygonize, overlay_regions, splinegonize, regions2step, plot_splinegons
if HAS_OCCT:
    from OCC.Display.SimpleGui import init_display
else:
    raise ImportError('PythonOCC is required for this script to run.')

# Use a sample image shipped with the code (https://stackoverflow.com/a/36476869/4892892)
script_dir = path.dirname(__file__)  # absolute directory the script is in
rel_image_path = 'data/1_labelimage.npy'
abs_image_path = path.join(script_dir, rel_image_path)
test_image = np.load(abs_image_path)

# Polygon representation of the label image
polygons = polygonize(test_image, connectivity=1, look_around=4, close=True)
# Plot the polygonized regions superposed on the original image
ax = overlay_regions(test_image, polygons)
plt.show()

# Splinegon representation of the label image
splinegons, _ = splinegonize(test_image, connectivity=1, look_around=4, degree_min=3,
                             degree_max=5, continuity='C2', tol=10)
# Plot the splinegonized regions
plot_splinegons(list(splinegons.values()))
# Write the geometry to a STEP file
data_dir = path.join(script_dir, 'data')
regions2step(list(splinegons.values()), path.join(data_dir, 'microstructure.stp'))
