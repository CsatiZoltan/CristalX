#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os.path as path

import numpy as np
import matplotlib.pyplot as plt

from grains import HAS_OCCT
from grains.geometry import polygonize, overlay_regions, fit_spline
if HAS_OCCT:
    from OCC.Display.SimpleGui import init_display
else:
    raise ImportError('PythonOCC is required for this script to run.')

# Create and plot a B-spline
random = np.random.random((10, 2))
spline = fit_spline(random, degree_max=11, degree_min=2, continuity='C0')
display, start_display, add_menu, add_function_to_menu = init_display()
display.DisplayShape(spline, update=True)

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
