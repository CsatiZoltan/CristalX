#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os.path as path

import numpy as np
from skimage.measure import regionprops

from grains.analysis import Analysis, plot_prop

# Use a sample image shipped with the code (https://stackoverflow.com/a/36476869/4892892)
script_dir = path.dirname(__file__)  # absolute directory the script is in
rel_image_path = 'data/1_labelimage.npy'
abs_image_path = path.join(script_dir, rel_image_path)

# Perform all image processing steps the class offers
scale = 1/29.55
label_image = np.load(abs_image_path)
GA = Analysis(label_image, interactive_mode=False)
GA.set_scale(scale)
GA.compute_properties()
GA.show_properties(gui=True)
GA.show_grains(grain_property='label')

properties = regionprops(GA.label_image, cache=False)
plot_prop(properties[187], scale)
