#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path as path

from grains.segmentation import Segmentation

# User-defined parameters
filter_window_size = 5
cluster_merging_threshold = 7

# Use a sample image shipped with the code (https://stackoverflow.com/a/36476869/4892892)
script_dir = path.dirname(__file__)  # absolute directory the script is in
rel_image_path = 'data/1_cropped.png'
abs_image_path = path.join(script_dir, rel_image_path)

# Perform all image processing steps the class offers
GS = Segmentation(abs_image_path, interactive_mode=False)
# Filtering helps avoiding over-segmentation
filtered = GS.filter_image(filter_window_size)
# Quick shift segmentation
segmented = GS.initial_segmentation(filtered)
# Reduce over-segmentation by merging clusters that are close enough
reduced = GS.merge_clusters(segmented, threshold=cluster_merging_threshold)
# The next lines are useful when one wants to edit the segmented regions manually. It is possible in
# ImagePy to modify the skeleton.
boundary = GS.find_grain_boundaries(reduced)
skeleton = GS.create_skeleton(boundary)
# We must reconstruct the segments from the skeleton
watershed = GS.watershed_segmentation(skeleton)
# Save the image as a png file
GS.save_image('1_labelimage.png', watershed, is_label_image=True)
# Save the final segmented image as a matrix for later processing
GS.save_array('1_labelimage', watershed)

