#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from os import path
from os.path import join

import numpy as np
from scipy.stats import kurtosis
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.morphology import dilation

from grains.simulation import nature_of_deformation
from grains.dic import DIC
from grains import HAS_TABLES


# Directory where we read from and write to
script_dir = os.path.dirname(__file__)  # absolute path to the directory the script is in
data_dir = join(script_dir, 'data')

################################################################################
#  Determine the ratio of intergranular and intragranular strain localization  #
################################################################################

# Crop data which belongs to the testing device (the indices were obtained using:
# `selected = plt.ginput(-1, timeout=-1))`
cropped_region = np.s_[800:2175, 860:3000]

# 1) Load the segmented image, representing the grains of the microstructure
microstructure = imread(join(data_dir, '2_labelimage.tiff'))
microstructure = microstructure[cropped_region]

# 2) Bands in the neighborhood of the interfaces
band_width = 3
microstructure = dilation(microstructure)

# 3) Process a series of full-field measurements (DIC)
dic_data = join(data_dir, 'FULLTEST_fields.hdf')
MEASUREMENT_COUNT = 382
PERCENTILE_COUNT = 10

# 4) User settings
show_plots = True
save_plots = False
time_instances = range(MEASUREMENT_COUNT)

# Analysis
matrix = np.zeros((PERCENTILE_COUNT+1, MEASUREMENT_COUNT))
kurt_band = []
kurt_bulk = []
visualize = show_plots or save_plots
if not path.isfile(dic_data):
    raise Exception('Measurement dataset {0} not found.'.format(dic_data))
if not HAS_TABLES:
    raise ImportError('The PyTables package is needed to load the dataset.')
import tables
measurement = tables.open_file(dic_data).root.res
for time_instance in time_instances:
    print('Time instance:', time_instance)
    # a) Load the displacement field
    displacement_field = measurement[time_instance, cropped_region[0], cropped_region[1], :]
    u = displacement_field[:, :, 0]
    v = displacement_field[:, :, 1]
    # b) Obtain the equivalent strain field
    dic = DIC(u, v)
    strain_tensor = dic.strain('Green-Lagrange')
    vonMises_strain = DIC.equivalent_strain(strain_tensor, 'von Mises')
    # c) Characterize the deformation in the current time step
    boundary_strain, bulk_strain, bands = nature_of_deformation(microstructure,
                                                                vonMises_strain,
                                                                band_width, visualize=visualize)
    if save_plots:
        plt.savefig(join(data_dir, 't_' + str(time_instance) + '.png'), dpi=500)
    if show_plots:
        plt.show(block=True)
    normalized_boundary_strain = boundary_strain[bands] / max(boundary_strain[bands])
    normalized_bulk_strain = bulk_strain[~bands] / max(bulk_strain[~bands])
    quantiles = np.percentile(normalized_boundary_strain, range(0, 101, PERCENTILE_COUNT))
    matrix[:, time_instance] = quantiles
    kurt_band.append(kurtosis(normalized_boundary_strain))
    kurt_bulk.append(kurtosis(normalized_bulk_strain))
