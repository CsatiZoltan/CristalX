# -*- coding: utf-8 -*-

"""
This module contains functions to provide geometrical descriptions of a label image.
Each labelled region is turned to a surface. This allows mesh generators to work on
the geometry instead of an image, thereby creating good quality meshes.
"""

import numpy as np
from skimage.segmentation import find_boundaries
from skimage.morphology import skeletonize
from skan import Skeleton


def build_skeleton(label_image, connectivity=1, detect_boundaries=True):
    """Builds skeleton connectivity of a label image.

    A single-pixel wide network is created, separating the labelled image regions. The resulting
    network contains information about how the regions are connected.

    Parameters
    ----------
    label_image : 2D ndarray with signed integer entries
        Label image, representing a segmented image.
    connectivity : {1,2}, optional
        A connectivity of 1 (default) means pixels sharing an edge will be considered neighbors.
        A connectivity of 2 means pixels sharing a corner will be considered neighbors.
    detect_boundaries : bool, optional
        When True, the image boundaries will be treated as part of the skeleton. This allows
        identifying boundary regions in the `skeleton2regions` function. The default is True.

    Returns
    -------
    skeleton_network : Skeleton
        Geometrical and topological information about the skeleton network of the input image.

    See Also
    --------
    skan.Skeleton

    """
    # 2D image, given as a numpy array is expected
    if type(label_image) != np.ndarray:
        raise Exception('Label image must be a numpy array (ndarray).')
    image_size = np.shape(label_image)
    if len(image_size) != 2:
        raise Exception('A 2D array is expected.')
    if not issubclass(label_image.dtype.type, np.signedinteger):
        raise Exception('Matrix entries must be positive integers.')

    # Surround the image with an outer region
    if detect_boundaries:
        label_image = np.pad(label_image, pad_width=1, mode='constant', constant_values=-1)

    # Find the boundaries of the label image and then extract its skeleton
    boundaries = find_boundaries(label_image, connectivity=connectivity)
    skeleton = skeletonize(boundaries)

    # Build the skeleton network using `skan`
    skeleton_network = Skeleton(skeleton, source_image=label_image, keep_images=True)
    return skeleton_network
