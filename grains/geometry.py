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


def polygon_orientation(polygon):
    """Determines whether a polygon is oriented clockwise or counterclockwise.

    Parameters
    ----------
    polygon : list
        Each element of the list denotes a vertex of the polygon and in turn is another list of two
        elements: the x and y coordinates of a vertex.

    Returns
    -------
    orientation : {'cw', 'ccw'}
        'cw': clockwise, 'ccw': counterclockwise orientation

    Notes
    -----
    The formula to determine the orientation is from https://stackoverflow.com/a/1165943/4892892.
    For simple polygons (polygons that admit a well-defined interior), a faster algorithm exits, see
    https://en.wikipedia.org/wiki/Curve_orientation#Orientation_of_a_simple_polygon.

    Examples
    --------
    >>> polygon = [[5, 0], [6, 4], [4, 5], [1, 5], [1, 0]]
    >>> polygon_orientation(polygon)
    'ccw'

    """
    n_vertex = len(polygon)
    edge_sum = 0
    for idx, vertex in enumerate(polygon):
        next_vertex = polygon[(idx + 1) % n_vertex]  # allow indexing past the last vertex
        edge_sum += (next_vertex[0] - vertex[0]) * (next_vertex[1] + vertex[1])
    if edge_sum > 0:
        orientation = 'cw'
    else:
        orientation = 'ccw'
    return orientation
