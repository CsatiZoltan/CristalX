# -*- coding: utf-8 -*-
"""
This module contains the Analysis class, responsible for the analysis of
segmented grain-based microstructures.
"""

from math import sqrt

import numpy as np
from scipy.spatial.distance import pdist
from skimage import io, color
from skimage.measure._find_contours import find_contours
from skimage.measure._marching_cubes_lewiner import marching_cubes_lewiner as marching_cubes
from skimage.measure import regionprops
import matplotlib.pyplot as plt
import pandas as pd


class Analysis:
    """Analysis of grain assemblies.

    Attributes
    ----------
    original_image : ndarray
        Matrix representing the initial, unprocessed image.
    save_location : str
        Directory where the processed images are saved
    """

    def __init__(self, label_image, interactive_mode=True):
        """Initialize the class with file paths and with some options

        Parameters
        ----------
        label_image : ndarray
            Segmented image.
        interactive_mode : bool, optional
            When True, images of each image manipulation step are plotted and
            details are shown in the console.
            Default is False.

        Returns
        -------
        None.
        """

        # Check inputs
        self.label_image = label_image
        self.__interactive_mode = interactive_mode
        self.scale = 1
        # Test whether it is a valid label image
        # Hack, delete it later (only needed for OOF)
        self.label_image = self.label_image + 1
        if not np.issubdtype(self.label_image.dtype, np.integer):
            raise Exception('Integer array is expected, received {0}.'.format(self.label_image.dtype))
        if not np.all(self.label_image > 0):
            raise Exception('Label image must contain positive integer entries')
        if np.ndim(self.label_image) != 2:
            raise Exception('Label image must be a matrix, received an array of dimension {0}.'.format(np.ndim(self.label_image)))
        # Optionally show the loaded image
        if self.__interactive_mode:
            io.imshow(color.label2rgb(self.label_image))
            io.show()
            print('Image successfully loaded.')
        # Initialize a dictionary holding the grain properties
        self.properties = {key: None for key in ['label', 'area', 'centroid',
                                                 'coordinate', 'diameter']}

    def set_scale(self, pixel_per_unit=1):
        """Defines a scale for performing computations in that unit.
        Image measures (area, diameter, etc.) are performed on a matrix
        corresponding to a label image. Therefore, the result of all the
        computations are obtained in pixel units. It is often of interest to
        access the results in physical units (mm, cm, inch, etc.). Manually
        converting pixels, pixel squares, etc. to pyhsical units, physical unit
        sqaures, etc. are tedious and error prone. Once the conversion between
        a pixel and a physical unit is given, all the subsequent calculations
        are performed in the desired physical unit.

        Parameters
        ----------
        pixel_per_unit : float or int or scalar ndarray, optional
            Number of pixels contained in a certain unit. The default is 1, in
            which case all measurements are performed in pixel units.

        Returns
        -------
        None.

        """
        correct_type = isinstance(pixel_per_unit, (int, float, np.ndarray))
        scalar = np.isscalar(pixel_per_unit)
        if not correct_type or not scalar:
            raise Exception('A scalar number expected.')
        self.scale = pixel_per_unit

    def compute_properties(self):
        """Determines relevant properties of the grains.
        The area of each grain is determined in the units previously given in
        the `set_scale` method.

        Parameters
        ----------
        window_size : int
            Size of the sampling window.
        image_matrix : 3D ndarray with size 3 in the third dimension, optional
            Input image to be filtered. If not given, the original image is used.

        Returns
        -------
        filtered_image : 3D ndarray with size 3 in the third dimension
            Filtered image, output of the median filter algorithm.
        """

        # Do not compute all the region properties (https://stackoverflow.com/a/36323763/4892892)
        properties = regionprops(self.label_image, cache=False)
        self.properties = dict(label=[prop.label for prop in properties],
                               area=[self.scale**2*prop.area for prop in properties],
                               centroid=[prop.centroid for prop in properties],
                               coordinate=[prop.coords for prop in properties],
                               equivalent_diameter=[self.scale*prop.equivalent_diameter for prop in properties],
                               feret_diameter=[self.scale*feret_diameter(prop) for prop in properties])
        if self.__interactive_mode:
            print('Grain properties determined.')

    def show_properties(self, gui=False):
        """Displays previously computed properties of the grains

        Parameters
        ----------
        gui : bool, optional
            If true, the grain properties are shown in a GUI. If false, they
            are printed to stdout. The default is False.
            The GUI requires the `dfgui` modul, which can be obtained from
            https://github.com/bluenote10/PandasDataFrameGUI

        Returns
        -------
        None.

        """
        table = pd.DataFrame(self.properties)
        # Print the table or show it in a GUI
        if gui:
            # PandasDataFrameGUI uses the wxwidgets backend of matplotlib
            import matplotlib
            matplotlib.use('WXAgg', force=True)
            import dfgui
            dfgui.show(table)
        else:
            pd.options.display.max_rows = None  # show all grains
            print(table)

    def show_grains(self, grain_property=None):
        """Display the grains, optionally with a property superposed.

        Parameters
        ----------
        grain_property : {None, 'area', 'centroid', 'coordinate',
                          'equivalent_diameter', 'feret_diameter', 'label'}
                         optional
            If not None, the selected property is shown on the grain as text.

        Returns
        -------
        None.

        """
        nlabel = len(np.unique(self.label_image))
        io.imshow(color.label2rgb(self.label_image, colors=np.random.random((nlabel, 3))))
        # Display information on the grains, if requested
        if grain_property is not None:
            grain_property = self.properties[grain_property]
            grain_centroid = self.properties["centroid"]
            for cent, prop in zip(grain_centroid, grain_property):
                cent = (cent[1], cent[0])  # image vs matrix indexing!
                plt.annotate('%.1f'%prop, cent, ha='center', va='center')
        plt.show()


def feret_diameter(prop):
    """Determines the maximum Feret diameter.

    Parameters
    ----------
    prop : RegionProperties
        Describes a labeled region.

    Returns
    -------
    max_feret_diameter : float
        Maximum Feret diameter of the region.

    See also
    --------
    skimage.measure.regionprops : Measure properties of labeled image regions

    Examples
    --------
    >>> import numpy as np
    >>> from skimage.measure import regionprops
    >>> image = np.ones((2,2), dtype=np.int8)
    >>> prop = regionprops(image)[0]
    >>> feret_diameter(prop)
    2.23606797749979

    """
    identity_convex_hull = np.pad(prop.convex_image, 2, mode='constant',
                                  constant_values=0)
    if prop._ndim == 2:
        coordinates = np.vstack(find_contours(identity_convex_hull, 0.5,
                                              fully_connected='high'))
    elif prop._ndim == 3:
        coordinates, _, _, _ = marching_cubes(identity_convex_hull, level=0.5)
    distances = pdist(coordinates, 'sqeuclidean')
    # TODO: allow computing the Feret diameter along a given direction (t,s)
    # For this, restrict the maximum distance to those pairwise distances whose
    # orientation (up to a given tolerance) is (t,s)
    max_feret_diameter = sqrt(np.max(distances))
    return max_feret_diameter


def plot_prop(prop, pixel_per_unit=1, show_axis=True):
    """Plots relevant region properties into a single figure.
    Four subfigures are created, giving the region's
        - image, its area and its center
        - filled image, its area
        - bounding box, its area
        - convex image, its area

    Parameters
    ----------
    prop : RegionProperties
        Describes a labeled region.
    pixel_per_unit : float or int, optional
            Number of pixels contained in a certain unit. The default is 1, in
            which case all measurements are performed in pixel units.

    Returns
    ------
    fig : matplotlib.figure.Figure
        The figure object is returned in case further manipulations are necessary.

    """
    scale = pixel_per_unit
    fig, ax = plt.subplots(2, 2, constrained_layout=True)
    # Use actual coordinates as appears in the label image
    extent = [prop.bbox[1], prop.bbox[3], prop.bbox[2], prop.bbox[0]]
    # Image, its area and its geometrical center
    ax[0, 0].imshow(prop.image, origin='upper', extent=extent)
    ax[0, 0].text(prop.centroid[1], prop.centroid[0], '+',
                  horizontalalignment='center', verticalalignment='center',
                  fontsize='x-large', fontweight='bold')
    ax[0, 0].set_title('Image, area: {0:.2f}'.format(scale**2*prop.area))
    # Image with holes filled, and its area
    ax[0, 1].imshow(prop.filled_image, origin='upper', extent=extent)
    ax[0, 1].set_title('Filled image, area: {0:.2f}'.format(scale**2*prop.filled_area))
    # Convex image, and its area
    ax[1, 0].imshow(prop.convex_image, origin='upper', extent=extent)
    ax[1, 0].set_title('Convex image, area: {0:.2f}'.format(scale**2*prop.convex_area))
    # Switch off axes if requested
    if not show_axis:
        for axis in ax.flat:
            axis.set_axis_off()
    plt.show()
    return fig
