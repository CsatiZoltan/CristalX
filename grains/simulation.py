#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions
---------
.. autosummary::
    :nosignatures:
    :toctree: functions/

    data_Pierre
    hallpetch_constants
    hallpetch
    hallpetch_plot
    change_domain
    nature_of_deformation

"""
from math import sqrt

import numpy as np
from scipy.stats import kurtosis
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.table import table
import matplotlib.gridspec as gridspec

from grains.analysis import label_image_skeleton, thicken_skeleton


def data_Pierre():
    """Yield stresses and average grain diameters from Pierre's thesis.

    Returns
    -------
    sigma_y : list of floats
        Yield stresses.
    d : list of floats
        Diameters of the grains.

    """
    sigma_y = [220, 75, 100]
    d = [0.03, 3, 0.5]
    return sigma_y, d


def hallpetch_constants(sigma_y, d):
    """Determines the two Hall-Petch constants.
    Given available measurements for the grains sizes and the yield stresses,
    the two constants in the Hall-Petch formula are computed.

    Parameters
    ----------
    sigma_y : list of floats
        Yield stresses.
    d : list of floats
        Diameters of the grains.

    Returns
    -------
    sigma_0 : float
        Starting stress for dislocation movement (material constant).
    k : float
        Strengthening coefficient (material constant).

    Notes
    -----
    If two data points are given in the inputs (corresponding to two measurements),
    the output parameters have unique values:
        k = (sigma_y[0] - sigma_y[1]) / (1/sqrt(d[0]) - 1/sqrt(d[1]))
        sigma_0 = sigma_y[0] - k/sqrt(d[0])
    If there are more than two measurements, the resulting linear system is
    overdetermined.
    In both cases, the outputs are determined using least squares fitting.

    """

    if not isinstance(sigma_y, list) or not isinstance(d, list):
        raise Exception("Inputs must be lists.")
    if len(sigma_y) != len(d):
        raise Exception("Inputs must have the same number of elements.")
    if len(d) < 2:
        raise Exception("At least two data points required for both inputs.")
    # Solve the (possibly overdetermined) linear system
    A = np.hstack([np.vstack(np.ones_like(d)), np.vstack(1/np.sqrt(d))])
    b = np.vstack(sigma_y)
    x = np.linalg.lstsq(A, b, rcond=None)[0]
    # Unpack numpy column vector
    sigma_0 = x[0][0]
    k = x[1][0]
    return sigma_0, k


def hallpetch(sigma_0, k, d):
    """Computes the yield stress from the Hall-Petch relation.

    Parameters
    ----------
    sigma_0 : float
        Starting stress for dislocation movement (material constant).
    k : float
        Strengthening coefficient (material constant).
    d : float or list of floats
        Diameter of the grain.

    Returns
    -------
    sigma_y : float or list of floats
        Yield stress.

    """
    if isinstance(d, list):
        sigma_y = [sigma_0 + k/sqrt(d_i) for d_i in d]
    else:
        sigma_y = sigma_0 + k/sqrt(d)
    return sigma_y


def hallpetch_plot(sigma_y, d, units=("MPa", "mm")):
    """Plots the Hall-Petch formula for given grain sizes and yield stresses.

    Parameters
    ----------
    sigma_y : ndarray
        Yield stresses.
    d : ndarray
        Grain diameters.
    units : 2-tuple of str, optional
        Units for the yield stress and the grain diameters.
        The default is ("MPa","mm").

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object is returned in case further manipulations are necessary.

    """

    fig, ax = plt.subplots()
    # Hall-Petch relation
    inv_sqrt_d = 1/np.sqrt(d)
    hp = ax.plot(inv_sqrt_d, sigma_y)
    # Yield stress for monocrystals
    sigma_0 = hallpetch_constants(sigma_y[0:2], d[0:2])[0]
    # Region of the plot that represents monocrystals
    d_max = np.amax(inv_sqrt_d)
    sigma_y_max = np.amax(sigma_y)
    monocrystals = Rectangle((0, sigma_0 - 0.25*(sigma_y_max - sigma_0)), d_max,
                             0.25*(sigma_y_max - sigma_0), angle=0,
                             facecolor="blue", alpha=0.5)
    ax.add_patch(monocrystals)
    ax.text(1/2*d_max, 1/8*(9*sigma_0 - sigma_y_max), "monocrystals",
            horizontalalignment="center", verticalalignment="center")
    # Axis title, limits and labels
    ax.set_xlim(0, d_max)
    ax.set_ylim(sigma_0 - 1/4*(sigma_y_max - sigma_0), sigma_y_max)
    xlabel = "$d^{{-\\frac{{1}}{{2}}}}$ " \
             "$\\left[\\mathrm{{{0}}}^{{-\\frac{{1}}{{2}}}}\\right]$".format(units[1])
    ylabel = "$\\sigma_y$ [{0}]".format(units[0])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Effect of the grain diameter on the yield strength")
    ax.legend(hp, ["Hall-Petch"])
    return fig


def change_domain(image, left, right, bottom, top, padding_value=np.nan):
    """Extends or crops the domain an image fills.

    The original image is extended, cropped or left unchanged in the left, right, bottom and top
    directions by padding the corresponding 2D or 3D array with a given value or removing
    existing elements. Non-integer image length is avoided by rounding up. This choice prevents
    ending up with an empty image during cropping or no added pixel during extension.

    Parameters
    ----------
    image : ndarray
        2D or 3D numpy array representing the original image.
    left, right : float
        When positive, the domain is extended, when negative, the domain is cropped by the given
        value relative to the image width.
    bottom, top : float
        When positive, the domain is extended, when negative, the domain is cropped by the given
        value relative to the image height.
    padding_value : float, optional
        Value for the added pixels. The default is numpy.nan.

    Returns
    -------
    changed_image : ndarray
        2D or 3D numpy array representing the modified domain.

    Examples
    --------
    Crop an image at the top, extended it at the bottom and on the left,
    and leave it unchanged on the right. Note the rounding for non-integer pixels.

    >>> import numpy as np
    >>> image = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
    >>> modified = change_domain(image, 0.5, 0, 1/3, -1/3, 0)
    >>> image  # no in-place modification, the original image is not overwritten
    array([[0.1, 0.2, 0.3],
           [0.4, 0.5, 0.6],
           [0.7, 0.8, 0.9]])
    >>> modified
    array([[0. , 0. , 0.4, 0.5, 0.6],
           [0. , 0. , 0.7, 0.8, 0.9],
           [0. , 0. , 0. , 0. , 0. ]])

    """
    # TODO: Negative resulting image size must not be allowed. It means that left + right > width
    #  and top + bottom>height must be fulfilled.
    height, width = np.shape(image)
    # number of pixels (entries in the matrix) to be added
    left = np.ceil(width*left).astype(np.int)
    right = np.ceil(width*right).astype(np.int)
    bottom = np.ceil(height*bottom).astype(np.int)
    top = np.ceil(height*top).astype(np.int)
    # Shrink the domain (crop the matrix)
    changed_image = image[abs(min(top, 0)):height - abs(min(bottom, 0)),
                          abs(min(left, 0)):width - abs(min(right, 0))]
    # Extend the domain (pad the matrix)
    pad_width = ((max(top, 0), max(bottom, 0)), (max(left, 0), max(right, 0)))
    changed_image = np.pad(changed_image, pad_width, constant_values=padding_value)
    return changed_image


def nature_of_deformation(microstructure, strain_field, interface_width=3, visualize=True):
    """Characterizes the intergranular/intragranular deformations.

    To decide whether the strain localization is intergranular (happens along grain boundaries,
    also called interfaces) or intragranular in a given microstructure, the strain field is
    projected on the microstructure. Here, by strain field we mean a scalar field, often called
    `equivalent strain` that is derived from a strain tensor.

    It is irrelevant for this function whether the strain field is obtained from a numerical
    simulation or from a (post-processed) full-field measurement. All what matters is that the
    strain field be available on a grid of the same size as the microstructure.

    The strain field is assumed to be localized on an interface if its neighborhood, with band
    width defined by the user, contains higher strain values than what is outside the band (i.e.
    the grain interiors). A too large band width identifies small grains to have boundary only,
    without any interior. This means that even if the strain field in reality localizes inside
    such small grains, the localization is classified as intergranular. However, even for
    extreme deformations, one should not expect that the strain localizes on an interface with a
    single-point width. Moreover, using a too small band width is susceptible to the exact position
    of the interfaces, which are extracted from the grain microstructure. A judicial balance needs
    to be achieved in practice.

    Parameters
    ----------
    microstructure : ndarray
        Labelled image corresponding to the segmented grains in a grain microstructure.
    strain_field : ndarray
        Discrete scalar field of the same size as the :code:`microstructure`.
    interface_width : int, optional
        Thickness of the band around the interfaces.
    visualize : bool, optional
        If True, three plots are created. Two of them show the deformation field within the bands
        and outside the bands. They are linked together, so when you pan or zoom on one,
        the other plot will follow. The third plot contains two histograms on top of each other,
        giving the frequency of the strain values within the bands and outside the bands.

    Returns
    -------
    boundary_strain : ndarray
        Copy of :code:`strain_field`, but values outside the band are set to NaN.
    bulk_strain : ndarray
        Copy of :code:`strain_field`, but values within the band are set to NaN.
    bands : ndarray
        Boolean array of the same size as :code:`strain_field`, with True values
        corresponding to the band.

    See Also
    --------
    :meth:`grains.dic.DIC.strain` : Computes a strain tensor from the displacement field.
    :meth:`grains.dic.DIC.equivalent_strain` : Extracts a scalar quantity from a strain tensor.
    :func:`matplotlib.pyplot.hist` : Plots a histogram.

    Notes
    -----
    1. From the modelling viewpoint, it is important to know whether the strain localizes to the
       grain boundaries or it is dominant within the grains as well. In the former case,
       simplifications in the models save computational time in the simulations.
    2. In dynamics, the evolution of the strain field is relevant. E.g. an initially
       intergranular deformation can turn into diffuse localization that occurs within the grains
       as well. In that case, a strain field must be obtained at each time step, and this
       function can be called for each such instance.

    Examples
    --------
    The following figure was created by this function with :code:`visualize` set to `True` and
    :code:`band_width` chosen to be 3.

    .. image:: ../images/intergranular-intragranular_deformations.*

    """
    # Bands around the interfaces
    bands = thicken_skeleton(label_image_skeleton(microstructure), interface_width)

    # Strain values in the bands and outside the bands
    boundary_strain = strain_field.copy()
    bulk_strain = strain_field.copy()
    boundary_strain[~bands] = np.nan
    bulk_strain[bands] = np.nan

    if visualize:
        # Set up figure for showing the results
        fig = plt.figure(tight_layout=True)
        gs = gridspec.GridSpec(2, 2)
        # Strain field localized onto the bands and onto the grain interiors
        ax_1 = fig.add_subplot(gs[0, 0])
        ax_1.imshow(boundary_strain)
        ax_1.set_title('Strain field in the bands around the grain boundaries')
        ax_2 = fig.add_subplot(gs[1, 0], sharex=ax_1, sharey=ax_1)
        ax_2.imshow(bulk_strain)
        ax_2.set_title('Strain field in the interior of the grains')
        # Histogram to show the strain value distributions
        ax_3 = fig.add_subplot(gs[0, 1])
        ax_3.set_title('Normalized strain distribution')
        normalized_boundary_strain = boundary_strain[bands] / max(boundary_strain[bands])
        normalized_bulk_strain = bulk_strain[~bands] / max(bulk_strain[~bands])
        _, _, h_boundary = ax_3.hist(boundary_strain[bands], bins='auto', alpha=0.5,  color='blue')
        _, _, h_bulk = ax_3.hist(bulk_strain[~bands], bins='auto', alpha=0.5, color='red')
        plt.legend([h_boundary, h_bulk], ['band', 'interior'])

        # Relevant statistical measures
        rows = ('mean', 'kurtosis')
        columns = ('band', 'interior')
        data = [[np.mean(normalized_boundary_strain), np.mean(normalized_bulk_strain)],
                [kurtosis(normalized_boundary_strain), kurtosis(normalized_bulk_strain)]]
        ax_4 = fig.add_subplot(gs[1, 1])
        ax_4.set_axis_off()
        data_table = table(ax_4, cellText=np.array(data).astype(str), rowLabels=rows,
                           colLabels=columns, loc='center')
        # data_table.auto_set_font_size(False)  # https://stackoverflow.com/a/15514091/4892892
        # data_table.set_fontsize(10)
    return boundary_strain, bulk_strain, bands
