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

"""

from math import sqrt

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


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
