#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 08:31:08 2020

@author: zoltan
"""

from math import sqrt

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def data_Pierre():
    """

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
