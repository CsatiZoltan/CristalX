"""
A trimmed version of the Gala project (https://github.com/janelia-flyem/gala)
with some additions (new function, added documentation for the existing ones).
Gala is licensed by the Janelia Farm License:
http://janelia-flyem.github.io/janelia_farm_license.html

Copyright 2012 HHMI. All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.  Redistributions in binary form must
reproduce the above copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided with the distribution.

Neither the name of HHMI nor the names of its contributors may be used to endorse or
promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


from numpy import unique, ones, minimum
from scipy.ndimage import grey_dilation, generate_binary_structure, minimum_filter


def imextendedmin(image, h, connectivity=1):
    """Extended-minima transform.
    The extended minima transform is the regional minima of the h-minima transform.
    The implementation follows the MATLAB function under the same name.

    Parameters
    ----------
    image : ndarray
        The input array on which to perform imextendedmin.
    h : float
        Any local minima shallower than this will be flattened.
    connectivity : int, optional
        Determines which elements are considered as neighbors of the central
        element. Elements up to a squared distance of `connectivity` from
        the center are considered neighbors. If connectivity=1, no diagonal
        elements are neighbors.

    Returns
    -------
    bool ndarray
        True at places of the extended minima.
    """
    return regional_minima(imhmin(image, h), connectivity)


def hminima(a, thresh):
    """Suppress all minima that are shallower than thresh.

    Parameters
    ----------
    a : array
        The input array on which to perform hminima.
    thresh : float
        Any local minima shallower than this will be flattened.

    Returns
    -------
    out : array
        A copy of the input array with shallow minima suppressed.
    """
    maxval = a.max()
    ainv = maxval-a
    return maxval - morphological_reconstruction(ainv-thresh, ainv)


imhmin = hminima


def morphological_reconstruction(marker, mask, connectivity=1):
    """Perform morphological reconstruction of the marker into the mask.

    See the Matlab image processing toolbox documentation for details:
    http://www.mathworks.com/help/toolbox/images/f18-16264.html
    """
    sel = generate_binary_structure(marker.ndim, connectivity)
    diff = True
    while diff:
        markernew = grey_dilation(marker, footprint=sel)
        markernew = minimum(markernew, mask)
        diff = (markernew-marker).max() > 0
        marker = markernew
    return marker


def regional_minima(a, connectivity=1):
    """Find the regional minima in an ndarray.
    As written in the MATLAB documentation of the imregionalmin function:
    "Regional minima are connected components of pixels with a constant
    intensity value, surrounded by pixels with a higher value."
    """
    values = unique(a)
    delta = (values - minimum_filter(values, footprint=ones(3)))[1:].min()
    marker = complement(a)
    mask = marker+delta
    return marker == morphological_reconstruction(marker, mask, connectivity)


def complement(a):
    return a.max()-a
