# -*- coding: utf-8 -*-

"""
This module provides utility functions that are useful for the **grains** package.
"""

import os
import zipfile


def compress(filename, level=9):
    """Creates a zip archive from a single file.

    Parameters
    ----------
    filename : str
        Name of the file to be compressed.
    level : int, optional
        Level of compression. Integers 0 through 9 are accepted. The default is 9.

    Returns
    -------
    None.

    See Also
    --------
    zipfile.ZipFile

    """
    name = os.path.splitext(filename)[0]
    output_file = name + '.zip'
    zip_options = dict(compression=zipfile.ZIP_DEFLATED, compresslevel=level)
    with zipfile.ZipFile(output_file, mode='w', **zip_options) as compressed:
        compressed.write(filename)


def decompress(filename):
    """Decompresses a zip archive.

    Parameters
    ----------
    filename : str
        Name of the zip archive.

    Returns
    -------
    data : dict
        The keys of the dictionary are the compressed file names (without extension),
        the corresponding values are their contents.

    See Also
    --------
    zipfile.ZipFile

    """
    data = {}
    with zipfile.ZipFile(filename, mode='r') as compressed:
        for file in compressed.namelist():
            with compressed.open(file, mode='r') as thisfile:
                name = os.path.splitext(thisfile.name)[0]
                data[name] = compressed.read(thisfile.name).decode()
    return data
