# -*- coding: utf-8 -*-

"""
This module provides utility functions that are useful for the **grains** package.
"""

import os
import zipfile


def duplicates(sequence):
    """Set of duplicate elements in a sequence.

    Parameters
    ----------
    sequence : sequence types (list, tuple, string, etc.)
        Sequence possibly containing repeating elements.

    Returns
    -------
    set
        Set of unique values.

    Notes
    -----
    Copied from https://stackoverflow.com/a/9836685/4892892

    Examples
    --------
    Note that the order of the elements in the resulting set does not matter.

    >>> a = [1, 2, 3, 2, 1, 5, 6, 5, 5, 5]  # list
    >>> duplicates(a)
    {1, 2, 5}
    >>> a = (1, 1, 0, -1, -1, 0)  # tuple
    >>> duplicates(a)
    {0, 1, -1}
    >>> a = 'abbcdkc'  # string
    >>> duplicates(a)
    {'c', 'b'}
    """
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in sequence if x in seen or seen_add(x))
    # turn the set into a list (as requested)
    return seen_twice


def toggle(lst):
    """Return True for False values and False for True values in a list.

    Parameters
    ----------
    lst : list
        An arbitrary list, possibly containing other lists.

    Returns
    -------
    list
        Element-wise logical not operator applied on the input list.

    Notes
    -----
    Solution taken from https://stackoverflow.com/a/51122372/4892892.

    Examples
    --------
    >>> toggle([True, False])
    [False, True]
    >>> toggle(['h', 0, 2.3, -2, 5, []])
    [False, True, False, False, False, True]

    """
    return list(map(lambda item: not item, lst))


def flatten_list(nested_list):
    """Merge a list of lists to a single list.

    Parameters
    ----------
    nested_list : list
        List containing other lists.

    Returns
    -------
    list
        Flattened list.

    Notes
    -----
    - Only a single level (i.e. list of lists) is handled, see the second example.
    - Several methods, such as list comprehension, monoid and loops, are proposed in
      https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists.
      Here, the list comprehension approach is used.

    Examples
    --------
    >>> nested_list = [['some'], ['items']]
    >>> flatten_list(nested_list)
    ['some', 'items']
    >>> multiply_nested_list = [[['item'], 'within', 'item']]
    >>> flatten_list(multiply_nested_list)
    [['item'], 'within', 'item']

    """
    return [item for sublist in nested_list for item in sublist]


def map_inplace(function, a_list):
    """Apply a function to each member of a list in-place.

    Parameters
    ----------
    function : function object
        Function to be applied to the entries of the list.
    a_list : list
        List.

    Notes
    -----
    A list comprehension or functional tools work on iterators, thereby not modifying the
    original container (https://stackoverflow.com/a/4148523/4892892). For in-place modification,
    the conventional for loop approach is used (https://stackoverflow.com/a/4148525/4892892).

    Examples
    --------
    >>> lst = ['2', 2]; func = lambda x: x*2
    >>> map_inplace(func, lst); lst
    ['22', 4]

    """
    for i in range(len(a_list)):
        a_list[i] = function(a_list[i])


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


def decompress(filename, path=None):
    """Decompresses the contents of a zip archive into the current directory.

    Parameters
    ----------
    filename : str
        Name of the zip archive.

    path : str, optional
        Directory to extract to. The default is the directory the function is called from.

    See Also
    --------
    zipfile.ZipFile

    """
    with zipfile.ZipFile(filename, mode='r') as compressed:
        compressed.extractall(path=path)


def decompress_inmemory(filename):
    """Decompresses the contents of a zip archive into a dictionary.

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
