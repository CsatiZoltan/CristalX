# -*- coding: utf-8 -*-

"""
This module provides *general* utility functions used by the **grains** package. The specific
helper functions reside in the proper module. For example, a function that works on a general
list goes here, but a computational geometry algorithm goes to the **geometry** module. The
functions in the **utils** module can be interesting for other projects too, partly because of
their general scope, and partly because of the few dependencies.
"""

import os
import zipfile

import numpy as np


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


def index_list(lst, indices):
    """Index a list by another list.

    Parameters
    ----------
    lst : list
        List to be indexed.
    indices : list
        Indices of the original list that will form the new list.

    Returns
    -------
    list
        Members of `lst`, selected by `indices`.

    Examples
    --------
    >>> index_list(['c', ['nested', 'list'], 13], [1, 2])
    [['nested', 'list'], 13]

    """
    return [lst[idx] for idx in indices]


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


def argsorted(sequence, reverse=False):
    """Return the indices that would sort a list or a tuple.

    Implementation is taken from https://stackoverflow.com/a/6979121/4892892.

    Parameters
    ----------
    sequence : list, tuple
        Input sequence in which the sorted indices to be found.
    reverse : bool
        If set to True, then the elements are sorted as if each comparison was reversed.

    Returns
    -------
    list
        List of indices that would sort the input list/tuple.

    See Also
    --------
    sorted
    numpy.argsort

    Examples
    --------
    >>> argsorted([2, 1.1, 1.1])
    [1, 2, 0]
    >>> argsorted([2, 1.1, 1.1], reverse=True)
    [0, 1, 2]
    >>> argsorted(())
    []

    """
    return sorted(range(len(sequence)), key=lambda k: sequence[k], reverse=reverse)


def map_inplace(function, __iterable):
    """Apply a function to each member of an iterable in-place.

    Parameters
    ----------
    function : function object
        Function to be applied to the entries of the iterable.
    __iterable : iterable
        Iterable.

    Notes
    -----
    Comprehensions or functional tools work on iterators, thereby not modifying the original
    container (https://stackoverflow.com/a/4148523/4892892). For in-place modification, the
    conventional for loop approach is used (https://stackoverflow.com/a/4148525/4892892).

    Examples
    --------
    >>> lst = ['2', 2]; func = lambda x: x*2
    >>> map_inplace(func, lst); lst
    ['22', 4]
    >>> lifespan = {'cat': 15, 'dog': 12}; die_early = lambda x: x/2
    >>> map_inplace(die_early, lifespan); lifespan
    {'cat': 7.5, 'dog': 6.0}

    """
    if type(__iterable) == dict:
        for key, value in __iterable.items():
            __iterable[key] = function(value)
    else:
        for i in range(len(__iterable)):
            __iterable[i] = function(__iterable[i])


def non_unique(array, axis=None):
    """Finds indices of non-unique elements in a 1D or 2D ndarray.

    Parameters
    ----------
    array : ndarray
        Array in which the non-unique elements are searched.
    axis : {None, 0, 1}, optional
        The axis to operate on. If None, `array` will be flattened. If an integer, the subarrays
        indexed by the given axis will be flattened and treated as the elements of a 1-D array with
        the dimension of the given axis. Object arrays or structured arrays that contain objects are
        not supported if the `axis` kwarg is used. The default is None.

    Returns
    -------
    nonunique_values : list
        Unique (individual, row or column) entries.
    nonunique_indices : list
        Each element of the list corresponds to non-unique elements, whose indices are given in
        a 1D numpy array.

    Examples
    --------
    In a 1D array, the repeated values and their indices are found by

    >>> val, idx = non_unique(np.array([1, -1, 0, -1, 2, 5, 0, -1]))
    >>> val
    [-1, 0]
    >>> idx
    [array([1, 3, 7]), array([2, 6])]

    In the matrix below, we can see that rows 0 and 2 are identical, as well as rows 1 and 4.

    >>> val, idx = non_unique(np.array([[1, 3], [2, 4], [1, 3], [-1, 0], [2, 4]]), axis=0)
    >>> val
    [array([1, 3]), array([2, 4])]
    >>> idx
    [array([0, 2]), array([1, 4])]

    By transposing the matrix above, the same holds for the columns.

    >>> val, idx = non_unique(np.array([[1, 2, 1, -1, 2], [3, 4, 3, 0, 4]]), axis=1)
    >>> val
    [array([1, 3]), array([2, 4])]
    >>> idx
    [array([0, 2]), array([1, 4])]

    If the dimensions along which to find the duplicates are not given, the input is flattened
    and the indexing happens in C-order (row-wise).

    >>> val, idx = non_unique(np.array([[1, 2, 1, -1, 2], [3, 4, 3, 0, 4]]))
    >>> val
    [1, 2, 3, 4]
    >>> idx
    [array([0, 2]), array([1, 4]), array([5, 7]), array([6, 9])]

    """
    if axis not in {None, 0, 1}:
        raise Exception('Parameter `axis` must be one of the following: None, 0, 1.')
    # Indices of the repeated elements
    _, idx, counts = np.unique(array, axis=axis, return_counts=True, return_inverse=True)
    distinct_indices = np.unique(idx)
    repeated_indices = distinct_indices[counts > 1]
    # Help in consistent indexing
    if axis is None:
        array = array.flatten()
    elif axis == 0:
        pass
    elif axis == 1:
        array = array.T
    # Find the repeated values and their positions at the same time
    nonunique_indices = []
    nonunique_values = []
    for i in repeated_indices:
        nonunique_indices.append(np.nonzero(idx == i)[0])
        nonunique_values.append(array[idx == i][0])
    return nonunique_values, nonunique_indices


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


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
