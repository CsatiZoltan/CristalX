# -*- coding: utf-8 -*-
"""
Extracting and processing meshes from .med files.
The functions were tested on the MEDCoupling API, version 9.4.0.

.. todo:: Support renumbering (https://docs.salome-platform.org/latest/dev/MEDCoupling/user/html/data_optimization.html).

Getting help:
-------------
* This module relies on the Python interface of MEDCoupling. `Click here <https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/index.html>`_ for the latest documentation.
* `User's manual <https://docs.salome-platform.org/latest/dev/MEDCoupling/user/html/index.html>`_ for the Python interface
* To know more about the MED file format, which is a specialization of HDF5, see the `documentation <https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/med-file.html>`_.
  For a discussion on the relation between the MED format and the APIS, see `this page
  <https://www.salome-platform.org/user-section/about/med>`_ and `that one
  <https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/library.html>`_.
* The definitions, such as `group`, used in this module are from the `development guide <https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/glossary.html>`_.
* A (mostly English) `tutorial <https://docs.salome-platform.org/latest/dev/MEDCoupling/tutorial/index.html>`_ for the Python interface to MEDCoupling is also useful.
  Particularly interesting are the `mesh manipulation examples <https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/medcouplingpyexamples.html#ExamplesMeshes>`_
* `Main page <https://salome-platform.org/user-section/documentation/current-release>`_ of the documentation
"""

import warnings

import numpy as np

from grains import HAS_MED
if HAS_MED:
    from MEDLoader import MEDFileData
    # from MEDCoupling import MEDCouplingUMesh, MEDCouplingMesh, MEDCouplingUMeshCellIterator
else:
    raise ImportError('This module requires the MEDCoupling module from Salome.')


def read_mesh(filename):
    """Reads a mesh file in .med format.
    Only one mesh, the first one, is supported. However, that mesh can contain groups.

    Parameters
    ----------
    filename : str
        Path to the mesh file.

    Returns
    -------
    MEDFileUMesh
        Represents an unstructured mesh. For details, see the manual on
        https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/classMEDCoupling_1_1MEDFileUMesh.html
    """
    # Read the mesh
    data = MEDFileData(filename)
    mesh = data.getMeshes()
    mesh_names = mesh.getMeshesNames()
    n_mesh = len(mesh_names)
    if n_mesh > 1:
        warnings.warn('{0} meshes detected, only the fist one is considered.'.format(n_mesh),
                      RuntimeWarning)
    mesh = mesh.getMeshWithName(mesh_names[0])
    if mesh.getMeshDimension() != 2:
        warnings.warn('Only 2D meshes are supported.', RuntimeWarning)
    return mesh


def get_nodes(mesh):
    """Obtains the nodes and the node groups of a mesh.

    Parameters
    ----------
    mesh : MEDFileUMesh
        Unstructured mesh.

    Returns
    -------
    nodes : ndarray
        2D numpy array with 2 columns, each row corresponding to a node, and the two columns
        giving the Cartesian coordinates of the nodes.
    node_groups : dict
        The keys in the dictionary are the node group names, while the values are list of integers,
        giving the nodes that belong to the particular group.

    See Also
    --------
    :func:`get_elements`,
    `getGroupArr <https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/classMEDCoupling_1_1MEDFileMesh.html#a4398c05f015e52d0d380eb39c6e4b942>`_

    """
    # Get the coordinates of all the nodes in the mesh
    nodes = mesh.getCoords().toNumPyArray()
    # Get the node groups
    group_names = mesh.getGroupsOnSpecifiedLev(1)
    node_groups = {}
    for group_name in group_names:
        node_groups[group_name] = mesh.getGroupArr(1, group_name).toNumPyArray()
    return nodes, node_groups


def get_elements(mesh, numbering='global'):
    """Obtains the elements for each group of a mesh.

    Elements of the same dimension as the mesh are collected (e.g. faces for a 2D mesh).

    .. todo:: put those elements that do not belong to any group into an automatically created group
    .. todo:: support ordering `elements` in alphabetical order
    .. todo:: implement the `'global'` strategy

    Parameters
    ----------
    mesh : MEDFileUMesh
        Unstructured mesh.
    numbering : {'global'}, optional
        Determines how to allocate element numbers in the mesh.
            'global': numbers the elements without taking into account which group they belong to.
            Use this strategy if you are not sure whether an element belongs to more than one group.
            'group': numbers the elements group-wise. This is much faster than the `'global'`
            strategy, but use this option if you are sure that the groups of the mesh do not
            contain common elements.
        The default is 'global'.

    Returns
    -------
    elements : ndarray
        Element-node connectivities in a 2D numpy array, in which each row corresponds to an
        element and the columns are the nodes of the elements. It is assumed that all the
        elements have the same number of nodes.
    element_groups : dict
        The keys in the dictionary are the element group names, while the values are list of
        integers, giving the elements that belong to the particular group.

    Warnings
    --------
    Currently, elements that do not fit into any groups are discarded.

    See Also
    --------
    get_nodes
    change_node_numbering

    Notes
    -----
    The element-node connectivities are read from the mesh. If you want to change the ordering
    of the nodes, use the :func:`change_node_numbering` function.

    Both this and the :func:`get_nodes` function relies on `getGroupsOnSpecifiedLev
    <https://docs.salome-platform.org/latest/dev/MEDCoupling/developer
    /classMEDCoupling_1_1MEDFileMesh.html#a2d59097b6d14b95c7d2aeee9f39b0438>`_ to obtain the groups
    based on a parameter, called `meshDimRelToMaxExt`. This parameter designates the relative
    dimension of the mesh entities whose IDs are required. If it is 1, it denotes the nodes. If
    0, entities of the same dimension as the mesh are meant (e.g. group of volumes for a 3D mesh,
    or group of faces for a 2D mesh). When -1, entities of spatial dimension immediately below
    that of the mesh are collected (e.g. group of faces for a 3D mesh, or group of edges for a
    2D mesh). For -2, entities of two dimensions below that of the mesh are fetched (e.g. group of
    edges for a 3D mesh).

    """
    group_names = mesh.getGroupsOnSpecifiedLev(0)
    n_element = mesh.getDistributionOfTypes(0)[1]
    elements = np.empty((n_element, 3), dtype=int)
    element_groups = {}
    if numbering == 'group':  # process the elements group-wise
        i = 0  # indexing elements in the whole mesh
        for group_name in group_names:
            local_mesh = mesh.getGroup(0, group_name)
            n_element_in_group = local_mesh.getNumberOfCells()
            local_elements = np.empty(n_element_in_group, dtype=int)
            for j, cell in enumerate(local_mesh):  # j: indexing elements within a group
                elements[i, :] = list(cell.getAllConn()[1:])
                local_elements[j] = i
                i += 1
            element_groups[group_name] = local_elements
        if i < n_element:  # some elements are not part of any group
            elements = elements[0:i]
            warnings.warn('Elements that do not belong to any of the groups were discarded.')
    elif numbering == 'global':
        raise ValueError('Strategy not yet implemented.')
    else:
        raise ValueError('Unknown strategy. Choose either "global" or "group".')
    return elements, element_groups


def change_node_numbering(elements, nodes, orientation='ccw'):
    """Changes element node numbering.

    Parameters
    ----------
    elements : ndarray
        Element-node connectivities in a 2D numpy array, in which each row corresponds to an
        element and the columns are the nodes of the elements. It is assumed that all the
        elements have the same number of nodes.
    nodes : ndarray
        2D numpy array with 2 columns, each row corresponding to a node, and the two columns
        giving the Cartesian coordinates of the nodes.
    orientation : {'ccw', 'cw'}
        Node numbering within an element, either `'ccw'` (counter-clockwise, default) or
        `'cw'` (clock-wise).

    Returns
    -------
    reordered_elements : ndarray
        Same format as the input `elements`, with the requested node ordering.

    Notes
    -----
    Supposed to be used with planar P1 or Q1 elements.

    Examples
    --------
    >>> import numpy as np
    >>> change_node_numbering(np.array([0, 1, 2]), np.array([[1, 1], [3, 5], [7,3]]))
    array([[2, 1, 0]])

    """
    if np.ndim(elements) == 1:  # single element, given as a vector
        elements = elements.reshape((1, np.size(elements)))  # reshape to a matrix
    if np.size(elements, 1) < 3:
        raise Exception('Elements must have at least 3 nodes.')
    reordered_elements = np.empty_like(elements, dtype=int)
    for i, element_nodes in enumerate(elements):
        coordinates = nodes[element_nodes]
        signed_area = _polygon_area(list(coordinates[:, 0]), list(coordinates[:, 1]))
        if (signed_area < 0 and orientation == 'ccw') or (signed_area > 0 and orientation == 'cw'):
            element_nodes = element_nodes[::-1]
        reordered_elements[i] = element_nodes
    return reordered_elements


def rotate_mesh(nodes, angle, point=(0, 0)):
    """Rotates a 2D mesh about a given point by a given angle.

    .. todo:: Also put it to the :mod:`grains.geometry` module (at least for rotating a single point).

    Parameters
    ----------
    nodes : ndarray
        2D numpy array with 2 columns, each row corresponding to a node, and the two columns
        giving the Cartesian coordinates of the nodes.
    angle : float
        Angle of rotation, in radians.
    point : list or tuple, optional
        Coordinates of the point about which the mesh is rotated. If not given, it is the origin.

    Returns
    -------
    rotated_nodes : ndarray
        Same format as the input `nodes`, with the requested rotation.

    Notes
    -----
    Rotating a point :math:`P` given by its coordinates in the global coordinate system as
    :math:`P(x,y)` around a point :math:`A(x,y)` by an angle :math:`\\alpha` is done as follows.

    1. The coordinates of :math:`P` in the local coordinate system, the origin of which is
    :math:`A`, is expressed as

    .. math::
        P(x',y') = P(x,y) - A(x,y).

    2. The rotation is performed in the local coordinate system as :math:`P'(x',y') = RP(x',y')`,
    where :math:`R` is the rotation matrix:

    .. math::
        R = \\begin{pmatrix}
            \\cos\\alpha & -\\sin\\alpha \\\\
            \\sin\\alpha & \\hphantom{-}\\cos\\alpha
            \\end{pmatrix}.

    3. The rotated point :math:`P'` is expressed in the original (global) coordinate system:

    .. math::
        P'(x,y) = P'(x',y') + A(x,y).

    """
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    rotated_nodes = np.empty_like(nodes)
    for i, node in enumerate(nodes):
        node_local = node - point
        rotated_node_local = rotation_matrix.dot(node_local)
        rotated_node = rotated_node_local + point
        rotated_nodes[i] = rotated_node
    return rotated_nodes


def write_inp(filename, nodes, elements, element_groups=None, node_groups=None):
    """Writes an unstructured mesh to an Abaqus .inp file.

    Parameters
    ----------
    filename : str
        Path to the mesh file.
    nodes : ndarray
        2D numpy array with 2 columns, each row corresponding to a node, and the two columns
        giving the Cartesian coordinates of the nodes.
    elements : ndarray
        Element-node connectivities in a 2D numpy array, in which each row corresponds to an
        element and the columns are the nodes of the elements. It is assumed that all the
        elements have the same number of nodes.
    element_groups : dict, optional
        The keys in the dictionary are the group names, while the values are list of integers,
        giving the elements that belong to the particular group.
    node_groups : dict, optional
        The keys in the dictionary are the node group names, while the values are list of integers,
        giving the nodes that belong to the particular group.


    Returns
    -------
    None

    See Also
    --------
    read_mesh
    get_nodes
    get_elements

    """
    # Node coordinates
    inp_file = '*NODE\n'
    for global_node_num, coordinate in enumerate(nodes, 1):
        # TODO: avoid reallocation by temporary strings as done below
        inp_file += str(global_node_num) + ', ' + str(list(coordinate))[1:-1] + '\n'

    # Element-node connectivities
    element_type = 'CPS3'  # 3-node linear plane stress element; can be changed in Abaqus
    inp_file += '*ELEMENT, TYPE=' + element_type + '\n'
    global_elem_num = 0
    for elem_nodes in elements:
        global_elem_num += 1
        # TODO: avoid reallocation by temporary strings as done below
        inp_file += str(global_elem_num) + ', ' + str(list(elem_nodes+1))[1:-1] + '\n'

    # Element groups
    element_sets = ''  # Grow a new string to avoid extensive reallocation
    for group_name, elements in element_groups.items():
        element_set = ''  # grow a new string to avoid extensive reallocation
        # Header
        element_set += '\n*ELSET, ELSET=' + group_name + '\n'
        # Elements belonging to this group
        for i, element in enumerate(elements+1, 1):
            element_set += str(element) + ', '
            if i % 16 == 0:  # Abaqus allows at most 16 entries per line
                element_set += '\n'
        element_sets += element_set
    inp_file += element_sets

    # Node groups
    node_sets = ''  # Grow a new string to avoid extensive reallocation
    for group_name, nodes in node_groups.items():
        node_set = ''  # grow a new string to avoid extensive reallocation
        # Header
        node_set += '\n*NSET, NSET=' + group_name + '\n'
        # Nodes belonging to this group
        for i, node in enumerate(nodes + 1, 1):
            node_set += str(node) + ', '
            if i % 16 == 0:  # Abaqus allows at most 16 entries per line
                node_set += '\n'
        node_sets += node_set
    inp_file += node_sets

    # Write to file
    with open(filename, 'w') as target:
        target.write(inp_file)


def element_area(element, element_nodes, node_coordinates):
    """Computes the area of an element.

    Parameters
    ----------
    element : int
        Element label.
    element_nodes : ndarray
        Element-node connectivities in a 2D numpy array, in which each row corresponds to an
        element and the columns are the nodes of the elements. It is assumed that all the
        elements have the same number of nodes.
    node_coordinates : ndarray
        2D numpy array with 2 columns, each row corresponding to a node, and the two columns
        giving the Cartesian coordinates of the nodes.

    Returns
    -------
    area : float
        Area of the element.

    See Also
    --------
    element_group_area

    """
    nodes = element_nodes[element]
    coordinates = node_coordinates[nodes]
    return _polygon_area(list(coordinates[:, 0]), list(coordinates[:, 1]))


def element_group_area(group, element_groups, element_nodes, node_coordinates):
    """Computes the area of an element group.

    Parameters
    ----------
    group : str
        Name of the element group
    element_groups : dict
        The keys in the dictionary are the group names, while the values are list of integers,
        giving the elements that belong to the particular group.
    element_nodes : ndarray
        Element-node connectivities in a 2D numpy array, in which each row corresponds to an
        element and the columns are the nodes of the elements. It is assumed that all the
        elements have the same number of nodes.
    node_coordinates : ndarray
        2D numpy array with 2 columns, each row corresponding to a node, and the two columns
        giving the Cartesian coordinates of the nodes.


    Returns
    -------
    area : float
        Area of the element group.

    See Also
    --------
    element_area

    """
    area = 0
    for element in element_groups[group]:
        area += element_area(element, element_nodes, node_coordinates)
    return area


def _polygon_area(x, y):
    """Computes the signed area of a non-self-intersecting, possibly concave, polygon.

    Directly taken from http://rosettacode.org/wiki/Shoelace_formula_for_polygonal_area#Python

    Parameters
    ----------
    x, y : list
        Coordinates of the consecutive vertices of the polygon.

    Returns
    -------
    float
        Area of the polygon.

    Warnings
    --------
    If numpy vectors are passed as inputs, the resulting area is incorrect! WHY?

    Notes
    -----
    The code is not optimized for speed and for numerical stability. Intended to be used to compute
    the area of finite element cells, in which case the numerical stability is not an issue
    (unless the cell is degenerate). As this function is called possibly as many times as the
    number of elements in the mesh, no input checking is performed.

    Examples
    --------
    >>> _polygon_area([0, 1, 1], [0, 0, 1])
    0.5

    """
    return 1/2*(sum(i*j for i, j in zip(x, y[1:]+y[:1])) -
                sum(i*j for i, j in zip(x[1:]+x[:1], y)))
