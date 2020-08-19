# -*- coding: utf-8 -*-
"""
This module implements computational geometry algorithms, needed for other modules.

All the examples assume that the modules `numpy` and `matplotlib.pyplot` were imported as `np`
and `plt`, respectively.

Classes
-------
.. autosummary::
    :nosignatures:

    Mesh
    TriMesh

Functions
---------
.. autosummary::

    _polygon_area

"""
from abc import ABC, abstractmethod

import numpy as np


class Mesh(ABC):
    """Data structure for a general mesh.

    This class by no means wants to provide intricate functionalities and does not strive to be
    efficient at all. Its purpose is to give some useful features the project relies on.
    The two main entities in the mesh are the `vertices` and the `cells`. They are expected to be
    passed by the user, so reading from various mesh files is not implemented. This keeps the
    class simple, requires few package dependencies and keeps the class focused as there are
    powerful tools to convert among mesh formats (see e.g. `meshio
    <https://github.com/nschloe/meshio>`_).

    Parameters
    ----------
    vertices : ndarray
        2D numpy array with 2 columns, each row corresponding to a vertex, and the two columns
        giving the Cartesian coordinates of the vertex.
    cells : ndarray
        Cell-vertex connectivities in a 2D numpy array, in which each row corresponds to a
        cell and the columns are the vertices of the cells. It is assumed that all the cells
        have the same number of vertices.

    """

    # @abstractmethod
    def __init__(self, vertices, cells):
        """
        .. todo::
            Error checking

        """
        self.vertices = vertices
        self.cells = cells
        self.vertex_sets = {}
        self.cell_sets = {}

    def create_cell_set(self, name, cells):
        """Forms a group from a set of cells.

        Parameters
        ----------
        name : str
            Name of the cell set.
        cells : list
            List of cells to be added to the set.

        Returns
        -------
        None

        """
        self.cell_sets[name] = cells

    def create_vertex_set(self, name, vertices):
        """Forms a group from a set of vertices.

        Parameters
        ----------
        name : str
            Name of the vertex set.
        vertices : list
            List of vertices to be added to the set.

        """
        self.vertex_sets[name] = vertices

    @staticmethod
    def _isvector(array):
        """Decides whether the input is a vector.

        Parameters
        ----------
        array : ndarray
            Numpy array to be checked.

        Returns
        -------
        bool
            True if the input is a 1D array or if it is a column or row vector. Otherwise, False.

        See Also
        --------
        _ismatrix

        """
        shape = np.shape(array)
        if len(shape) == 1:
            return True
        elif len(shape) == 2 and shape[0] == 1 or shape[1] == 1:
            return True
        else:
            return False

    @staticmethod
    def _ismatrix(array):
        """Decides whether the input is a matrix.

        Parameters
        ----------
        array : ndarray
            Numpy array to be checked.

        Returns
        -------
        bool
            True if the input is a 2D array. Otherwise, False.

        See Also
        --------
        _isvector

        """
        shape = np.shape(array)
        return len(shape) == 2


class TriMesh(Mesh):
    """Unstructured triangular mesh.

    Vertices and cells are both stored as numpy arrays. This makes the simple mesh manipulations
    easy and provides interoperability with the whole scientific Python stack.

    Parameters
    ----------
    vertices : ndarray
        2D numpy array with 2 columns, each row corresponding to a vertex, and the two columns
        giving the Cartesian coordinates of the vertex.
    cells : ndarray
        Cell-vertex connectivities in a 2D numpy array, in which each row corresponds to a
        cell and the columns are the vertices of the cells. It is assumed that all the cells
        have the same number of vertices.

    """

    def __init__(self, vertices, cells):
        super().__init__(vertices, cells)

    def change_vertex_numbering(self, orientation, inplace=False):
        """Changes cell vertex numbering.

        Parameters
        ----------
        orientation : {'ccw', 'cw'}
            Vertex numbering within a cell, either `'ccw'` (counter-clockwise, default) or
            `'cw'` (clock-wise).
        inplace : bool, optional
            If True, the vertex ordering is updated in the mesh object. The default is False.

        Returns
        -------
        reordered_cells : ndarray
            Same format as the :code:`cells` member variable, with the requested vertex ordering.

        Notes
        -----
        Supposed to be used with planar P1 or Q1 finite s.

        Examples
        --------
        >>> mesh = TriMesh(np.array([[1, 1], [3, 5], [7,3]]), np.array([0, 1, 2]))
        >>> mesh.change_vertex_numbering('ccw')
        array([[2, 1, 0]])

        """
        if np.ndim(self.cells) == 1:  # single cell, given as a vector
            self.cells = self.cells.reshape((1, np.size(self.cells)))  # reshape to a matrix
        if np.size(self.cells, 1) < 3:
            raise Exception('Cells must have at least 3 vertices.')
        reordered_cells = np.empty_like(self.cells, dtype=int)
        for i, cell_vertices in enumerate(self.cells):
            coordinates = self.vertices[cell_vertices]
            signed_area = _polygon_area(list(coordinates[:, 0]), list(coordinates[:, 1]))
            if (signed_area < 0 and orientation == 'ccw') or \
               (signed_area > 0 and orientation == 'cw'):
                cell_vertices = cell_vertices[::-1]
            reordered_cells[i] = cell_vertices
        if inplace:
            self.cells = reordered_cells
        return reordered_cells

    def scale(self, factor, inplace=False):
        """Scales the geometry by modifying the coordinates of the vertices.

        Parameters
        ----------
        factor : float
            Each vertex coordinate is multiplied by this non-negative number.
        inplace : bool, optional
            If True, the vertex positions are updated in the mesh object. The default is False.

        Returns
        -------
        None.

        """
        scaled_vertices = factor*self.vertices
        if inplace:
            self.vertices = scaled_vertices
        return scaled_vertices

    def rotate(self, angle, point=(0, 0), inplace=False):
        """Rotates a 2D mesh about a given point by a given angle.

        Parameters
        ----------
        angle : float
            Angle of rotation, in radians.
        point : list or tuple, optional
            Coordinates of the point about which the mesh is rotated. If not given, it is the origin.
        inplace : bool, optional
            If True, the vertex positions are updated in the mesh object. The default is False.

        Returns
        -------
        rotated_vertices : ndarray
            Same format as the :code:`vertices` member variable, with the requested rotation.

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
        rotation_matrix = np.array(
            [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        rotated_vertices = np.empty_like(self.vertices)
        for i, vertex in enumerate(self.vertices):
            vertex_local = vertex - point
            rotated_vertex_local = rotation_matrix.dot(vertex_local)
            rotated_vertex = rotated_vertex_local + point
            rotated_vertices[i] = rotated_vertex
        if inplace:
            self.vertices = rotated_vertices
        return rotated_vertices

    def cell_area(self, cell):
        """Computes the area of a cell.

        Parameters
        ----------
        cell : int
            Cell label.

        Returns
        -------
        area : float
            Area of the cell.

        See Also
        --------
        cell_set_area

        """
        vertices = self.cells[cell]
        coordinates = self.vertices[vertices]
        return _polygon_area(list(coordinates[:, 0]), list(coordinates[:, 1]))

    def cell_set_area(self, cell_set):
        """Computes the area of a cell set.

        Parameters
        ----------
        cell_set : str
            Name of the cell set.

        Returns
        -------
        area : float
            Area of the cell set.

        See Also
        --------
        cell_area

        """
        area = 0
        for cell in self.cell_sets[cell_set]:
            area += self.cell_area(cell)
        return area

    def write_inp(self, filename):
        """Writes the mesh to an Abaqus .inp file.

        Parameters
        ----------
        filename : str
            Path to the mesh file.

        Returns
        -------
        None

        """
        # Node coordinates
        inp_file = '*NODE\n'
        for global_node_num, coordinate in enumerate(self.vertices, 1):
            # TODO: avoid reallocation by temporary strings as done below
            inp_file += str(global_node_num) + ', ' + str(list(coordinate))[1:-1] + '\n'

        # Cell-vertex connectivities
        element_type = 'CPS3'  # 3-node linear plane stress element; can be changed in Abaqus
        inp_file += '*ELEMENT, TYPE=' + element_type + '\n'
        global_elem_num = 0
        for elem_vertices in self.cells:
            global_elem_num += 1
            # TODO: avoid reallocation by temporary strings as done below
            inp_file += str(global_elem_num) + ', ' + str(list(elem_vertices+1))[1:-1] + '\n'

        # Element groups
        element_sets = ''  # Grow a new string to avoid extensive reallocation
        for group_name, elements in self.cell_sets.items():
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
        for group_name, nodes in self.vertex_sets.items():
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
    number of cells in the mesh, no input checking is performed.

    Examples
    --------
    >>> _polygon_area([0, 1, 1], [0, 0, 1])
    0.5

    """
    return 1/2*(sum(i*j for i, j in zip(x, y[1:]+y[:1])) -
                sum(i*j for i, j in zip(x[1:]+x[:1], y)))
