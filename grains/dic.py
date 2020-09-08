# -*- coding: utf-8 -*-
"""
This module is used to process field values obtained by digital image correlation (DIC).
Plotting, comparison with numerical solutions, etc. are implemented.

Classes
-------
.. autosummary::
    :nosignatures:

    DIC

"""
import numpy as np
import matplotlib.pyplot as plt

from grains.geometry import TriMesh


class DIC:
    """
    Performs operations on digital image correlation data.

    Notes
    -----
    Throughout the class methods, the following terms are used. When not precised, `image` refers to
    the DIC data plotted as an image. We can perceive the image as a `pixel grid`, in which the
    vertices are the centres of the pixels. An image coordinate system :math:`(X,Y)` can be defined
    on this grid such that the vertices have integer coordinates, the origin :math:`A` is at the
    top left hand corner, and the :math:`Y`-axis points towards the right, while the :math:`X`-axis
    points to the bottom. When the experimental result is compared with a numerical solution
    (assumed to be available at the `nodes` of a `mesh`), one needs to map values defined on the
    DIC grid onto the nodes of the mesh, and vice versa. The mesh exists on the physical domain,
    for which a physical coordinate system :math:`(x,y)` is attached to. Hence, there is a need
    to express :math:`(X,Y)` in terms of :math:`(x,y)`.
    Let :math:`A(a_x,a_y)` be the origin of the :math:`(X,Y)` coordinate system given in terms of
    the :math:`x,y` coordinates, and :math:`s[\\mathrm{pixel/mm}]` is the scale for the DIC image.
    The coordinate transformation is then given by

    .. math::
        \\begin{gather}
           x = a_x + \\frac{Y}{s} \\\\
           y = a_y - \\frac{X}{s}
        \\end{gather}

    The DIC grid in the physical coordinate system :math:`(x,y)` is called `physical grid`.

    It should be noted that the special alignment of :math:`(x,y)` with respect to :math:`(X,Y)`
    is not a major constraint. In practical applications (when a Cartesian coordinate system is
    used), the orientation of :math:`(x,y)` in this way is a convention.

    """

    def __init__(self, u, v):
        self.u = u
        self.v = v
        self.origin = None
        self.scale = None

    def plot_displacement(self, component, ax=None):
        """Plots a component of the displacement field.

        Parameters
        ----------
        component : {1, 2}
            Component to be plotted, where 1 corresponds to the first, 2 corresponds to the
            second component of the displacement field.
        ax : matplotlib.axes.Axes, optional
            The `Axes` instance the grid resides in. The default is None, in which case a new
            `Axes` within a new figure is created.

        Returns
        -------
        None

        See Also
        --------
        plot_strain

        """
        if component == 1:
            field = self.u
        elif component == 2:
            field = self.v
        else:
            raise ValueError('Component must be either 1 or 2.')
        if not ax:
            ax = plt.figure().add_subplot()
        ax.set_aspect('equal')
        ax.imshow(field)

    def plot_strain(self, component, minval=None, maxval=None, legend=True):
        """Plots a component of the infinitesimal strain tensor.

        The partial derivatives of the displacement field are computed with numerical
        differentiation.

        Parameters
        ----------
        component : tuple, {(1,1), (1,2), (2,1), (2,2)}
            Component to be plotted,
            where

                - (1,1) denotes :math:`\\varepsilon_{11}`
                - (1,2) and (2,1) denote :math:`\\varepsilon_{12} = \\varepsilon_{21}`
                - (2,2) denotes :math:`\\varepsilon_{22}`

                for the infinitesimal strain tensor

            .. math::
                \\varepsilon = \\begin{pmatrix}
                                   \\varepsilon_{11} & \\varepsilon_{12} \\\\
                                   \\varepsilon_{21} & \\varepsilon_{22}
                               \\end{pmatrix}

        Returns
        -------
        None

        See Also
        --------
        plot_displacement,
        :func:`numpy.gradient`

        """
        if component == (1, 1):
            strain = np.gradient(self.u)[1]
            label = r'$\varepsilon_{xx}$'
        elif component == (2, 2):
            strain = np.gradient(self.u)[0]
            label = r'$\varepsilon_{yy}$'
        elif component in {(1, 2), (2, 1)}:
            strain = 0.5 * (np.gradient(self.u)[0] + np.gradient(self.v)[1])
            label = r'$\varepsilon_{xy}$'
        else:
            raise ValueError('Incorrect strain tensor component.')
        plt.matshow(strain, interpolation='none', vmin=minval, vmax=maxval)
        if legend:
            cbar = plt.colorbar(orientation='horizontal', format='%.2f', aspect=100,
                                label=label)
            cbar.ax.set_xlabel(label, fontsize=30)

    def set_transformation(self, origin, pixels_per_physicalunit):
        """Sets the transformation rule between the pixel and the physical coordinate systems.

        To determine the position and the size of the DIC image in the physical space,
        a transformation rule needs to be given, as described in the `Notes` section of the
        :class:`~grains.dic.DIC` class.

        Parameters
        ----------
        origin : tuple
            2-tuple of float, the coordinates of the origin of the DIC grid (upper left corner)
            in the physical coordinate system.
        pixels_per_physicalunit : float

        Returns
        -------
        None

        """
        self.origin = origin
        self.scale = pixels_per_physicalunit

    def plot_pixelgrid(self, ax=None):
        """Plots the DIC grid in the image coordinate system.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The `Axes` instance the grid resides in. The default is None, in which case a new
            `Axes` within a new figure is created.

        Returns
        -------
        ax : matplotlib.axes.Axes

        See Also
        --------
        plot_physicalgrid

        Notes
        -----
        For a DIC image with :math:`m \\times n` number of pixels, the number of grid lines is
        :math:`(m + 1) + (n + 1)`. Plotting all these lines would not only slow down the program,
        it would also make the grid practically indistinguishable from a filled rectangle.
        Therefore, for high resolution images, only grid lines around the boundary of the image
        are plotted. The target resolution above which this strategy is used can be given in the
        class constructor.


        """
        # Coordinate values spanning the pixel grid
        image_height, image_width = np.shape(self.u)
        X = np.arange(0, image_height)
        Y = np.arange(0, image_width)
        # Draw the grid
        X_grid, Y_grid = np.meshgrid(X, Y)
        if not ax:
            ax = plt.figure().add_subplot()
        ax.set_aspect('equal')
        ax.plot(Y_grid, X_grid, 'k.', markersize=1)
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        return ax

    def plot_physicalgrid(self, ax=None):
        """Plots the DIC grid in the physical coordinate system.

        This method is only available once the relation between the image coordinate system and
        the physical coordinate system has been set up by the :meth:`set_transformation` method.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The `Axes` instance the grid resides in. The default is None, in which case a new
            `Axes` within a new figure is created.

        Returns
        -------
        axes : matplotlib.axes.Axes

        See Also
        --------
        plot_pixelgrid

        Notes
        -----
        See in the :meth:`plot_pixelgrid` method.

        """
        # Coordinate values spanning the pixel grid
        image_height, image_width = np.shape(self.u)
        X = np.arange(0, image_height)
        Y = np.arange(0, image_width)
        # Coordinate values spanning the physical grid
        x = self.origin[0] + Y/self.scale
        y = self.origin[1] - np.flip(X)/self.scale
        # Draw the grid
        x_grid, y_grid = np.meshgrid(x, y)
        if not ax:
            ax = plt.figure().add_subplot()
        ax.set_aspect('equal')
        ax.plot(x_grid, y_grid, 'k.', markersize=1)
        return ax

    def plot_superimposedmesh(self, mesh, *args, **kwargs):
        """Plots the DIC grid with a mesh superimposed on it.

        Since the finite element mesh represents a physical domain, the DIC grid is plotted in
        the physical coordinate system, which is determined by the :meth:`set_transformation`
        method.

        .. note::
            Maybe no need to couple this module with the `geometry` module just for the
            sake of this function. It is probably easier to simply pass the nodes and the
            elements, and call `triplot`.

        Parameters
        ----------
        mesh : Mesh
            A subclass of :class:`grains.geometry_new.Mesh` to be plotted along with the physical
            DIC grid.
        args : optional
            The same as for :meth:`matplotlib.axes.Axes.plot`. The settings given here influence
            the mesh plotting.
        kwargs : optional
            The same as for :meth:`matplotlib.axes.Axes.plot`. The settings given here influence
            the mesh plotting.

        Returns
        -------
        ax : matplotlib.axes.Axes
            The Axes object on which the plot is drawn.

        See Also
        --------
        :meth:`grains.geometry_new.TriMesh.plot`

        Examples
        --------
        Let us create a grid and a random mesh.

        >>> x_grid, y_grid = np.mgrid[-1:1:100j, 1:2:50j]
        >>> exact_solution = lambda x, y: 1 - (x + y**2) * np.sign(x)
        >>> grid = DIC(np.random.rand(*np.shape(x_grid.T)), np.random.rand(*np.shape(x_grid.T)))
        >>> grid.set_transformation((-1, 2), 50)
        >>> n_nodes = 100  # modify this to see how good the interpolated solution is
        >>> msh = TriMesh(*TriMesh.sample_mesh(1, n_nodes))
        >>> grid.plot_superimposedmesh(msh, linewidth=3, markerfacecolor='b')
        >>> plt.show()

        """
        ax = plt.figure().add_subplot()
        self.plot_physicalgrid(ax)
        mesh.plot(ax=ax, *args, **kwargs)
        return ax
