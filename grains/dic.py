# -*- coding: utf-8 -*-
"""
This module is used to process field values obtained by digital image correlation (DIC).
Plotting, comparison with numerical solutions, etc. are implemented.

Classes
-------
.. autosummary::
    :nosignatures:
    :toctree: classes/

    DIC

Functions
---------
.. autosummary::
    :toctree: functions/

    plot_strain

"""
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from deprecation import deprecated

from grains import __version__
from grains.geometry import TriMesh


class DIC:
    """Performs operations on digital image correlation data.

    Parameters
    ----------
    u, v : ndarray
        The two components of the displacement field, discretized on a rectangular grid.

    Raises
    ------
    ValueError
        If the displacement components are not discretized on the same grid or if the grid size
        is not at least 3-by-3. This latter requirement is not a problem in practice (which
        camera is not capable of shooting photos in 9 pixels resolution?), while it allows the
        evaluation of the numerical derivatives with second order accuracy on the boundaries too.

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
        if np.shape(u) <= (3, 3) or np.shape(u) <= (3, 3):
            raise ValueError('Displacement field must be given on a grid with minimum size 3-by-3.')
        if np.size(u) != np.size(v):
            raise ValueError('Displacement components must have the same size.')
        self.u = u
        self.v = v
        self.origin = None
        self.scale = None

    def strain(self, strain_measure):
        r"""Computes the strain tensor from the displacement vector.

        The symmetric strain tensor :math:`\boldsymbol{\varepsilon}` with components

        .. math::
            \boldsymbol{\varepsilon} = \begin{pmatrix}
                                             \varepsilon_{11} & \varepsilon_{12} \\
                                             \varepsilon_{21} & \varepsilon_{22}
                                         \end{pmatrix}

        is computed for the given strain measure. The partial derivatives of the displacement
        field :math:`(u,v)`, available on an :math:`m \times n` grid, are computed with numerical
        differentiation with second order accuracy.

        Parameters
        ----------
        strain_measure : {'infinitesimal', 'Green-Lagrange'}
            One of the following strain measures.

            - 'infinitesimal'

                .. math::
                    \varepsilon_{11} = \frac{\partial u}{\partial x} \\
                    \varepsilon_{12} = \frac{1}{2} \left( \frac{\partial u}{\partial y} +
                    \frac{\partial v}{\partial x} \right) \\
                    \varepsilon_{22} = \frac{\partial v}{\partial y} \\

            - 'Green-Lagrange'

                .. math::
                    \varepsilon_{11} = \frac{1}{2} \left( 2\frac{\partial u}{\partial x} +
                    \left(\frac{\partial u}{\partial x}\right)^2  + \left(\frac{\partial v}{
                    \partial x}\right)^2 \right) \\
                    \varepsilon_{12} = \frac{1}{2} \left( \frac{\partial u}{\partial y} +
                    \frac{\partial v}{\partial x} + \frac{\partial u}{\partial x} \frac{\partial
                    u}{\partial y} + \frac{\partial v}{\partial x} \frac{\partial v}{\partial
                    y} \right) \\
                    \varepsilon_{22} = \frac{1}{2} \left( 2\frac{\partial v}{\partial y} +
                    \left(\frac{\partial u}{\partial y}\right)^2  + \left(\frac{\partial v}{
                    \partial y}\right)^2 \right) \\

        Returns
        -------
        strain_tensor : ndarray
            Components of the strain tensor as an :math:`m \times n \times 3` array. The first
            two dimensions correspond to the grid points they are determined at, the third dimension
            gives the components of the tensor in the following order: :math:`\varepsilon_{11},
            \varepsilon_{12}, \varepsilon_{22}`.

        See Also
        --------
        equivalent_strain
        plot_strain
        :func:`numpy.gradient`

        Notes
        -----
        The derivatives of the displacement field are computed in the the physical coordinate
        system.

        Examples
        --------
        Consider a rectangular body with a linear displacement function in horizontal
        direction and no displacement vertically. The deformation of the body is such that
        only the horizontal strain is nonzero.

        >>> u = np.array([[0, 1e-3, 2e-3, 3e-3], [0, 1e-3, 2e-3, 3e-3], [0, 1e-3, 2e-3, 3e-3]])
        >>> v = np.zeros((3,4))
        >>> d = DIC(u, v)
        >>> E = d.strain('infinitesimal')
        >>> E[:, :, 0]
        array([[0.001, 0.001, 0.001, 0.001],
               [0.001, 0.001, 0.001, 0.001],
               [0.001, 0.001, 0.001, 0.001]])
        >>> np.allclose(E[:, :, 2], 0)
        True
        >>> np.allclose(E[:, :, 1], 0)
        True

        The Green-Lagrange strain tensor is slightly different. Comparing with the infinitesimal
        strain tensor shows how good the assumption of small deformations is.

        >>> E_gl = d.strain('Green-Lagrange')
        >>> E_gl[:, :, 0] - E[:, :, 0]
        array([[5.e-07, 5.e-07, 5.e-07, 5.e-07],
               [5.e-07, 5.e-07, 5.e-07, 5.e-07],
               [5.e-07, 5.e-07, 5.e-07, 5.e-07]])

        The other two strain components remain zero.

        >>> np.allclose(E_gl[:, :, 1:2], 0)
        True

        """
        strain_tensor = np.zeros((*self.u.shape, 3))
        dudx = np.gradient(self.u, edge_order=2)[1]
        dudy = np.gradient(self.u, edge_order=2)[0]
        dvdx = np.gradient(self.v, edge_order=2)[1]
        dvdy = np.gradient(self.v, edge_order=2)[0]
        if strain_measure == 'infinitesimal':
            strain_tensor[:, :, 0] = dudx
            strain_tensor[:, :, 1] = 0.5*(dudy + dvdx)
            strain_tensor[:, :, 2] = dvdy
        elif strain_measure == 'Green-Lagrange':
            strain_tensor[:, :, 0] = 0.5*(2*dudx + dudx**2 + dvdx**2)
            strain_tensor[:, :, 1] = 0.5*(dudy + dvdx + dudx*dudy + dvdx*dvdy)
            strain_tensor[:, :, 2] = 0.5*(2*dvdy + dudy**2 + dvdy**2)
        else:
            raise Exception('Strain measure {0} is not implemented.'.format(strain_measure))
        return strain_tensor

    @staticmethod
    def equivalent_strain(strain_tensor, strain_measure):
        r"""Computes a scalar quantity from the strain tensor.

        Parameters
        ----------
        strain_tensor : ndarray
            Components of the strain tensor :math:`\boldsymbol{\varepsilon}` as an
            :math:`m \times n \times 3` array. The first two dimensions correspond to the
            grid points they are determined at, the third dimension gives the components
            of the tensor in the following order: :math:`\varepsilon_{11}, \varepsilon_{12},
            \varepsilon_{22}`.
        strain_measure : {'von Mises'}
            One of the following strain measures.

            - 'von Mises'

                .. math::
                    \varepsilon_M := \sqrt{\frac{2}{3} \boldsymbol{\varepsilon}^{\mathrm{dev}} :
                    \boldsymbol{\varepsilon}^{\mathrm{dev}}} = \sqrt{\frac{2}{3} \left(
                    \boldsymbol{\varepsilon}:\boldsymbol{\varepsilon} - \frac{1}{3}(\mathrm{tr}
                    \boldsymbol{\varepsilon})^2 \right)}

                where :math:`\boldsymbol{\varepsilon}^{\mathrm{dev}}` is the deviatoric part of
                the strain tensor. This can further be simplified (see the notes below) to

                .. math::
                    \varepsilon_M = \frac{2}{3}\sqrt{\varepsilon_{11}^2 + \varepsilon_{22}^2 +
                    3\varepsilon_{12}^2 - \varepsilon_{11}\varepsilon_{22}}

        Returns
        -------
        ndarray
            Equivalent strain at the grid points, given as an :math:`m \times n` numpy array.

        See Also
        --------
        strain

        Notes
        -----
        Under plane stress conditions, :math:`\varepsilon_{33}` is not zero, but it is computed as

        .. math::
            \varepsilon_{33} = -\frac{\nu}{E}(\sigma_{11} + \sigma_{22})

        Since the stresses are not known from the DIC, one must settle with plane strain
        conditions where :math:`\varepsilon_{33} = 0`. This is what we follow in the :class:`DIC`
        class. We remark that in stereo-DIC multiple cameras are used, which allows the
        measurement of :math:`\varepsilon_{33}`.

        """
        eps_11, eps_12, eps_22 = strain_tensor.transpose((2, 0, 1))  # unpacking works in 1st dim
        if strain_measure == 'von Mises':
            return 2/3*np.sqrt(eps_11**2 + eps_22**2 + 3*eps_12**2 - eps_11*eps_22)

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

    @deprecated(deprecated_in='1.1.0', current_version=__version__, removed_in='1.3.0',
                details='This method will be modified to perform the plotting only. The '
                        'computation of the various strain measures will be done in the '
                        ':meth:`strain` method. The function that replaces this function '
                        'will keep the same name and is currently implemented outside this '
                        'class.')
    def plot_strain(self, component, minval=None, maxval=None, colorbar=True):
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
        minval, maxval : float, optional
            Set the color scaling for the image by fixing the values that map to the colormap
            color limits. If :code:`minval` is not provided, the limit is determined from the
            minimum value of the data. Similarly for :code:`maxval`.
        colorbar : bool, optional
            If True, a horizontal colorbar is placed below the strain map. The default is True.

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
        if colorbar:
            cbar = plt.colorbar(orientation='horizontal', format='%.2f', aspect=100,
                                label=label)
            cbar.ax.set_xlabel(label, fontsize=30)
        else:
            plt.gca().set_xlabel(label, fontsize=30)

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

    def project_onto_mesh(self, nodes, method='linear'):
        """Project the experimental displacement field onto a mesh.

        Parameters
        ----------
        nodes : ndarray
            2D numpy array with 2 columns, each row corresponding to a node of the mesh,
            and the two columns giving the Cartesian coordinates of the nodes.

        Returns
        -------
        u_nodes, v_nodes : ndarray
            1D numpy arrays (vectors), the components of the projected displacement field at the
            given :code:`nodes`.
        method : {'nearest', 'linear', 'cubic'}, optional
            Type of the interpolation. The default is linear.

        See Also
        --------
        project_onto_grid

        Examples
        --------
        Suppose we have measured field data, obtained by DIC. The experimental field values are
        obtained as an image. The pixel centers of this image form a grid, each grid point having
        an associated scalar value, the value of the measured field. In this example, we pretend
        that the measured field is given as a function. The DIC grid and the measured field values
        can be seen in the upper left and lower left figures, respectively.

        .. plot::

            >>> from grains.dic import DIC
            >>> from grains.geometry import TriMesh
            >>> x_grid, y_grid = np.mgrid[-1:1:100j, 1:2:50j]
            >>> exact_solution = lambda x, y: 1 - (x + y**2) * np.sign(x)
            >>> grid = DIC(exact_solution(x_grid, y_grid).T, exact_solution(x_grid,y_grid).T)

            The finite element solution is obtained at the nodes (upper right figure) of the mesh.
            Here, we assumed that the nodes are sampled from a uniformly random distribution on
            :math:`[-1,1] \\times [1, 2]`.

            >>> n_nodes = 100  # modify this to see how good the interpolated solution is
            >>> mesh = TriMesh(*TriMesh.sample_mesh(1, n_nodes))

            For the sake of this example, we also assume that the nodal values are "exact", i.e. they
            are the function values at the nodes. In reality, of course, this will not be the case,
            but this allows us to investigate the effect moving from the FE mesh to the DIC grid.

            >>> grid.set_transformation((-1, 2), 50)
            >>> fe_values = exact_solution(mesh.vertices[:, 0], mesh.vertices[:, 1])
            >>> interpolated = grid.project_onto_mesh(mesh.vertices, 'nearest')

            The FE solution available at the nodes are interpolated at the DIC grid points, as shown
            in the lower right figure. Interpolation with continuous functions cannot resolve well the
            discontinuity present in the "exact solution", i.e. in the measurement data. A discontinuous
            manufactured solution was intentionally prepared to illustrate this.

            >>> ax = []
            >>> ax.append(plt.subplot(221))
            >>> grid.plot_physicalgrid(ax[0])
            >>> plt.title('DIC grid')
            >>> ax.append(plt.subplot(222))
            >>> mesh.plot('k.', ax=ax[1], markersize=1)
            >>> plt.title('FE nodes')
            >>> ax.append(plt.subplot(223))
            >>> plt.imshow(exact_solution(x_grid, y_grid).T, extent=(-1, 1, 1, 2), origin='lower',
            ...            vmin=-4, vmax = 5)
            >>> plt.title('Exact solution on the DIC grid')
            >>> ax.append(plt.subplot(224))
            >>> interpolated[0][np.isnan(interpolated[0])] = 0
            >>> mesh.associate_field(interpolated[0])
            >>> mesh.plot_field(0, show_mesh=False, ax=ax[3])
            >>> plt.title('FE solution interpolated at the DIC grid')
            >>> plt.show()

        You are invited to modify this example to simulate a finer mesh by increasing
        :code:`n_nodes`. Also try different interpolation techniques.

        """
        # Coordinate values spanning the pixel grid
        image_height, image_width = np.shape(self.u)
        X = np.arange(0, image_height)
        Y = np.arange(0, image_width)
        # Coordinate values spanning the physical grid
        x = self.origin[0] + Y/self.scale
        y = self.origin[1] - np.flip(X)/self.scale
        # Interpolate the grid values at the nodes
        xx, yy = np.meshgrid(x, y)
        u_nodes = griddata((np.reshape(xx, (len(x) * len(y),)), np.reshape(yy, (len(x) * len(y),))),
                           np.reshape(self.u, (len(x) * len(y),)), nodes, method=method)
        v_nodes = griddata((np.reshape(xx, (len(x) * len(y),)), np.reshape(yy, (len(x) * len(y),))),
                           np.reshape(self.v, (len(x) * len(y),)), nodes, method=method)
        return u_nodes, v_nodes

    def project_onto_grid(self, nodes, nodal_values, method='linear'):
        """Project a numerically computed field onto the DIC grid.

        .. todo::
            Use the example code of this method to create the plotting method of this class.
            When done, rewrite the example with the methods.

        Parameters
        ----------
        nodes : ndarray
            2D numpy array with 2 columns, each row corresponding to a node of the mesh,
            and the two columns giving the Cartesian coordinates of the nodes.
        nodal_values : ndarray
            1D numpy array (vector), the numerically computed field at the nodes of the mesh.
        method : {'nearest', 'linear', 'cubic'}, optional
            Type of the interpolation. The default is linear.

        Returns
        -------
        ndarray
            2D numpy array, the projected field values at the physical DIC grid points.

        See Also
        --------
        project_onto_mesh,
        :func:`scipy.interpolate.griddata`

        Examples
        --------
        The DIC grid and the nodal finite element values (interpolated in-between, which is an
        errenous way to visualize discontinuous functions) can be seen in the upper left and lower
        left figures, respectively.

        .. plot::

            >>> from grains.dic import DIC
            >>> from grains.geometry import TriMesh
            >>> x_grid, y_grid = np.mgrid[-1:1:100j, 1:2:50j]
            >>> exact_solution = lambda x, y: 1 - (x + y**2) * np.sign(x)
            >>> grid = DIC(exact_solution(x_grid, y_grid).T, exact_solution(x_grid,y_grid).T)
            >>> grid.set_transformation((-1, 2), 50)

            The finite element solution is obtained at the nodes (upper right figure) of the mesh.
            Here, we assumed that the nodes are sampled from a uniformly random distribution on
            :math:`[-1,1] \\times [1, 2]`.

            >>> n_nodes = 100  # modify this to see how good the interpolated solution is

            For the sake of this example, we also assume that the nodal values are "exact", i.e.
            they admit the function values at the nodes. In reality, of course, this will
            not be the case, but this allows us to investigate the effect moving from the FE mesh
            to the DIC grid.

            >>> mesh = TriMesh(*TriMesh.sample_mesh(1, n_nodes))
            >>> fe_values = exact_solution(mesh.vertices[:, 0], mesh.vertices[:, 1])
            >>> mesh.associate_field(fe_values)
            >>> interpolated = grid.project_onto_grid(mesh.vertices, fe_values, method='linear')

            The FE field available at the nodes are interpolated at the DIC grid points, as shown
            in the lower right figure. Note that no extrapolation is performed, so values are not
            computed at the grid points lying outside the convex hull of the finite element nodes.

            >>> ax = []
            >>> ax.append(plt.subplot(221))
            >>> plt.plot(x_grid, y_grid, 'k.', ms=1)
            >>> plt.title('DIC grid')
            >>> ax.append(plt.subplot(222))
            >>> plt.plot(mesh.vertices[:, 0], mesh.vertices[:, 1], 'k.', ms=1)
            >>> plt.title('FE nodes')
            >>> ax.append(plt.subplot(223))
            >>> mesh.plot_field(0, ax=ax[-1])
            >>> plt.title('FE field')
            >>> ax.append(plt.subplot(224))
            >>> plt.imshow(interpolated, extent=(-1, 1, 1, 2), origin='lower', vmin=-4, vmax = 5)
            >>> plt.title('FE field interpolated at the DIC grid')
            >>> plt.tight_layout()
            >>> for a in ax:
            ...     a.set_aspect('equal', 'box')
            >>> plt.show()

        You are invited to modify this example to simulate a finer mesh by increasing
        :code:`n_nodes`. Also try different interpolation techniques.

        """
        # Coordinate values spanning the pixel grid
        image_height, image_width = np.shape(self.u)
        X = np.arange(0, image_height)
        Y = np.arange(0, image_width)
        # Coordinate values spanning the physical grid
        x = self.origin[0] + Y/self.scale
        y = self.origin[1] - np.flip(X)/self.scale
        # Interpolate the nodal values at the grid points
        xx, yy = np.meshgrid(x, y)
        return griddata(nodes, nodal_values, (xx, yy), method=method)

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
        ax : matplotlib.axes.Axes

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
            A subclass of :class:`grains.geometry.Mesh` to be plotted along with the physical
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
        :meth:`grains.geometry.TriMesh.plot`

        Examples
        --------
        Let us create a grid and a random mesh (the same as in the `Examples` section of the
        :meth:`project_onto_mesh`).

        .. plot::

            >>> from grains.dic import DIC
            >>> from grains.geometry import TriMesh
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


@deprecated(deprecated_in='1.1.0', current_version=__version__, removed_in='1.3.0',
            details='This function will be a static method of the :code:`DIC` class. '
                    'See also the deprecation warning of :meth:`DIC.plot_strain`.')
def plot_strain(strain, minval=None, maxval=None, colorbar=True, label=''):
    r"""Plots a scalar strain field.

    Parameters
    ----------
    strain : ndarray
        Scalar strain field sampled on an m-by-n grid, given as a 2D numpy array.
    minval, maxval : float, optional
        Set the color scaling for the image by fixing the values that map to the colormap
        color limits. If :code:`minval` is not provided, the limit is determined from the
        minimum value of the data. Similarly for :code:`maxval`.
    colorbar : bool, optional
        If True, a horizontal colorbar is placed below the strain map. The default is True.
    label : str, optional
        Label describing the strain field. LaTeX code is also accepted, e.g.
        `r'$\varepsilon_{yy}$'`. If not given, no text is displayed.

    Returns
    -------
    None

    See Also
    --------
    :meth:`DIC.strain`

    """
    plt.matshow(strain, interpolation='spline36', vmin=minval, vmax=maxval, alpha=1)
    if colorbar:
        cbar = plt.colorbar(orientation='horizontal', format='%.2f', aspect=100,
                            label=label)
        cbar.ax.set_xlabel(label, fontsize=30)
    else:
        plt.gca().set_xlabel(label, fontsize=30)
