# -*- coding: utf-8 -*-

"""
This module contains functions to provide geometrical descriptions of a label image.
Each labelled region is turned to a surface. This allows mesh generators to work on
the geometry instead of an image, thereby creating good quality meshes.

The following functions require PythonOCC version 0.18 to be installed to allow spline
manipulations:

- `fit_spline`
- `branches2spline`

"""

import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections
from skimage.segmentation import find_boundaries
from skimage.morphology import skeletonize
from skimage.color import label2rgb
from skan import Skeleton, summarize

from grains import HAS_OCCT

if HAS_OCCT:
    from OCC.gp import gp_Pnt
    from OCC.TColgp import TColgp_Array1OfPnt
    from OCC.GeomAPI import GeomAPI_PointsToBSpline
else:
    warnings.warn('PythonOCC is not available. Some functions cannot be used.', ImportWarning)

from .utils import toggle, index_list, non_unique


def build_skeleton(label_image, connectivity=1, detect_boundaries=True):
    """Builds skeleton connectivity of a label image.

    A single-pixel wide network is created, separating the labelled image regions. The resulting
    network contains information about how the regions are connected.

    Parameters
    ----------
    label_image : 2D ndarray with signed integer entries
        Label image, representing a segmented image.
    connectivity : {1,2}, optional
        A connectivity of 1 (default) means pixels sharing an edge will be considered neighbors.
        A connectivity of 2 means pixels sharing a corner will be considered neighbors.
    detect_boundaries : bool, optional
        When True, the image boundaries will be treated as part of the skeleton. This allows
        identifying boundary regions in the `skeleton2regions` function. The default is True.

    Returns
    -------
    skeleton_network : Skeleton
        Geometrical and topological information about the skeleton network of the input image.

    See Also
    --------
    skan.Skeleton

    """
    # 2D image, given as a numpy array is expected
    if type(label_image) != np.ndarray:
        raise Exception('Label image must be a numpy array (ndarray).')
    image_size = np.shape(label_image)
    if len(image_size) != 2:
        raise Exception('A 2D array is expected.')
    if not issubclass(label_image.dtype.type, np.signedinteger):
        raise Exception('Matrix entries must be positive integers.')

    # Surround the image with an outer region
    if detect_boundaries:
        label_image = np.pad(label_image, pad_width=1, mode='constant', constant_values=-1)

    # Find the boundaries of the label image and then extract its skeleton
    boundaries = find_boundaries(label_image, connectivity=connectivity)
    skeleton = skeletonize(boundaries)

    # Build the skeleton network using `skan`
    skeleton_network = Skeleton(skeleton, source_image=label_image, keep_images=True)
    return skeleton_network


def skeleton2regions(skeleton_network, look_around=2):
    """Determines the regions bounded by a skeleton network.

    This function can be perceived as an intermediate step between a skeleton network and
    completely geometrical representation of the regions. That is, it keeps the key topological
    information required to create a fully geometrical description, but it also contains the
    coordinates of the region boundaries. The outputs of this function can be used to build
    different region representations.

    Parameters
    ----------
    skeleton_network : Skeleton
        Geometrical and topological information about the skeleton network of a label image.
    look_around : int, optional
        A junction is considered part of a region if it is at most `look_around` pixel far from it.
        The default is 2. For further details, see the Notes below.

    Returns
    -------
    region_branches : dict
        For each region it contains the branch indices that bound that region.
    branch_coordinates : list
        Coordinates of the points on each branch.
    branch_regions : dict
        For each region it contains the neighboring regions.
        This auxiliary data is not essential as it can be restored from `region_branches`.
        However, it is computed as temporary data needed for `region_branches`.

    See Also
    --------
    build_skeleton
    overlay_regions

    Notes
    -----
    Although the algorithms were created to require minimum user intervention, some parameters
    must be fine-tuned so as to achieve an optimal result in identifying the regions. Visualization
    plays an important role in it. Full automation is either not possible or would require a huge
    computational cost. The shortcomings of the algorithms in this function are the following:

    - It is assumed that only branches that connect two junctions form the boundary of a region.
      This rule ensures that end points are not taken into account. However, this assumption also
      rules out the identification of region being contained in another region as the embedded
      region would be described by an isolated cycle.
    - The recognition of which branches form a region is based on the premise that a junction
      belongs to a region if its n-pixel neighbourhood contains a pixel from that region.
      Ideally, n=1 would be used, meaning that the single-pixel width skeleton is located at most
      1 pixel afar from the regions it lies among. This is true but the junctions of the skeleton
      can be farther than 1 pixel from a region. Hence, `n` has to be a parameter of our model.
      Increasing `n` helps in including junctions (and hence the connecting branches to it) to
      regions, which actually belong there. On the other hand, if `n` is too large, junctions
      that do not belong to the region are also included. Currently, we recommend trying
      different parameters `n`, plot the reconstructed regions over the label image using the
      `overlay_regions` function, and see how good the result is. As a heuristic, start with `n=2`.
    - There are configurations in which two junctions are part of a region but those two
      junctions are connected by more than one branches (typically two). The question is: which
      branch to choose as a boundary part of the region? The answer is: the one that is entirely
      inside the region. Testing it is probably not easy or is costly, therefore, we rely on a
      heuristic argument: the branch with the shortest length is chosen.

    """
    if not isinstance(skeleton_network, Skeleton):
        raise Exception('Skeleton object is expected.')
    # Extract branch-junction connectivities and the coordinates of the junctions
    S = skeleton_network
    skeleton_data = summarize(S)
    mask = skeleton_data['branch-type'] == 2  # only junction-to-junction connections create regions
    endpoints_src = skeleton_data['node-id-src'][mask].to_numpy()
    endpoints_dst = skeleton_data['node-id-dst'][mask].to_numpy()
    image_size = np.shape(S.source_image)
    branch_junctions = np.transpose(np.vstack((endpoints_src, endpoints_dst)))
    junctions = np.unique([endpoints_src, endpoints_dst])
    junction_coordinates = S.coordinates[junctions, :]

    # Find which regions are incident to a junction
    junction_regions = {key: None for key in junctions}
    region_junctions = {}
    # TODO: Simplify the for-loop by using e.g. enumerate or
    #  https://discuss.codecademy.com/t/loop-two-variables-simultaneously-in-python-3/261808/2
    for i in range(len(junctions)):
        # Snap junction to the nearest image coordinate
        junction_coord = np.round(junction_coordinates[i, :]).astype(np.uint32)
        # Look-around for the neighboring pixels (be careful on the image boundaries)
        neighbor_idx = np.s_[
                       max(junction_coord[0] - look_around, 0):
                       min(junction_coord[0] + look_around+1, image_size[0]),
                       max(junction_coord[1] - look_around, 0):
                       min(junction_coord[1] + look_around+1, image_size[1])]
        neighbors = S.source_image[neighbor_idx]
        neighboring_regions = np.unique(neighbors)
        # Save junction-region and the region-junction connectivities
        # TODO: perhaps no need for the region-junction connectivities
        junction_regions[junctions[i]] = neighboring_regions
        for region in neighboring_regions:
            if region not in region_junctions:
                region_junctions[region] = [junctions[i]]
            else:
                region_junctions[region].append(junctions[i])

    # Determine which regions neighbor a branch
    branch_regions = {}
    for i, branch in enumerate(branch_junctions):
        neighboring_regions = np.intersect1d(junction_regions[branch[0]],
                                             junction_regions[branch[1]])
        branch_regions[i] = neighboring_regions

    # For each region, find the branches that bound it
    region_branches = {}
    for branch, regions in branch_regions.items():
        for region in regions:
            if region not in region_branches:
                region_branches[region] = [branch]
            else:
                region_branches[region].append(branch)

    # More than one branch can connect two junctions. In that case, leave only one.
    branch_lengths = S.path_lengths()
    for i in region_branches.keys():
        branches = region_branches[i]
        junctions = region_junctions[i]
        if len(branches) > len(junctions):
            # Branches that connect the same two junctions (i.e. multiple edges in the graph)
            # branch_junctions[region_branches[i]]
            _, id = non_unique(branch_junctions[branches], 0)
            # Only the shortest branch bounds this region
            edges_to_remove = []
            for multiple_edges in id:
                global_multiple_edges = index_list(branches, multiple_edges)
                lengths = branch_lengths[global_multiple_edges]
                idx_min_length = np.argmin(lengths)
                other_edges = np.setdiff1d(range(len(multiple_edges)), idx_min_length)
                edges_to_remove.append(multiple_edges[other_edges])
                # Remove the current region among the regions that connect to the removed edges
                for branch in index_list(global_multiple_edges, other_edges):
                    regions = branch_regions[branch]
                    kept_regions = regions[i != regions]
                    branch_regions[branch] = kept_regions
            correct_branches = np.delete(branches, edges_to_remove).tolist()
            region_branches[i] = correct_branches

    branch_coordinates = [S.path_coordinates(i) for i in range(S.n_paths) if mask[i]]
    # Return outputs
    return region_branches, branch_coordinates, branch_regions


def polygon_orientation(polygon):
    """Determines whether a polygon is oriented clockwise or counterclockwise.

    Parameters
    ----------
    polygon : list
        Each element of the list denotes a vertex of the polygon and in turn is another list of two
        elements: the x and y coordinates of a vertex.

    Returns
    -------
    orientation : {'cw', 'ccw'}
        'cw': clockwise, 'ccw': counterclockwise orientation

    Notes
    -----
    The formula to determine the orientation is from https://stackoverflow.com/a/1165943/4892892.
    For simple polygons (polygons that admit a well-defined interior), a faster algorithm exits, see
    https://en.wikipedia.org/wiki/Curve_orientation#Orientation_of_a_simple_polygon.

    Examples
    --------
    >>> polygon = [[5, 0], [6, 4], [4, 5], [1, 5], [1, 0]]
    >>> polygon_orientation(polygon)
    'ccw'

    """
    n_vertex = len(polygon)
    edge_sum = 0
    for idx, vertex in enumerate(polygon):
        next_vertex = polygon[(idx + 1) % n_vertex]  # allow indexing past the last vertex
        edge_sum += (next_vertex[0] - vertex[0]) * (next_vertex[1] + vertex[1])
    if edge_sum > 0:
        orientation = 'cw'
    else:
        orientation = 'ccw'
    return orientation


def branches2boundary(branches, orientation='ccw'):
    """Joins connecting branches so that they form the boundary of a surface.

    This function assumes that you already know that a set of branches form the boundary of a
    surface, but you want to know the ordering. Clockwise or counterclockwise orientation is
    supported. A branch is given by its two end points but it can also contain intermediate points
    in between, as in the general case. The points on the branches are assumed to be ordered. If
    certain branches are not used in forming the boundary of a surface, they are excluded from the
    output lists.

    Parameters
    ----------
    branches : list
        Each element of the list gives N>=2 points on the branch, ordered from one end point to the
        other. If N=2, the two end points are meant. The points are provided as an Nx2 ndarray,
        the first column giving the x, the second column giving the y coordinates of the points.
    orientation : {'cw', 'ccw'}, optional
        Clockwise ('cw') or counterclockwise ('ccw') orientation of the surface boundary.
        The default is 'ccw'.

    Returns
    -------
    order : list
        Order of the branches so that they form the boundary of a surface.
    redirected : list
        A list of bool with True value if the orientation of the corresponding branch had to be
        swapped to form the surface boundary.
    polygon : ndarray
        The polygon formed by connecting the points of the branches along the boundary. It is given
        given as an Mx2 ndarray, where M is the number of unique points on the boundary (i.e. only
        one end point is kept for two connecting branches).
        This auxiliary data is not essential as it can be restored from the set of branches,
        and their ordering. However, it is computed as temporary data needed for determining the
        orientation of the boundary.

    Examples
    --------
    >>> import numpy as np
    >>> branches = [np.array([[1, 1], [1.5, 2], [2, 3]]), np.array([[1, 1], [-1, 2]]),
    ... np.array([[1.5, -3], [2, 3]]), np.array([[1.5, -3], [-1, 2]])]
    >>> order, redirected, polygon = branches2boundary(branches, orientation='cw')
    >>> order
    [0, 2, 3, 1]
    >>> redirected
    [False, True, True, False]
    >>> polygon
    array([[ 1. ,  1. ],
           [ 1.5,  2. ],
           [ 2. ,  3. ],
           [ 1.5, -3. ],
           [-1. ,  2. ]])

    """
    # The path is the collection of consecutively added branches. When there are no more branches
    # to add, the path becomes the boundary of a surface.
    # First, we identify the branches that form the boundary.
    n_branch = len(branches)
    redirected = [False for i in range(n_branch)]
    # Start with an arbitrary branch, the first one
    last_branch = 0  # index of the lastly added branch to the path
    order = [0]
    first_vertex = branches[last_branch][0, :]
    last_vertex = branches[last_branch][-1, :]
    # Visit all vertices of the would-be surface and find the chain of connecting branches
    for vertex in range(n_branch - 1):
        if np.allclose(last_vertex, first_vertex):  # one complete cycle is finished
            break
        # Search for branches connecting to the last vertex of the path
        for branch_index, branch in enumerate(branches):
            if branch_index == last_branch:  # exclude the last branch of the path
                continue
            # Check if the first or second end point of the branch connects to the last vertex
            first_endpoint = branch[0, :]
            second_endpoint = branch[-1, :]
            vertex_connects_to_first_endpoint = np.allclose(first_endpoint, last_vertex)
            vertex_connects_to_second_endpoint = np.allclose(second_endpoint, last_vertex)
            if vertex_connects_to_first_endpoint or vertex_connects_to_second_endpoint:
                order.append(branch_index)
                last_branch = branch_index
            if vertex_connects_to_first_endpoint:
                last_vertex = second_endpoint
                break
            elif vertex_connects_to_second_endpoint:
                last_vertex = first_endpoint
                redirected[branch_index] = True  # change the orientation of the connecting branch
                break
            # TODO: add checks against edge cases

    # Polygonize the boundary branches so as to determine the orientation in the next step
    polygon = []
    for branch_index in order:
        branch = branches[branch_index].copy()
        if redirected[branch_index]:
            branch = np.flipud(branch)
        polygon.append(branch[0:-1])
    polygon = np.vstack(polygon)

    # Handle the orientation of the polygon
    if polygon_orientation(polygon) != orientation:
        order.reverse()
        redirected = toggle(redirected)
        polygon = np.flipud(polygon)

    return order, redirected, polygon


def polygonize(label_image, connectivity=1, detect_boundaries=True, look_around=2,
               orientation='ccw', close=False):
    """Polygon representation of a label image.

    Parameters
    ----------
    label_image : 2D ndarray with signed integer entries
        Label image, representing a segmented image.
    connectivity : {1,2}, optional
        A connectivity of 1 (default) means pixels sharing an edge will be considered neighbors.
        A connectivity of 2 means pixels sharing a corner will be considered neighbors.
    detect_boundaries : bool, optional
        When True, the image boundaries will be treated as part of the skeleton. This allows
        identifying boundary regions in the `skeleton2regions` function. The default is True.
    look_around : int, optional
        A junction is considered part of a region if it is at most `look_around` pixel far from it.
        The default is 2. For further details, see the Notes below.
    orientation : {'cw', 'ccw'}, optional
        Clockwise ('cw') or counterclockwise ('ccw') orientation of the polygons.
        The default is 'ccw'.
    close : bool, optional
        When True, one vertex in the polygons is repeated to indicate that the polygons are
        indeed closed. The default is False.

    Returns
    -------
    polygons : dict
        The keys in the dictionary correspond to the labels of the input image, while the values
        are ndarray objects with two columns, the x and y coordinates of the polygons.

    See Also
    --------
    skeleton2regions
    segments2polygon

    Examples
    --------
    >>> test_image = np.array([
    ...   [1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ...   [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]],
    ...  dtype=np.int8)
    >>> polygons = polygonize(test_image, connectivity=1)

    """
    polygons = {}
    # Build the skeleton network from the label image
    S = build_skeleton(label_image, connectivity, detect_boundaries)
    # Identify the regions from the skeleton
    region_branches, branch_coordinates, _ = skeleton2regions(S, look_around)
    # Represent each region as a polygon
    for region, branches in region_branches.items():
        if region == -1:  # artificial outer region
            continue
        points_on_boundary = index_list(branch_coordinates, branches)
        _, _, poly = branches2boundary(points_on_boundary, orientation)
        if close:
            poly = np.vstack((poly, poly[0, :]))
        polygons[region] = poly
    # overlay_skeleton_networkx(S.graph, S.coordinates, image=label_image)
    plt.show()
    return polygons


def branches2splines(branches, degree_min=3, degree_max=8, continuity='C2', tol=1e-3):
    """Lays splines on branches.

    Parameters
    ----------
    branches : list
        Each element of the list gives N>=2 points on the branch, ordered from one end point to the
        other. If N=2, the two end points are meant. The points are provided as an Nx2 ndarray,
        the first column giving the x, the second column giving the y coordinates of the points.

    Other Parameters
    ----------------
    degree_min, degree_max, continuity, tol
        See the :py:func:`fit_spline` function.

    Returns
    -------
    list
        Each member of the list is a `Geom_BSplineCurve` object, the B-spline approximation of
        the corresponding input branch.
        For details on the resulting spline, see the OpenCASCADE documentation:
        https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_geom___b_spline_curve.html

    See Also
    --------
    fit_spline

    """
    return [fit_spline(branch, degree_min, degree_max, continuity, tol) for branch in branches]


def fit_spline(points, degree_min=3, degree_max=8, continuity='C2', tol=1e-3):
    """Approximates a set of points on the plane with a B-spline.

    Parameters
    ----------
    points : ndarray
        2D ndarray with N>=2 rows and 2 columns, representing the points on the plane, ordered
        from one end point to the other. The first column gives the x, the second column gives the
        the y coordinates of the points.
    degree_min : int, optional
        Minimum degree of the spline. The default is 3.
    degree_max : int, optional
        Maximum degree of the spline. The default is 8.
    continuity : {'C0', 'G1', 'C1', 'G2', 'C2', 'C3', 'CN'}, optional
        The continuity of the spline will be at least `continuity`. The default is 'C2'.
        For their meanings, consult with
        https://www.opencascade.com/doc/occt-7.4.0/refman/html/_geom_abs___shape_8hxx.html
    tol : float, optional
        The distance from the points to the spline will be lower than `tol`. The default is 1e-3.

    Returns
    -------
    spline : Geom_BSplineCurve
        For details on the resulting spline, see the OpenCASCADE documentation:
        https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_geom___b_spline_curve.html

    Raises
    ------
    ValueError
        If the minimum degree of the B-spline to be constructed is greater than its maximum degree.

    """
    if degree_min > degree_max:
        raise ValueError('Minimum degree must not be lower than the maximum degree.')
    # Create points from the input array that OCCT understands
    n_point = np.size(points, 0)
    pts = TColgp_Array1OfPnt(0, n_point-1)
    for i in range(n_point):
        pts.SetValue(i, gp_Pnt(float(points[i][0]), float(points[i][1]), 0))
    # Construct the spline
    continuity = _spline_continuity_enum(continuity)
    spline = GeomAPI_PointsToBSpline(pts, degree_min, degree_max, continuity, tol).Curve()
    return spline


def plot_polygon(vertices, **kwargs):
    """Plots a polygon.

    Parameters
    ----------
    vertices : ndarray
        2D ndarray of size Nx2, with each row designating a vertex and the two columns
        giving the x and y coordinates of the vertices, respectively.
    **kwargs : Line2D properties, optional
        Keyword arguments accepted by matplotlib.pyplot.plot

    Returns
    -------
    None

    See Also
    --------
    matplotlib.pyplot.plot

    Examples
    --------
    >>> plot_polygon(np.array([[1, 1], [2, 3], [1.5, -3], [-1, 2]]), marker='o');  plt.show()

    """
    # Close the polygon (repeat one vertex) if not yet closed
    first_vertex = vertices[0, :]
    last_vertex = vertices[-1, :]
    closed = np.allclose(first_vertex, last_vertex)
    if not closed:
        vertices = np.vstack((vertices, first_vertex))
    plt.plot(vertices[:, 0], vertices[:, 1], **kwargs)


def overlay_regions(label_image, polygons, axes=None):
    """Plots a label image, and overlays polygonal regions over it.

    Parameters
    ----------
    label_image : 2D ndarray with signed integer entries
        Label image, representing a segmented image.
    polygons : dict
        The keys in the dictionary correspond to the labels of the input image, while the values
        are ndarray objects with two columns, the x and y coordinates of the polygons.
        This format is respected by the output of the `polygonize` function.
    axes : matplotlib.axes.Axes, optional
        An Axes object on which to draw. If None, a new one is created.

    Returns
    -------
    axes : matplotlib.axes.Axes
        The Axes object on which the plot is drawn.

    See Also
    --------
    polygonize
    matplotlib.collections.LineCollection

    """
    # TODO: support an option to give the number of identified regions as a title
    if axes is None:
        _, axes = plt.subplots()
    # Extract the polygons from the dictionary and convert them to a list to ease the plotting.
    # At the same time, swap the x and y coordinates so that the polygons are expressed in the
    # same coordinate system as the label image.
    polygons = [polygons[i][:, ::-1] for i in polygons.keys()]
    # Plot the polygons efficiently
    axes.add_collection(collections.LineCollection(polygons, colors='black'))
    # Plot the label image, with each color corresponding to a different region
    random_colors = np.random.random((len(polygons), 3))
    plt.imshow(label2rgb(label_image, colors=random_colors))
    # Axis settings
    axes.set_aspect('equal')
    axes.set_axis_off()
    plt.tight_layout()
    return axes


def _spline_continuity_enum(continuity):
    """Enumeration value corresponding to the continuity of a B-spline.

    Parameters
    ----------
    continuity : {'C0', 'G1', 'C1', 'G2', 'C2', 'C3', 'CN'}
        Continuity of the B-spline.

    Returns
    -------
    enum : int
        Integer value corresponding to continuity key in the OpenCASCADE API:
        https://www.opencascade.com/doc/occt-7.4.0/refman/html/_geom_abs___shape_8hxx.html

    """
    if continuity == 'C0':
        enum = 0
    elif continuity == 'G1':
        enum = 1
    elif continuity == 'C1':
        enum = 2
    elif continuity == 'G2':
        enum = 3
    elif continuity == 'C2':
        enum = 4
    elif continuity == 'C3':
        enum = 5
    elif continuity == 'CN':
        enum = 6
    else:
        raise ValueError('Choose one of {"C0", "G1", "C1", "G2", "C2", "C3", "CN"}.')
    return enum


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
