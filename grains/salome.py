# -*- coding: utf-8 -*-
"""
The documentation generated from this file is available on
https://cristalx.readthedocs.io/en/latest/api.html#module-grains.salome.

The aim of this :ref:`module <salome module>` is to manage geometry and mesh operations on
two-dimensional microstructures that tessellate a domain. This has the important consequence that
the domain is assumed to consist of non-overlapping shapes that cover it. After generating a
`conforming` mesh, each element belongs to a single group and no element exists outside the
group. More precisely, elements that do not belong to any group are not taken into account by the
:class:`Mesh` class. They can still be accessed by the functions of Salome, but the use of such
lower-level abstractions is against the philosophy of encapsulation this module provides.

The non-goals of this module include everything not related to the tessellation nature of 2D
microstructures. The emphasis is on readability and not on speed.

This file can be used either as a module or as a script.

.. _salome module:

Using as a module
-----------------
Developers, implementing new features, will use `salome.py` as a Python module. Although this
module is part of the :mod:`.grains` package in the `CristalX project
<https://github.com/CsatiZoltan/CristalX>`_, it does not rely on the other modules of
:mod:`grains`. This ensures that it can be used standalone and :ref:`run as a script <salome
script>` in the Salome environment. It only uses language constructs available in Python 3.6,
and external packages shipped with Salome 9.4.0 onwards. To enable debugging, code completion
and other useful development methods, consult with the documentation on ...
Most classes contain a protected variable that holds the underlying Salome object. Unless you
debug, it is not necessary to directly deal with Salome objects programmatically.

.. _salome script:

Using as a script
-----------------
When `salome.py` is run as a script, the contents in the :code:`if __name__ == "__main__":` block
is executed. Edit it to suit your needs. Similarly to the case when :ref:`used as a module
<salome module>`, the script can only be run from Salome's own Python interpreter: either from
the shell or from the GUI. To run it from the shell (including the GUI's built-in Python command
prompt), type

.. code-block:: python

    exec(open("<path_to_CristalX>/grains/salome.py", "rb").read())

You can also execute the script from the GUI by clicking on :menuselection:`File --> Load Script`...


Classes
-------
.. autosummary::
    :nosignatures:
    :toctree: classes/

    Geometry
    Face
    Edge
    Interface
    Mesh
    FaceMesh
    InterfaceMesh
    CohesiveZone
    GUI

"""
from enum import Enum
from functools import reduce

import numpy as np
import matplotlib.path as mpltPath

try:
    import salome
    # kernel-related modules
    from salome.kernel.studyedit import getStudyEditor
    salome.salome_init()
    import salome_notebook

    # geometry-related packages and modules
    import GEOM
    from salome.geom import geomBuilder
    import SALOMEDS

    # mesh-related packages and modules
    import SMESH
    from salome.smesh import smeshBuilder
except:
    raise ImportError('This module requires the salome module from Salome.')


# Initializations
geompy = geomBuilder.New()
smesh = smeshBuilder.New()
studyEditor = getStudyEditor()
notebook = salome_notebook.NoteBook()


# Classes and functions of this module
# TODO: Put the `show` method of Edge and Face into a common superclass.
class Face:
    """Closed part of a plane.

    A :class:`Face` knows about the :class:`Edge` objects that bound it.

    Parameters
    ----------
    face : GEOM_Object of shape type 'FACE'
        The main Salome object wrapped by this class.
    name : str
        Name of the face.

    See Also
    --------
    `GEOM_Object <https://docs.salome-platform.org/latest/gui/GEOM/geompy_doc
    /interfaceGEOM_1_1GEOM__Object.html>`_,
    `Shape type <https://docs.salome-platform.org/latest/gui/GEOM/geompy_doc/namespaceGEOM.html
    #a82a00e336c65dad4cc04b65563b26eb5>`_

    """
    def __init__(self, face, name):
        self._face = face
        self.name = name
        self.edges = []


class Edge:
    """A shape corresponding to a curve, and bounded by a vertex at each extremity.

    An :class:`Edge` knows about the :class:`Face` it is part of, and the faces neighboring it.

    Parameters
    ----------
    edge : GEOM_Object of shape type 'EDGE'
        The main Salome object wrapped by this class.
    name : str
        Name of the edge.

    See Also
    --------
    `GEOM_Object <https://docs.salome-platform.org/latest/gui/GEOM/geompy_doc
    /interfaceGEOM_1_1GEOM__Object.html>`_,
    `Shape type <https://docs.salome-platform.org/latest/gui/GEOM/geompy_doc/namespaceGEOM.html
    #a82a00e336c65dad4cc04b65563b26eb5>`_

    """
    def __init__(self, edge, name):
        self._edge = edge
        self.name = name
        self.incident_face = None
        self.neighboring_faces = None

    def length(self):
        """Length of the edge.

        Returns
        -------
        float
            Length of the edge.

        See Also
        --------
        `geomBuilder.BasicProperties <https://docs.salome-platform.org/latest/gui/GEOM/geompy_doc
        /group__l2__measure.html#ga6d60abd33031977af29b8036d001bf8b>`_

        """
        return geompy.BasicProperties(self._edge)[0]


class Interface:
    """An edge between two faces.

    Similar to an :class:`Edge`, but two neighboring :class:`Face` objects share a common
    :class:`Interface`. An :class:`Interface` knows about the two faces that it separates.

    Parameters
    ----------
    edge : GEOM_Object of shape type 'EDGE'
        The main Salome object wrapped by this class.
    name : str
        Name of the interface.
    neighboring_faces : list of Face
        The two neighboring faces.

    See Also
    --------
    Edge, Face

    """
    def __init__(self, edge, name, neighboring_faces):
        self._edge = edge
        self.name = name
        self.neighboring_faces = neighboring_faces

    def length(self):
        """Length of the interface.

        Returns
        -------
        float
            Length of the interface.

        See Also
        --------
        `geomBuilder.BasicProperties <https://docs.salome-platform.org/latest/gui/GEOM/geompy_doc
        /group__l2__measure.html#ga6d60abd33031977af29b8036d001bf8b>`_

        """
        return geompy.BasicProperties(self._edge)[0]


class Geometry:
    """Represents the geometrical entities of a two-dimensional tessellated domain.

    A :class:`Geometry` object knows about the faces that tessellate the domain, and about the
    edges and interfaces that separate the faces.

    Parameters
    ----------
    name : str, optional
        Name of the geometry.

    See Also
    --------
    Face, Edge, Interface

    """
    def __init__(self, name='microstructure'):
        self.name = name
        self.geometry = None
        self.faces = []
        self.edges = []
        self.interfaces = []

    def load(self, step_file):
        """Loads the geometry from a STEP file.

        Parameters
        ----------
        step_file : str
            The STEP file containing the geometry.

        Returns
        -------
        None.

        """
        self.geometry = geompy.ImportSTEP(step_file, True, True)
        geompy.addToStudy(self.geometry, self.name)

    def extract_faces(self):
        """Decomposes the geometry into faces.

        Returns
        -------
        None.

        See Also
        --------
        extract_edges

        """
        faces = geompy.ExtractShapes(self.geometry, geompy.ShapeType["FACE"], False)
        for face_num, face in enumerate(faces, 1):
            face_name = 'Face_'+str(face_num)
            face_object = Face(face, name=face_name)
            self.faces.append(face_object)

    def extract_edges(self):
        """Decomposes each face into edges.

        This method must be called after the :meth:`~Geometry.extract_faces` method, otherwise it
        has no effect.

        Returns
        -------
        None.

        See Also
        --------
        extract_faces

        """
        i_edge = 0  # edge counter
        for face in self.faces:
            face_name = face.name
            geompy.addToStudyInFather(self.geometry, face._face, face_name)
            edges = geompy.ExtractShapes(face._face, geompy.ShapeType["EDGE"], False)
            for edge in edges:
                i_edge += 1
                edge_name = 'Edge_'+str(i_edge)
                edge_object = Edge(edge, edge_name)
                self.edges.append(edge_object)
                face.edges.append(edge_object)  # face-edge mapping
                edge_object.incident_face = face  # edge-face mapping
                geompy.addToStudyInFather(face._face, edge, edge_name)

    def create_interfaces(self):
        """Constructs unique interfaces that separate the faces.

        Based on the edges (obtained by exploding the mesh), interfaces are created.
        Interfaces are unique separators of two neighboring faces. In other words,

        - if the edge is a boundary edge, no interface is created,
        - two neighboring faces have two overlapping edges, of which one is defined
          to be an interface.

        It is assumed that the edges and faces of the geometry have already been obtained.

        Returns
        -------
        None.

        See Also
        --------
        extract_faces, extract_edges

        """
        # Neighboring faces to each edge
        for i_int, overlapping_edges in enumerate(self._find_overlapping_edges(), 1):
            # The interface is defined on that Salome Edge object that will hold a mesh
            edge_with_a_better_name = self._has_smaller_ID(overlapping_edges)
            interface_name = 'Interface_' + str(i_int)
            interface_object = Interface(edge_with_a_better_name._edge, name=interface_name,
                                         neighboring_faces=[overlapping_edges[0].incident_face,
                                                            overlapping_edges[1].incident_face])
            self.interfaces.append(interface_object)
            # Display the newly created interface and its neighboring faces in Salome's GUI tree
            geompy.addToStudyInFather(self.geometry, interface_object._edge, interface_name)
            for face in interface_object.neighboring_faces:
                geompy.addToStudyInFather(interface_object._edge, face._face, face.name)

    def _find_overlapping_edges(self):
        """Finds edges that are on top of each other.

        Overlapping edges have the same length. Although its converse is not true in general,
        we will assume so. This part of the algorithm (i.e. deciding the overlapping edges) can
        later be refined.

        Returns
        -------
        overlapping_edges : list
            Each member of the list contains a list of (supposedly) two Edge objects.

        See Also
        --------
        create_interfaces, Edge.length

        """
        # Length of the edges
        edge_lengths = {}
        for face in self.faces:
            for edge in face.edges:
                edge_lengths[edge] = edge.length()
        # Invert the edge - edge length mapping (https://stackoverflow.com/a/20672375/4892892)
        rev_edge_lengths = {}
        for key, value in edge_lengths.items():
            rev_edge_lengths.setdefault(value, set()).add(key)
        # Edges on the boundary cannot overlap as they are not between two faces
        overlapping_edges = list(rev_edge_lengths.values())
        overlapping_edges = [list(edges) for edges in overlapping_edges if len(edges) == 2]
        return overlapping_edges

    @staticmethod
    def _has_smaller_ID(edges):
        """Selects the edge with a smaller ID.

        When generating mesh with the NETGEN plugin of Salome, if several edges overlap, only the
        edge with the smallest ID holds a mesh. The purpose of this function is to find the edge
        with the smaller ID out of two overlapping edges.

        Parameters
        ----------
        edges : list of Edge
            List of two Edge objects.

        Returns
        -------
        Edge
            Either the first or the second element of the input list, depending on which of them
            has a smaller ID.

        See Also
        --------
        create_interfaces

        Notes
        -----
        For a detailed discussion on this highly important issue, see the
        `corresponding forum thread <https://www.salome-platform.org/forum/forum_11/79066443>`_.

        """
        # ID of Salome's Face object the edge is incident to
        face1_id = salome.myStudy.FindObject(edges[0].incident_face.name).GetID()
        face1_registered_position = int(face1_id[face1_id.rfind(':') + 1:])
        face2_id = salome.myStudy.FindObject(edges[1].incident_face.name).GetID()
        face2_registered_position = int(face2_id[face2_id.rfind(':') + 1:])
        if face1_registered_position < face2_registered_position:
            return edges[0]
        else:
            return edges[1]


class Mesh:
    """Performs mesh manipulations on a tessellated geometry.

    Parameters
    ----------
    geometry : Geometry
        Geometry object on which the mesh exists.
    name : str, optional
        Name of the mesh.

    See Also
    --------
    Geometry, FaceMesh, InterfaceMesh

    """

    class ElementType(Enum):
        """Subset of the element types recognized by Salome.

        This enumeration is for convenience. Only those elements of Salome are considered that
        are relevant for the :class:`Mesh` class.

        See Also
        --------
        SMESH.ElementType

        """
        ALL = SMESH.ALL
        NODE = SMESH.NODE
        EDGE = SMESH.EDGE
        FACE = SMESH.FACE

    def __init__(self, geometry, name='Mesh'):
        self._geometry = geometry
        self._mesh = smesh.Mesh(geometry.geometry, name)
        self.face_meshes = []
        self.interface_meshes = []

    def generate(self):
        """Generates a mesh on the geometry.

        .. todo::
            Do not hardcode values and explain the need for consistent orientation.

        Returns
        -------

        """
        NETGEN_1D_2D = self._mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
        NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
        NETGEN_2D_Parameters_1.SetMaxSize(246.598)
        NETGEN_2D_Parameters_1.SetMinSize(0.000230144)
        NETGEN_2D_Parameters_1.SetSecondOrder(0)
        NETGEN_2D_Parameters_1.SetOptimize(1)
        NETGEN_2D_Parameters_1.SetFineness(2)
        NETGEN_2D_Parameters_1.SetChordalError(-1)
        NETGEN_2D_Parameters_1.SetChordalErrorEnabled(0)
        NETGEN_2D_Parameters_1.SetUseSurfaceCurvature(1)
        NETGEN_2D_Parameters_1.SetFuseEdges(1)
        NETGEN_2D_Parameters_1.SetWorstElemMeasure(0)
        NETGEN_2D_Parameters_1.SetQuadAllowed(0)
        NETGEN_2D_Parameters_1.SetCheckChartBoundary(0)
        isDone = self._mesh.Compute()  # TODO: do sth in case of error
        number_of_faces_reoriented = self._mesh.Reorient2D(self._mesh, [0, 0, 1], [0, 0, 0])
        if not self._mesh.IsCoherentOrientation2D():
            raise Exception('This is embarrassing!')

    def obtain_face_meshes(self):
        """Retrieves the elements of the mesh on each face.

        Returns
        -------
        None.

        """
        for face in self._geometry.faces:
            mesh_on_face = self._mesh.GroupOnGeom(face._face, name=face.name, typ=SMESH.FACE)
            face_mesh_object = FaceMesh(mesh_on_face, name=face.name, on_face=face)
            self.face_meshes.append(face_mesh_object)

    def obtain_interface_meshes(self):
        """Obtains the 1D interfacial mesh for each interface.

        Returns
        -------
        None.

        """
        for i, interface in enumerate(self._geometry.interfaces, 1):
            mesh_on_edge = self._mesh.GroupOnGeom(interface._edge, typ=SMESH.EDGE)
            if mesh_on_edge.IsEmpty():
                raise Exception('')
            mesh_on_edge.SetName(interface.name)
            interface_mesh_object = InterfaceMesh(mesh_on_edge, name=interface.name,
                                                  on_interface=interface)
            self.interface_meshes.append(interface_mesh_object)
            # geompy.addToStudyInFather(self._geometry.geometry, interface._edge, interface.name)
            # Remove the published interfaces from the study
            # n_interface = len(self._geometry.interfaces)
            # for i in range(1, n_interface + 1):
            #     edge_to_unpublish = salome.myStudy.FindObject('I_' + str(i))
            #     studyEditor.removeItem(edge_to_unpublish)
            # self._geometry.interfaces = []
            # i_int = 0  # interface (unique edge) counter
            # for edge in self._geometry.edges:
            #     # if len(faces) == 1:
            #     #    continue
            #     mesh_on_edge = self._mesh.GroupOnGeom(edge._edge, typ=SMESH.EDGE)
            #     if mesh_on_edge.IsEmpty():  # the edge does not contain a 1D mesh
            #         self._mesh.RemoveGroup(mesh_on_edge)
            #     else:  # the edge is an interface as it holds a 1D mesh
            #         i_int += 1
            #         interface_name = 'Interface_' + str(i_int)
            #         interface_object = Interface(edge._edge, name=interface_name)
            #         interface_object.neighboring_faces = edge.neighboring_faces
            #         self._geometry.interfaces.append(interface_object)
            #         mesh_on_edge.SetName(interface_name)
            #         interface_mesh_object = InterfaceMesh(mesh_on_edge, name=interface_name)
            #         self.interface_meshes.append(interface_mesh_object)
            #         geompy.addToStudyInFather(self._geometry.geometry, edge._edge, interface_name)
            #         # It is useful to show the incident faces to an interface
            #         for face in edge.neighboring_faces:
            #             geompy.addToStudyInFather(edge._edge, face._face, face.name)
            # # Remove the published edges from the study
            # n_edge = len(self._geometry.edges)
            # for i in range(1, n_edge + 1):
            #     edge_to_unpublish = salome.myStudy.FindObject('Edge_' + str(i))
            #     studyEditor.removeItem(edge_to_unpublish)

    def incident_elements(self, edge, element_type=None):
        """Searches for elements incident to an edge.

        Parameters
        ----------
        edge : list of int
            An edge of an element, given by its two nodes.
        element_type : Mesh.ElementType, optional
            Perform the search for the given element type only.

        Returns
        -------
        list of int
            Element IDs that are incident to the given edge.

        See Also
        --------
        smeshBuilder.Mesh.GetNodeInverseElements

        """
        if not element_type:
            element_type = Mesh.ElementType.ALL
        elements_to_first_node = self._mesh.GetNodeInverseElements(edge[0], element_type.value)
        elements_to_second_node = self._mesh.GetNodeInverseElements(edge[1], element_type.value)
        return list(set(elements_to_first_node).intersection(elements_to_second_node))
        # An alternative, slower, implementation. The function takes an element set as an
        # optional argument to restrict the search. Otherwise, the whole mesh is searched.
        # incident_elements = []
        # possible_elements = element_set.elements() if element_set else self._mesh.GetNodesId()
        # for elem in possible_elements:
        #     if not set(edge).difference(self._mesh.GetElemNodes(elem)):
        #         incident_elements.append(elem)
        # return incident_elements

    def incident_face_mesh(self, interface_mesh):
        """Face meshes incident to an interface mesh.

        .. todo::
            Use this method in _affected_elements as well.

        Parameters
        ----------
        interface_mesh : InterfaceMesh
            Interface mesh for which the connecting face meshes are sought.

        Returns
        -------
        face_mesh : list of FaceMesh
            Face meshes incident to an interface mesh.

        Raises
        ------
        Exception
            If no face mesh is incident to the interface mesh.

        """
        neighboring_faces = interface_mesh.on_interface.neighboring_faces
        # neighboring_face_meshes = [face_mesh for face_mesh in self.face_meshes if
        #                            face_mesh.on_face in neighboring_faces]
        neighboring_face_meshes = []
        for face in neighboring_faces:
            for face_mesh in self.face_meshes:
                if face_mesh.on_face is face:
                    neighboring_face_meshes.append(face_mesh)
        an_interface_node = interface_mesh.nodes()[0]
        incident_face_mesh = [face_mesh for face_mesh in neighboring_face_meshes if
                              an_interface_node in face_mesh.nodes()]
        if not incident_face_mesh:
            Exception('No incident face mesh found for the given interface mesh.')
        return incident_face_mesh

    def generate_element_nodes(self, elements):
        """Nodes of selected elements, returned one at a time.

        Parameters
        ----------
        elements : iterable
            Element IDs.

        Yields
        ------
        list of int
            The first entry of the list is the element ID, the remaining entries are the node IDs
            of the element.

        """
        for element in elements:
            nodes = self._mesh.GetElemNodes(element)
            yield [element, *nodes]

    def one_ring(self, node, definition='connecting'):
        """Elements around a node.

        Parameters
        ----------
        node : int
            Node ID for which the one-ring is searched.
        definition : {'connecting', 'surrounding'}, optional
            What you mean by neighboring elements. See the notes below.

        Returns
        -------
        list
            List of integers (element IDs).

        Notes
        -----
        One should make a distinction between elements `connecting` to a node and elements
        `surrounding` a node. For a mesh with no overlapping nodes, the two definitions give the
        same elements. However, if multiple nodes are located at the same geometrical point, it can
        happen that the incident elements are not connected to the same node.

        """
        if definition == 'connecting':
            return self._mesh.GetElementsByNodes([node], elemType=SMESH.FACE)
        elif definition == 'surrounding':
            coord = self._mesh.GetNodeXYZ(node)
            return self._mesh.FindElementsByPoint(*coord, elementType=SMESH.FACE)
        else:
            raise ValueError('Choose between "connecting" and "surrounding".')

    def point_in_element(self, element, point):
        """Checks whether a point is in an element.

        This method is implemented for 2D meshes only.

        Parameters
        ----------
        element : int
            Element of the mesh.
        point : tuple of float
            Point coordinates (x,y).

        Returns
        -------
        bool
            True if the given element contains the given point.

        Raises
        ------
        Exception
            If the mesh is not two-dimensional.
        ValueError
            If the element does not exist in the mesh.

        Notes
        -----
        This method calls an efficient `matplotlib` function to determine whether a point is in
        a polygon. For alternative implementations, see `this discussion
        <https://stackoverflow.com/questions/36399381/
        whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python>`_.

        See Also
        --------
        matplotlib.path.Path.contains_point

        """
        if self._mesh.MeshDimension() != 2:
            raise Exception('This method is applicable ony 2D meshes only.')
        nodes = self._mesh.GetElemNodes(element)
        if nodes is None:
            raise ValueError('Element {0} does not exist in the mesh.'.format(element))
        element_vertices = [self._mesh.GetNodeXYZ(node)[:2] for node in nodes]
        path = mpltPath.Path(element_vertices)
        return path.contains_point(point)

    def element_edge_normal(self, element, edge):
        """Outward-pointing unit normal to an element edge.

        The edge is assumed to be planar.

        Parameters
        ----------
        element : int
            ID of the element.
        edge : list of int
            Edge of the element for which the normal is to be found. The edge is given by its
            two nodes.

        Returns
        -------
        normal : tuple of float
            Outward-pointing unit normal.

        See Also
        --------
        point_in_element

        """
        # Tangent and unit normal vectors on the edge
        N1 = np.array(self._mesh.GetNodeXYZ(edge[0])[:2])
        N2 = np.array(self._mesh.GetNodeXYZ(edge[1])[:2])
        tangent = N2 - N1
        edge_length = np.sqrt(tangent.dot(tangent))
        tangent = tangent/edge_length
        normal = np.array([tangent[1], -tangent[0]])
        # Shift the midpoint of the edge in the direction of the normal vector
        midpoint = 1/2*(N1 + N2)
        perturbation = 1e-10  # TODO: use a non-absolute measure
        shifted_midpoint = midpoint + perturbation*normal
        # Determine the outward-pointing normal
        normal = -normal if self.point_in_element(element, shifted_midpoint) else normal
        normal = tuple(normal)
        return normal


class FaceMesh:
    """Mesh on a face, part of the whole mesh.

    Parameters
    ----------
    face_mesh : smeshBuilder.Mesh.GroupOnGeom
        The main Salome object wrapped by this class.
    name : str
        Name of the face mesh.
    on_face : Face
        Geometrical face on which this mesh exists.

    """
    def __init__(self, face_mesh, name, on_face):
        self._face_mesh = face_mesh
        self.name = name
        self.on_face = on_face

    def nodes(self):
        """Retrieves the nodes of the face mesh.

        Returns
        -------
        list of int
            Nodes belonging to the face mesh.

        See Also
        --------
        :func:`SMESH.SMESH_GroupBase.GetNodeIDs`

        """
        return self._face_mesh.GetNodeIDs()

    def elements(self):
        """Retrieves the elements of the face mesh.

        Returns
        -------
        list of int
            Elements belonging to the face mesh.

        See Also
        --------
        :func:`SMESH.SMESH_IDSource.GetIDs`

        """
        return self._face_mesh.GetIDs()


class InterfaceMesh:
    """Mesh on an interface, part of the whole mesh.

    Parameters
    ----------
    interface_mesh : smeshBuilder.Mesh.GroupOnGeom
        The main Salome object wrapped by this class.
    name : str
        Name of the interface mesh.
    on_interface : Interface
        Interface on which this mesh exists.

    """

    def __init__(self, interface_mesh, name, on_interface):
        self._interface_mesh = interface_mesh
        self.name = name
        self.on_interface = on_interface

    def endpoint_nodes(self):
        """Nodes at the extremities of the interface mesh.

        Returns
        -------
        ep_nodes : list of int
            Nodes of the end points of the interface on which the interface mesh is defined. If
            the interface is open, it has two end points. When closed, the two end points
            coincide and instead of the two coinciding nodes, a single node is returned.

        """
        # Interface on which the interface mesh exists
        interface = self.on_interface
        # End points of the interface
        ep_1 = geompy.PointCoordinates(geompy.GetFirstVertex(interface._edge))
        ep_2 = geompy.PointCoordinates(geompy.GetLastVertex(interface._edge))
        # Find the end point nodes
        ep_nodes = []
        for node in self.nodes():
            node_coord = tuple(self._interface_mesh.GetMesh().GetNodeXYZ(node))
            if node_coord in [ep_1, ep_2]:
                ep_nodes.append(node)
        return ep_nodes

    def elements_by_nodes(self, nodes):
        """Connecting elements to given nodes.

        Parameters
        ----------
        nodes : list of int
            Nodes for which we want to find the connecting elements.

        Returns
        -------
        list of int
            Elements that are incident to the given nodes.

        See Also
        --------
        smeshBuilder.Mesh.GetElementsByNodes

        """
        connecting_elements = self._interface_mesh.GetMesh().GetElementsByNodes(nodes, SMESH.EDGE)
        own_elements = self.elements()
        return list(set(connecting_elements).intersection(own_elements))

    def nodes(self):
        """Retrieves the nodes of the interface mesh.

        Returns
        -------
        list of int
            Nodes belonging to the interface mesh.

        See Also
        --------
        :func:`SMESH.SMESH_GroupBase.GetNodeIDs`

        """
        return self._interface_mesh.GetNodeIDs()

    def elements(self):
        """Retrieves the elements of the interface mesh.

        Returns
        -------
        list of int
            Elements belonging to the interface mesh.

        See Also
        --------
        :func:`SMESH.SMESH_IDSource.GetIDs`

        """
        return self._interface_mesh.GetIDs()


class CohesiveZone:
    """Constructs zero-thickness elements along the interfaces.

    Parameters
    ----------
    mesh : Mesh
        Mesh into which the cohesive elements will be inserted.

    """
    def __init__(self, mesh):
        # TODO: `mesh` should be public data.
        self._mesh = mesh
        self.cohesive_elements = {}
        self._cohesive_elements_created = False

    def decouple_faces(self):
        """Decouples the face meshes along the interfaces.

        The algorithm consists of two main steps. First, new interface meshes are created that
        overlap with the existing ones and contain independent nodes and interface elements. In
        the same step, the incident face mesh nodes are updated to reflect the topological
        changes. However, in this method, extra nodes are introduced at the junctions,
        leading to a kinematic inconsistency. Therefore, the extraneous interface mesh nodes are
        renumbered in the second step of the algorithm.

        Returns
        -------
        None.

        See Also
        --------
        create_cohesive_elements

        """
        self._enrich_interfaces()
        self._correct_junction_nodes()

    def create_cohesive_elements(self):
        """Creates zero-thickness quadrilateral elements along the interfaces.

        It is necessary that the mesh has already been decoupled along the interfaces by the
        :meth:`decouple_faces` method. That method introduced duplicated nodes and edge elements
        along the interfaces. The purpose of this method is to tie each interface (edge) element
        to its corresponding duplicate in order to form a four-noded zero-thickness element,
        referred to as cohesive element. The bottom edge of the new cohesive element corresponds
        to the original edge element, while its top edge is formed by the duplicated interface
        edge element.

        Returns
        -------
        cohesive_elements : list
            List of nodes that form the cohesive elements. The node numbering follows the
            `node ordering of Salome <https://docs.salome-platform.org/latest/gui/SMESH/
            connectivity.html#connectivity-page>`_, which is the same as the `node ordering in
            Abaqus <https://abaqus-docs.mit.edu/2017/English/SIMACAEELMRefMap/
            simaelm-r-cohesive2d.htm#simaelm-r-cohesive2d-t-nodedef1>`_.

        See Also
        --------
        decouple_faces, _generate_cohesive_element, smeshBuilder.Mesh.AddFace

        """
        if self._cohesive_elements_created:
            print('Cohesive elements already created.')
            return self.cohesive_elements
        # Form the cohesive elements from the overlapping interface meshes
        n_int = len(self._mesh.interface_meshes) // 2
        for i_int in range(n_int):
            original_interface_mesh = self._mesh.interface_meshes[i_int]
            duplicated_interface_mesh = self._mesh.interface_meshes[i_int + n_int]
            interface_name = 'Interface ' + str(i_int)
            self.cohesive_elements[interface_name] = []
            for orig_elem, dupl_elem in zip(original_interface_mesh.elements(),
                                            duplicated_interface_mesh.elements()):
                cohesive_element = self._generate_cohesive_element(orig_elem, dupl_elem)
                elem = self._mesh._mesh.AddFace(cohesive_element)
                self.cohesive_elements[interface_name].append(elem)
        self._cohesive_elements_created = True
        return self.cohesive_elements

    def _enrich_interfaces(self):
        """Inserts new interface elements and nodes into the mesh.

        Although Salome has built-in functionality for duplicating nodes and creating elements,
        even accessible from the GUI with `Modification -> Transformation -> Duplicate Nodes
        and/or Elements`, it does not work with multiple intersecting interfaces or for closed
        interfaces. The reason is that the first step of the two-step procedure Salome performs
        fails in such situations. Therefore, this method uses a modified algorithm for the first
        step, and then calls the second step. These steps are the following:

        1. Find the elements (called affected elements) in the mesh whose node numbers need to be
        changed due to the topological changes in the mesh caused by the introduction of new nodes.

        2. The affected elements are fed to an existing function in Salome, which returns the 1D
        elements it creates from the duplicated nodes. The new interface mesh is stored in the
        :class:`CohesiveZone` object.

        Returns
        -------
        None.

        See Also
        --------
        decouple_faces, _affected_elements, smeshBuilder.Mesh.DoubleNodeElemGroups,
        smeshBuilder.Mesh.MakeGroupByIds

        """
        new_interface_mesh_objects = []
        for i, interface_mesh in enumerate(self._mesh.interface_meshes):
            # The node numbers of the following face elements will be modified
            affected_elements = self._affected_elements(interface_mesh)
            # The "DoubleNodeElemGroups" function of Salome requires us to create a group from
            # the face elements
            a = self._mesh._mesh.MakeGroupByIds(str(i+1), SMESH.FACE, list(affected_elements))
            # Create a new interface mesh; Salome will construct a new mesh group
            new_interface_mesh = self._mesh._mesh.DoubleNodeElemGroups(
                [interface_mesh._interface_mesh], [], [a], True, False)
            if new_interface_mesh is None:  # sometimes Salome cannot perform this operation; WHY?
                raise Exception('Interface {0} could not be enriched.'.format(i+1))
            interface_name = 'Interface_' + str(i+1) + '_doubled'
            new_interface_mesh.SetName(interface_name)
            new_interface_mesh_objects.append(InterfaceMesh(new_interface_mesh, name=interface_name,
                                                            on_interface=interface_mesh.on_interface))
        self._mesh.interface_meshes.extend(new_interface_mesh_objects)

    def _affected_elements(self, interface_mesh):
        """Face elements whose nodes must be renumbered when duplicating an interface mesh.

        Parameters
        ----------
        interface_mesh : InterfaceMesh
            The original interface mesh that will be duplicated.

        Returns
        -------
        elements : set
            Elements of the mesh that require node renumbering.

        See Also
        --------
        _enrich_interfaces

        """
        # Elements connecting to the interface mesh
        yield_one_ring = (self._mesh.one_ring(node, definition='connecting') for node in
                          interface_mesh.nodes())
        neighboring_face_elements = reduce(set().union, yield_one_ring)
        # First neighboring face mesh to this interface mesh
        face_name = interface_mesh.on_interface.neighboring_faces[0].name
        incident_face_mesh = None
        for face_mesh in self._mesh.face_meshes:
            if face_mesh.name == face_name:
                incident_face_mesh = face_mesh
                break
        if not incident_face_mesh:
            raise Exception('Could not find a face incident to this interface mesh.')
        # Only those elements are considered that are part of the face the interface is incident to
        elements = neighboring_face_elements.intersection(incident_face_mesh.elements())
        return elements

    def _correct_junction_nodes(self):
        """Post-processing to handle inconsistent interface nodes at the junctions.

        The interface-wise creation of new edge elements in :meth:`_affected_elements` may result
        in edge element nodes that do not connect the opposite face element nodes on the two
        sides of the interface. This function checks the edge element nodes at the junctions and
        renumbers them so that they hold the same label as the face element nodes they connect to.

        Returns
        -------
        None.

        See Also
        --------
        decouple_faces, _affected_elements, smeshBuilder.Mesh.FindCoincidentNodesOnPart,
        smeshBuilder.Mesh.ChangeElemNodes

        """
        # Check if all interfaces could be enriched in the previous step
        n_interface = len(self._mesh._geometry.interfaces)
        n_interface_mesh = len(self._mesh.interface_meshes)
        if 2*n_interface != n_interface_mesh:
            raise ValueError('Not all interfaces could be enriched.')
        #
        for i, interface in enumerate(self._mesh._geometry.interfaces):
            interface_mesh_original = self._mesh.interface_meshes[i]
            interface_mesh_duplicated = self._mesh.interface_meshes[i+n_interface]
            interface_meshes = [interface_mesh_original, interface_mesh_duplicated]
            face1 = interface.neighboring_faces[0]
            face2 = interface.neighboring_faces[1]
            face_mesh1 = self._mesh.face_meshes[self._mesh._geometry.faces.index(face1)]
            face_mesh2 = self._mesh.face_meshes[self._mesh._geometry.faces.index(face2)]
            face_meshes = [face_mesh2, face_mesh1]
            for interface_mesh, face_mesh in zip(interface_meshes, face_meshes):
                endpoint_nodes = interface_mesh.endpoint_nodes()
                if len(endpoint_nodes) == 1:  # closed interface
                    continue  # closed interfaces are not influenced by other interfaces
                else:  # open interface
                    for node in endpoint_nodes:  # possible inconsistency at the junctions
                        node_collection = face_mesh.nodes()
                        node_collection.append(node)
                        coincident_nodes = self._mesh._mesh.FindCoincidentNodesOnPart(
                            node_collection, 1e-5)
                        if not coincident_nodes:  # if coincident_nodes = [], nodes are consistent
                            continue
                        else:
                            coincident_nodes = coincident_nodes[0]
                            if len(coincident_nodes) != 2:
                                raise Exception(
                                    'len(coincident_nodes) = {0}'.format(len(coincident_nodes)))
                            coincident_nodes = set(coincident_nodes)
                            connecting_face_node = list(coincident_nodes.difference([node]))[0]
                            # Find the interface element that holds the junction node
                            interface_elem = interface_mesh.elements_by_nodes([node])
                            if len(interface_elem) != 1:
                                raise ValueError('One element must contain the junction node.')
                            # Replace that node with the one on the face mesh, determined above
                            interface_elem = interface_elem[0]
                            old_nodes = self._mesh._mesh.GetElemNodes(interface_elem)
                            new_nodes = old_nodes
                            new_nodes[old_nodes.index(node)] = connecting_face_node
                            self._mesh._mesh.ChangeElemNodes(interface_elem, new_nodes)

    def _generate_cohesive_element(self, bottom_element, top_element):
        """Creates a zero-thickness quadrilateral element.

        Parameters
        ----------
        bottom_element : int
            Edge element that will form the bottom edge of the cohesive element.
        top_element: int
            Edge element that will form the top edge of the cohesive element. It is assumed that
            the top element geometrically overlaps with the bottom element.

        Returns
        -------
        list of int
            The four nodes of the cohesive element, numbered counter-clockwise. The node
            ordering adheres to the `node numbering in Abaqus <https://abaqus-docs.mit.edu/2017/
            English/SIMACAEELMRefMap/simaelm-r-cohesive2d.htm#simaelm-r-cohesive2d-t-nodedef1>`_.

        See Also
        --------
        :meth:`Mesh.incident_elements`, :meth:`Mesh.element_edge_normal`,
        smeshBuilder.Mesh.GetElemNodes, smeshBuilder.Mesh.GetNodeXYZ

        """
        # Construct a cohesive element from the interface element pair
        bottom_nodes = self._mesh._mesh.GetElemNodes(bottom_element)
        top_nodes = self._mesh._mesh.GetElemNodes(top_element)
        # Face element the interface element is attached to
        face_elem = self._mesh.incident_elements(bottom_nodes, Mesh.ElementType.FACE)
        # Outward-pointing normal vector to the face element edge
        normal = np.array(self._mesh.element_edge_normal(face_elem[0], bottom_nodes))
        # Translate a node on the top element in the direction of the normal vector
        P2 = np.array(self._mesh._mesh.GetNodeXYZ(top_nodes[1])[:2])
        perturbation = 1e-10  # TODO: use a non-absolute measure
        P2_shifted = P2 + perturbation*normal
        # Ensure that the cohesive element has CCW orientation
        N1 = np.array(self._mesh._mesh.GetNodeXYZ(bottom_nodes[0])[:2])
        N2 = np.array(self._mesh._mesh.GetNodeXYZ(bottom_nodes[1])[:2])
        cross_product = np.cross(N2 - N1, P2_shifted - N2)
        if cross_product > 0:
            return [bottom_nodes[0], bottom_nodes[1], top_nodes[1], top_nodes[0]]
        else:
            return [bottom_nodes[1], bottom_nodes[0], top_nodes[0], top_nodes[1]]


class GUI:
    """Using GUI-related functionalities in Salome.

    Notes
    -----
    A part of Salome's GUI is exposed to Python. To get an idea of what is available, see
    https://docs.salome-platform.org/latest/gui/GUI/text_user_interface.html

    """
    component_map = {Geometry: 'GEOM', Face: 'GEOM', Edge: 'GEOM', Interface: 'GEOM',
                     Mesh: 'SMESH', FaceMesh: 'SMESH', InterfaceMesh: 'SMESH'}

    @staticmethod
    def update_object_browser():
        """Refreshes Salome's object browser.

        Only makes sense if executed with the GUI enabled.

        Returns
        -------
        None

        See Also
        --------
        has_desktop

        """
        if GUI.has_desktop:
            salome.sg.updateObjBrowser()

    @staticmethod
    def has_desktop():
        """Indicates if the Salome GUI is running.

        Returns
        -------
        bool
            True if Salome's GUI is available, False otherwise.

        """
        return salome.sg.hasDesktop()

    class SalomeNoDesktop(Exception):
        """Raised when Salome is run without desktop, but a desktop functionality is invoked."""
        pass

    @staticmethod
    def assert_salome_desktop():
        """Checks if Salome's GUI is available, and raises an exception if it is not.

        This function acts as a helper function when relying on Salome's GUI.

        Raises
        ------
        SalomeNoDesktop
            If Salome's GUI is not available.

        """
        if not GUI.has_desktop():
            raise GUI.SalomeNoDesktop

    @classmethod
    def show(cls, obj, show_only=False):
        """Shows objects in Salome's GUI.

        .. todo:: Support list of objects.

        Parameters
        ----------
        obj : iterable
            The object(s) to be shown in Salome. Objects of the following classes are supported:
            Geometry, Face, Edge, Interface, Mesh, FaceMesh, InterfaceMesh.
        show_only : bool, optional
            If True, the other objects are hidden. The default value is False.

        Returns
        -------
        None

        Raises
        ------
        SalomeNoDesktop
            If Salome's GUI is not available.
        TypeError
            If :code:`obj` is not an object that can be displayed.
        ValueError
            If the given object does not exist in the Salome study.

        Examples
        --------
        For example, you can display an interface mesh and a face mesh by calling

        .. code-block:: python

                GUI.show([interface_mesh, face_mesh])

        where :code:`interface_mesh` and :code:`face_mesh` are :class:`InterfaceMesh` and
        :class:`FaceMesh` objects respectively. This way of using the `show` method provides
        great flexibility as different types of objects can be handled at the same time.

        """
        GUI.assert_salome_desktop()
        component = cls._get_component(obj)
        if not component:
            raise TypeError('Only objects of the following classes are supported: {0}.'.format(
                            [supported_class.__name__ for supported_class in cls.component_map]))
        study_objects = salome.myStudy.FindObjectByName(obj.name, component)
        if len(study_objects) == 0:
            raise ValueError('No object found by name "{0}".'.format(obj.name))
        if len(study_objects) > 1:
            print('More than one object found by the name "{0}". Showing the first one.'.format(
                  obj.name))
        object_id = study_objects[0].GetID()
        if show_only:
            salome.sg.DisplayOnly(object_id)
        else:
            salome.sg.Display(object_id)
        salome.sg.UpdateView()

    @classmethod
    def view(cls, view='top'):
        """Sets the viewpoint.

        Parameters
        ----------
        view : {'front', 'back', 'top', 'bottom', 'left', 'right'}, optional
            Position from which the scene is viewed. The default is 'top'.

        Returns
        -------
        None

        """
        cls.assert_salome_desktop()
        if view == 'front':
            salome.sg.ViewFront()
        elif view == 'back':
            salome.sg.ViewBack()
        elif view == 'top':
            salome.sg.ViewTop()
        elif view == 'bottom':
            salome.sg.ViewBottom()
        elif view == 'left':
            salome.sg.ViewLeft()
        elif view == 'right':
            salome.sg.ViewRight()
        else:
            raise ValueError('Unsupported view.')

    @classmethod
    def _get_component(cls, obj):
        """Determines the component of an object.

        This function maps a class of this module to the Salome module the class uses. For instance,
        class :class:`Face` is mapped to 'GEOM'.

        Parameters
        ----------
        obj
            Any object for which the component name is looked for.

        Returns
        -------
        str or None
            The name of the component the object belongs to. If an object of an unsupported class
            is given, None is returned. For the list of supported classes, see the
            :code:`component_map` member of :class:`GUI` class.

        See Also
        --------
        show

        """
        return cls.component_map.get(type(obj))


# The following code is executed when called from within Salome's Python interpreter
if __name__ == "__main__":

    # Read the microstructure and decompose it into primitives
    microstructure = Geometry()
    microstructure.load('/home/zoltan/University/Lille/Work/Code/CristalX/micro.step')
    microstructure.extract_faces()
    microstructure.extract_edges()
    microstructure.create_interfaces()

    # Generate the mesh
    mesh = Mesh(microstructure)
    mesh.generate()

    # Fetch the mesh on the faces and on the interfaces
    mesh.obtain_face_meshes()
    mesh.obtain_interface_meshes()

    # Create cohesive elements
    cohesive_zone = CohesiveZone(mesh)
    cohesive_zone.decouple_faces()
    cohesive_zone.create_cohesive_elements()

    # Refresh the GUI to make the objects appear
    GUI.update_object_browser()
    GUI.view('top')
