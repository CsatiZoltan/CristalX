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

"""
from enum import Enum

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
