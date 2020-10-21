# Geometry and mesh processing in Salome



The *cad* module, as shown in the [Algorithms](algorithms.md#from-image-to-geometry) section, allows to approximate a segmented image with planar spline surfaces (splinegons) and those surfaces can be written to a STEP file. Now, we demonstrate on the sample microstructure how to repair the geometry and generate a conforming mesh on the splinegons. All the manipulations are done in Salome 9.4.0.



## Geometry

1. Import the geometry
   To import a STEP file (obtained by writing the splinegons into a STEP file), select the *Geometry* module and then *File -> Import -> STEP*.

2. Repair the geometry

   For the sample microstructure, one grain has not been identified, i.e. it contains a hole. To fill the hole, we perform the following steps (always select the object that was created the last time): 
   1. *Repair -> Close Contour* and select the (in our sample) two boundary curves of the unidentified grain.
   2. *Repair -> Suppress Holes* and select the two boundary curves you selected in the previous step. Salome informs you that a face will be created in place of the hole.
   3. *Repair -> Limit tolerance* and set the tolerance to 1e-3.
   4. *Repair -> Sewing* and set 1e-2 for the tolerance.

3. Extract the grains
   We will need the mesh on each grain, therefore, we extract the grain faces by "exploding" the microstructure using *New Entity -> Explode*. In the dialog box, select the sewn geometry (end result of the geometry repairing workflow above) as *Main Object* and *Face* as *Sub-shapes Type*. Salome properly identifies the 250 sub-shapes, i.e. the grains.

4. Fix the "artificial" grains

   The sample microstructure describes the central part of a tensile specimen. However, to adhere to the Saint-Venant principle, the loading is exerted further from this central zone. Hence, we created long enough regions in the direction of the prescribed load by attaching two rectangular faces to the central region. These faces are not precise rectangles because the boundary of the microstructure does not consist of straight segments. Nevertheless, from now on, we will use the term rectangle to describe the extended region. To achieve perfect matching between the rectangle and the boundary of the microstructure, one side of the rectangle must contain the same lines as the boundary of the microstructure. The other three sides will be straight line segments.

   1. We want to extract the boundary lines of the microstructure. To do that, we explode the boundary grains into edges with *New Entity -> Explode*. The *Main Object* in the dialog box is a selected boundary grain, and the *Sub-shape Type* is *Edge*. Do this with all the boundary grains.
   2. Create the three sides of the rectangle
      1. Create the four vertices of the rectangle with *New Entry -> Basic -> Point*. Two of these points are the extremities of the microstructure, the other two are given such that
         -  the rectangle's three sides become parallel to the coordinate axes (so that boundary conditions are easily prescribed)
         -  the length of the longer side of the rectangle is long enough compared to the microstructure (we applied a ratio of 5)
         
      2. Connect these vertices with three lines, which forms the three sides of the rectangle: *New Entity -> Build -> Edge*.
   3. Create the rectangle by connecting the three line segments constructed in the previous step with the boundary lines of the fourth side, obtained by exploding the boundary grains: *New Entity -> Build -> Wire*.
   4. Delete the original "artificial" grain as it will be replaced by the new one.
   5. Create the surface enclosed by the rectangle: *New Entity -> Build -> Face*, and select the wire created before. Perform steps 1-5 for the other rectangle (on the other side of the microstructure) as well.
   6. In order to create a conforming mesh on the whole domain (i.e. on the microstructure and on the two rectangles), we need to create a compound surface: *New Entity -> Build -> Compound* and select the grains plus the two rectangles.
   7. It is not enough to create the mesh, we also need to know which elements belong to which grains. To allow this, we explode the compound surface: select the compound object created in the previous step as the *Main Object* in *New Entity -> Explode*. The *Sub-shape Type* is *Face*, as before.
   



The geometry is now impeccable, let us start meshing.



## Mesh

Select the *Mesh* module.

1. Create the mesh

   1. Choose *Mesh -> Create Mesh* from the menu. The geometry on which the mesh will be created is the compound surface.

   2. The mesh type in this study is *Triangular* and the algorithm is *NETGEN 1D-2D*. Choose *NETGEN 2D Parameters* as a hypothesis and experiment with the settings to suit your needs. Accept the changes.

   3. Choose *Mesh -> Compute* to generate the mesh with the chosen settings.

   If you want to alter the mesh, change the hypothesis and compute the new mesh.

2. Obtain the sub-meshes for the grains

   The segmented image contains regions that are not part of the recrystallized region, but they belong to the homogeneous region. We want to handle them in the same way as the two artificial grains. We do not need to merge the surfaces, just put them in the same element group. Therefore, we create

   -  a common element group for the elements lying in the homogeneous region
   -  one element group for each grain in the heterogeneous (recrystallized) region

   using *Mesh -> Create Groups from Geometry* and selecting all the 250 faces of the geometry. To form a single group for the heterogeneous region, merge the corresponding element groups in *Mesh -> Union Groups*. Once done, delete the groups that were used in forming the union.

3. We will explicitly need to prescribe boundary conditions on the left and right sides of the rectangular domain. For this, the nodes on those sides are collected with *Mesh -> Create Group* and set the *Elements Type* to *Node*. Select the nodes on the left edge and give the node set a name. Repeat it for the nodes on the right.

3. Export the mesh
   Click on *File -> Export -> MED file*. The sub-meshes (element groups, node groups) are also exported.



The mesh is ready for further processing. We have another [guide](med) that discusses how to handle the MED file from within Python.