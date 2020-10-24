# Processing a .med file



After exporting the mesh from Salome to a MED file, we may want to perform certain operations on it. The [MEDCoupling](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/index.html) tool of Salome provides C++ and Python APIs for this purpose. However, that requires the user to

-  have Salome installed as those APIs are available from the Salome kernel
-  get to know the API

Moreover, it can happen that some mesh processing functionalities they may want to use does not exist. Since meshes consisting of cells of the same type (e.g. triangles) can be represented as homogeneous and contiguous arrays, converting the mesh from MED to numpy arrays seems a reasonable choice. This is what our *med* module does: it provides a thin wrapper around *MEDCoupling* to extract the mesh and the defined groups (cell and vertex groups) from the MED file and convert them to numpy arrays. This way, the user who deals with numerical modelling can implement their mesh processing algorithms based on numpy arrays, which is fast and straightforward. Furthermore, the person who performs the CAD operations and has Salome installed, can use our *med* module to export the mesh to numpy arrays so that the numerical analyst can directly work on it without having to have Salome installed and without any knowledge on the *MEDCoupling* API.

If you want to know more about the implementation details, read the documentation for the [med module](med_module).



## Using the module

To use our *med* module, access to the *MEDLoader* module from Salome is required. In the following, we assume that Salome has been installed and added to the path, so the command `salome` is available. If you 
- want to use the Python REPL:
  
   ```bash
   $ salome shell -- python
   Python 3.6.5 (default, Dec 16 2019, 16:42:15) 
   [GCC 7.3.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   ```
   
   Now, you have access to the *MEDCoupling* Python API. E.g.
   ```python
   >>> from MEDLoader import MEDFileData
   ```
   
- want to use the PyCharm IDE

   ```bash
   $ salome shell "path_to_PyCharm/bin/pycharm.sh"
   ```
   This will start PyCharm. In the IDE (PyCharm in this example), set the interpreter to that of Salome's Python. For me, it is located at `BINARIES-UB18.04/Python/bin/python3`, where `UB18.04` refers to the fact that I [download](https://www.salome-platform.org/downloads/current-version)ed a pre-compiled Salome binaries for Ubuntu 18.04.

To learn more about the `salome` command, read the [manual](https://docs.salome-platform.org/latest/tui/KERNEL/salome_command.html).