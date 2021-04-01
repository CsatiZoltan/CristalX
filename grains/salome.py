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

"""
