# -*- coding: utf-8 -*-
import warnings
warnings.simplefilter('module', ImportWarning)  # print import warnings once in a module
# warnings.simplefilter('always', RuntimeWarning)

__version__ = '1.1.0'

try:
    import OCC
    HAS_OCCT = True
except ImportError:
    HAS_OCCT = False

try:
    import MEDLoader
    HAS_MED = True
except ImportError:
    HAS_MED = False

try:
    import tables
    HAS_TABLES = True
except ImportError:
    HAS_TABLES = False

try:
    import pyinstrument
    HAS_PYINSTRUMENT = True
except ImportError:
    HAS_PYINSTRUMENT = False
