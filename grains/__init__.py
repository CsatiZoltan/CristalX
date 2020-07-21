# -*- coding: utf-8 -*-
import warnings
warnings.simplefilter('module', ImportWarning)  # print import warnings once in a module
warnings.simplefilter('always', RuntimeWarning)

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

