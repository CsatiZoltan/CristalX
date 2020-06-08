# -*- coding: utf-8 -*-
import warnings
warnings.simplefilter('module', ImportWarning)  # print import warnings once in a module

try:
    import OCC
    HAS_OCCT = True
except ImportError:
    HAS_OCCT = False
