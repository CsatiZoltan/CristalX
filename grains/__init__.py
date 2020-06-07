# -*- coding: utf-8 -*-

try:
    import OCC
    HAS_OCCT = True
except ImportError:
    HAS_OCCT = False
