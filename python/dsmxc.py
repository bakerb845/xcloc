#!/usr/bin/env python3

#Purpose: Python interface to the utilities that compute the diffraction
#         stack image of the correlograms.
#Copyright: Ben Baker distributed under the MIT license.

import sys 
from ctypes import cdll
from ctypes import c_int
from ctypes import c_bool
from ctypes import c_float
from ctypes import c_double
from ctypes import byref
from ctypes import POINTER
from numpy import zeros
from numpy import array
from numpy import amax
from numpy import amin
from numpy import reshape
from numpy import float32
from numpy import float64
from numpy import sqrt
from numpy import ascontiguousarray 
from numpy import linspace
from math import pi  
from xclocTypes import xclocTypes as xctypes

class dsmxc:
    ##
    # @defgroup dsmxc Diffraction Stack Migration of Correlograms
    # @brief This class computes the diffraction stack migration image of the
    #        correlograms.
    # @ingroup pyxcloc
    def __init__(self, xcloc_lib):
        xcloc_lib.xcloc_dsmxc_initialize.argtypes = (c_int, # ntables
                                                     c_int, # ngrd
                                                     c_int, # nxcPairs
                                                     c_int, # nptsInXCs
                                                     c_double, # dt
                                                     POINTER(c_int), # xcPairs
                                                     POINTER(c_int))
        xcloc_lib.xcloc_dsmxc_finalize.argtypes = None
        xcloc_lib.xcloc_dsmxc_getImage64f.argtypes = (c_int,
                                                      POINTER(c_double), POINTER(c_int))
        xcloc_lib.xcloc_dsmxc_getImage32f.argtypes = (c_int,
                                                      POINTER(c_float), POINTER(c_int))
        xcloc_lib.xcloc_dsmxc_setCorrelograms64f.argtpyes = (c_int,
                                                             c_int,
                                                             c_int,
                                                             POINTER(c_double),
                                                             POINTER(c_int))
        xcloc_lib.xcloc_dsmxc_setCorrelograms32f.argtpyes = (c_int,
                                                             c_int,
                                                             c_int,
                                                             POINTER(c_float),
                                                             POINTER(c_int))
        xcloc_lib.xcloc_dsmxc_setTable64fF.argtypes = (c_int, c_int,
                                                       POINTER(c_double), POINTER(c_int))
        xcloc_lib.xcloc_dsmxc_setTable32fF.argtypes = (c_int, c_int,
                                                       POINTER(c_float), POINTER(c_int))
        xcloc_lib.xcloc_dsmxc_compute.argtypes = [POINTER(c_int)]

        self.lib = xcloc_lib
        return
