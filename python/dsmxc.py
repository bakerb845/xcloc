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
from numpy import unique
from numpy import float64
from numpy import sqrt
from numpy import ascontiguousarray 
from numpy import linspace
from math import pi  
from xclocTypes import xclocTypes as xctypes

class dsmxc:
    ##
    # @defgroup pydsmxc Diffraction Stack Migration of Correlograms
    # @brief This class computes the diffraction stack migration image of the
    #        correlograms.
    # @ingroup pyxcloc
    def __init__(self, xcloc_lib):
        xcloc_lib.xcloc_dsmxc_initialize.argtypes = (c_int, # ngrd
                                                     c_int, # nxcPairs
                                                     c_int, # nptsInXCs
                                                     c_double, # dt
                                                     POINTER(c_int), # xcPairs
                                                     c_int, # verbosity
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
        xcloc_lib.xcloc_dsmxc_setTable64f.argtypes = (c_int, c_int,
                                                      POINTER(c_double), POINTER(c_int))
        xcloc_lib.xcloc_dsmxc_setTable32f.argtypes = (c_int, c_int,
                                                      POINTER(c_float), POINTER(c_int))
        xcloc_lib.xcloc_dsmxc_compute.argtypes = [POINTER(c_int)]

        self.lib = xcloc_lib
        self.linit = False
        return

    def initialize(self, ngrd, nptsInXCs, xcPairs, dt,
                   verbose=xctypes.XCLOC_PRINT_WARNINGS):
        """!
        @param[in] ngrd       Number of points in migration grid.
        @param[in] nptsInXCs  Number of points in correlograms.  This must be
                              positive and an odd number.
        @param[in] xcPairs    This is a [nxcPairs x 2] matrix that defines
                              the signal indices comprising a correlation
                              pair.
        @param[in] dt         Sampling period (seconds) of correlograms. 
        @param[in] verbose    Controls verbosity.
        @retval ierr          Error flag where 0 indicates success.
        @ingroup dsmxc 
        """
        fname = '%s::%s'%(self.__class__.__name__, self.initialize.__name__)
        if (ngrd < 1):
            print("%s: ngrd=%d must be positive"%(fname, ngrd))
            return -1
        if (nptsInXCs < 1 or nptsInXCs%2 == 0):
            print("%s: nptsInXCs=%d must be positive and odd"%(fname, nptsInXCs))
            return -1
        if (dt <= 0.0):
            print("%s: sampling period=%f must be positive"%(fname, dt))
        nxcPairs = xcPairs.shape[0]
        xcPairs = xcPairs.flatten(order='C') 
        xcPairs = ascontiguousarray(xcPairs)
        xcPairsPtr = xcPairsPtr.ctypes.data_as(POINTER(c_int)) 
        ierr = c_int(1)
        self.lib.xcloc_dsmxc_initialize(ngrd, nxcPairs, nptsInXCs, dt, xcPairsPtr,
                                        verbose, ierr)
        ierr = ierr.value
        if (ierr != 0):
            print("%s: Failed to initialize dsmxc"%fname)
            return ierr 
        self.linit = True
        return ierr 

    def finalize(self):
        """!
        Releases memory on the diffraction stack migration module.
        @ingroup dsmxc
        """
        self.lib.xcloc_dsmxc_finalize()
        self.linit = False
        return

    def __exit__(self):
        self.finalize()
        return
