#!/usr/bin/env python3
"""
Purpose: Python interface to the generic xcloc utilities.
Copyright: Ben Baker distributed under the MIT license.
"""
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

class utils:
    """
    Generic utilities to help with using the library.
    """
    def __init__(self, xcloc_lib):
        xcloc_lib.xcloc_utils_computeDefaultXCTable.argtypes = (c_bool, # ldoAutoCorrs
                                                                c_int,  # nsignals
                                                                c_int,  # nwork
                                                                c_int,  # numbering 
                                                                POINTER(c_int), # nxcs
                                                                POINTER(c_int), # xcPairs
                                                                POINTER(c_int))
        self.lib = xcloc_lib

    def computeDefaultXCTable(self, nsignals, ldoAutoCorrs = False):
        """
        Computes a default cross-correlation table pair.
 
        Inputs
        ------
        nsignals : int
            Number of signals.

        ldoAutoCorrs : bool
            If true then compute auto-correlations.
            Otherwise, only compute cross correlations.

        Returns
        -------
        xcPairs : matrix 
            On successful exit this is a [nxcs x 2] table of cross correlation
            pairs.  Each row is a correlation pair between signal (column 1) and
            signal (column 2).  The signal indices are C numbered.  
            Otherwise, this is None.
        """
        fname = '%s::%s'%(self.__class__.__name__, self.computeDefaultXCTable.__name__)
        ierr = c_int(1)
        nxcs = c_int(1)
        # Space query
        nwork =-1 
        self.lib.xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                                   xctypes.XCLOC_FORTRAN_NUMBERING,
                                                   byref(nxcs), None,
                                                   byref(ierr)) 
        if (ierr.value != 0): 
            print("%s: Error computing workspace size"%fname)
            return None
        nwork = 2*nxcs.value
        xcPairs = ascontiguousarray(zeros(nwork, dtype=c_int))
        xcPairsPtr = xcPairs.ctypes.data_as(POINTER(c_int))
        self.lib.xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                                   xctypes.XCLOC_FORTRAN_NUMBERING,
                                                   byref(nxcs), xcPairsPtr,
                                                   byref(ierr))
        if (ierr.value != 0): 
            print("%s: Error computing default table!"%fname)
            return None
        nxcs = nxcs.value
        xcPairs = xcPairs.reshape([nxcs, 2], order='C') - 1 # C numbering
        return xcPairs
