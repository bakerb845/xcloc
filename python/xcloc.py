#!/usr/bin/env python3
"""
Python interface to the xcloc interferometric location software.

Copyright: Ben Baker distributed under the MIT license.
"""
from ctypes import cdll
from ctypes import c_int
from ctypes import c_double
from ctypes import byref
from ctypes import POINTER
from numpy import zeros
from numpy import array
from numpy import reshape
from numpy import float64
from numpy import sqrt
from numpy import ascontiguousarray 
from numpy import linspace
import os
#import mpi4py

(
   XCLOC_HIGH_ACCURACY, # High accuracy vector math calculations
   XCLOC_LOW_ACCURACY,  # Low accuracy vector math calculations
   XCLOC_EP_ACCURACY    # Really fast/inaccurate vector math calculations
) = map(c_int, range(3))

(
   XCLOC_SINGLE_PRECISION, # Single precision
   XCLOC_DOUBLE_PRECISION, # Double precision
) = map(c_int, range(2))

class xcloc:
    def __init__(self,
                 xcloc_path=os.environ['LD_LIBRARY_PATH'].split(os.pathsep),
                 xcloc_library='libxcloc_shared.so'):
 
        self.fdxc = fdxc(xcloc_path, xcloc_library)
    def __exit__(self):
        self.fdxc.finalize()
        return


class fdxc:
    """
    This class computes cross-correlations in via the Fourier transform. 
    """
    def __init__(self,
                 xcloc_path=os.environ['LD_LIBRARY_PATH'].split(os.pathsep),
                 xcloc_library='libxcloc_shared.so'):
        lfound = False
        for path in xcloc_path:
            xcloc_path = os.path.join(path, xcloc_library)
            if (os.path.isfile(xcloc_path)):
                lfound = True
                break
        if (lfound):
            xcloc_lib = cdll.LoadLibrary(xcloc_path) 
        else:
            print("Couldn't find libxcloc")
            return
        xcloc_lib.xcloc_initializeF.argtypes = None 
        xcloc_lib.xcloc_finalizeF.argtypes = None

        xcloc_lib.xcloc_fdxc_initialize.argtypes = (c_int,
                                                    c_int,
                                                    c_int,
                                                    c_int,
                                                    c_int,
                                                    c_int,
                                                    POINTER(c_int))
        xcloc_lib.xcloc_fdxc_computePhaseCorrelograms.argtypes = [POINTER(c_int)]
        xcloc_lib.xcloc_fdxc_computeCrossCorrelograms.argtypes = [POINTER(c_int)]
        xcloc_lib.xcloc_fdxc_finalize.argtypes = None 
        self.lib = xcloc_lib

    def __exit__(self):
        self.finalize()
        return

    def initialize(self, npts, nsignals, nptsPad=None,
                   verbose=0,
                   precision=XCLOC_SINGLE_PRECISION,
                   accuracy=XCLOC_HIGH_ACCURACY):
        """
        Initializes the Fourier domain based cross-correlation class.

        Inputs
        ------
        npts : int
            Number of points in input signals.
        nsignals : int
            Number of signals that will be input.

        Optional Inputs
        ---------------
        nptsPad : int
            This is a tuning parameter for avoiding pathologic DFT lengths.
            It is related to the length of the cross-correlograms by:
              xcLen = 2*nptsPad - 1
        verbose : int
            Controls the verbosity of the Fortran module.  0 is quiet.
        precision : int
            Controls the precision (single or double) in the Fortran modules.
        accuracy: int
            Controls the accuracy of some of the vectorized calculations.
            For single precision it is recommended to use high accuracy.
            It is always recommended to avoided extended precision.
 
        Returns
        -------
        ierr : int
            0 indicates success
        """
        self.lib.xcloc_initializeF(2)
        # Input checks 
        if (npts < 1):
            print("npts=%d must be positive"%npts)
            return -1
        if (nsignals < 2):
            print("nsignals=%d must exceed 2"%nsignals)
            return -1
        if (precision != XCLOC_SINGLE_PRECISION and
            precision != XCLOC_DOUBLE_PRECISION):
            print("precision=%d must be 0 or 1"%precision)
            return -1
        # If nptsPad is set then check that it makes sense; otherwise set it
        if (nptsPad != None):
            if (nptsPad < npts):
                print("nptsPad must be greater than npts", npts, nptsPad)
                return -1
        else:
            nptsPad = npts
        # Fire up the library
        ierr = c_int(1)
        self.lib.xcloc_fdxc_initialize(npts, nsignals, nptsPad, verbose,
                                       precision, accuracy, byref(ierr))
        return ierr 

    def finalize(self):
        """
        Finalizes the frequency domain cross-correlation library.
        """
        self.lib.xcloc_fdxc_finalize()
        self.lib.xcloc_finalizeF()
        return

 

if __name__ == "__main__":
    xcloc = xcloc(xcloc_path=['/home/bakerb25/C/xcloc/lib'])
    xcloc.fdxc.initialize(100, 10)
    xcloc.fdxc.finalize()
