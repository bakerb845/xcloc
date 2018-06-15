#!/usr/bin/env python3
"""
Python interface to the xcloc interferometric location software.

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
import os
#import mpi4py

class xcloc:
    (
       XCLOC_HIGH_ACCURACY, # High accuracy vector math calculations
       XCLOC_LOW_ACCURACY,  # Low accuracy vector math calculations
       XCLOC_EP_ACCURACY    # Really fast/inaccurate vector math calculations
    ) = map(int, range(3))

    (
       XCLOC_SINGLE_PRECISION, # Single precision
       XCLOC_DOUBLE_PRECISION  # Double precision
    ) = map(int, range(2))

    (
       XCLOC_C_NUMBERING,       # C numbering
       XCLOC_FORTRAN_NUMBERING  # Fortran numbering
    ) = map(int, range(2))

    def __init__(self,
                 xcloc_path=os.environ['LD_LIBRARY_PATH'].split(os.pathsep),
                 xcloc_library='libxcloc_shared.so'):
        # Load the library
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
        ##################################################################################
        #                                     XCLOC                                      #
        ##################################################################################
        xcloc_lib.xcloc_initializeF.argtypes = None 
        xcloc_lib.xcloc_finalizeF.argtypes = None

        # Hook up the modules
        self.utils = utils(xcloc_lib)
        self.fdxc = fdxc(xcloc_lib, self.utils)
        self.dsmxc = dsmxc(xcloc_lib)
        self.lib = xcloc_lib

    def __exit__(self):
        self.fdxc.finalize()
        return

class dsmxc:
    """
    This class computes the diffraction stack migration image of the
    correlograms.
    """
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
        ierr = c_int(1)
        nxcs = c_int(1)
        # Space query
        nwork =-1 
        self.lib.xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                                   xcloc.XCLOC_FORTRAN_NUMBERING,
                                                   byref(nxcs), None,
                                                   byref(ierr)) 
        if (ierr.value != 0): 
            print("Error computing workspace size")
            return None
        nwork = 2*nxcs.value
        xcPairs = ascontiguousarray(zeros(nwork, dtype=c_int))
        xcPairsPtr = xcPairs.ctypes.data_as(POINTER(c_int))
        self.lib.xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                                   xcloc.XCLOC_FORTRAN_NUMBERING,
                                                   byref(nxcs), xcPairsPtr,
                                                   byref(ierr))
        if (ierr.value != 0): 
            print("Error computing default table!")
            return None
        nxcs = nxcs.value
        xcPairs = xcPairs.reshape([nxcs, 2], order='C') - 1 # C numbering
        return xcPairs


class fdxc:
    """
    This class computes cross-correlations in via the Fourier transform. 
    """
    def __init__(self, xcloc_lib, utils):
        xcloc_lib.xcloc_fdxc_initialize.argtypes = (c_int, # npts
                                                    c_int, # nsignals
                                                    c_int, # nptsPad
                                                    c_int, # nxcs
                                                    POINTER(c_int), # xcPairs
                                                    c_int, # verbose
                                                    c_int, # precision
                                                    c_int, # accuracy
                                                    POINTER(c_int))
        xcloc_lib.xcloc_fdxc_setSignals64f.argtypes = (c_int,
                                                       c_int,
                                                       c_int,
                                                       POINTER(c_double),
                                                       POINTER(c_int))
        xcloc_lib.xcloc_fdxc_setSignals32f.argtypes = (c_int,
                                                       c_int,
                                                       c_int,
                                                       POINTER(c_float),
                                                       POINTER(c_int))
        xcloc_lib.xcloc_fdxc_getNumberOfSignals.argtypes = (POINTER(c_int),
                                                            POINTER(c_int))
        xcloc_lib.xcloc_fdxc_getNumberOfCorrelograms.argtypes = (POINTER(c_int),
                                                                 POINTER(c_int))
        xcloc_lib.xcloc_fdxc_getCorrelogramLength.argtypes = (POINTER(c_int),
                                                              POINTER(c_int))
        xcloc_lib.xcloc_fdxc_getPrecision.argtypes = [POINTER(c_int)]
        xcloc_lib.xcloc_fdxc_getCorrelograms64f.argtypes = (c_int,
                                                            c_int,
                                                            POINTER(c_double),
                                                            POINTER(c_int))
        xcloc_lib.xcloc_fdxc_getCorrelograms32f.argtypes = (c_int,
                                                            c_int,
                                                            POINTER(c_float),
                                                            POINTER(c_int))
        xcloc_lib.xcloc_fdxc_setSignal32fF.argtypes = (c_int,
                                                       c_int, 
                                                       POINTER(c_float),
                                                       POINTER(c_int))
        xcloc_lib.xcloc_fdxc_computePhaseCorrelograms.argtypes = [POINTER(c_int)]
        xcloc_lib.xcloc_fdxc_computeCrossCorrelograms.argtypes = [POINTER(c_int)]
        xcloc_lib.xcloc_fdxc_finalize.argtypes = None
        self.lib = xcloc_lib
        self.utils = utils

    def __exit__(self):
        self.finalize()
        return

    def initialize(self, npts, nsignals, nptsPad=None,
                   xcPairs=None,
                   verbose=0,
                   precision=xcloc.XCLOC_SINGLE_PRECISION,
                   accuracy=xcloc.XCLOC_HIGH_ACCURACY):
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
        #self.lib.xcloc_initializeF(2)
        # Input checks 
        if (npts < 1):
            print("npts=%d must be positive"%npts)
            return -1
        if (nsignals < 2):
            print("nsignals=%d must exceed 2"%nsignals)
            return -1
        if (precision != xcloc.XCLOC_SINGLE_PRECISION and
            precision != xcloc.XCLOC_DOUBLE_PRECISION):
            print("precision=%d must be 0 or 1"%precision)
            return -1
        # If nptsPad is set then check that it makes sense; otherwise set it
        if (not nptsPad is None):
            if (nptsPad < npts):
                print("nptsPad must be greater than npts", npts, nptsPad)
                return -1
        else:
            nptsPad = npts
        if (xcPairs is None):
            xcPairs = self.utils.computeDefaultXCTable(nsignals)
            xcPairs = xcPairs + 1
        if (amin(xcPairs) == 0):
            print("xcPairs must be Fortran indexed")
            return -1
        nxcs = xcPairs.shape[0]
        xcPairs = xcPairs.flatten(order='C')
        xcPairsPtr = xcPairs.ctypes.data_as(POINTER(c_int))
        # Fire up the library
        ierr = c_int(1)
        self.lib.xcloc_fdxc_initialize(npts, nsignals, nptsPad,
                                       nxcs, xcPairsPtr,
                                       verbose,
                                       precision, accuracy, byref(ierr))
        return ierr 

    def __getCorrelograms__(self):
        nxcs = c_int(1) 
        nptsInXCs = c_int(1)
        ierr = c_int(1)
        self.lib.xcloc_fdxc_getNumberOfCorrelograms(nxcs, byref(ierr))
        if (ierr.value != 0):
            print("Error getting number of xcs")
            return None
        nxcs = nxcs.value
        self.lib.xcloc_fdxc_getCorrelogramLength(nptsInXCs, byref(ierr))
        if (ierr.value != 0):
            print("Error getting xc lengths")
            return None
        nptsInXCs = nptsInXCs.value
        precision = c_int(1)
        self.lib.xcloc_fdxc_getPrecision(precision)
        precision = precision.value
        if (precision == xcloc.XCLOC_SINGLE_PRECISION):
            xcs = ascontiguousarray(zeros(nxcs*nptsInXCs, dtype=c_float))
            xcsPtr = xcs.ctypes.data_as(POINTER(c_float))
            self.lib.xcloc_fdxc_getCorrelograms32f(nptsInXCs, nxcs,
                                                   xcsPtr, byref(ierr))
            if (ierr.value != 0):
                print("Error getting float correlograms")
                return None
        elif (precision == xcloc.XCLOC_DOUBLE_PRECISION):
            xcs = ascontiguousarray(zeros(nxcs*nptsInXCs, dtype=c_double))
            xcsPtr = xcs.ctypes.data_as(POINTER(c_double))
            self.lib.xcloc_fdxc_getCorrelograms64f(nptsInXCs, nxcs,
                                                   xcsPtr, byref(ierr))
            if (ierr.value != 0):
                print("Error getting double correlograms")
                return None
        else:
            print("Unknown precision", precision)
            return None
        xcs = xcs.reshape([nxcs, nptsInXCs], order='C')
        return xcs

    def computePhaseCorrelograms(self):
        """
        Computes the phase correlograms.
 
        Returns
        -------
        xcs : matrix
           On successful exit this contains the phase correlograms in a matrix
           with shape [nxcs x nptsInXCs].
           On failure this is None.
        """
        ierr = c_int(1)
        self.lib.xcloc_fdxc_computePhaseCorrelograms(ierr)
        ierr = ierr.value
        if (ierr != 0):
            print("Error computing phase correlograms")
            return None
        xcs = self.__getCorrelograms__()
        return xcs

    def computeCrossCorrelograms(self):
        """
        Computes the cross correlograms.

        Returns
        -------
        xcs : matrix
           On successful exit this contains the cross correlograms in a matrix
           with shape [nxcs x nptsInXCs].
           On failure this is None.
        """
        ierr = c_int(1)
        self.lib.xcloc_fdxc_computeCrossCorrelograms(ierr)
        ierr = ierr.value
        if (ierr != 0):
            print("Error computing cross correlograms")
            return None
        xcs = self.__getCorrelograms__()
        return xcs

    def setSignals(self, signals):
        """
        Sets the matrix of signals to correlate.

        Inputs
        ------
        signals : matrix 
            A matrix of input data with dimension [nsignals x npts].

        Returns
        -------
        ierr : int
            0 indicates success. 
        """
        ierr = c_int(1)
        nsignals = signals.shape[0]
        npts = signals.shape[1]
        lds = npts # Leading dimension of fortran array
        signals = reshape(signals, signals.size, order='C')
        if (signals.dtype == float32):
            signals = ascontiguousarray(signals, float32)
            signalPtr = signals.ctypes.data_as(POINTER(c_float)) 
            self.lib.xcloc_fdxc_setSignals32f(lds, npts, nsignals, signalPtr,
                                              byref(ierr))
            if (ierr.value != 0):
                print("Failed to call setSignals32f")
        elif (signals.dtype == float64):
            signals = ascontiguousarray(signals, float64)
            signalPtr = signals.ctypes.data_as(POINTER(c_double))
            self.lib.xcloc_fdxc_setSignals64f(lds, npts, nsignals, signalPtr,
                                              byref(ierr))
            if (ierr.value != 0):
                print("Failed to call setSignals64f")
        else:
            print("Precision must be float32 or float64")
        ierr = ierr.value
        return ierr

    def finalize(self):
        """
        Finalizes the frequency domain cross-correlation library.
        """
        self.lib.xcloc_fdxc_finalize()
        #self.lib.xcloc_finalizeF()
        return

def unit_test(xcloc_path):
    xcloc = xcloc(xcloc_path=[xcloc_path])
    npts = 10     # Length of signal
    nsignals = 5  # Number of signals
    # Initialize - default is single precision + high accuracy
    ierr = xcloc.fdxc.initialize(npts, nsignals)
    if (ierr != 0):
        sys.exit("initialization failed")
    # Set some signals to correlate
    signals = zeros([nsignals, npts])
    signals[0,:] = [ 1, 2, 3, 4, 5, -1,  0, 0,  0,  0]
    signals[1,:] = [ 3, 2, 1, 3, 2,  1,  0, 0,  0,  0]
    signals[2,:] = [-1, 2,-1,-2, 1, -1,  0, 0,  0,  0]
    signals[3,:] = [ 4,-2, 3,-1, 5, -1,  0, 0,  0,  0]
    signals[4,:] = [ 0,-3,-4, 5,-2,  3,  0, 0,  0,  0]
    ierr = xcloc.fdxc.setSignals(signals)
    if (ierr != 0):
        sys.exit("failed to set signals")
    # Compute the correlograms
    xcs = xcloc.fdxc.computeCrossCorrelograms()
    if xcs is None:
        print("failed to compute correlograms")
    print(xcs)
    # Clean up
    xcloc.fdxc.finalize()
    return 

if __name__ == "__main__":
    xcloc = xcloc(xcloc_path=['/home/bakerb25/C/xcloc/lib'])
    npts = 10     # Length of signal
    nsignals = 5  # Number of signals
    # Initialize - default is single precision + high accuracy
    xcloc.fdxc.initialize(npts, nsignals, precision=xcloc.XCLOC_DOUBLE_PRECISION)
    # Set some signals to correlate
    signals = zeros([nsignals, npts])
    signals[0,:] = [ 1, 2, 3, 4, 5, -1,  0, 0,  0,  0]
    signals[1,:] = [ 3, 2, 1, 3, 2,  1,  0, 0,  0,  0]
    signals[2,:] = [-1, 2,-1,-2, 1, -1,  0, 0,  0,  0]
    signals[3,:] = [ 4,-2, 3,-1, 5, -1,  0, 0,  0,  0]
    signals[4,:] = [ 0,-3,-4, 5,-2,  3,  0, 0,  0,  0]
    xcloc.fdxc.setSignals(signals)
    # Compute the correlograms
    xcs = xcloc.fdxc.computeCrossCorrelograms()
    print(xcs) 
    # Clean up
    xcloc.fdxc.finalize()
    
