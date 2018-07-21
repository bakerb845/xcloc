#!/usr/bin/env python3
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
import os

class fdxc:
    ##
    # @defgroup pyfdxc Frequency Domain Cross-Correlation
    # @brief Python interface to the FFT-based correlation computation routines.
    # @ingroup pyxcloc
    # @copyright Ben Baker distributed under the MIT license.
    def __init__(self, xcloc_lib, utils):
        xcloc_lib.xcloc_fdxc_initialize.argtypes = (c_int, # npts
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
                   verbose=xctypes.XCLOC_PRINT_WARNINGS,
                   precision=xctypes.XCLOC_SINGLE_PRECISION,
                   accuracy=xctypes.XCLOC_HIGH_ACCURACY):
        """!
        @brief Initializes the Fourier domain based cross-correlation class.

        @param npts       Number of points in input signals.
        @param nsignals   Number of signals that will be input.

        @param nptsPad    This is a tuning parameter for avoiding pathologic DFT
                          lengths.  It is related to the length of the
                          cross-correlograms by: xcLen = 2*nptsPad - 1
        @param verbose    Controls the verbosity of the Fortran module. 
                          0 is quiet.
        @param precision  Controls the precision (single or double) in the
                          Fortran modules.
        @param accuracy   Controls the accuracy of some of the vectorized
                          calculations.  For single precision it is
                          recommended to use high accuracy.  It is always
                          recommended to avoided extended precision.
 
        @retval ierr 0 indicates success
        @ingroup pyfdxc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.initialize.__name__)
        # Input checks 
        if (npts < 1):
            print("%s: npts=%d must be positive"%(fname, npts))
            return -1
        if (nsignals < 2):
            print("%s: nsignals=%d must exceed 2"%(fname, nsignals))
            return -1
        if (precision != xctypes.XCLOC_SINGLE_PRECISION and
            precision != xctypes.XCLOC_DOUBLE_PRECISION):
            print("%s: precision=%d must be 0 or 1"%(fname, precision))
            return -1
        # If nptsPad is set then check that it makes sense; otherwise set it
        if (not nptsPad is None):
            if (nptsPad < npts):
                print("%s: nptsPad must be greater than npts"%(fname, npts, nptsPad))
                return -1
        else:
            nptsPad = npts
        if (xcPairs is None):
            xcPairs = self.utils.computeDefaultXCTable(nsignals)
            xcPairs = xcPairs + 1
        if (amin(xcPairs) == 0): 
            print("%s: xcPairs must be Fortran indexed"%fname)
            return -1
        nxcs = xcPairs.shape[0]
        xcPairs = xcPairs.flatten(order='C')
        xcPairsPtr = xcPairs.ctypes.data_as(POINTER(c_int))
        # Fire up the library
        ierr = c_int(1)
        self.lib.xcloc_fdxc_initialize(npts, nptsPad,
                                       nxcs, xcPairsPtr,
                                       verbose,
                                       precision, accuracy, byref(ierr))
        ierr = ierr.value
        if (ierr != 0):
            print("%s: Failed to initialize cross-correlator"%fname)
        return ierr 

    def __getCorrelograms__(self):
        fname = '%s::%s'%(self.__class__.__name__, self.__getCorrelograms__.__name__)
        nxcs = c_int(1) 
        nptsInXCs = c_int(1)
        ierr = c_int(1)
        self.lib.xcloc_fdxc_getNumberOfCorrelograms(nxcs, byref(ierr))
        if (ierr.value != 0):
            print("%s: Error getting number of xcs"%fname)
            return None
        nxcs = nxcs.value
        self.lib.xcloc_fdxc_getCorrelogramLength(nptsInXCs, byref(ierr))
        if (ierr.value != 0):
            print("%s: Error getting xc lengths"%fname)
            return None
        nptsInXCs = nptsInXCs.value
        precision = c_int(1)
        self.lib.xcloc_fdxc_getPrecision(precision)
        precision = precision.value
        if (precision == xctypes.XCLOC_SINGLE_PRECISION):
            xcs = ascontiguousarray(zeros(nxcs*nptsInXCs, dtype=c_float))
            xcsPtr = xcs.ctypes.data_as(POINTER(c_float))
            self.lib.xcloc_fdxc_getCorrelograms32f(nptsInXCs, nxcs,
                                                   xcsPtr, byref(ierr))
            if (ierr.value != 0):
                print("%s: Error getting float correlograms"%fname)
                return None
        elif (precision == xctypes.XCLOC_DOUBLE_PRECISION):
            xcs = ascontiguousarray(zeros(nxcs*nptsInXCs, dtype=c_double))
            xcsPtr = xcs.ctypes.data_as(POINTER(c_double))
            self.lib.xcloc_fdxc_getCorrelograms64f(nptsInXCs, nxcs,
                                                   xcsPtr, byref(ierr))
            if (ierr.value != 0):
                print("%s: Error getting double correlograms"%fname)
                return None
        else:
            print("%s: Unknown precision %d"%(fname, precision))
            return None
        xcs = xcs.reshape([nxcs, nptsInXCs], order='C')
        return xcs

    def computePhaseCorrelograms(self):
        """!
        @brief Computes the phase correlograms.

        @retval xcs  On successful exit this contains the phase correlograms
                     in a matrix with shape [nxcs x nptsInXCs].
        @retval xcs  On failure this is None.

        @ingroup pyfdxc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.computePhaseCorrelograms.__name__)
        ierr = c_int(1)
        self.lib.xcloc_fdxc_computePhaseCorrelograms(ierr)
        ierr = ierr.value
        if (ierr != 0):
            print("%s: Error computing phase correlograms"%fname)
            return None
        xcs = self.__getCorrelograms__()
        return xcs

    def computeCrossCorrelograms(self):
        """!
        @brief Computes the cross correlograms.

        @retval xcs  On successful exit this contains the cross correlograms
                     in a matrix with shape [nxcs x nptsInXCs].
        @retval xcs  On failure this is None.
        @ingroup pyfdxc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.computeCrossCorrelograms.__name__)
        ierr = c_int(1)
        self.lib.xcloc_fdxc_computeCrossCorrelograms(ierr)
        ierr = ierr.value
        if (ierr != 0):
            print("%s: Error computing cross correlograms"%fname)
            return None
        xcs = self.__getCorrelograms__()
        return xcs

    def setSignals(self, signals):
        """!
        @brief Sets the matrix of signals to correlate.

        @param signals  A matrix of input data with dimension [nsignals x npts].
        @retval ierr    0 indicates success. 
        @ingroup pyfdxc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.setSignals.__name__)
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
                print("%s: Failed to call setSignals32f"%fname)
        elif (signals.dtype == float64):
            signals = ascontiguousarray(signals, float64)
            signalPtr = signals.ctypes.data_as(POINTER(c_double))
            self.lib.xcloc_fdxc_setSignals64f(lds, npts, nsignals, signalPtr,
                                              byref(ierr))
            if (ierr.value != 0):
                print("%s: Failed to call setSignals64f"%fname)
        else:
            print("%s: Precision must be float32 or float64"%fname)
            return -1
        ierr = ierr.value
        return ierr

    def finalize(self):
        """!
        @brief Finalizes the frequency domain cross-correlation library.
        @ingroup pyfdxc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.finalize.__name__)
        self.lib.xcloc_fdxc_finalize()
        #self.lib.xcloc_finalizeF()
        return
