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

class spxc:
    ##
    # @defgroup spxc Signals Processing of Correlograms
    # @brief This class helps filter the cross-correlograms prior to migration.
    #        To mitigate deconstructive interface in the sidelobes due to an
    #        imperfect velocity model it is useful to either compute an envelope
    #        or RMS Filter of the correlograms prior to migration.
    # @ingroup pyxcloc
    # @copyright Ben Baker distributed under the MIT license.
    def __init__(self, xcloc_lib):
        xcloc_lib.xcloc_spxc_initialize.argtypes = (c_int, # number of filter coefficients
                                                    c_int, # filter type
                                                    POINTER(c_int))
        xcloc_lib.xcloc_spxc_filterXCsOutOfPlace64f.argtypes = (c_int, #ldxc
                                                                c_int, #nptsInXCs
                                                                c_int, #nxcs
                                                                POINTER(c_double), # xcs
                                                                POINTER(c_double), # xcsFilt
                                                                POINTER(c_int))
        xcloc_lib.xcloc_spxc_filterXCsOutOfPlace32f.argtypes = (c_int, #ldxc
                                                                c_int, #nptsInXCs
                                                                c_int, #nxcs
                                                                POINTER(c_float), # xcs
                                                                POINTER(c_float), # xcsFilt
                                                                POINTER(c_int))
        xcloc_lib.xcloc_spxc_filterXCsInPlace64f.argtypes = (c_int, #ldxc
                                                             c_int, #nptsInXCs
                                                             c_int, #nxcs
                                                             POINTER(c_double), # xcs
                                                             POINTER(c_int))
        xcloc_lib.xcloc_spxc_filterXCsInPlace32f.argtypes = (c_int, #ldxc
                                                             c_int, #nptsInXCs
                                                             c_int, #nxcs
                                                             POINTER(c_float), # xcs
                                                             POINTER(c_int))
        xcloc_lib.xcloc_spxc_finalize.argtypes = None
        self.lib = xcloc_lib
        return

    def __exit__(self):
        self.finalize()
        return

    def initialize(self,
                   n=301,
                   ftype=xctypes.XCLOC_SPXC_ENVELOPE_FILTER):
        """!
        @brief Initializes the module module that will filter the correlograms.
        @param[in] n       Number of filtering coefficients.  This must be an
                           odd number.
        @param[in] ftype   The filtering type.  This can be an envelope,
                           or an RMS filter.
        @retval ierr   0 indicates success.
        @ingroup spxc
        """ 
        fname = '%s::%s'%(self.__class__.__name__, self.initialize.__name__)
        if (ftype == xctypes.XCLOC_SPXC_ENVELOPE_FILTER or
            ftype == xctypes.XCLOC_SPXC_RMS_FILTER):
            if (n%2 != 1):
                print("%s: Adding a point to the filter length"%fname)
                n = n + 1
        else:
            if (ftype != xctypes.XCLOC_SPXC_DONOT_FILTER):
                print("%s: Invalid filter type"%fname)
                return -1
        ierr = c_int(1)
        self.lib.xcloc_spxc_initialize(n, ftype, ierr)
        ierr = ierr.value
        if (ierr != 0):
            print("%s: Failed to initialize spxc"%fname)
        return ierr

    def filter(self, signals):
        """!
        @brief Filters many signals using an RMS averaging or envelope.
        @param[in] signals  A matrix of input data with dimension
                            [nsignals x npts].
        @retval  sigProc  A matrix of envelopes or RMS filtered input
                          signals.  This has dimension [nsignals x npts].
        @retval  sigProc  On failure this is None.
        @ingroup spxc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.filter.__name__)
        ierr = c_int(1)
        nsignals = signals.shape[0]
        npts = signals.shape[1]
        lds = npts # Leading dimension of fortran array
        signals = reshape(signals, signals.size, order='C')
        sigProc = zeros([nsignals, npts], dtype=signals.dtype)
        if (signals.dtype == float32):
            signals = ascontiguousarray(signals, float32)
            signalPtr = signals.ctypes.data_as(POINTER(c_float)) 
            sigProc = ascontiguousarray(sigProc, float32)
            sigProcPtr = sigProc.ctypes.data_as(POINTER(c_float))
            self.lib.xcloc_spxc_filterXCsOutOfPlace32f(lds, npts, nsignals,
                                                       signalPtr, sigProcPtr,
                                                       byref(ierr))
            if (ierr.value != 0): 
                print("%s: Failed to filter float signals"%fname)
        elif (signals.dtype == float64):
            signals = ascontiguousarray(signals, float64)
            signalPtr = signals.ctypes.data_as(POINTER(c_double))
            sigProc = ascontiguousarray(sigProc, float64)
            sigProcPtr = sigProc.ctypes.data_as(POINTER(c_double))
            self.lib.xcloc_spxc_filterXCsOutOfPlace64f(lds, npts, nsignals,
                                                       signalPtr, sigProcPtr,
                                                       byref(ierr))
            if (ierr.value != 0): 
                print("%s: Failed to filter double signals"%fname)
        else:
            print("%s: Precision must be float32 or float64"%fname)
            return None
        sigProc = sigProc.reshape([nsignals, npts], order='C')
        return sigProc
        

    def finalize(self):
        """!
        Finalizes the filtering library.
        @ingroup spxc
        """
        self.lib.xcloc_spxc_finalize()
        return
