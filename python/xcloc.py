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
from math import pi 
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

    (
       XCLOC_SPXC_DONOT_FILTER,    # No filtering
       XCLOC_SPXC_ENVELOPE_FILTER, # Envelope of correlograms
       XCLOC_SPXC_RMS_FILTER       # RMS filter of correlograms
    ) = map(int, range(3))

    def __init__(self,
                 xcloc_path=os.environ['LD_LIBRARY_PATH'].split(os.pathsep),
                 xcloc_library='libxcloc_shared.so'):
        fname = '%s::%s'%(self.__class__.__name__, self.__init__.__name__)
        # Load the library
        lfound = False
        for path in xcloc_path:
            xcloc_lib = os.path.join(path, xcloc_library)
            if (os.path.isfile(xcloc_lib)):
                lfound = True
                break
        if (lfound):
            xcloc_lib = cdll.LoadLibrary(xcloc_lib)
        else:
            print("%s: Couldn't find libxcloc"%fname)
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
        self.spxc = spxc(xcloc_lib)
        self.lib = xcloc_lib

    def __exit__(self):
        self.fdxc.finalize()
        return
####################################################################################################
#                                              DSM XC                                              #
####################################################################################################
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
####################################################################################################
#                                            Utilities                                             #
####################################################################################################
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
                                                   xcloc.XCLOC_FORTRAN_NUMBERING,
                                                   byref(nxcs), None,
                                                   byref(ierr)) 
        if (ierr.value != 0): 
            print("%s: Error computing workspace size"%fname)
            return None
        nwork = 2*nxcs.value
        xcPairs = ascontiguousarray(zeros(nwork, dtype=c_int))
        xcPairsPtr = xcPairs.ctypes.data_as(POINTER(c_int))
        self.lib.xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                                   xcloc.XCLOC_FORTRAN_NUMBERING,
                                                   byref(nxcs), xcPairsPtr,
                                                   byref(ierr))
        if (ierr.value != 0): 
            print("%s: Error computing default table!"%fname)
            return None
        nxcs = nxcs.value
        xcPairs = xcPairs.reshape([nxcs, 2], order='C') - 1 # C numbering
        return xcPairs

####################################################################################################
#                                            FDXC                                                  #
####################################################################################################
class fdxc:
    """
    This class computes cross-correlations in via the Fourier transform. 
    """
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
        fname = '%s::%s'%(self.__class__.__name__, self.initialize.__name__)
        # Input checks 
        if (npts < 1):
            print("%s: npts=%d must be positive"%(fname, npts))
            return -1
        if (nsignals < 2):
            print("%s: nsignals=%d must exceed 2"%(fname, nsignals))
            return -1
        if (precision != xcloc.XCLOC_SINGLE_PRECISION and
            precision != xcloc.XCLOC_DOUBLE_PRECISION):
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
        if (precision == xcloc.XCLOC_SINGLE_PRECISION):
            xcs = ascontiguousarray(zeros(nxcs*nptsInXCs, dtype=c_float))
            xcsPtr = xcs.ctypes.data_as(POINTER(c_float))
            self.lib.xcloc_fdxc_getCorrelograms32f(nptsInXCs, nxcs,
                                                   xcsPtr, byref(ierr))
            if (ierr.value != 0):
                print("%s: Error getting float correlograms"%fname)
                return None
        elif (precision == xcloc.XCLOC_DOUBLE_PRECISION):
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
        """
        Computes the phase correlograms.
 
        Returns
        -------
        xcs : matrix
           On successful exit this contains the phase correlograms in a matrix
           with shape [nxcs x nptsInXCs].
           On failure this is None.
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
        """
        Computes the cross correlograms.

        Returns
        -------
        xcs : matrix
           On successful exit this contains the cross correlograms in a matrix
           with shape [nxcs x nptsInXCs].
           On failure this is None.
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
        ierr = ierr.value
        return ierr

    def finalize(self):
        """
        Finalizes the frequency domain cross-correlation library.
        """
        fname = '%s::%s'%(self.__class__.__name__, self.finalize.__name__)
        self.lib.xcloc_fdxc_finalize()
        #self.lib.xcloc_finalizeF()
        return

####################################################################################################
#                                    Signals Processing for XCs                                    #
####################################################################################################

class spxc:
    """
    This class helps filter the cross-correlograms prior to migration.
    To mitigate deconstructive interface in the sidelobes due to an
    imperfect velocity model it is useful to either compute an envelope
    or RMS Filter of the correlograms prior to migration.
    """
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
                   ftype=xcloc.XCLOC_SPXC_ENVELOPE_FILTER):
        """
        Initializes the module module that will filter the correlograms.

        Inputs
        ------
        n : int
           Number of filtering coefficients.  This must be an odd number.
        ftype : int
           The filtering type.  This can be none, an envelope, or an RMS filter.
           If it is none, then the value of n does not matter.

        Returns
        -------
        ierr : int
           0 indicates success.
        """ 
        fname = '%s::%s'%(self.__class__.__name__, self.initialize.__name__)
        if (ftype == xcloc.XCLOC_SPXC_ENVELOPE_FILTER or
            ftype == xcloc.XCLOC_SPXC_RMS_FILTER):
            if (n%2 != 1):
                print("%s: Adding a point to the filter length"%fname)
                n = n + 1
        else:
            if (ftype != xcloc.XCLOC_SPXC_DONOT_FILTER):
                print("%s: Invalid filter type"%fname)
                return -1
        ierr = c_int(1)
        self.lib.xcloc_spxc_initialize(n, ftype, ierr)
        ierr = ierr.value
        if (ierr != 0):
            print("%s: Failed to initialize spxc"%fname)
        return ierr

    def compute(self, signals):
        """
        Computes the RMS or envelope of many signals.

        Inputs
        ------
        signals : matrix 
            A matrix of input data with dimension [nsignals x npts].

        Returns
        -------
        sigProc : matrix
            A matrix of envelopes or RMS filtered input signals.  This has
            dimension [nsignals x npts]. 
            On failure this is None.
        """
        fname = '%s::%s'%(self.__class__.__name__, self.compute.__name__)
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
        """
        Finalizes the filtering library.
        """
        self.lib.xcloc_spxc_finalize()
        return

####################################################################################################
#                                         Unit tests                                               #
####################################################################################################
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
    # Compute signal envelopes
    from numpy import sin
    from numpy import exp
    import matplotlib.pyplot as plt
    t = linspace(0, 3, 3001)
    #nsignals = 600
    smat = zeros([nsignals, len(t)])#, dtype=float32)
    for i in range(nsignals):
        smat[i,:] = sin(2.*pi*7*t)*exp(-t/2)
    xcloc.spxc.initialize()
    env = xcloc.spxc.compute(smat)
    #plt.plot(smat[2,:])
    #plt.plot(env[2,:])
    #plt.show()
    xcloc.spxc.finalize()
    # Clean up
    xcloc.fdxc.finalize()
    
