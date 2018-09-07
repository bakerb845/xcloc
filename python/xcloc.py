#!/usr/bin/env python3

import sys
from ctypes import cdll
from ctypes import c_int
from ctypes import c_int64
from ctypes import c_bool
from ctypes import c_float
from ctypes import c_double
from ctypes import byref
from ctypes import POINTER
from numpy import zeros
from numpy import array
from numpy import unique
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
import types as xctypes
from dsmxc import dsmxc
from fdxc import fdxc
from spxc import spxc
from utils import utils
from xclocTypes import xclocTypes as xctypes 
#import mpi4py

##  
# @defgroup pyxcloc pyxcloc
# @brief Interface to the xcloc library.
# @copyright Ben Baker distributed under the MIT license.
class xcloc:
    ##
    # @defgroup pyxcloc_xcloc xcloc
    # @brief Interface to the xcloc library.
    # @copyright Ben Baker distributed under the MIT license.
    # @ingroup pyxcloc 
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
        xcloc_lib.xcloc_initialize.argtypes = [c_int, # npts
                                               c_int, # nptsPad
                                               c_int, # nxcs
                                               c_int, # xcSignalType2migrate
                                               c_double, # dt
                                               c_int, # ngrd
                                               c_int, # nfcoeffs
                                               c_int, # ftype
                                               POINTER(c_int), # xcPairs
                                               c_int, # verbose
                                               c_int, # prec
                                               c_int, # accuracy
                                               POINTER(c_int)]
        xcloc_lib.xcloc_initialize.restype = None
        xcloc_lib.xcloc_signalToTableIndex.argtypes = [c_int, # is
                                                       POINTER(c_int), # it
                                                       POINTER(c_int)] 
        xcloc_lib.xcloc_signalToTableIndex.restype = None
        xcloc_lib.xcloc_setTable64f.argtypes = [c_int, # tableNumber
                                                c_int, # ngrd
                                                POINTER(c_double), # table
                                                POINTER(c_int)]
        xcloc_lib.xcloc_setTable64f.restype = None
        xcloc_lib.xcloc_setTable32f.argtypes = [c_int, # tableNumber
                                                c_int, # ngrd
                                                POINTER(c_float), # table
                                                POINTER(c_int)]
        xcloc_lib.xcloc_setTable32f.restype = None
        xcloc_lib.xcloc_setSignals64f.argtypes = [c_int, # ldx
                                                  c_int, # npts
                                                  c_int, # nsignals
                                                  POINTER(c_double), # signals
                                                  POINTER(c_int)]
        xcloc_lib.xcloc_setSignals64f.restype = None
        xcloc_lib.xcloc_getImageMax.argtypes = [POINTER(c_int), # Max index
                                                POINTER(c_float), # Max value
                                                POINTER(c_int)]
        xcloc_lib.xcloc_getImageMax.restype = None
        xcloc_lib.xcloc_setSignals32f.argtypes = [c_int, # ldx
                                                  c_int, # npts
                                                  c_int, # nsignals
                                                  POINTER(c_float), # signals
                                                  POINTER(c_int)]
        xcloc_lib.xcloc_setSignals32f.restype = None
        xcloc_lib.xcloc_getImage64f.argtypes = [c_int, #nwork
                                                POINTER(c_double), # image
                                                POINTER(c_int)]
        xcloc_lib.xcloc_getImage64f.restype = None
        xcloc_lib.xcloc_getImage32f.argtypes = [c_int, #nwork
                                                POINTER(c_float), # image
                                                POINTER(c_int)]
        xcloc_lib.xcloc_getImage32f.restype = None
        xcloc_lib.xcloc_getNumberOfGridPointsInImage.argtypes = [POINTER(c_int), # ngrd
                                                                 POINTER(c_int)]
        xcloc_lib.xcloc_getNumberOfGridPointsInImage.restype = None
        xcloc_lib.xcloc_compute.argtypes = [POINTER(c_int)]
        xcloc_lib.xcloc_compute.restype = None
        xcloc_lib.xcloc_finalize.argtypes = None
        xcloc_lib.xcloc_finalize.restype = None
        ##################################################################################
        #                                       xclocMPI                                 #
        ##################################################################################
        try:
            xcloc_lib.xclocMPI_initialize.argtypes = [c_int64, # MPI_Fint
                                                      c_int, # root
                                                      c_int, # dsmGroupsize
                                                      c_int, # npts
                                                      c_int, # nptsPad
                                                      c_int, # nxcs
                                                      c_int, # signal2migrate
                                                      c_double, # dt
                                                      c_int, # ngrd
                                                      c_int, # nfCoeffs
                                                      c_int, # filtertype
                                                      POINTER(c_int), # xcPairs
                                                      c_int,  # verbose
                                                      c_int,  # precision
                                                      c_int,  # accuracy
                                                      POINTER(c_int)]
            xcloc_lib.xclocMPI_initialize.restype = None
            xcloc_lib.xclocMPI_setXCTypeToMigrate.argtypes = [c_int, # root
                                                              c_int, # xcTypeToMigrate
                                                              POINTER(c_int)]
            xcloc_lib.xclocMPI_setXCTypeToMigrate.restype = None
            xcloc_lib.xclocMPI_getGlobalRootProcessID.argtypes = [c_int, # root
                                                                  POINTER(c_int)]
            xcloc_lib.xclocMPI_getGlobalRootProcessID.restype = None
            xcloc_lib.xclocMPI_finalize.argtypes = None
            xcloc_lib.xclocMPI_finalize.restype = None
            self.lhaveMPI = True
            print("yes")
        except:
            self.lhaveMPI = False

        # Hook up the modules
        self.utils = utils(xcloc_lib)
        self.fdxc = fdxc(xcloc_lib, self.utils)
        self.dsmxc = dsmxc(xcloc_lib)
        self.spxc = spxc(xcloc_lib)
        self.lib = xcloc_lib
        self.linit = False
        self.lhaveImage = False
        self.lhaveXCs = False

    def __exit__(self):
        self.fdxc.finalize()
        self.dsmxc.finalize()
        self.spxc.finalize()
        self.finalize()
        return

    def initialize(self, dt, npts, ngrd,
                   nsignals=None,
                   nptsPad=None,
                   xcPairs=None,
                   nfCoeff = 301,
                   ftype = xctypes.XCLOC_SPXC_ENVELOPE_FILTER,
                   xcSignalToMigrate=xctypes.XCLOC_MIGRATE_PHASE_XCS,
                   verbose=xctypes.XCLOC_PRINT_WARNINGS,
                   precision=xctypes.XCLOC_SINGLE_PRECISION,
                   accuracy=xctypes.XCLOC_HIGH_ACCURACY):
        """!
        @brief Initializes the xcloc module.
        @param dt                 Sampling period of signals (seconds).
        @param npts               Number of points in signals.
        @param nsignals           Number of signals.
        @param ngrd               Number of grid points in travel times fields.
        @param nptsPad            The correlogram length is 2*npts-1.  If this
                                  is a large semi-prime number then performance
                                  can degrade in the Fourier Transform.  nptsPad
                                  is an override that must be greater than or
                                  equal to npts to mitigate this. 
        @param ntCoeff            If the correlograms are to processed then this
                                  is the number of FIR coefficients.  This must
                                  be odd and positive.
        @param ftype              Typically the travel time model is inadequate.
                                  To mitigate destructive interference at the
                                  true source location it is advisable to, in 
                                  some way, smooth the correlograms and ensure
                                  that they are positive.  Strategies for this
                                  include computing the RMS or envelope of the
                                  correlograms.
        @param xcSignalToMigrate  This defines the signal to type to be
                                  migrated.  Typically the phase correlograms
                                  are migrated but it is also possible to
                                  migrate the cross correlograms.
        @param verbose            Controls the verbosity of the module.
        @param precision          Controls the floating precision model used
                                  in the calculations.
        @param accuracy           Controls the accuracy of the vectorized
                                  calculations.
        @retval 0 indicates success.
        @ingroup pyxcloc_xcloc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.initialize.__name__)
        # Input checks 
        if (npts < 1): 
            print("%s: npts=%d must be positive"%(fname, npts))
            return -1
        if (not (xcPairs is None)):
            nsignals = len(unique(xcPairs.flatten()))
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
        # Number of grid points in migration image must be positive
        if (ngrd < 1):
            print("%s: No grid points"%fname)
            return -1
        # Check signal migration type
        if (xcSignalToMigrate < 0 or xcSignalToMigrate > 1):
            print("%s: Invalid type of signal to migrate"%fname)
            return -1
        # Check the xc filtering options
        if (ftype < 0 or ftype > 2):
            print("%s: Invalid filtering type"%fname) 
            return -1
        if (ftype > 0):
            if (nfCoeff < 1):
                print("%s: No filtering coefficents"%fname)
                return -1
            if (nfCoeff%2 == 0):
                print("%s: Adding one filtering coefficient"%fname)
                nfCoeff = nfCoeff + 1

        nxcs = xcPairs.shape[0]
        xcPairs = xcPairs.flatten(order='C')
        xcPairsPtr = xcPairs.ctypes.data_as(POINTER(c_int))

        ierr = c_int(1) 
        self.lib.xcloc_initialize(npts, nptsPad, nxcs,
                                  xcSignalToMigrate,
                                  dt, ngrd, 
                                  nfCoeff, ftype,
                                  xcPairsPtr,
                                  verbose, precision, accuracy, ierr)
        ierr = ierr.value
        if (ierr != 0):
            print("%s: Failed to initialize xcloc"%fname)
            return -1
        self.linit = True
        return 0

    def compute(self):
        """!
        @brief Now that the travel time tables and signals are set this computes
               the diffraction stack image of the correlograms.
        @retval 0 indicates success.
        @ingroup pyxcloc_xcloc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.compute.__name__) 
        self.lhaveImage = False
        self.lhaveXCs = False
        if (not self.linit):
            print("%s: Module not initialized"%fname)
            return -1
        ierr = c_int(1)
        self.lib.xcloc_compute(ierr)
        ierr = ierr.value
        if (ierr != 0):
            print("%s: Error in xcloc calculation"%fname)
        self.lhaveImage = True
        self.lhaveXCs = True
        return ierr

    def setTable(self, signalNumber, table,
                 numbering=xctypes.XCLOC_C_NUMBERING):
        """!
        @brief Sets the travel time table corresponding to the signalNumber'th
               signal.
        @param signalNumber  Signal number.  The numbering either starts 0 or 1
                             depending on numbering. 
        @param table         The travel times from the signal'th receiver 
                             to all points in the medium.  This is an array
                             of dimension [ngrd].
        @param numbering     Defines the numbering.
        @retval 0 indicates success.
        @ingroup pyxcloc_xcloc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.setTable.__name__)
        if (not self.linit):
            print("%s: Module not initialized"%fname)
        i = signalNumber
        if (numbering == xctypes.XCLOC_C_NUMBERING): 
            i = i + 1
        it = c_int(1)
        ierr = c_int(1)
        self.lib.xcloc_signalToTableIndex(i, it, ierr) 
        it = it.value
        ngrd = table.size
        if (ierr.value != 0):
            print("%s: Could not find signal %d"%(fname, signalNumber)) 
            return -1
        ngrd = table.size
        table = table.flatten(order='C') 
        if (table.dtype == float32):
            table = ascontiguousarray(table, float32)
            tablePtr = table.ctypes.data_as(POINTER(c_float))
            self.lib.xcloc_setTable32f(it, ngrd, tablePtr, ierr)
            if (ierr.value != 0):
                print("%s: Failed to call setTable32f"%fname)
        elif (table.dtype == float64):
            table = ascontiguousarray(table, float64)
            tablePtr = table.ctypes.data_as(POINTER(c_double))
            self.lib.xcloc_setTable64f(it, ngrd, tablePtr, ierr)
            if (ierr.value != 0):
                print("%s: Failed to call setTable64f"%fname)
        else:
            print("%s: Precision must be float32 or float64"%fname)
            return -1
        ierr = ierr.value
        return ierr
         
 
    def setSignals(self, signals):
        """!
        @brief Sets the matrix of signals to correlate.

        @param signals  A matrix of input data with dimension [nsignals x npts].
        @retval 0 indicates success. 
        @ingroup pyxcloc_xcloc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.setSignals.__name__)
        if (not self.linit): 
            sys.exit("%s: xcloc is not initialized"%fname)
            return -1
        ierr = c_int(1)
        nsignals = signals.shape[0]
        npts = signals.shape[1]
        lds = npts # Leading dimension of fortran array
        signals = reshape(signals, signals.size, order='C')
        if (signals.dtype == float32):
            signals = ascontiguousarray(signals, float32)
            signalPtr = signals.ctypes.data_as(POINTER(c_float))
            self.lib.xcloc_setSignals32f(lds, npts, nsignals, signalPtr,
                                         byref(ierr))
            if (ierr.value != 0):
                print("%s: Failed to call setSignals32f"%fname)
        elif (signals.dtype == float64):
            signals = ascontiguousarray(signals, float64)
            signalPtr = signals.ctypes.data_as(POINTER(c_double))
            self.lib.xcloc_setSignals64f(lds, npts, nsignals, signalPtr,
                                         byref(ierr))
            if (ierr.value != 0):
                print("%s: Failed to call setSignals64f"%fname)
        else:
            print("%s: Precision must be float32 or float64"%fname)
            return -1
        ierr = ierr.value
        return ierr
   
    def getImage(self, dtype=float64):
        """!
        @brief Gets the DSM image.
        @retval An array containing the DSM image.
        @retval None indicates a failure.
        @ingroup pyxcloc_xcloc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.getImage.__name__)
        if (not self.lhaveImage):
            print("%s: DSM image not yet computed"%fname)
            return None
        ngrd = c_int(1)
        ierr = c_int(1)
        self.lib.xcloc_getNumberOfGridPointsInImage(ngrd, ierr)
        ngrd = ngrd.value
        if (ierr.value != 0):
            print("%s: Failed to get image size"%fname)
            return None
        if (dtype == float64):
            grd = ascontiguousarray(zeros(ngrd, dtype=c_double))
            grdPtr = grd.ctypes.data_as(POINTER(c_double)) 
            self.lib.xcloc_getImage64f(ngrd, grdPtr, ierr)
        else:
            grd = ascontiguousarray(zeros(ngrd, dtype=c_float))
            grdPtr = grd.ctypes.data_as(POINTER(c_float))
            self.lib.xcloc_getImage32f(ngrd, grdPtr, ierr)
        if (ierr.value != 0):
            print("%s: Failed to get image"%fname)
            return None
        return grd

    def getImageMax(self):
        """!
         @brief Gets the maximum value and index of the DSM image.
                The max index is most likely to correspond to the event
                location.
         @retval The max index (C numbered) and max value of the DSM image.
         @retval None, None indicates a failure.
         @ingroup pyxcloc_xcloc
        """
        fname = '%s::%s'%(self.__class__.__name__, self.getImageMax.__name__)
        if (not self.lhaveImage):
            print("%s: DSM image not yet computed"%fname)
            return None, None
        maxIndex = c_int(1)
        maxValue = c_float(1)
        ierr = c_int(1)
        self.lib.xcloc_getImageMax(maxIndex, maxValue, ierr)
        if (ierr.value != 0):
            print("%s: Failed to get max value and index of DSM image"%fname)
            return None, None
        imax = maxIndex.value - 1 # Fortran to C
        vmax = maxValue.value
        return imax, vmax

    def finalize(self):
        """!
        @brief Releases memory on the xcloc module.
        @ingroup pyxcloc_xcloc
        """
        self.lib.xcloc_finalize()
        self.linit = False
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
    xcloc.fdxc.initialize(npts, nsignals, precision=xctypes.XCLOC_DOUBLE_PRECISION)
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
    env = xcloc.spxc.filter(smat)
    #plt.plot(smat[2,:])
    #plt.plot(env[2,:])
    #plt.show()
    xcloc.spxc.finalize()
    # Clean up
    xcloc.fdxc.finalize()
    
