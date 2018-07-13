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
import os
import types as xctypes
from dsmxc import dsmxc
from fdxc import fdxc
from spxc import spxc
from utils import utils
from xclocTypes import xclocTypes as xctypes 
#import mpi4py

class xcloc:
    ##
    # @class pyxcloc pyxcloc 
    # @defgroup pyxcloc pyxcloc
    # @brief Interface to the xcloc library.
    # @copyright Ben Baker distributed under the MIT license.
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
        self.dsmxc.finalize()
        self.spxc.finalize()
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
    
