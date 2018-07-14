#!/usr/bin/env python3

class xclocTypes:
    ##  
    # @defgroup xclocTypes Types
    # @brief Defines types in the xcloc module. 
    # @ingroup pyxcloc
    # @copyright Ben Baker distributed under the MIT license.
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

    (
    XCLOC_PRINT_NONE,     # Print nothing.
    XCLOC_PRINT_ERRORS,   # Print errors.
    XCLOC_PRINT_WARNINGS, # Print errors and warnings.
    XCLOC_PRINT_INFO,     # Print errors, warnings, and general information.
    XCLOC_PRINT_DEBUG     # Print everything.
    ) = map(int, range(-1,4))
