!> @brief Holds constants for Fortran modules.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_CONSTANTS
      USE ISO_C_BINDING
      !> Single precision.
      INTEGER(C_INT), PARAMETER :: XCLOC_SINGLE_PRECISION = 0
      !> Double precision.
      INTEGER(C_INT), PARAMETER :: XCLOC_DOUBLE_PRECISION = 1
      !> High accuracy MKL vector math calculations
      INTEGER(C_INT), PARAMETER :: XCLOC_HIGH_ACCURACY = 0
      !> Low accuracy MKL vector math calculations
      INTEGER(C_INT), PARAMETER :: XCLOC_LOW_ACCURACY = 1
      !> Really fast and inaccurate MKL vector math calculations.
      INTEGER(C_INT), PARAMETER :: XCLOC_EP_ACCURACY = 2
      !> Verbosity - print nothing.
      INTEGER(C_INT), PARAMETER :: XCLOC_PRINT_NONE =-1
      !> Verbosity - print errors only.
      INTEGER(C_INT), PARAMETER :: XCLOC_PRINT_ERRORS = 0
      !> Verbosity - print errors and warnings.
      INTEGER(C_INT), PARAMETER :: XCLOC_PRINT_WARNINGS = 1
      !> Verbosity - print errors, warnings, info.
      INTEGER(C_INT), PARAMETER :: XCLOC_PRINT_INFO = 2
      !> Verbosity - print everything.
      INTEGER(C_INT), PARAMETER :: XCLOC_PRINT_DEBUG = 3
      !> Do not do filter cross-correlgrams.
      INTEGER(C_INT), PARAMETER :: XCLOC_SPXC_DONOT_FILTER = 0
      !> Compute envelope of cross-correlograms.
      INTEGER(C_INT), PARAMETER :: XCLOC_SPXC_ENVELOPE_FILTER = 1
      !> Compute RMS of cross-correlograms.
      INTEGER(C_INT), PARAMETER :: XCLOC_SPXC_RMS_FILTER = 2
      !> Single complex 0.0
      COMPLEX(C_FLOAT_COMPLEX), PARAMETER :: czero = CMPLX(0.0, 0.0)
      !> Double complex 0.0
      COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: zzero = DCMPLX(0.d0, 0.d0)
      !> Default block size for migration loop.
      INTEGER(C_INT), PARAMETER :: XCLOC_DEFAULT_BLOCKSIZE = 1024
END MODULE
