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
      !> Single complex 0.0
      COMPLEX(C_FLOAT_COMPLEX), PARAMETER :: czero = CMPLX(0.0, 0.0)
      !> Double complex 0.0
      COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: zzero = DCMPLX(0.d0, 0.d0)
END MODULE
