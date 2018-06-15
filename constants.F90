!> @brief Holds constants for Fortran modules.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_CONSTANTS
      USE ISO_C_BINDING
      !> Single precision.
      INTEGER(C_INT), PARAMETER :: XCLOC_SINGLE_PRECISION = 0
      !> Double precision.
      INTEGER(C_INT), PARAMETER :: XCLOC_DOUBLE_PRECISION = 1
      !> C numbering.
      INTEGER(C_INT), PARAMETER :: XCLOC_C_NUMBERING = 0 
      !> Fortran numbering.
      INTEGER(C_INT), PARAMETER :: XCLOC_FORTRAN_NUMBERING = 1
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
      !----------------------------------------------------------------------------------!
      !                             Public/Private Subroutines                           !
      !----------------------------------------------------------------------------------!
      PUBLIC :: xcloc_constants_isValidPrecision
      PUBLIC :: xcloc_constants_isValidAccuracy
      CONTAINS
!========================================================================================!
!                                         Begin the Code                                 !
!========================================================================================!
!>    @brief Determines if the precision is supported.
!>    @param[in] prec   Precision to test.
!>    @result True if the precision is supported.
!>    @result False if the precision is not supported.
      LOGICAL(C_BOOL) FUNCTION xcloc_constants_isValidPrecision(prec) &
      RESULT(isValid) BIND(C, NAME='xcloc_constants_isValidPrecision')
      INTEGER(C_INT), VALUE, INTENT(IN) :: prec
      isValid = .TRUE.
      IF (prec /= XCLOC_SINGLE_PRECISION .AND. prec /= XCLOC_DOUBLE_PRECISION) THEN
         WRITE(*,900) prec, XCLOC_SINGLE_PRECISION, XCLOC_DOUBLE_PRECISION
         isValid = .FALSE.
      ENDIF
  900 FORMAT("xcloc_constants_isValidPrecision: prec=", I4, " must be single=", I2, &
             "or double=", I2)
      RETURN
      END
!>    @brief Determines if the accuracy is supported.
!>    @param[in] accuracy   Accuracy to test.
!>    @result True if the accuracy is supported.
!>    @result False if the accuracy is not supported.
      LOGICAL(C_BOOL) FUNCTION xcloc_constants_isValidAccuracy(acc) &
      RESULT(isValid) BIND(C, NAME='xcloc_constants_isValidAccuracy')
      INTEGER(C_INT), VALUE, INTENT(IN) :: acc 
      isValid = .TRUE.
      IF (acc /= XCLOC_HIGH_ACCURACY .AND. &
          acc /= XCLOC_LOW_ACCURACY  .AND. &
          acc /= XCLOC_EP_ACCURACY) THEN
         WRITE(*,900) acc, XCLOC_HIGH_ACCURACY, XCLOC_LOW_ACCURACY, XCLOC_EP_ACCURACY
         isValid = .FALSE.
      ENDIF
  900 FORMAT("xcloc_constants_isValidAccuracy: accuracy=", I4, " must be highe=", I2, &
             "or low=", I2, "or extended precision=", I2)
      RETURN
      END
!>    @brief Determines if the numbering is supported.
!>    @param[in] num   Numbering to test.
!>    @result True if the numbering is supported.
!>    @result False if the numbering is not supported.
      LOGICAL(C_BOOL) FUNCTION xcloc_constants_isValidNumbering(num) &
      RESULT(isValid) BIND(C, NAME='xcloc_constants_isValidNumbering')
      INTEGER(C_INT), VALUE, INTENT(IN) :: num
      isValid = .TRUE.
      IF (num /= XCLOC_C_NUMBERING .AND. num /= XCLOC_FORTRAN_NUMBERING) THEN
         WRITE(*,900) num, XCLOC_C_NUMBERING, XCLOC_FORTRAN_NUMBERING
         isValid = .FALSE.
      ENDIF
  900 FORMAT("xcloc_constants_isValidNumbering: num=", I4, " must be C=", I2, &
             "or Fortran=", I2) 
      RETURN
      END
END MODULE
