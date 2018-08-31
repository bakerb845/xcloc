!> @defgroup constants Constants
!> @ingroup xcloc
!> @ingroup xcloc_mpi
!> @brief Holds constants for Fortran modules.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_CONSTANTS
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      IMPLICIT NONE
      !> @ingroup constants
      !> Single precision.
      INTEGER(C_INT), PARAMETER :: XCLOC_SINGLE_PRECISION = 0
      !> @ingroup constants
      !> Double precision.
      INTEGER(C_INT), PARAMETER :: XCLOC_DOUBLE_PRECISION = 1
      !> @ingroup constants
      !> C numbering.
      INTEGER(C_INT), PARAMETER :: XCLOC_C_NUMBERING = 0 
      !> @ingroup constants
      !> Fortran numbering.
      INTEGER(C_INT), PARAMETER :: XCLOC_FORTRAN_NUMBERING = 1
      !> @ingroup constants
      !> High accuracy MKL vector math calculations
      INTEGER(C_INT), PARAMETER :: XCLOC_HIGH_ACCURACY = 0
      !> @ingroup constants
      !> Low accuracy MKL vector math calculations
      INTEGER(C_INT), PARAMETER :: XCLOC_LOW_ACCURACY = 1
      !> @ingroup constants
      !> Really fast and inaccurate MKL vector math calculations.
      INTEGER(C_INT), PARAMETER :: XCLOC_EP_ACCURACY = 2
      !> @ingroup constants
      !> Verbosity - print nothing.
      INTEGER(C_INT), PARAMETER :: XCLOC_PRINT_NONE =-1
      !> @ingroup constants
      !> Verbosity - print errors only.
      INTEGER(C_INT), PARAMETER :: XCLOC_PRINT_ERRORS = 0
      !> @ingroup constants
      !> Verbosity - print errors and warnings.
      INTEGER(C_INT), PARAMETER :: XCLOC_PRINT_WARNINGS = 1
      !> @ingroup constants
      !> Verbosity - print errors, warnings, info.
      INTEGER(C_INT), PARAMETER :: XCLOC_PRINT_INFO = 2
      !> @ingroup constants
      !> Verbosity - print everything.
      INTEGER(C_INT), PARAMETER :: XCLOC_PRINT_DEBUG = 3
      !> @ingroup constants
      !> Do not do filter cross-correlgrams.
      INTEGER(C_INT), PARAMETER :: XCLOC_SPXC_DONOT_FILTER = 0
      !> @ingroup constants
      !> Compute envelope of cross-correlograms.
      INTEGER(C_INT), PARAMETER :: XCLOC_SPXC_ENVELOPE_FILTER = 1
      !> @ingroup constants
      !> Compute RMS of cross-correlograms.
      INTEGER(C_INT), PARAMETER :: XCLOC_SPXC_RMS_FILTER = 2
      !> @ingroup constants
      !> Compute the diffraction stack image of phase correlograms.
      INTEGER(C_INT), PARAMETER :: XCLOC_MIGRATE_PHASE_XCS = 0
      !> @ingroup constants
      !> Compute the diffraction stack image of the cross correlograms.
      INTEGER(C_INT), PARAMETER :: XCLOC_MIGRATE_XCS = 1
      !> @ingroup constants
      !> Single complex 0.0
      COMPLEX(C_FLOAT_COMPLEX), PARAMETER :: czero = CMPLX(0.0, 0.0, KIND=4)
      !> @ingroup constants
      !> Double complex 0.0
      COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: zzero = CMPLX(0.d0, 0.d0, KIND=8)
      !> @ingroup constants
      !> Default block size for migration loop.
      INTEGER(C_INT), PARAMETER :: XCLOC_DEFAULT_BLOCKSIZE = 1024
      !----------------------------------------------------------------------------------!
      !                             Public/Private Subroutines                           !
      !----------------------------------------------------------------------------------!
      PUBLIC :: xcloc_constants_isValidPrecision
      PUBLIC :: xcloc_constants_isValidAccuracy
      PUBLIC :: xcloc_constants_isValidSignalToMigrate
      PUBLIC :: xcloc_constants_isValidNumbering
      PUBLIC :: xcloc_constants_isValidFilteringType
      CONTAINS
!========================================================================================!
!                                         Begin the Code                                 !
!========================================================================================!
!>    @brief Determines if the precision is supported.
!>    @param[in] prec   Precision to test.
!>    @result True if the precision is supported.
!>    @result False if the precision is not supported.
!>    @ingroup constants
      LOGICAL(C_BOOL) FUNCTION xcloc_constants_isValidPrecision(prec) &
      RESULT(isValid) BIND(C, NAME='xcloc_constants_isValidPrecision')
      INTEGER(C_INT), VALUE, INTENT(IN) :: prec
      isValid = .TRUE.
      IF (prec /= XCLOC_SINGLE_PRECISION .AND. prec /= XCLOC_DOUBLE_PRECISION) THEN
         WRITE(ERROR_UNIT,900) prec, XCLOC_SINGLE_PRECISION, XCLOC_DOUBLE_PRECISION
         isValid = .FALSE.
      ENDIF
  900 FORMAT("xcloc_constants_isValidPrecision: prec=", I4, " must be single=", I2, &
             "or double=", I2)
      RETURN
      END
!>    @brief Determines if the accuracy is supported.
!>    @param[in] acc   Accuracy to test.
!>    @result True if the accuracy is supported.
!>    @result False if the accuracy is not supported.
!>    @ingroup constants
      LOGICAL(C_BOOL) FUNCTION xcloc_constants_isValidAccuracy(acc) &
      RESULT(isValid) BIND(C, NAME='xcloc_constants_isValidAccuracy')
      INTEGER(C_INT), VALUE, INTENT(IN) :: acc 
      isValid = .TRUE.
      IF (acc /= XCLOC_HIGH_ACCURACY .AND. &
          acc /= XCLOC_LOW_ACCURACY  .AND. &
          acc /= XCLOC_EP_ACCURACY) THEN
         WRITE(ERROR_UNIT,900) acc, XCLOC_HIGH_ACCURACY, &
                               XCLOC_LOW_ACCURACY, XCLOC_EP_ACCURACY
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
!>    @ingroup constants
      LOGICAL(C_BOOL) FUNCTION xcloc_constants_isValidNumbering(num) &
      RESULT(isValid) BIND(C, NAME='xcloc_constants_isValidNumbering')
      INTEGER(C_INT), VALUE, INTENT(IN) :: num
      isValid = .TRUE.
      IF (num /= XCLOC_C_NUMBERING .AND. num /= XCLOC_FORTRAN_NUMBERING) THEN
         WRITE(ERROR_UNIT,900) num, XCLOC_C_NUMBERING, XCLOC_FORTRAN_NUMBERING
         isValid = .FALSE.
      ENDIF
  900 FORMAT("xcloc_constants_isValidNumbering: num=", I4, " must be C=", I2, &
             "or Fortran=", I2) 
      RETURN
      END
!>    @brief Determines if the signal to migrate is supported.
!>    @param[in] s2m   Sigmnal to migrate.
!>    @result True if the type of signal to migrate is supported.
!>    @result False if the type of signal to migrate is not supported.
!>    @ingroup constants
      LOGICAL(C_BOOL) FUNCTION xcloc_constants_isValidSignalToMigrate(s2m) &
      RESULT(isValid) BIND(C, NAME='xcloc_constants_isValidSignalToMigrate')
      INTEGER(C_INT), VALUE, INTENT(IN) :: s2m
      isValid = .TRUE.
      IF (s2m /= XCLOC_MIGRATE_PHASE_XCS .AND. s2m /= XCLOC_MIGRATE_XCS) THEN
         WRITE(ERROR_UNIT,900) s2m, XCLOC_MIGRATE_PHASE_XCS, XCLOC_MIGRATE_XCS
         isValid = .FALSE.
      ENDIF
  900 FORMAT("xcloc_constants_isValidSignalToMigrate: s2m=", I4, " must be =", I2, &
             "or =", I2) 
      RETURN
      END
!>    @brief Determines if the filtering type is supported.
!>    @param[in] ftype
!>    @result True if the type of filtering is supported.
!>    @result False if the  type of filtering is not suported.
!>    @ingroup constants
      LOGICAL(C_BOOL) FUNCTION xcloc_constants_isValidFilteringType(ftype) &
      RESULT(isValid) BIND(C, NAME='xcloc_constants_isValidFilteringType')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ftype
      isValid = .TRUE.
      IF (ftype /= XCLOC_SPXC_DONOT_FILTER .AND. &
          ftype /= XCLOC_SPXC_ENVELOPE_FILTER .AND. &
          ftype /= XCLOC_SPXC_RMS_FILTER) THEN
         WRITE(ERROR_UNIT,900) ftype, XCLOC_SPXC_DONOT_FILTER, &
                               XCLOC_SPXC_ENVELOPE_FILTER, XCLOC_SPXC_RMS_FILTER
         isValid = .FALSE.
      ENDIF
  900 FORMAT("xcloc_constants_isValidFilteringType: ftype=", I4, &
             " must be envelope=", I2, &
             " or rms=", I2, &
             " or none=", I2)
      RETURN
      END
END MODULE
