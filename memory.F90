!> @defgroup memory Memory Utilities
!> @ingroup xcloc
!> @brief Memory handling routines.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_MEMORY
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
      IMPLICIT NONE
      PUBLIC :: xcloc_memory_padLength
      CONTAINS
      !==================================================================================!
      !                                      Begin the Code                              !
      !==================================================================================!
      !> @ingroup memory
      INTEGER(C_INT) FUNCTION xcloc_memory_padLength(alignment, sizeof_dataType, n) &
      RESULT(padLength) &
      BIND(C, NAME='xcloc_meomry_padLength')
      USE ISO_C_BINDING
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: alignment, sizeof_dataType
      INTEGER(C_INT), VALUE, INTENT(IN) :: n 
      INTEGER(C_INT) xmod 
      padLength = 0  
      xmod = MOD(n*INT(sizeof_dataType), INT(alignment))
      IF (xmod /= 0) padLength = (INT(alignment) - xmod)/INT(sizeof_dataType)
      padLength = n + padLength
      RETURN
      END FUNCTION
END MODULE
