!> @brief Memory handling routines.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_MEMORY
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
      PUBLIC :: xcloc_memory_padLength
      CONTAINS
      !==================================================================================!
      !                                      Begin the Code                              !
      !==================================================================================!
      INTEGER(C_INT) FUNCTION xcloc_memory_padLength(alignment, sizeof_dataType, n) &
      RESULT(padLength) &
      BIND(C, NAME='xcloc_meomry_padLength')
      USE ISO_C_BINDING
      IMPLICIT NONE 
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: alignment, sizeof_dataType
      INTEGER(C_INT), VALUE, INTENT(IN) :: n 
      INTEGER(C_INT) xmod 
      padLength = 0  
      xmod = MOD(n*sizeof_dataType, INT(alignment))
      IF (xmod /= 0) padLength = (INT(alignment) - xmod)/sizeof_dataType
      padLength = n + padLength
      RETURN
      END FUNCTION
END MODULE
