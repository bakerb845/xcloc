!> @brief Some Fortran interfaces to IPP functions.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_IPPS
  INTERFACE
      !> @brief Converts a double array to a float array.
      !> @param[out] pSrc  Double precision source array.
      !> @param[in] pDst   Float precision copy of pSrc.
      !> @param[in] len    Length of arrays.
      INTEGER(C_INT) FUNCTION ippsConvert_64f32f(pSrc, pDst, len) &
      BIND(C, NAME='ippsConvert_64f32f')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: len
      REAL(C_DOUBLE), INTENT(IN) :: pSrc(len)
      REAL(C_FLOAT), INTENT(OUT) :: pDst(len)
      END FUNCTION
      !> @brief Converts a float array to a double array.
      !> @param[out] pSrc  Float precision source array.
      !> @param[in] pDst   Double precision copy of pSrc.
      !> @param[in] len    Length of arrays.
      INTEGER(C_INT) FUNCTION ippsConvert_32f64f(pSrc, pDst, len) &
      BIND(C, NAME='ippsConvert_64f32f')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: len 
      REAL(C_FLOAT), INTENT(IN) :: pSrc(len)
      REAL(C_DOUBLE), INTENT(OUT) :: pDst(len)
      END FUNCTION

   END INTERFACE
END MODULE
