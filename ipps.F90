!> @brief Some Fortran interfaces to IPP functions.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_IPPS
  INTEGER, PARAMETER :: ipp32f = 13
  INTEGER, PARAMETER :: ipp64f = 19
  INTEGER, PARAMETER :: ippStsNoErr = 0
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
      !> @brief Gets the single-rate FIR filter state size.
      !> @param[in] tapsLen   Number of taps.
      !> @param[out] pSpecSize  Size of the internal constant specification structure.
      !> @param[out] pBufSize   Size of the work buffer required for FIR filtering.
      !> @result 0 indicates success.
      INTEGER(C_INT) FUNCTION ippsFIRSRGetSize_finter64f(tapsLen, pSpecSize, pBufSize) &
      BIND(C, NAME='ippsFIRSRGetSize_finter64f')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: tapsLen
      INTEGER(C_INT), INTENT(OUT) :: pSpecSize, pBufSize
      END FUNCTION
      !> @brief Gets the single-rate FIR filter state size.
      !> @param[in] tapsLen   Number of taps.
      !> @param[out] pSpecSize  Size of the internal constant specification structure.
      !> @param[out] pBufSize   Size of the work buffer required for FIR filtering.
      !> @result 0 indicates success.
      INTEGER(C_INT) FUNCTION ippsFIRSRGetSize_finter32f(tapsLen, pSpecSize, pBufSize) &
      BIND(C, NAME='ippsFIRSRGetSize_finter32f')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: tapsLen
      INTEGER(C_INT), INTENT(OUT) :: pSpecSize, pBufSize
      END FUNCTION
      !> @brief Comptues a Kaiser window
      !> @param[in,out] kaiser  On input contains the points at which to compute the
      !>                        Kaiser window.
      !> @param[in,out] kaiser  On exit contains this is the Kaiser window.
      !> @param[in] n           Number of elements in n.
      !> @param[in] alpha       Controls shape of Kaiser window. 
      !> @result 0 indicates success.
      INTEGER(C_INT) FUNCTION ippsWinKaiser_64f_I(pSrcDst, n, alpha) &
      BIND(C, NAME='ippsWinKaiser_64f_I')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: n 
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: alpha
      REAL(C_DOUBLE), INTENT(INOUT) :: pSrcDst(n) 
      END FUNCTION
      !> @brief Comptues a Kaiser window
      !> @param[in,out] kaiser  On input contains the points at which to compute the
      !>                        Kaiser window.
      !> @param[in,out] kaiser  On exit contains this is the Kaiser window.
      !> @param[in] n           Number of elements in n.
      !> @param[in] alpha       Controls shape of Kaiser window. 
      !> @result 0 indicates success.
      INTEGER(C_INT) FUNCTION ippsWinKaiser_32f_I(pSrcDst, n, alpha) &
      BIND(C, NAME='ippsWinKaiser_32f_I')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: n
      REAL(C_FLOAT), VALUE, INTENT(IN) :: alpha
      REAL(C_FLOAT), INTENT(INOUT) :: pSrcDst(n)
      END FUNCTION

   END INTERFACE
END MODULE
