!> Generic utilities to simplify using the library.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_UTILS
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
      IMPLICIT NONE

      PUBLIC :: xcloc_utils_computeDefaultXCTable
      CONTAINS
!========================================================================================!
!                                     Begin the Code                                     !
!========================================================================================!
!>    @brief Computes the default cross-correlation table pairs.  The pairs can be viewed
!>           as a symmetric matrix.
!>
!>    @param[in] ldoAutoCorrs  If true then the auto-correlations are to be performed in
!>                             addition to the cross-correlations; i.e., the main diagonal
!>                             of the cross-correlation matrix is computed.
!>    @param[in] ldoAutoCorrs  If false then only the cross-correlations are computed.
!>    @param[in] nsignals      Number of signals.
!>    @param[in] nwork         The size of xcPairs which should be equal to 2*nxcs.
!>    @param[in] nwork         If nwork is negative then this is a space query and
!>                             xcPairs will not be accessed.
!>    @param[in] numbering     If numbering is XCLOC_FORTRAN_NUMBERING then the signal
!>                             indices begin at 1.
!>    @param[in] numbering     If numbering is XCLOC_C_NUMBERING then the signal indices
!>                             begin at 0.
!>
!>    @param[out] nxcs         Number of cross-correlations.
!>    @param[out] xcPairs      If nwork > 0 then this contains the pairs such that the
!>                             i'th signal is to be correlated with the j'th signal.
!>                             This has dimension [nwork] but only the first 2*nxcs 
!>                             indices will be set.
!>    @param[out] ierr         0 indicates success.
!>
      SUBROUTINE xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, &
                                                   nwork, numbering,       &
                                                   nxcs, xcPairs, ierr)    &
      BIND(C, NAME='xcloc_utils_computeDefaultXCTable')
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: ldoAutoCorrs
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsignals, numbering, nwork
      INTEGER(C_INT), INTENT(OUT) :: nxcs, xcPairs(nwork), ierr
      INTEGER i, j, indx
      ierr = 0
      nxcs = 0
      IF (nsignals < 1) THEN
         WRITE(*,900)
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.ldoAutoCorrs .AND. nsignals < 2) THEN
         WRITE(*,901)
         ierr = 1
         RETURN
      ENDIF
      ! Only the superdiagonal
      IF (.NOT.ldoAutoCorrs) THEN
         nxcs = (nsignals*(nsignals - 1))/2
         IF (nwork < 0) RETURN ! space inquiry
         xcPairs(:) = 0
         IF (nwork < 2*nxcs) THEN
            WRITE(*,905) 2*nxcs
            ierr = 1
            RETURN
         ENDIF
         DO i=1,nsignals
            DO j=i+1,nsignals
               indx = nsignals*(i - 1) - (i*(i+1))/2 + j
               xcPairs(2*(indx-1)+1) = i
               xcPairs(2*(indx-1)+2) = j
            ENDDO
         ENDDO
      ! Upper triangle (including diagonal)
      ELSE
         nxcs = (nsignals*(nsignals + 1))/2
         IF (nwork < 0) RETURN ! space inquiry
         xcPairs(:) = 0
         IF (nwork < 2*nxcs) THEN
            WRITE(*,905) 2*nxcs
            ierr = 1
            RETURN
         ENDIF
         DO i=1,nsignals
            DO j=i,nsignals
               indx = nsignals*(i - 1) - ((i-1)*i)/2 + j
               xcPairs(2*(indx-1)+1) = i
               xcPairs(2*(indx-1)+2) = j
            ENDDO
         ENDDO
      ENDIF
      IF (numbering == XCLOC_C_NUMBERING) xcPairs(:) = xcPairs(:) - 1
  900 FORMAT('xcloc_utils_computeDefaultXCTable: No signals')
  901 FORMAT('xcloc_utils_computeDefaultXCTable: At least 2 signals required')
  905 FORMAT('xcloc_utils_computeDefaultXCTable: Error - nwork must be >=', I5)
      RETURN
      END
END MODULE
