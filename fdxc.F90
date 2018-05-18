!> @brief Computes the cross-correlograms via the Fourier transform.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_FDXC
      USE ISO_C_BINDING

      !----------------------------------------------------------------------------------!
      !                                 Private Variables                                !
      !----------------------------------------------------------------------------------!
      !> This is a map from the ixc'th correlation to the signals pairs to be computed.
      !> It is column major matrix with dimension [2 x nxcs_].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: xcPairs_(:) 
      !> Length of the input signals.
      INTEGER, PRIVATE, SAVE :: npts_ = 0
      !> The number of signals to correlate. 
      INTEGER, PRIVATE, SAVE :: nsignals_ = 0
      !> The length of the signals to transform.  This can mitigate the pathologic
      !> case where the signal transform lengths are large (semi)prime numbers which
      !> make the DFT very expensive.
      INTEGER, PRIVATE, SAVE :: nptsPad_ = 0
      !> The number of cross-correlations.
      INTEGER, PRIVATE, SAVE :: nxcs_ = 0
      !> Controls verbosity where 0 will report on errors only.
      INTEGER, PRIVATE, SAVE :: verbose_ = 0
      !> If true then the cross-correlation pairs table (xcPairs) is set.
      LOGICAL, PRIVATE, SAVE :: lhaveTable_ = .FALSE.
      !----------------------------------------------------------------------------------!
      !                            Public/Private Subroutines                            !
      !----------------------------------------------------------------------------------!
      PUBLIC :: xcloc_fdxc_initialize
      PUBLIC :: xcloc_fdxc_finalize
      PUBLIC :: xcloc_fdxc_setXCTableF
      PUBLIC :: xcloc_fdxc_computeDefaultXCTableF
      PRIVATE :: initializeFFTW
      CONTAINS
!========================================================================================!
!                                     Begin the Code                                     !
!========================================================================================!
!>    @brief Initializes the Fourier transforms for the nsignals each of length npts.
!>
!>    @param[in] npts      Number of points in each input signal.
!>    @param[in] nsignals  Number of signals which will be set on the module.
!>    @param[in] nptsPad   A tuning parameter to mitigate DFT lengths that could 
!>                         potentially be large semi-prime numbers. 
!>    @param[in] verbose   Controls the verbosity of the module.  0 is quiet.
!>
!>    @param[out] ierr     0 indicates success.
!>
      SUBROUTINE xcloc_fdxc_initialize(npts, nsignals, nptsPad, verbose, ierr) &
      BIND(C, NAME='xcloc_fdxc_initialize')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: npts, nptsPad, nsignals, verbose
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      CALL xcloc_fdxc_finalize()
      IF (npts < 1 .OR. nsignals < 2 .OR. nptsPad < npts) THEN
         IF (npts < 1) WRITE(*,905) npts
         IF (nsignals < 2) WRITE(*,906) nsignals
         IF (nptsPad < npts) WRITE(*,907) nptsPad, npts
         ierr = 1
         RETURN
      ENDIF
      ! Set the input variables
      npts_     = npts
      nptsPad_  = nptsPad
      nsignals_ = nsignals
      verbose_  = verbose
      ! Format statements
  905 FORMAT("xcloc_fdxc_initialize: npts=", I8, "must be positive")
  906 FORMAT("xcloc_fdxc_initialize: nsignals=", I8, "must be at least 2")
  907 FORMAT("xcloc_fdxc_initialize: nptsPad=", I8, "must be greater than nts=", I8)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases the memory in the module and reset the variables.
!>
!>    @copyright Ben Baker distributed under the MIT license.
      SUBROUTINE xcloc_fdxc_finalize( ) & 
      BIND(C, NAME='xcloc_fdxc_finalize')
      IMPLICIT NONE
      IF (ALLOCATED(xcPairs_)) DEALLOCATE(xcPairs_)
      npts_ = 0
      nptsPad_ = 0
      nsignals_ = 0
      verbose_ = 0
      lhaveTable_ = .FALSE.
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the cross-correlation table pairs.
!>    @param[in] nxcs     Number of cross-correlations.
!>    @param[in] xcPairs  This is a [2 x nxcs] matrix in column major format where
!>                        the indices, (2*(ixc-1)+1, 2*(ixc-1)+2), map to the 
!>                        (i,j)'th signal pair comprising a correlation.
!>    @param[out] ierr    0 indicates success 
!>
      SUBROUTINE xcloc_fdxc_setXCTableF(nxcs, xcPairs, ierr) &
      BIND(C, NAME='xcloc_fdxc_setXCTableF')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: nxcs
      INTEGER(C_INT), INTENT(IN) :: xcPairs(2*nxcs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 1
      lhaveTable_ = .FALSE.
      IF (nxcs < 1) THEN
         WRITE(*, 900) nxcs
         RETURN
      ENDIF
      IF (MINVAL(xcPairs) < 1 .OR. MAXVAL(xcPairs) > nsignals_) THEN
         IF (MINVAL(xcPairs) < 1) WRITE(*,901) MINVAL(xcPairs)
         IF (MAXVAL(xcPairs) > nsignals_) WRITE(*,902) MAXVAL(xcPairs), nsignals_
         RETURN
      ENDIF
      ! Allocate space and copy table
      ierr = 0
      IF (ALLOCATED(xcPairs_)) DEALLOCATE(xcPairs_)
      ALLOCATE(xcPairs_(2*nxcs))
      xcPairs_(1:2*nxcs) = xcPairs(1:2*nxcs) 
      lhaveTable_ = .TRUE.
      ! Format statements
  900 FORMAT('xcloc_fdxc_setXCTableF: Error nxcs must be positive', I5)
  901 FORMAT('xcloc_fdxc_setXCTableF: minval(xcPairs)', I4, 'must be positive')
  902 FORMAT('xcloc_fdxc_setXCTableF: minval(xcPairs)', I4, 'cannot exceed', I4)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the default cross-correlation table pairs.  The pairs can be viewed
!>           as a symmetric matrix.
!>
!>    @param[in] ldoAutoCorrs  If true then the auto-correlations are to be performed in
!>                             addition to the cross-correlations; i.e., the main diagonal
!>                             of the cross-correlation matrix is computed.
!>    @param[in] ldoAutoCorrs  If false then only the cross-correlations are computed.
!>    @param[in] nwork         The size of xcPairs which should be equal to 2*nxcs.
!>    @param[in] nwork         If nwork is negative then this is a space query and
!>                             xcPairs will not be accessed.
!>
!>    @param[out] nxcs         Number of cross-correlations.
!>    @param[out] xcPairs      If nwork > 0 then this contains the pairs such that the
!>                             i'th signal is to be correlated with the j'th signal.
!>                             This has dimension [nwork] but only the first 2*nxcs 
!>                             indices will be set.
!>    @param[out] ierr         0 indicates success.
!> 
      SUBROUTINE xcloc_fdxc_computeDefaultXCTableF(ldoAutoCorrs, nwork,   &
                                                    nxcs, xcPairs, ierr) &
      BIND(C, NAME='xcloc_fdxc_computeDefaultXCTableF')
      IMPLICIT NONE
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: ldoAutoCorrs
      INTEGER(C_INT), VALUE, INTENT(IN) :: nwork
      INTEGER(C_INT), INTENT(OUT) :: nxcs, xcPairs(nwork), ierr
      INTEGER i, j, indx
      ierr = 0
      nxcs = 0
      IF (nsignals_ < 1) THEN
         WRITE(*,*) 'xcloc_fdxc_computeDefaultXCTableF: No signals!'
         ierr = 1
         RETURN
      ENDIF
      ! Only the superdiagonal
      IF (.NOT.ldoAutoCorrs) THEN
         nxcs = (nsignals_*(nsignals_ - 1))/2
         IF (nwork < 0) RETURN ! space inquiry
         xcPairs(:) = 0
         IF (nwork < 2*nxcs) THEN
            WRITE(905,*) 2*nxcs
            ierr = 1
            RETURN
         ENDIF
         DO i=1,nsignals_
            DO j=i+1,nsignals_
               indx = nsignals_*(i - 1) - (i*(i+1))/2 + j;
               xcPairs(2*(indx-1)  ) = i
               xcPairs(2*(indx-1)+1) = j
            ENDDO
         ENDDO
      ! Upper triangle (including diagonal)
      ELSE
         nxcs = (nsignals_*(nsignals_ + 1))/2
         IF (nwork < 0) RETURN ! space inquiry
         xcPairs(:) = 0
         IF (nwork < 2*nxcs) THEN
            WRITE(905,*) 2*nxcs
            ierr = 1
            RETURN
         ENDIF
         DO i=1,nsignals_
            DO j=1,nsignals_
               indx = nsignals_*(i - 1) - ((i-1)*i)/2 + j
               xcPairs(2*(indx-1)  ) = i
               xcPairs(2*(indx-1)+1) = j
            ENDDO
         ENDDO
      ENDIF
 905  FORMAT('xcloc_fdxc_computeDefaultXCTableF: Error - nwork must be >=', I5)
      RETURN
      END
!========================================================================================!
      SUBROUTINE initializeFFTW(ierr)
      INCLUDE 'fftw/fftw3.f03'
      INTEGER, INTENT(OUT) :: ierr
      ierr = 0
      RETURN
      END

END MODULE
