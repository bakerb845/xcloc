!> @brief Computes the cross-correlograms via the Fourier transform.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_FDXC
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
      USE XCLOC_IPPS
      IMPLICIT NONE
      !----------------------------------------------------------------------------------!
      !                                 Private Variables                                !
      !----------------------------------------------------------------------------------!
      !> This is a map from the ixc'th correlation to the signals pairs to be computed.
      !> It is column major matrix with dimension [2 x nxcs_].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: xcPairs_(:) 
      !> Length of the input signals.
      INTEGER, PRIVATE, SAVE :: npts_ = 0
      !> Length of the time domain cross-correlograms.
      INTEGER, PRIVATE, SAVE :: nptsInXCs_ = 0
      !> Length of the Fourier transforms.
      INTEGER, PRIVATE, SAVE :: nptsInFTs_ = 0 
      !> The number of signals to correlate. 
      INTEGER, PRIVATE, SAVE :: nsignals_ = 0
      !> The padding in the data.  This will be greater than or equal to nptsInXCs_.
      INTEGER, PRIVATE, SAVE :: dataOffset_ = 0
      !> The padding in the FTs.  This will be greater than or equal to nptsInFTs_. 
      INTEGER, PRIVATE, SAVE :: ftOffset_ = 0
      !> Precision of module.
      INTEGER, PRIVATE, SAVE :: precision_ = XCLOC_SINGLE_PRECISION
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
      !> Holds the cross-correlograms.  This is an array of dimension 
      !> [dataOffset_ x nxcs_].
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, SAVE :: xcs64f_(:)
      !> Holds the cross-correlograms.  This is an array of dimension 
      !> [dataOffset_ x nxcs_].
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, SAVE :: xcs32f_(:)
      !> Holds the input signals to transform.  This is an array of dimension 
      !> [dataOffset_ x nsignals_].
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, SAVE :: inputSignals64f_(:)
      !> Holds the input signals to transform.  This is an array of dimension
      !> [dataOffset_ x nsignals].
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, SAVE :: inputSignals32f_(:)
      !> Holds the Fourier transformed input signals.  This is an array of dimension
      !> [ftOffset_ x nsignals_].  
      COMPLEX(C_FLOAT_COMPLEX),  PRIVATE, ALLOCATABLE, SAVE :: inputFTs32f_(:)
      !> Holds the Fourier transformed input signals.  This is an array of dimension
      !> [ftOffset_ x nsignals_].
      COMPLEX(C_DOUBLE_COMPLEX), PRIVATE, ALLOCATABLE, SAVE :: inputFTs64f_(:)
      !> Holds the Fourier transformed cross-correlograms.  This is an array
      !> of dimension [ftOffset x nxcs_].
      COMPLEX(C_FLOAT_COMPLEX),  TARGET, PRIVATE, ALLOCATABLE, SAVE :: xcFTs32f_(:)
      !> Holds the Fourier transformed cross-correlograms.  This is an array
      !> of dimension [ftOffset x nxcs_].
      COMPLEX(C_DOUBLE_COMPLEX), TARGET, PRIVATE, ALLOCATABLE, SAVE :: xcFTs64f_(:)
      !> Plan for transforming from time domain to frequency domain.
      TYPE(C_PTR), PRIVATE, SAVE :: forwardPlan_
      !> Plan for transforming from frequency domain to time domain correlograms.
      TYPE(C_PTR), PRIVATE, SAVE :: inversePlan_
      !> Accuracy of the MKL computations.
      INTEGER(KIND=8), PRIVATE :: accuracy_
      !> Flag indicating FFTw has been initialized.
      LOGICAL, PRIVATE, SAVE :: linitFFTw = .FALSE.
      INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: alignment = 64
      INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: sizeof_float = 4
      INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: sizeof_double = 8
      INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: sizeof_float_complex = 8
      INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: sizeof_double_complex = 16
      !----------------------------------------------------------------------------------!
      !                            Public/Private Subroutines                            !
      !----------------------------------------------------------------------------------!
      PUBLIC :: xcloc_fdxc_initialize
      PUBLIC :: xcloc_fdxc_finalize
      PUBLIC :: xcloc_fdxc_setXCTableF
      PUBLIC :: xcloc_fdxc_computeDefaultXCTableF
      PUBLIC :: xcloc_fdxc_setSignal32fF
      PUBLIC :: xcloc_fdxc_setSignals32fF
      PUBLIC :: xcloc_fdxc_computePhaseCorrelograms
      PUBLIC :: xcloc_fdxc_computeCrossCorrelograms

      PRIVATE :: xcloc_fdxc_computeFDCorrelations
      PRIVATE :: xcloc_fdxc_initializeFFTW
      PRIVATE :: xcloc_fdxc_forwardTransform
      PRIVATE :: xcloc_fdxc_inverseTransform
      PRIVATE :: padLength
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
!>    @param[in] prec      Controls the precision of the module.  
!>    @param[in] accuracy  Controls the accuracy of the vector calculations in MKL.
!>
!>    @param[out] ierr     0 indicates success.
!>
      SUBROUTINE xcloc_fdxc_initialize(npts, nsignals, nptsPad,       &
                                       verbose, prec, accuracy, ierr) &
      BIND(C, NAME='xcloc_fdxc_initialize')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: accuracy, npts, nptsPad, nsignals, &
                                           verbose, prec
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
      IF (prec /= XCLOC_SINGLE_PRECISION .AND. prec /= XCLOC_DOUBLE_PRECISION) THEN
         WRITE(*,908) prec
         ierr = 1
      ENDIF 
      CALL xcloc_fdxc_setAccuracy(accuracy, ierr)
      ! Set the input variables
      npts_       = npts
      nptsPad_    = nptsPad
      nsignals_   = nsignals
      verbose_    = verbose
      nptsInXCs_  = 2*nptsPad_ - 1   ! Length of the cross-correlations
      nptsInFTs_  = nptsInXCs_/2 + 1 ! Number of points in the fourier transforms
      dataOffset_ = padLength(alignment, sizeof_float,         nptsInXCs_)
      ftOffset_   = padLength(alignment, sizeof_float_complex, nptsInFTs_)
      ! Format statements
  905 FORMAT("xcloc_fdxc_initialize: npts=", I8, "must be positive")
  906 FORMAT("xcloc_fdxc_initialize: nsignals=", I8, "must be at least 2")
  907 FORMAT("xcloc_fdxc_initialize: nptsPad=", I8, "must be greater than nts=", I8)
  908 FORMAT("xcloc_fdxc_initialize: prec=", I4, "must be single=",I2, "or double=", I2)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the accuracy of some of the underlying MKL vectorized math functions. 
!>    @param[in] accuracy   XCLOC_HIGH_ACCURACY will use full precision.
!>    @param[in] accuracy   XCLOC_LOW_ACCURACY will discard the last 2 bits.
!>    @param[in] accuracy   XCLOC_EP_ACCURACY will use about 1/2 precision.
!>    @param[out] ierr      0 indicates success.
      SUBROUTINE xcloc_fdxc_setAccuracy(accuracy, ierr) &
      BIND(C, NAME='xcloc_fdxc_setAccuracy')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'mkl_vml.f90'
      INTEGER(C_INT), VALUE, INTENT(IN) :: accuracy
      INTEGER(KIND=4), INTENT(OUT) :: ierr
      IF (accuracy == XCLOC_HIGH_ACCURACY) THEN
         accuracy_ = vmlsetmode(IOR(VML_HA, VML_ERRMODE_ERRNO)) !, VML_ERRMODE_STDERR))
      ELSEIF (accuracy == XCLOC_LOW_ACCURACY) THEN
         accuracy_ = vmlsetmode(IOR(VML_LA, VML_ERRMODE_ERRNO)) !VML_ERRMODE_STDERR))
      ELSEIF (accuracy == XCLOC_EP_ACCURACY) THEN
         accuracy_ = vmlsetmode(IOR(VML_EP, VML_ERRMODE_ERRNO)) !VML_ERRMODE_STDERR))
      ELSE
         WRITE(*,905) accuracy
         ierr = 1
         RETURN
      ENDIF
      ierr = vmlgeterrstatus() 
      IF (ierr /= VML_STATUS_OK) THEN
         WRITE(*,900)
      ELSE
         ierr = 0 
      ENDIF
      accuracy_ = vmlgetmode() ! Get the accuracy however MKL defines it
  900 FORMAT('xcloc_fdxc_setAccuracy: Error setting accuracy')
  905 FORMAT('xcloc_fdxc_setAccuracy: Invalid accuracy', I4)
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
      INCLUDE 'fftw/fftw3.f03'
      IF (ALLOCATED(xcPairs_))         DEALLOCATE(xcPairs_)
      IF (ALLOCATED(inputSignals32f_)) DEALLOCATE(inputSignals32f_)
      IF (ALLOCATED(inputSignals64f_)) DEALLOCATE(inputSignals64f_)
      IF (ALLOCATED(xcs32f_))          DEALLOCATE(xcs32f_)
      IF (ALLOCATED(xcs64f_))          DEALLOCATE(xcs64f_)
      IF (ALLOCATED(inputFTs32f_))     DEALLOCATE(inputFTs32f_)
      IF (ALLOCATED(inputFTs64f_))     DEALLOCATE(inputFTs64f_)
      IF (ALLOCATED(xcFTs32f_))        DEALLOCATE(xcFTs32f_)
      IF (ALLOCATED(xcFTs64f_))        DEALLOCATE(xcFTs64f_) 
      IF (lhaveTable_) THEN
         CALL FFTWF_DESTROY_PLAN(forwardPlan_)
         CALL FFTWF_DESTROY_PLAN(inversePlan_)
      ENDIF
      npts_ = 0
      nptsPad_ = 0
      nsignals_ = 0
      verbose_ = 0
      nptsInXCs_ = 0
      nptsInFTs_ = 0
      dataOffset_ = 0
      ftOffset_ = 0
      linitFFTw = .FALSE.
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
      nxcs_ = nxcs
      IF (ALLOCATED(xcPairs_)) DEALLOCATE(xcPairs_)
      ALLOCATE(xcPairs_(2*nxcs_))
      xcPairs_(1:2*nxcs_) = xcPairs(1:2*nxcs)
      lhaveTable_ = .TRUE.
      ! Allocate rest of space
      CALL xcloc_fdxc_initializeFFTW(ierr)
      IF (ierr /= 0) WRITE(*,903)
      ! Format statements
  900 FORMAT('xcloc_fdxc_setXCTableF: Error nxcs must be positive', I5)
  901 FORMAT('xcloc_fdxc_setXCTableF: minval(xcPairs)', I4, 'must be positive')
  902 FORMAT('xcloc_fdxc_setXCTableF: minval(xcPairs)', I4, 'cannot exceed', I4)
  903 FORMAT('xcloc_fdxc_setXCTableF: Error initializing FFTs')
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
               indx = nsignals_*(i - 1) - (i*(i+1))/2 + j
               xcPairs(2*(indx-1)+1) = i
               xcPairs(2*(indx-1)+2) = j
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
            DO j=i,nsignals_
               indx = nsignals_*(i - 1) - ((i-1)*i)/2 + j
               xcPairs(2*(indx-1)+1) = i
               xcPairs(2*(indx-1)+2) = j
            ENDDO
         ENDDO
      ENDIF
 905  FORMAT('xcloc_fdxc_computeDefaultXCTableF: Error - nwork must be >=', I5)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Convenience routine to set all the signals on the module.
!>    @param[in] ldx       Leading dimension of x.  This cannot be less than npts_.
!>    @param[in] npts      Number of points in each signal.
!>    @param[in] nsignals  Number of signals.
!>    @param[in] x         Signals to set.  This is an [ldx x nsignal] matrix.
!>                         with leading dimension ldx. 
!>    @param[out] ierr     0 indicates success.
      SUBROUTINE xcloc_fdxc_setSignals32fF(ldx, npts, nsignals, x, ierr)  &
      BIND(C, NAME='xcloc_fdxc_setSignals32fF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldx, npts, nsignals
      REAL(C_FLOAT), INTENT(IN) :: x(ldx*nsignals)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i, ix
      ierr = 0
      IF (ldx < npts .OR. npts /= npts_ .OR. nsignals /= nsignals_) THEN
         IF (ldx < npts) WRITE(*,900) ldx, npts
         IF (npts /= npts_) WRITE(*,901) npts_
         IF (nsignals /= nsignals_) WRITE(*,902) nsignals_
         ierr = 1
         RETURN
      ENDIF
      DO i=1,nsignals_
         ix = (i - 1)*ldx + 1
         CALL xcloc_fdxc_setSignal32fF(i, npts, x(ix), ierr)
         IF (ierr /= 0) THEN
            WRITE(*,910) i
            RETURN
         ENDIF
      ENDDO
  900 FORMAT('xcloc_fdxc_setSignals32fF: Error ldx=', I6, '<', 'npts=', I6)
  901 FORMAT('xcloc_fdxc_setSignals32fF: Error expecting npts=', I6)
  902 FORMAT('xcloc_fdxc_setSignals32fF: Error expecting nsignals=', I6)
  910 FORMAT('xcloc_fdxc_setSignals32fF: Error setting signal index', I4)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE xcloc_fdxc_getCorrelograms32f(ldxc, nxcs, xcs, ierr)  &
      BIND(C, NAME='xcloc_fdxc_getCorrelograms32f')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldxc, nxcs
      REAL(C_FLOAT), INTENT(OUT) :: xcs(ldxc*nxcs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i1, ixc
      ierr = 0
      IF (ldxc < nptsInXCs_) THEN
         ierr = 1
         RETURN
      ENDIF
      IF (nxcs < nxcs_) THEN
         ierr = 1
         RETURN
      ENDIF
      DO ixc=1,nxcs_
         i1 = (ixc - 1)*ldxc + 1
         CALL xcloc_fdxc_getCorrelogram32fF(ixc, ldxc, xcs(i1), ierr)
      ENDDO
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the corrNumber'th cross-correlation on the module.
!>    @param[in] corrNumber  Correlation number to get.  This must be in the rnage
!>                           [1, nxcs_].
!>    @param[in] lwork       Space allocated to xc.
!>    @param[out] xc         The corrNumber'th cross-correlation.  This has dimension
!>                           [lwork] however only the first nptsInXCs_ points are
!>                           accessed.
!>    @param[out] ierr       0 indicates success.
      SUBROUTINE xcloc_fdxc_getCorrelogram32fF(corrNumber, lwork, xc, ierr) &
      BIND(C, NAME='xcloc_fdxc_getCorrelogram32fF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: corrNumber, lwork
      REAL(C_FLOAT), INTENT(OUT) :: xc(*)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i1
      ierr = 0
      IF (lwork < nptsInXCs_) THEN
         WRITE(*,900) lwork, nptsInXCs_
         ierr = 1
         RETURN
      ENDIF
      IF (corrNumber < 1 .OR. corrNumber > nxcs_) THEN
         WRITE(*,905) corrNumber, nxcs_,  nxcs_
         ierr = 1
         RETURN
      ENDIF
      i1 = (corrNumber - 1)*dataOffset_ + 1
      IF (precision_ == XCLOC_SINGLE_PRECISION) THEN
         xc(1:nptsInXCs_) = xcs32f_(i1:i1+nptsInXCs_-1)
      ELSE
         ierr = ippsConvert_64f32f(xcs64f_(i1), xc, nptsInXCs_)
      ENDIF
  900 FORMAT('xcloc_fdxc_getCorrelogram32fF: lworC=', I6, 'must be at least', I6)
  905 FORMAT('xcloc_fdxc_getCorrelogram32fF: corrNumber',I6, 'must be in range [1',I6,']')
      RETURN 
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the signalNumber'th input signal on the module.
!>    @param[in] signalNumber  Number of signal to set.  This must be in the range
!>                             [1, nsignals_].
!>    @param[in] npts          Number of points in the signal.  This cannot exceed npts_.
!>    @param[in] x             Signal to set.  This is an array of dimension [npts].
!>    @param[out] ierr         0 indicates success.
      SUBROUTINE xcloc_fdxc_setSignal32fF(signalNumber, npts, x, ierr) &
      BIND(C, NAME='xcloc_fdxc_setSignal32fF')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: signalNumber, npts
      REAL(C_FLOAT), INTENT(IN) :: x(npts)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i1, i2
      ierr = 1
      IF (npts < 1) RETURN ! nothing to do
      IF (npts > npts_) THEN
         WRITE(*,900) npts_
         RETURN
      ENDIF
      IF (signalNumber < 1 .OR. signalNumber > nsignals_) THEN
         WRITE(*,905) signalNumber
         RETURN
      ENDIF
      ierr = 0 
      i1 = (signalNumber - 1)*dataOffset_ + 1 
      i2 = signalNumber*dataOffset_
      IF (precision_ == XCLOC_SINGLE_PRECISION) THEN
         inputSignals32f_(i1:i1+npts-1) = x(1:npts)
      ELSE
         ierr = ippsConvert_32f64f(x, inputSignals64f_(i1), npts)
      ENDIF
  900 FORMAT('xcloc_fdxc_setSignal32fF: Error expecting npts=', I5)
  905 FORMAT('xcloc_fdxc_setSignal32fF: Error signalNumber must be in range [1,',I4,']')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the phase correlograms.
!>    @param[out] ierr  0 indicates success.
      SUBROUTINE xcloc_fdxc_computePhaseCorrelograms(ierr) &
      BIND(C, NAME='xcloc_fdxc_computePhaseCorrelograms')
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(OUT) :: ierr
      LOGICAL(C_BOOL), PARAMETER :: lphaseCorr = .TRUE.
      ! Convert input signals to time domain
      CALL xcloc_fdxc_forwardTransform(ierr)
      IF (ierr /= 0) THEN
         WRITE(*,900)  
         RETURN
      ENDIF
      ! Compute the phase correlograms 
      CALL xcloc_fdxc_computeFDCorrelations(lphasecorr, ierr)
      ! Inverse transform back to the time domain  
      CALL xcloc_fdxc_inverseTransform(ierr)
      IF (ierr /= 0) THEN
         WRITE(*,910)
         RETURN
      ENDIF
  900 FORMAT('xcloc_fdxc_computePhaseCorrelograms: Error computing forward transforms') 
  910 FORMAT('xcloc_fdxc_computePhaseCorrelograms: Error computing inverse transforms')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the cross-correlograms.
!>    @param[out] ierr  0 indicates success.
      SUBROUTINE xcloc_fdxc_computeCrossCorrelograms(ierr) &
      BIND(C, NAME='xcloc_fdxc_computeCrossCorrelograms')
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(OUT) :: ierr
      LOGICAL(C_BOOL), PARAMETER :: lphaseCorr = .FALSE.
      ! Convert input signals to time domain
      CALL xcloc_fdxc_forwardTransform(ierr)
      IF (ierr /= 0) THEN
         WRITE(*,900)  
         RETURN
      ENDIF
      ! Compute the phase correlograms 
      CALL xcloc_fdxc_computeFDCorrelations(lphasecorr, ierr)
      ! Inverse transform back to the time domain  
      CALL xcloc_fdxc_inverseTransform(ierr)
      IF (ierr /= 0) THEN
         WRITE(*,910)
         RETURN
      ENDIF
  900 FORMAT('xcloc_fdxc_computeCrossCorrelograms: Error computing forward transforms')
  910 FORMAT('xcloc_fdxc_computeCrossCorrelograms: Error computing inverse transforms')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the (phase) correlations in the frequency domain.
!>    @param[in] lphaseCorr  If true then compute the phase correlations.
      SUBROUTINE xcloc_fdxc_computeFDCorrelations(lphaseCorr, ierr)
      IMPLICIT NONE
      INCLUDE 'mkl_vml.f90'
      LOGICAL(C_BOOL), INTENT(IN) :: lphaseCorr
      INTEGER, INTENT(OUT) :: ierr
      INTEGER i, indx, ixc, iw, j, jndx, kndx, n
      COMPLEX(C_DOUBLE_COMPLEX), CONTIGUOUS, POINTER :: xcPtr64z(:)
      COMPLEX(C_FLOAT_COMPLEX), CONTIGUOUS, POINTER :: xcPtr32c(:)
      REAL(C_DOUBLE), ALLOCATABLE :: mag64f_(:)
      REAL(C_FLOAT), ALLOCATABLE :: mag32f_(:)
      REAL(C_FLOAT), PARAMETER :: tol32 = EPSILON(1.0)*10.0
      REAL(C_DOUBLE), PARAMETER :: tol64 = EPSILON(1.d0)*10.d0 
      ierr = 0
      IF (precision_ == XCLOC_SINGLE_PRECISION) THEN
!        !$OMP PARALLEL DEFAULT(NONE) &
!        !$OMP SHARED(accuracy_, ftOffset_, inputFTs32f_, lphaseCorr, nxcs_) &
!        !$OMP SHARED(nptsInFTs_, xcFTs32f_, xcPairs_) &
!        !$OMP PRIVATE(i, indx, iw, ixc, j, jndx, kndx, mag32f_, n) &
!        !$OMP FIRSTPRIVATE(xcPtr32c)
         n = nptsInFTs_
         IF (lphaseCorr) THEN
            ALLOCATE(mag32f_(n))
            NULLIFY(xcPtr32c)
         ENDIF
         ! Loop on the number of cross-correlations
!        !$OMP DO
         DO ixc=1,nxcs_
            i = xcPairs_(2*(ixc-1)+1)
            j = xcPairs_(2*(ixc-1)+2)
            indx = (i - 1)*ftOffset_ + 1
            jndx = (j - 1)*ftOffset_ + 1
            kndx = (ixc - 1)*ftOffset_ + 1
            ! Compute u1*conj(u2)
            CALL vmcMulByConj(n, inputFTs32f_(indx), inputFTs32f_(jndx), &
                              xcFTs32f_(kndx), accuracy_)
            ! Normalize by the magnitude at each frequency
            IF (lphaseCorr) THEN
               CALL vmcAbs(n, xcFTs32f_(kndx), mag32f_, accuracy_)
               ! Avoid division by zero
               mag32f_(1:n) = MAX(mag32f_(1:n), tol32)
               ! Safely normalize
               xcPtr32c => xcFTs32f_(kndx:kndx+n-1)
               !$OMP DO SIMD ALIGNED(xcPtr32c, mag32f_: 16)
               DO iw=1,n
                  xcPtr32c(iw) = xcPtr32c(iw)/mag32f_(iw)
                  !xcFTs32f_(kndx-1+iw) = xcFTs32f_(kndx-1+iw)/mag32f_(iw)
               ENDDO
               NULLIFY(xcPtr32c)
            ENDIF
         ENDDO ! Loop on number of cross-correlations
         IF (lphaseCorr) DEALLOCATE(mag32f_)
!        !$OMP END PARALLEL
      ELSE
         !$OMP PARALLEL DEFAULT(NONE) &
         !$OMP SHARED(accuracy_, ftOffset_, inputFTs64f_, lphaseCorr, nxcs_) &
         !$OMP SHARED(nptsInFTs_, xcFTs64f_, xcPairs_) &
         !$OMP PRIVATE(i, indx, iw, ixc, j, jndx, kndx, mag64f_, n, xcPtr64z)
         n = nptsInFTs_
         IF (lphaseCorr) THEN
            NULLIFY(xcPtr64z)
            ALLOCATE(mag64f_(n))
         ENDIF
         ! Loop on the number of cross-correlations
         !$OMP DO
         DO ixc=1,nxcs_
            i = xcPairs_(2*(ixc-1)+1)
            j = xcPairs_(2*(ixc-1)+2)
            indx = (i - 1)*ftOffset_ + 1 
            jndx = (j - 1)*ftOffset_ + 1 
            kndx = (ixc - 1)*ftOffset_ + 1 
            ! Compute u1*conj(u2)
            CALL vmzMulByConj(n, inputFTs64f_(indx), inputFTs64f_(jndx), &
                              xcFTs64f_(kndx), accuracy_)
            ! Normalize by the magnitude at each frequency
            IF (lphaseCorr) THEN
               CALL vmzAbs(n, xcFTs64f_(kndx), mag64f_, accuracy_)
               ! Avoid division by zero
               mag64f_(1:n) = MAX(mag64f_(1:n), tol64)
               ! Safely normalize
               xcPtr64z => xcFTs64f_(kndx:kndx+n-1) 
               !$OMP DO SIMD ALIGNED(xcPtr64z, mag64f_: 16)
               DO iw=1,n
                  xcPtr64z(iw) = xcPtr64z(iw)/mag64f_(iw)
               ENDDO
               NULLIFY(xcPtr64z)
            ENDIF
         ENDDO ! Loop on number of cross-correlations
         IF (lphaseCorr) DEALLOCATE(mag64f_)
         !$OMP END PARALLEL
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Forward transforms the input signals.
!>    @param[out] ierr  0 indicates success.
      SUBROUTINE xcloc_fdxc_forwardTransform(ierr) 
      IMPLICIT NONE
      INCLUDE 'fftw/fftw3.f03'
      INTEGER, INTENT(OUT) :: ierr
      ierr = 0
      IF (nsignals_ == 0) RETURN ! nothing to do
      IF (.NOT.linitFFTw) THEN
         WRITE(*,900) 
         ierr = 1
         RETURN
      ENDIF
      IF (precision_ == XCLOC_SINGLE_PRECISION) THEN
         CALL FFTWF_EXECUTE_DFT_R2C(forwardPlan_, inputSignals32f_, inputFTs32f_)
      ELSE
         CALL FFTW_EXECUTE_DFT_R2C(forwardPlan_, inputSignals64f_, inputFTs64f_)
      ENDIF
  900 FORMAT('xcloc_fdxc_forwardTransform: Transforms not yet initialized')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Brings the frequency domain correlations back to the time domain.
!>    @param[out] ierr  0 indicates success.
      SUBROUTINE xcloc_fdxc_inverseTransform(ierr)
      IMPLICIT NONE
      INCLUDE 'fftw/fftw3.f03'
      INTEGER, INTENT(OUT) :: ierr
      REAL(C_FLOAT), ALLOCATABLE :: work32(:)
      REAL(C_FLOAT) scal32
      INTEGER indx, jndx, ixc, ncopy1, ncopy2
      ierr = 0 
      IF (nsignals_ == 0) RETURN ! nothing to do
      IF (.NOT.linitFFTw) THEN
         WRITE(*,900)
         ierr = 1
         RETURN
      ENDIF
      IF (precision_ == XCLOC_SINGLE_PRECISION) THEN
         scal32 = 1.0/FLOAT(nptsInXCs_)
         CALL FFTWF_EXECUTE_DFT_C2R(inversePlan_, xcFTs32f_, xcs32f_)
         ALLOCATE(work32(ftOffset_)); work32(:) = 0.0
         ncopy1 = nptsInXCs_/2
         ncopy2 = ncopy1 + 1
         DO ixc=1,nxcs_
            indx = (ixc - 1)*dataOffset_ + 1
            jndx = (ixc - 1)*dataOffset_ + ncopy2
            work32(ncopy1+1:nptsInXCs_)     = xcs32f_(indx:indx+ncopy1)
            work32(1:ncopy1) = xcs32f_(jndx+1:jndx+ncopy1)
            xcs32f_(indx:indx+nptsInXCs_-1) = work32(1:nptsInXCs_)*scal32

            !indx = (ixc - 1)*dataOffset_ + ncopy1 + 1
            !jndx = (ixc - 1)*dataOffset_ + 1
            !work32(1:ncopy2) = xcs32f_(jndx:jndx+ncopy2-1)/scal32
            !xcs32f_(jndx:jndx+ncopy1-1) = xcs32f_(indx+1:indx+ncopy1)/scal32
            !xcs32f_(indx:indx+ncopy2-1) = work32(1:ncopy2) 
         ENDDO
         DEALLOCATE(work32)
      ELSE
         CALL FFTW_EXECUTE_DFT_C2R(inversePlan_, xcFTs64f_, xcs64f_)
      ENDIF
  900 FORMAT('xcloc_fdxc_inverseTransform: Transforms not yet initialized')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Initializes the FFTw descriptors. 
!>    @result 0 indicates success.
      SUBROUTINE xcloc_fdxc_initializeFFTW(ierr)
      IMPLICIT NONE
      INCLUDE 'fftw/fftw3.f03'
      INTEGER, INTENT(OUT) :: ierr
      INTEGER(C_INT) nf(1), ni(1), inembed(1), onembed(1)
      INTEGER(C_INT), PARAMETER :: rank = 1    ! Computing multiple 1D transforms
      INTEGER(C_INT), PARAMETER :: istride = 1 ! Distance elements in input column
      INTEGER(C_INT), PARAMETER :: ostride = 1 ! Distance elements in output column
      ierr = 0
      linitFFTw = .FALSE.
      IF (nsignals_ < 1 .OR. nptsInXCs_ < 1) RETURN ! Not enough info to initialize
      IF (verbose_ > 0) WRITE(*,*) 'xcloc_fdxc_initializeFFTW: Initializing FFTs...'
      IF (ALLOCATED(inputSignals32f_)) DEALLOCATE(inputSignals32f_)
      IF (ALLOCATED(inputFTs32f_))     DEALLOCATE(inputFTs32f_)
      IF (ALLOCATED(xcFTs32f_))        DEALLOCATE(xcFTs32f_)
      IF (ALLOCATED(xcs32f_))          DEALLOCATE(xcs32f_)
      IF (ALLOCATED(inputSignals64f_)) DEALLOCATE(inputSignals64f_)
      IF (ALLOCATED(inputFTs64f_))     DEALLOCATE(inputFTs64f_)
      IF (ALLOCATED(xcFTs64f_))        DEALLOCATE(xcFTs64f_)
      IF (ALLOCATED(xcs64f_))          DEALLOCATE(xcs64f_)
      inembed(1) = 0
      onembed(1) = 0
      nf(1) = nptsInXCs_ ! Each signal has length nptsInXCs_
      ni(1) = nptsInXCs_ ! Each signal has length nptsInXCs_
      IF (precision_ == XCLOC_SINGLE_PRECISION) THEN
         ALLOCATE(inputSignals32f_(dataOffset_*nsignals_)); inputSignals32f_(:) = 0.0
         ALLOCATE(inputFTs32f_(ftOffset_*nsignals_)); inputFTs32f_(:) = czero
         ALLOCATE(xcFTs32f_(ftOffset_*nxcs_)); xcFTs32f_(:) = czero
         ALLOCATE(xcs32f_(dataOffset_*nxcs_)); xcs32f_(:) = 0.0
         forwardPlan_ = FFTWF_PLAN_MANY_DFT_R2C(rank, nf, nsignals_,       &
                                                inputSignals32f_, inembed, &
                                                istride, dataOffset_,      &
                                                inputFTs32f_, onembed,     &
                                                ostride, ftOffset_,        &
                                                FFTW_PATIENT)
!        inversePlan_ = FFTWF_PLAN_MANY_DFT_R2C(rank, ni, nxcs_,      &
!                                               xcs32f_, inembed,     &
!                                               istride, dataOffset_, &
!                                               xcFTs32f_, onembed,   &
!                                               ostride, ftOffset_,   & 
!                                               FFTW_PATIENT)
         inversePlan_ = FFTWF_PLAN_MANY_DFT_C2R(rank, ni, nxcs_,        &
                                                xcFTs32f_, inembed,     &
                                                istride, ftOffset_,     &
                                                xcs32f_, onembed,       &
                                                ostride, dataOffset_,   & 
                                                FFTW_PATIENT)
      ELSE
         ALLOCATE(inputSignals64f_(dataOffset_*nsignals_)); inputSignals64f_(:) = 0.d0
         ALLOCATE(inputFTs64f_(ftOffset_*nsignals_)); inputFTs64f_(:) = zzero
         ALLOCATE(xcFTs64f_(ftOffset_*nxcs_)); xcFTs64f_(:) = zzero
         ALLOCATE(xcs64f_(dataOffset_*nxcs_)); xcs64f_(:) = 0.d0
         forwardPlan_ = FFTW_PLAN_MANY_DFT_R2C(rank, nf, nsignals_,       &
                                               inputSignals64f_, inembed, &
                                               istride, dataOffset_,      &
                                               inputFTs64f_, onembed,     &
                                               ostride, ftOffset_,        &
                                               FFTW_PATIENT)
         inversePlan_ = FFTW_PLAN_MANY_DFT_C2R(rank, ni, nxcs_,        &
                                               xcFTs64f_, inembed,     &
                                               istride, ftOffset_,     &
                                               xcs64f_, onembed,       &
                                               ostride, dataOffset_,   &
                                               FFTW_PATIENT)
      ENDIF
      linitFFTw = .TRUE.
      RETURN
      END

      INTEGER(C_INT) FUNCTION padLength(alignment, sizeof_dataType, n) &
      BIND(C, NAME='padLength')
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
