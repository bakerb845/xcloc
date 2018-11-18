!> @defgroup spxc Cross Correlogram Signals Processing
!> @ingroup xcloc
!> @ingroup xcloc_mpi
!> @brief Performs some signals processing on the cross-correlograms.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_SPXC
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
      USE XCLOC_UTILS
      USE XCLOC_IPPS
      IMPLICIT NONE
      !> @ingroup spxc
      !> Sparse real Hilbert transform coefficients.
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: sparseHfiltR64f_
      !> @ingroup spxc
      !> Sparse imaginary Hilbert transform coefficients.
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: sparseHfiltI64f_
      !> @ingroup spxc
      !> Sparse real Hilbert transform coefficients for single precision.
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: sparseHfiltR32f_
      !> @ingroup spxc
      !> Sparse imaginary Hilbert transform coefficients for single precision.
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: sparseHfiltI32f_
      !> @ingroup spxc
      !> Real Hilbert transformer coefficients.
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: hfiltR64f_
      !> @ingroup spxc
      !> Imaginary Hilbert transformer coefficients.
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: hfiltI64f_
      !> @ingroup spxc
      !> Real Hilbert transformer coefficients for single precision.
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: hfiltR32f_
      !> @ingroup spxc
      !> Imaginary Hilbert transformer coefficients for single precision.
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: hfiltI32f_
      !> @ingroup spxc
      !> FIR averaging filter coefficients.
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: rms64f_
      !> @ingroup spxc
      !> FIR averaging filter for single precision. 
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: rms32f_
      !> @ingroup spxc
      !> Non-zero indices (C indexed) of real FIR coefficients.
      INTEGER(C_INT), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: nzReIndices_
      !> @ingroup spxc
      !> Non-zero indices (C indexed) of imaginary FIR coefficients.
      INTEGER(C_INT), PRIVATE, ALLOCATABLE, DIMENSION(:), SAVE :: nzImIndices_
      !!> Accuracy.
      !INTEGER(C_INT), PRIVATE, SAVE :: accuracy_ = XCLOC_HIGH_ACCURACY
      !> @ingroup spxc
      !> Filtering type.
      INTEGER(C_INT), PRIVATE, SAVE :: ftype_ = XCLOC_SPXC_DONOT_FILTER
      !> @ingroup spxc
      !> Number of filter coefficients.
      INTEGER(C_INT), PRIVATE, SAVE :: nCoeffs_ = 0
      !> @ingroup spxc
      !> Number of real FIR coefficients (including zeros).
      INTEGER(C_INT), PRIVATE, SAVE :: nReCoeffs_ = 0
      !> @ingroup spxc
      !> Number of imaginary FIR coefficients (including zeros).
      INTEGER(C_INT), PRIVATE, SAVE :: nImCoeffs_ = 0
      !> @ingroup spxc
      !> Filter order.
      INTEGER(C_INT), PRIVATE, SAVE :: nOrder_ = 0
      !> @ingroup spxc
      !> Number of non-zero real coefficients.
      INTEGER(C_INT), PRIVATE, SAVE :: nnzReCoeffs_ = 0
      !> @ingroup spxc
      !> Number of non-zero imaginary coefficients.
      INTEGER(C_INT), PRIVATE, SAVE :: nnzImCoeffs_ = 0
      !> @ingroup spxc
      !> If true then this is a Type III FIR filter and is sparse.
      LOGICAL, PRIVATE, SAVE :: lisTypeIII_ = .TRUE.
      !> @ingroup spxc
      !> \f$ \pi \f$
      DOUBLE PRECISION, PARAMETER, PRIVATE :: pi = 3.14159265358979323846d0

      !----------------------------------------------------------------------------------!
      !                             Public/Private Subroutines                           !
      !----------------------------------------------------------------------------------!
      PUBLIC :: xcloc_spxc_initialize
      PUBLIC :: xcloc_spxc_finalize
      PUBLIC :: xcloc_spxc_filterXCsInPlace64f
      PUBLIC :: xcloc_spxc_filterXCsInPlace32f
      PUBLIC :: xcloc_spxc_filterXCsOutOfPlace64f
      PUBLIC :: xcloc_spxc_filterXCsOutOfPlace32f

      PRIVATE :: xcloc_spxc_rmsFilter_design
      PRIVATE :: xcloc_spxc_envelope_design
      PRIVATE :: xcloc_spxc_envelope_apply64f
      PRIVATE :: xcloc_spxc_envelope_apply32f
      PRIVATE :: xcloc_spxc_rmsFilter_apply64f
      PRIVATE :: xcloc_spxc_rmsFilter_apply32f
      PRIVATE :: sinc
      CONTAINS
!========================================================================================!
!                                       Begin the Code                                   !
!========================================================================================!
!>    Initializes the processing module.
!>    @param[in] n         The number of filter coefficients.  If the processing is
!>                         to be done then this must be odd.
!>    @param[in] ftype     XCLOC_SPXC_DONOT_FILTER will not filter correlograms.
!>    @param[in] ftype     XCLOC_SPXC_ENVELOPE_FILTER will apply envelope to correlograms.
!>    @param[in] ftype     XCLOC_SPXC_RMS_FILTER will comptue RMS of correlograms.
!>    @param[out] ierr     0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_initialize(n, ftype, ierr) &
      BIND(C, NAME='xcloc_spxc_initialize')
      INTEGER(C_INT), VALUE, INTENT(IN) :: n, ftype
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      CALL xcloc_spxc_finalize()
      IF (.NOT. xcloc_constants_isValidFilteringType(ftype)) THEN
         WRITE(ERROR_UNIT,905) ftype
         ierr = 1
         RETURN
      ENDIF
      IF (ftype == XCLOC_SPXC_DONOT_FILTER) THEN
         ftype_ = ftype 
      ELSE
         ! Check inputs
         IF (n < 1) THEN
            WRITE(ERROR_UNIT,900)
            ierr = 1
            RETURN
         ENDIF
         ! Nothing to do
         IF (n == 1) THEN
            ftype_ = XCLOC_SPXC_DONOT_FILTER
            WRITE(ERROR_UNIT,901) 
            RETURN
         ENDIF
         ! Require FIR filters have an odd number of coefficients for easy group delay
         nCoeffs_ = n
         IF (MOD(n, 2) == 0) THEN
            WRITE(ERROR_UNIT,904)
            nCoeffs_ = n + 1
         ENDIF
         nOrder_ = nCoeffs_ - 1 
         ! Design the filters
         IF (ftype == XCLOC_SPXC_ENVELOPE_FILTER) THEN
            CALL xcloc_spxc_envelope_design(ierr)
            IF (ierr /= 0) THEN
               WRITE(ERROR_UNIT,906)
               RETURN
            ENDIF
            ftype_ = XCLOC_SPXC_ENVELOPE_FILTER
         ELSEIF (ftype == XCLOC_SPXC_RMS_FILTER) THEN
            CALL xcloc_spxc_rmsFilter_design(ierr)
            IF (ierr /= 0) THEN
               WRITE(ERROR_UNIT,907)
               RETURN
            ENDIF
            ftype_ = XCLOC_SPXC_RMS_FILTER
         ELSE
            WRITE(ERROR_UNIT,905) ftype
            ierr = 1
            RETURN
         ENDIF
      ENDIF
  900 FORMAT('xcloc_spxc_initialize: Number of filter coefficients must be posisitve')
  901 FORMAT('xcloc_spxc_initialize: nTaps=1; skipping filtering')
  904 FORMAT('xcloc_spxc_initialize: Adding one filter coefficient')
  905 FORMAT('xcloc_spxc_initialize: Invalid filter type ', I0)
  906 FORMAT('xcloc_spxc_initialize: Error designing Hilbert transformer')
  907 FORMAT('xcloc_spxc_initialize: Error designing RMS filter')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Filters the correlograms out of place.
!>    @param[in] ldxc       Leading dimension of correlograms.  This must be at least
!>                          nptsInXCs.
!>    @param[in] nptsInXCs  Number of points in correlograms. 
!>    @param[in] nxcs       Number of correlograms.
!>    @param[in] xcs        Correlograms to filter.  This is an [ldxc x nxcs] matrix in
!>                          column major format.
!>    @param[out] xcsFilt   Filtered correlograms.  This is an [ldxc x nxcs] matrix in
!>                          column major format.
!>    @param[out] ierr      0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_filterXCsOutOfPlace64f(ldxc, nptsInXCs, nxcs, &
                                                   xcs, xcsFilt, ierr)    &
      BIND(C, NAME='xcloc_spxc_filterXCsOutOfPlace64f')
      INTEGER(C_INT), VALUE, INTENT(IN) ::  ldxc, nxcs, nptsInXCs
      REAL(C_DOUBLE), INTENT(IN) :: xcs(1:ldxc*nxcs)
      REAL(C_DOUBLE), INTENT(OUT) :: xcsFilt(1:ldxc*nxcs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      xcsFilt(1:ldxc*nxcs) = xcs(1:ldxc*nxcs)
      CALL xcloc_spxc_filterXCsInPlace64f(ldxc, nptsInXCs, nxcs, xcsFilt, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
      ENDIF
  900 FORMAT('xcloc_spxc_filterXCsOutOfPlace64f: In place filtering failed')
      RETURN
      END
!>    @brief Filters the correlograms out of place.
!>    @param[in] ldxc       Leading dimension of correlograms.  This must be at least
!>                          nptsInXCs.
!>    @param[in] nptsInXCs  Number of points in correlograms. 
!>    @param[in] nxcs       Number of correlograms.
!>    @param[in] xcs        Correlograms to filter.  This is an [ldxc x nxcs] matrix in
!>                          column major format.
!>    @param[out] xcsFilt   Filtered correlograms.  This is an [ldxc x nxcs] matrix in
!>                          column major format.
!>    @param[out] ierr      0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_filterXCsOutOfPlace32f(ldxc, nptsInXCs, nxcs, &
                                                   xcs, xcsFilt, ierr)    &
      BIND(C, NAME='xcloc_spxc_filterXCsOutOfPlace32f')
      INTEGER(C_INT), VALUE, INTENT(IN) ::  ldxc, nxcs, nptsInXCs
      REAL(C_FLOAT), INTENT(IN) :: xcs(1:ldxc*nxcs)
      REAL(C_FLOAT), INTENT(OUT) :: xcsFilt(1:ldxc*nxcs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      xcsFilt(1:ldxc*nxcs) = xcs(1:ldxc*nxcs)
      CALL xcloc_spxc_filterXCsInPlace32f(ldxc, nptsInXCs, nxcs, xcsFilt, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1 
      ENDIF
  900 FORMAT('xcloc_spxc_filterXCsOutOfPlace32f: In place filtering failed')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Filters the correlograms in place.
!>    @param[in] ldxc       Leading dimension of correlograms.  This must be at least
!>                          nptsInXCs.
!>    @param[in] nptsInXCs  Number of points in correlograms. 
!>    @param[in] nxcs       Number of correlograms.
!>    @param[in,out] xcs    Correlograms to filter.  This is an [ldxc x nxcs] matrix in
!>                          column major format.
!>    @param[in,out] xcs    On exit these are the filtered correlograms.
!>    @param[out] ierr      0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_filterXCsInPlace64f(ldxc, nptsInXCs, nxcs, &
                                                xcs, ierr)             &
      BIND(C, NAME='xcloc_spxc_filterXCsInPlace64f')
      INTEGER(C_INT), VALUE, INTENT(IN) ::  ldxc, nxcs, nptsInXCs
      REAL(C_DOUBLE), INTENT(INOUT) :: xcs(1:ldxc*nxcs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER lds, npts, nsignals
      ierr = 0 
      lds = ldxc
      nsignals = nxcs
      npts = nptsInXCs
      IF (ldxc <  nptsInXCs .OR. nxcs < 1 .OR. nptsInXCs < 1) THEN
         IF (ldxc < nptsInXCs) WRITE(ERROR_UNIT,900)
         IF (nxcs < 1) WRITE(ERROR_UNIT,901)
         IF (nptsInXCs < 1) WRITE(ERROR_UNIT,902)
         ierr = 1 
         RETURN
      ENDIF
      IF (ftype_ == XCLOC_SPXC_ENVELOPE_FILTER) THEN
         CALL xcloc_spxc_envelope_apply64f(lds, npts, nsignals, xcs, ierr)
         IF (ierr /= 0) WRITE(ERROR_UNIT,903)
      ELSEIF (ftype_ == XCLOC_SPXC_RMS_FILTER) THEN
         CALL xcloc_spxc_rmsFilter_apply64f(lds, npts, nsignals, xcs, ierr)
         IF (ierr /= 0) WRITE(ERROR_UNIT,904)
      ELSE
         IF (ftype_ /= XCLOC_SPXC_DONOT_FILTER) THEN
            WRITE(ERROR_UNIT,905)
            ierr = 1
         ENDIF
         RETURN
      ENDIF
  900 FORMAT('xcloc_spxc_filterXCsInPlace64f: Error ldxc < nptsInXCs')
  901 FORMAT('xcloc_spxc_filterXCsInPlace64f: No correlograms')
  902 FORMAT('xcloc_spxc_filterXCsInPlace64f: No points in signals')
  903 FORMAT('xcloc_spxc_filterXCsInPlace64f: Envelope filter failed')
  904 FORMAT('xcloc_spxc_filterXCsInPlace64f: RMS filter failed') 
  905 FORMAT('xcloc_spxc_filterXCsInPlace64f: Internal error - unknown filter')
      RETURN
      END
!>    @brief Filters the correlograms in place.
!>    @param[in] ldxc       Leading dimension of correlograms.  This must be at least
!>                          nptsInXCs.
!>    @param[in] nptsInXCs  Number of points in correlograms. 
!>    @param[in] nxcs       Number of correlograms.
!>    @param[in,out] xcs    Correlograms to filter.  This is an [ldxc x nxcs] matrix in
!>                          column major format.
!>    @param[in,out] xcs    On exit these are the filtered correlograms.
!>    @param[out] ierr      0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_filterXCsInPlace32f(ldxc, nptsInXCs, nxcs, &
                                                xcs, ierr)             &
      BIND(C, NAME='xcloc_spxc_filterXCsInPlace32f')
      INTEGER(C_INT), VALUE, INTENT(IN) ::  ldxc, nxcs, nptsInXCs
      REAL(C_FLOAT), INTENT(INOUT) :: xcs(1:ldxc*nxcs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER lds, npts, nsignals
      ierr = 0 
      lds = ldxc
      nsignals = nxcs
      npts = nptsInXCs
      IF (ldxc <  nptsInXCs .OR. nxcs < 1 .OR. nptsInXCs < 1) THEN
         IF (ldxc < nptsInXCs) WRITE(ERROR_UNIT,900)
         IF (nxcs < 1) WRITE(ERROR_UNIT,901)
         IF (nptsInXCs < 1) WRITE(ERROR_UNIT,902)
         ierr = 1 
         RETURN
      ENDIF
      IF (ftype_ == XCLOC_SPXC_ENVELOPE_FILTER) THEN
         CALL xcloc_spxc_envelope_apply32f(lds, npts, nsignals, xcs, ierr)
         IF (ierr /= 0) WRITE(ERROR_UNIT,903)
      ELSEIF (ftype_ == XCLOC_SPXC_RMS_FILTER) THEN
         CALL xcloc_spxc_rmsFilter_apply32f(lds, npts, nsignals, xcs, ierr)
         IF (ierr /= 0) WRITE(ERROR_UNIT,904)
      ELSE
         IF (ftype_ /= XCLOC_SPXC_DONOT_FILTER) THEN
            WRITE(ERROR_UNIT,905)
            ierr = 1
         ENDIF
         RETURN
      ENDIF
  900 FORMAT('xcloc_spxc_filterXCsInPlace32f: Error ldxc < nptsInXCs')
  901 FORMAT('xcloc_spxc_filterXCsInPlace32f: No correlograms')
  902 FORMAT('xcloc_spxc_filterXCsInPlace32f: No points in signals')
  903 FORMAT('xcloc_spxc_filterXCsInPlace32f: Envelope filter failed')
  904 FORMAT('xcloc_spxc_filterXCsInPlace32f: RMS filter failed') 
  905 FORMAT('xcloc_spxc_filterXCsInPlace32f: Internal error - unknown filter')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases memory and resets variables on the XC signals processing module.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_finalize( ) &
      BIND(C, NAME='xcloc_spxc_finalize')
      IF (ALLOCATED(rms64f_)) DEALLOCATE(rms64f_)
      IF (ALLOCATED(rms32f_)) DEALLOCATE(rms32f_)
      IF (ALLOCATED(hfiltR64f_)) DEALLOCATE(hfiltR64f_) 
      IF (ALLOCATED(hfiltR32f_)) DEALLOCATE(hfiltR32f_)
      IF (ALLOCATED(hfiltI64f_)) DEALLOCATE(hfiltI64f_)
      IF (ALLOCATED(hfiltI32f_)) DEALLOCATE(hfiltI32f_)
      IF (ALLOCATED(nzReIndices_)) DEALLOCATE(nzReIndices_)
      IF (ALLOCATED(nzImIndices_)) DEALLOCATE(nzImIndices_)
      IF (ALLOCATED(sparseHfiltR64f_)) DEALLOCATE(sparseHfiltR64f_)
      IF (ALLOCATED(sparseHfiltI64f_)) DEALLOCATE(sparseHfiltI64f_)
      IF (ALLOCATED(sparseHfiltR32f_)) DEALLOCATE(sparseHfiltR32f_)
      IF (ALLOCATED(sparseHfiltI32f_)) DEALLOCATE(sparseHfiltI32f_)
      lisTypeIII_ = .TRUE.
      nnzReCoeffs_ = 0
      nnzImCoeffs_ = 0
      nReCoeffs_ = 0
      nImCoeffs_ = 0
      nCoeffs_ = 0
      nOrder_ = 0
      ftype_ = XCLOC_SPXC_DONOT_FILTER
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Designs the averaging coefficients for the RMS filter.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_rmsFilter_design(ierr)
      INTEGER, INTENT(OUT) :: ierr
      ierr = 0
      ALLOCATE(rms64f_(nCoeffs_))
      ALLOCATE(rms32f_(nCoeffs_))
      rms64f_(1:nCoeffs_) = 1.d0/DBLE(nCoeffs_)
      rms32f_(1:nCoeffs_) = SNGL(rms64f_(1:nCoeffs_))
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Designs the FIR Hilbert transformer.
!>    @param[out] ierr 0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_envelope_design(ierr)
      INTEGER, INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: kaiser(:), sinct(:), t(:)
      DOUBLE PRECISION di, dn, gain, xfact
      INTEGER i
      DOUBLE PRECISION, PARAMETER :: beta = 8.d0
      DOUBLE PRECISION, PARAMETER :: fc = 1.d0
      DOUBLE PRECISION, PARAMETER :: fc2 = 0.5d0*fc
      ierr = 0
      lisTypeIII_ = .TRUE.
      IF (MOD(nOrder_, 2) /= 0) THEN
         lisTypeIII_ = .FALSE.
         WRITE(ERROR_UNIT,900) 
         ierr = 1
         RETURN
      ENDIF 
      ! Create the Kaiser window
      ALLOCATE(kaiser(nCoeffs_)); kaiser(:) = 1.d0
      xfact = 2.d0/DBLE(nCoeffs_ - 1)*beta
      ierr = ippsWinKaiser_64f_I(kaiser, nCoeffs_, xfact) ! only good to 7 digits?
      IF (ierr /= ippStsNoErr) THEN
         WRITE(ERROR_UNIT,905)
         ierr = 1
         RETURN
      ENDIF
      ! Create the ideal lowpass filter
      ALLOCATE(t(nCoeffs_))
      dn = DBLE(nCoeffs_)
      DO i=1,nCoeffs_
         di = DBLE(i - 1)
         t(i) = fc2*(0.5d0*(1.d0 - dn) + di) 
      ENDDO
      ALLOCATE(sinct(nCoeffs_))
      CALL sinc(nCoeffs_, t, sinct)
      ! Create the complex valued FIR coefficients with 12.66 of O & S:
      ! hfilt = sinc(t)*exp(i*pi*t) is what matlab uses.  However,
      ! Oppenheim has sin*sinc in 12.67.  i'll stick with matlab for consistency.
      ALLOCATE(hfiltR64f_(nCoeffs_))
      ALLOCATE(hfiltI64f_(nCoeffs_))
      nReCoeffs_ = nCoeffs_
      nImCoeffs_ = nCoeffs_
      DO i=1,nCoeffs_
         hfiltR64f_(i) = kaiser(i)*sinct(i)*DCOS(pi*t(i))
         hfiltI64f_(i) = kaiser(i)*sinct(i)*DSIN(pi*t(i))
      ENDDO
      gain = SUM(hfiltR64f_)
      hfiltR64f_(:) = hfiltR64f_(:)/gain
      hfiltI64f_(:) = hfiltI64f_(:)/gain
      ! Fix the Type IV case
      IF (lisTypeIII_) THEN
         DO i=1,nCoeffs_
            hfiltR64f_(i) = 0.d0
            IF (i == (nCoeffs_+1)/2) hfiltR64f_(i) = 1.d0
            IF (MOD(i, 2) == 1)      hfiltI64f_(i) = 0.d0
         ENDDO
         ALLOCATE(hfiltR32f_(nCoeffs_)); hfiltR32f_(:) = SNGL(hfiltR64f_(:))
         ALLOCATE(hfiltI32f_(nCoeffs_)); hfiltI32f_(:) = SNGL(hfiltI64f_(:))
         ! Sparsify coefficients
         nnzReCoeffs_ = 1
         nnzImCoeffs_ = (nCoeffs_ + 1)/2 - 1
         ALLOCATE(nzReIndices_(nnzReCoeffs_))
         ALLOCATE(nzImIndices_(nnzImCoeffs_))
         nzReIndices_(1) = (nCoeffs_+1)/2
         DO i=1,nnzImCoeffs_
            nzImIndices_(i) = 2*i
         ENDDO
         ALLOCATE(sparseHfiltR64f_(nnzReCoeffs_))
         ALLOCATE(sparseHfiltI64f_(nnzImCoeffs_))
         ALLOCATE(sparseHfiltR32f_(nnzReCoeffs_))
         ALLOCATE(sparseHfiltI32f_(nnzImCoeffs_))
         sparseHfiltR64f_(:) = hfiltR64f_(nzReIndices_(:))
         sparseHfiltI64f_(:) = hfiltI64f_(nzImIndices_(:))
         sparseHfiltR32f_(:) = SNGL(sparseHfiltR64f_(:))
         sparseHfiltI32f_(:) = SNGL(sparseHfiltI64f_(:))
         nzReIndices_(:) = nzReIndices_(:) - 1 ! C index
         nzImIndices_(:) = nzImIndices_(:) - 1 ! C index
      ELSE
         WRITE(ERROR_UNIT,900) 
         ierr = 1
         RETURN
      ENDIF
!do i=1,nCoeffs_
!print *, hfiltR64f_(i), hfiltI64f_(i)
!enddo
!print *, nzImIndices_
!print *, nnzReCoeffs_, nnzImCoeffs_
!print *, nzReIndices_
!print *, nzImIndices_
!print *, SNGL(hfiltR64f_)
!print *, Sngl(HFILTI64f_)
!print *, sparseHfiltR64f_
!print *, sparseHfiltI64f_ 
  900 FORMAT('xcloc_spxc_envelope_design: Type IV filter not implemented')
  905 FORMAT('xcloc_spxc_envelope_design: Failed making kaiser window')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the upper envelope of the data.
!>    @param[in] lds       Leading dimension of data.  This must be at least npts.
!>    @param[in] npts      Number of points in each signal.
!>    @param[in] nsignals  Number of signals to filter.
!>    @param[in,out] x     On input this is the column major data matrix with dimension
!>                         [lds x nsignals].
!>    @param[in,out] x     On exit the input is replaced by the the upper envelope of
!>                         the data.
!>    @param[out] ierr     0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_envelope_apply64f(lds, npts, nsignals, x, ierr)
      USE ISO_C_BINDING
      INTEGER, INTENT(IN) :: nsignals, lds, npts
      REAL(C_DOUBLE), INTENT(INOUT) :: x(lds*nsignals)
      INTEGER, INTENT(OUT) :: ierr
      INTERFACE
         INTEGER(C_INT) &
         FUNCTION xcloc_firFilter_envelope64f(lds, npts, nsignals, &
                                              nReCoeffs, reCoeffs, &
                                              nImCoeffs, imCoeffs, &
                                              x) &
         BIND(C, NAME='xcloc_firFilter_envelope64f')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), VALUE, INTENT(IN) :: nsignals, lds, npts,  &
                                              nReCoeffs, nImCoeffs
         REAL(C_DOUBLE), INTENT(IN) :: reCoeffs(nReCoeffs), imCoeffs(nImCoeffs)
         REAL(C_DOUBLE), INTENT(INOUT) :: x(lds*npts)
         END FUNCTION
      END INTERFACE
      ierr = xcloc_firFilter_envelope64f(lds, npts, nsignals, &
                                         nReCoeffs_, hfiltR64f_, &
                                         nImCoeffs_, hfiltI64f_, &
                                         x)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
      ENDIF
  900 FORMAT('xcloc_spxc_envelope_apply64f: Failed to compute envelope')
      RETURN
      END
!>    @brief Computes the upper envelope of the data
!>    @param[in] lds       Leading dimension of data.  This must be at least npts.
!>    @param[in] npts      Number of points in each signal.
!>    @param[in] nsignals  Number of signals to filter.
!>    @param[in,out] x     On input this is the column major data matrix with dimension
!>                         [lds x nsignals].
!>    @param[in,out] x     On exit the input is replaced by the the upper envelope of
!>                         the data.
!>    @param[out] ierr     0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_envelope_apply32f(lds, npts, nsignals, x, ierr)
      USE ISO_C_BINDING
      INTEGER, INTENT(IN) :: nsignals, lds, npts
      REAL(C_FLOAT), INTENT(INOUT) :: x(lds*nsignals)
      INTEGER, INTENT(OUT) :: ierr
      INTERFACE
         INTEGER(C_INT) &
         FUNCTION xcloc_firFilter_envelope32f(lds, npts, nsignals, &
                                              nnzReCoeffs, nzReCoeffs, reCoeffs, &
                                              nnzImCoeffs, nzImCoeffs, imCoeffs, &
                                              x) &
         BIND(C, NAME='xcloc_firFilter_envelope32f')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), VALUE, INTENT(IN) :: nsignals, lds, npts,  &
                                              nnzReCoeffs, nnzImCoeffs
         INTEGER(C_INT), INTENT(IN) :: nzReCoeffs(nnzReCoeffs), nzImCoeffs(nnzImCoeffs)
         REAL(C_FLOAT), INTENT(IN) :: reCoeffs(nnzReCoeffs), imCoeffs(nnzImCoeffs)
         REAL(C_FLOAT), INTENT(INOUT) :: x(lds*npts)
         END FUNCTION
      END INTERFACE
      ierr = xcloc_firFilter_envelope32f(lds, npts, nsignals, &
                                         nnzReCoeffs_, nzReIndices_, sparseHfiltR32f_, &
                                         nnzImCoeffs_, nzImIndices_, sparseHfiltI32f_, &
                                         x)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1 
      ENDIF
  900 FORMAT('xcloc_spxc_envelope_apply32f: Failed to compute envelope')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Applies the RMS filter to the data.
!>    @param[in] lds       Leading dimension of data.  This must be at least npts.
!>    @param[in] npts      Number of points in each signal.
!>    @param[in] nsignals  Number of signals to filter.
!>    @param[in,out] x     On input this is the column major data matrix with dimension
!>                         [lds x nsignals].
!>    @param[in,out] x     On exit the input is replaced by the the RMS filtered data.
!>    @param[out] ierr     0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_rmsFilter_apply64f(lds, npts, nsignals, x, ierr)
      USE ISO_C_BINDING
      INTEGER, INTENT(IN) :: nsignals, lds, npts
      REAL(C_DOUBLE), INTENT(INOUT) :: x(lds*nsignals)
      INTEGER, INTENT(OUT) :: ierr
      INTERFACE
         INTEGER(C_INT) &
         FUNCTION xcloc_firFilter_rmsFilter64f(lds, npts, nsignals, &
                                               tapsLen, taps, x)    &   
         BIND(C, NAME='xcloc_firFilter_rmsFilter64f')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), VALUE, INTENT(IN) :: nsignals, lds, npts, tapsLen
         REAL(C_DOUBLE), INTENT(IN) :: taps(tapsLen)
         REAL(C_DOUBLE), INTENT(INOUT) :: x(lds*npts)
         END FUNCTION
      END INTERFACE
      ierr = xcloc_firFilter_rmsFilter64f(lds, npts, nsignals, nCoeffs_, &
                                          rms64f_, x)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
      ENDIF
  900 FORMAT('xcloc_spxc_rmsFilter_apply64f: Failed to apply RMS filter')
      RETURN
      END
!>    @brief Applies the RMS filter to the data.
!>    @param[in] lds       Leading dimension of data.  This must be at least npts.
!>    @param[in] npts      Number of points in each signal.
!>    @param[in] nsignals  Number of signals to filter.
!>    @param[in,out] x     On input this is the column major data matrix with dimension
!>                         [lds x nsignals].
!>    @param[in,out] x     On exit the input is replaced by the the RMS filtered data.
!>    @param[out] ierr     0 indicates success.
!>    @ingroup spxc
      SUBROUTINE xcloc_spxc_rmsFilter_apply32f(lds, npts, nsignals, x, ierr)
      USE ISO_C_BINDING
      INTEGER, INTENT(IN) :: nsignals, lds, npts
      REAL(C_FLOAT), INTENT(INOUT) :: x(lds*nsignals)
      INTEGER, INTENT(OUT) :: ierr
      INTERFACE
         INTEGER(C_INT) &
         FUNCTION xcloc_firFilter_rmsFilter32f(lds, npts, nsignals, &
                                               tapsLen, taps, x)    &
         BIND(C, NAME='xcloc_firFilter_rmsFilter32f')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), VALUE, INTENT(IN) :: nsignals, lds, npts, tapsLen
         REAL(C_FLOAT), INTENT(IN) :: taps(tapsLen)
         REAL(C_FLOAT), INTENT(INOUT) :: x(lds*npts)
         END FUNCTION
      END INTERFACE
      ierr = xcloc_firFilter_rmsFilter32f(lds, npts, nsignals, nCoeffs_, &
                                          rms32f_, x)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
      ENDIF
  900 FORMAT('xcloc_spxc_rmsFilter_apply32f: Failed to apply RMS filter')
      RETURN
      END

!========================================================================================!
!                              Private functions                                         !
!========================================================================================!
!>    @brief Computes \f$ \mathrm{sinc}(x) = \frac{\sin(x)}{\pi x} \f$.
!>    @param[in] n    Length of vector x.
!>    @param[in] x    Dependent variable at which to compute sinc.
!>    @param[out] sc  sinc(x).
!>    @ingroup spxc
      SUBROUTINE SINC(n, x, sc)
      INTEGER, INTENT(IN) :: n
      DOUBLE PRECISION, INTENT(IN) :: x(n)
      DOUBLE PRECISION, INTENT(OUT) :: sc(n)
      DOUBLE PRECISION pix
      INTEGER i
      DO i=1,n
         pix = pi*x(i)
         sc(i) = 1.d0 ! sinc will evaluate to 1 if x is 0
         IF (x(i) /= 0.d0) sc(i) = DSIN(pix)/pix
      ENDDO
      RETURN
      END
END MODULE
