!> @brief Performs some signals processing on the cross-correlograms.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_SPXC
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
      USE XCLOC_UTILS

      !> Real Hilbert transformer coefficients.
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, SAVE :: hfiltR64f_(:)
      !> Imaginary Hilbert transformer coefficients.
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, SAVE :: hfiltI64f_(:)
      !> Real Hilbert transformer coefficients for single precision.
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, SAVE :: hfiltR32f_(:)
      !> Imaginary Hilbert transformer coefficients for single precision.
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, SAVE :: hfiltI32f_(:)
      !> FIR averaging filter coefficients.
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, SAVE :: rms64f_(:)
      !> FIR averaging filter for single precision. 
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, SAVE :: rms32f_(:)
      !> Precision (single or double precision).
      INTEGER(C_INT), PRIVATE, SAVE :: precision_ = XCLOC_SINGLE_PRECISION
      !> Accuracy.
      INTEGER(C_INT), PRIVATE, SAVE :: accuracy_ = XCLOC_HIGH_ACCURACY
      !> Filtering type.
      INTEGER(C_INT), PRIVATE, SAVE :: ftype_ = XCLOC_SPXC_DONOT_FILTER
      !> Number of filter coefficients.
      INTEGER(C_INT), PRIVATE, SAVE :: nCoeffs_ = 0
      !> Filter order.
      INTEGER(C_INT), PRIVATE, SAVE :: norder_ = 0

      PUBLIC :: xcloc_spxc_initialize
      PUBLIC :: xcloc_spxc_finalize

      PRIVATE :: xcloc_spxc_rmsFilter_design
      PRIVATE :: xcloc_spxc_envelope_design
      PRIVATE :: xcloc_spxc_rmsFilter_apply64f
      PRIVATE :: xcloc_spxc_rmsFilter_apply32f
      CONTAINS
      !==================================================================================!
      !                                     Begin the Code                               !
      !==================================================================================!
      !> Initializes the processing module.
      !> @param[in] n         The number of filter coefficients.  If the processing is
      !>                      to be done then this must be odd.
      !> @param[in] ftype     XCLOC_SPXC_DONOT_FILTER will not filter correlograms.
      !> @param[in] ftype     XCLOC_SPXC_ENVELOPE_FILTER will apply envelope to
      !>                      correlograms.
      !> @param[in] ftype     XCLOC_SPXC_RMS_FILTER will comptue RMS of correlograms.
      !> @param[in] prec      XCLOC_SINGLE_PRECISION will create single precision filters.
      !> @param[in] prec      XCLOC_DOUBLE_PRECISION will create double precision filters.
      !> @param[in] accuracy  Controls the accuracy during the envelope computation.
      !> @param[out] ierr     0 indicates success.
      SUBROUTINE xcloc_spxc_initialize(n, ftype, prec, accuracy, ierr)
      INTEGER(C_INT), VALUE, INTENT(IN) :: n, ftype, prec, accuracy
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (ftype == XCLOC_SPXC_DONOT_FILTER) THEN
         ftype_ = ftype 
      ELSE
         ! Check inputs
         IF (n < 1) THEN
            WRITE(*,900)
            ierr = 1
            RETURN
         ENDIF
         ! Nothing to do
         IF (n == 1) THEN
            ftype_ = XCLOC_SPXC_DONOT_FILTER
            WRITE(*,901) 
            RETURN
         ENDIF
         IF (.NOT.xcloc_constants_isValidPrecision(prec)) THEN
            WRITE(*,902)
            ierr = 1
            RETURN
         ENDIF
         IF (.NOT.xcloc_constants_isValidAccuracy(accuracy)) THEN
            WRITE(*,903)
            ierr = 1
            RETURN
         ENDIF
         ! Require FIR filters have an odd number of coefficients for easy group delay
         precision_ = prec
         nCoeffs_ = n
         IF (MOD(n, 2) == 0) THEN
            WRITE(*,904)
            nCoeffs_ = n + 1
         ENDIF
          
         IF (ftype == XCLOC_SPXC_ENVELOPE_FILTER) THEN

         ELSEIF (ftype == XCLOC_SPXC_RMS_FILTER) THEN

         ELSE
            WRITE(*,905) ftype
            ierr = 1
            RETURN
         ENDIF
      ENDIF
  900 FORMAT('xcloc_spxc_initialize: Number of filter coefficients must be posisitve')
  901 FORMAT('xcloc_spxc_initialize: nTaps=1; skipping filtering')
  902 FORMAT('xcloc_spxc_initialize: Invalid precision')
  903 FORMAT('xcloc_spxc_initialize: Invalid accuracy')
  904 FORMAT('xcloc_spxc_initialize: Adding one filter coefficient')
  905 FORMAT('xcloc_spxc_initialize: Invalid filter type', I4)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases memory and resets variables on the XC signals processing module.
      SUBROUTINE xcloc_spxc_finalize( )
      IF (ALLOCATED(rms64f_))    DEALLOCATE(rms64f_)
      IF (ALLOCATED(rms32f_))    DEALLOCATE(rms32f_)
      IF (ALLOCATED(hfiltR64f_)) DEALLOCATE(hfiltR64f_) 
      IF (ALLOCATED(hfiltR32f_)) DEALLOCATE(hfiltR32f_)
      IF (ALLOCATED(hfiltI64f_)) DEALLOCATE(hfiltI64f_)
      IF (ALLOCATED(hfiltI32f_)) DEALLOCATE(hfiltI32f_)
      nCoeffs_ = 0
      ftype_ = XCLOC_SPXC_DONOT_FILTER
      precision_ = XCLOC_SINGLE_PRECISION
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Designs the averaging coefficients for the RMS filter.
!>    @param[out] ierr  0 indicates success.
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
      !> @brief Designs the hilbert transformer.
      !> @param[out] ierr 0 indicates success.
      SUBROUTINE xcloc_spxc_envelope_design(ierr)
      INTEGER, INTENT(OUT) :: ierr
      ierr = 0
      RETURN
      END

      !> @brief Applies the envelope.
      SUBROUTINE xcloc_spxc_envelope_apply()

      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Applies the RMS filter to the data.
!>    @param[in] nsignals  Number of signals to filter.
!>    @param[in] lds       Leading dimension of data.  This must be at least npts.
!>    @param[in] npts      Number of points in each signal.
!>    @param[in,out] x     On input this is the column major data matrix with dimension
!>                         [lds x nsignals].
!>    @param[in,out] x     On exit the input is replaced by the the RMS filtered data.
!>    @param[out] ierr     0 indicates success.
      SUBROUTINE xcloc_spxc_rmsFilter_apply64f(nsignals, lds, npts, x, ierr)
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nsignals, lds, npts
      REAL(C_DOUBLE), INTENT(INOUT) :: x(lds*nsignals)
      INTEGER, INTENT(OUT) :: ierr
      INTERFACE
         INTEGER(C_INT) &
         FUNCTION xcloc_firFilter_rmsFilter64f(nsignals, lds, npts, &
                                               tapsLen, taps, x)    &   
         BIND(C, NAME='xcloc_firFilter_rmsFilter64f')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), VALUE, INTENT(IN) :: nsignals, lds, npts, tapsLen
         REAL(C_DOUBLE), INTENT(IN) :: taps(tapsLen)
         REAL(C_DOUBLE), INTENT(INOUT) :: x(lds*npts)
         END FUNCTION
      END INTERFACE
      ierr = xcloc_firFilter_rmsFilter64f(nsignals, lds, npts, nCoeffs_, &
                                          rms64f_, x)
      IF (ierr /= 0) THEN
         WRITE(*,900)
         ierr = 1
      ENDIF
  900 FORMAT('xcloc_spxc_rmsFilter_apply64f: Failed to apply RMS filter')
      RETURN
      END

!>    @brief Applies the RMS filter to the data.
!>    @param[in] nsignals  Number of signals to filter.
!>    @param[in] lds       Leading dimension of data.  This must be at least npts.
!>    @param[in] npts      Number of points in each signal.
!>    @param[in,out] x     On input this is the column major data matrix with dimension
!>                         [lds x nsignals].
!>    @param[in,out] x     On exit the input is replaced by the the RMS filtered data.
!>    @param[out] ierr     0 indicates success.
      SUBROUTINE xcloc_spxc_rmsFilter_apply32f(nsignals, lds, npts, x, ierr)
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nsignals, lds, npts
      REAL(C_FLOAT), INTENT(INOUT) :: x(lds*nsignals)
      INTEGER, INTENT(OUT) :: ierr
      INTERFACE
         INTEGER(C_INT) &
         FUNCTION xcloc_firFilter_rmsFilter32f(nsignals, lds, npts, &
                                               tapsLen, taps, x)    &
         BIND(C, NAME='xcloc_firFilter_rmsFilter32f')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), VALUE, INTENT(IN) :: nsignals, lds, npts, tapsLen
         REAL(C_FLOAT), INTENT(IN) :: taps(tapsLen)
         REAL(C_FLOAT), INTENT(INOUT) :: x(lds*npts)
         END FUNCTION
      END INTERFACE
      ierr = xcloc_firFilter_rmsFilter32f(nsignals, lds, npts, nCoeffs_, &
                                          rms32f_, x)
      IF (ierr /= 0) THEN
         WRITE(*,900)
         ierr = 1
      ENDIF
  900 FORMAT('xcloc_spxc_rmsFilter_apply32f: Failed to apply RMS filter')
      RETURN
      END


END MODULE
