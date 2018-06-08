!> @brief Performs some signals processing on the cross-correlograms.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_SPXC
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS

      !> Real Hilbert transformer coefficients.
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, SAVE :: hfiltR64f_(:)
      !> Imaginary Hilbert transformer coefficients.
      REAL(C_DOUBLE), PRIVATE, ALLOCATABLE, SAVE :: hfiltI64f_(:)
      !> Real Hilbert transformer coefficients for single precision.
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, SAVE :: hfiltR32f_(:)
      !> Imaginary Hilbert transformer coefficients for single precision.
      REAL(C_FLOAT), PRIVATE, ALLOCATABLE, SAVE :: hfiltI32f_(:)
      !> Accuracy.
      INTEGER(C_INT), PRIVATE, SAVE :: accuracy_ = XCLOC_HIGH_ACCURACY
      !> Filtering type.
      INTEGER(C_INT), PRIVATE, SAVE :: ftype_ = XCLOC_SPXC_DONOT_FILTER
      !> Filter order.
      INTEGER(C_INT), PRIVATE, SAVE :: norder_ = 0
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
      RETURN
      END

      SUBROUTINE xcloc_spxc_rmsFilter_design( )

      RETURN
      END
      !> @brief Designs the hilbert transformer.
      SUBROUTINE xcloc_spxc_envelope_design(n )

      RETURN
      END
      !> @brief Applies the envelope.
      SUBROUTINE xcloc_spxc_envelope_apply( )

      RETURN
      END


END MODULE
