!> @defgroup xcloc xcloc
!> @brief The xcloc serial library.
!> @author Ben Baker 
!> @copyright Ben Baker distributed under the MIT license.

!> @defgroup xcloc_xcloc Cross-Correlation Based Event Location 
!> @brief Cross-correlation event location utilities.
!> @ingroup xcloc
!> @copyright Ben Baker distributed under the MIT license. 
MODULE XCLOC
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      USE XCLOC_FDXC
      USE XCLOC_DSMXC
      USE XCLOC_SPXC
      USE XCLOC_CONSTANTS
      IMPLICIT NONE

      !> @ingroup xcloc
      !> The number of cross-correlograms.
      INTEGER, PRIVATE, SAVE :: nptsInXCs_ = 0
      !> @ingroup xcloc
      !> The number of cross-correlations pairs.
      INTEGER, PRIVATE, SAVE :: nxcs_ = 0 
      !> @ingroup xcloc
      INTEGER, PRIVATE, SAVE :: ngrd_ = 0
      !> Signal to migrate.
      INTEGER, PRIVATE, SAVE :: signalMigrationType_ = XCLOC_MIGRATE_PHASE_XCS
      !> Controls verbosity.
      INTEGER, PRIVATE, SAVE :: verbose_ = XCLOC_PRINT_WARNINGS
      !> The precision of the module.
      INTEGER, PRIVATE, SAVE :: precision_ = XCLOC_SINGLE_PRECISION
      !> A flag indicating whether or not xcloc has been initialized.
      LOGICAL, PRIVATE, SAVE :: linit_ = .FALSE.
      !> Determines whether or not the signals have been set.
      LOGICAL, PRIVATE, SAVE :: lhaveSignals_ = .FALSE.
      !> A flag determining if the travel time tables have been set.
      LOGICAL(C_BOOL), PRIVATE, SAVE :: lhaveAllTables_ = .FALSE.
      !> A flag determining if the migration image is computed.
      LOGICAL, PRIVATE, SAVE :: lhaveImage_ = .FALSE.
      !> A flag indicating the cross-correlograms are computed.
      LOGICAL, PRIVATE, SAVE :: lhaveXCs_ = .FALSE.

      PUBLIC :: xcloc_initialize
      PUBLIC :: xcloc_finalize
      PUBLIC :: xcloc_signalToTableIndex
      PUBLIC :: xcloc_setTable64f
      PUBLIC :: xcloc_setTable32f
      PUBLIC :: xcloc_setSignals64f
      PUBLIC :: xcloc_setSignals32f
      PUBLIC :: xcloc_setSignalToMigrate
      PUBLIC :: xcloc_setXCFilter
      PUBLIC :: xcloc_getPrecision
      PUBLIC :: xcloc_compute
      CONTAINS
!========================================================================================!
!                                     Begin the Code                                     !
!========================================================================================! 
!>    @brief Initializes xcloc.
!>    @param[in] npts      Number of points in input signals.
!>    @param[in] nptsPad   A tuning parameter to mitigate DFT lengths that could
!>                         potentially be large semi-prime numbers.
!>    @param[in] nxcs      Number of cross-correlograms.
!>    @param[in] s2m       Signal type to migrate.  This can be XCLOC_MIGRATE_PHASE_XCS
!>                         to migrate phase correlograms or XCLOC_MIGRATE_XCS
!>                         to migrate cross-correlograms.
!>    @param[in] dt        Sampling period (seconds) of signals.
!>    @param[in] ngrd      Number of points in migration grid.
!>    @param[in] nfcoeffs  The number of filter coefficients.  If the processing is
!>                         to be done then this must be odd.
!>    @param[in] ftype     XCLOC_SPXC_DONOT_FILTER will not filter correlograms.
!>    @param[in] ftype     XCLOC_SPXC_ENVELOPE_FILTER will apply envelope to correlograms.
!>    @param[in] ftype     XCLOC_SPXC_RMS_FILTER will comptue RMS of correlograms.
!>    @param[in] xcPairs   This is a [2 x nxcs] matrix in column major format where 
!>                         the indices, (2*(ixc-1)+1, 2*(ixc-1)+2), map to the (i,j)'th
!>                         signal pair comprising a correlation.
!>    @param[in] verbose   Controls the verbosity of the module.  0 is quiet.
!>    @param[in] prec      Controls the precision of the module.  
!>    @param[in] accuracy  Controls the accuracy of the vector calculations in MKL.
!>    @param[out] ierr     0 indicates success.
!>    @ingroup xcloc_xcloc
      SUBROUTINE xcloc_initialize(npts, nptsPad, nxcs,           &
                                  s2m, dt, ngrd,                 &
                                  nfcoeffs, ftype,               &
                                  xcPairs,                       &
                                  verbose, prec, accuracy, ierr) &
      BIND(C, NAME='xcloc_initialize')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: npts, nptsPad, nxcs, ngrd, nfcoeffs, ftype, &
                                           s2m, verbose, prec, accuracy
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt
      INTEGER(C_INT), INTENT(IN) :: xcPairs(2*nxcs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER nptsInXCs
      ierr = 0
      CALL xcloc_finalize()
      CALL xcloc_setSignalToMigrate(s2m, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,895)
         CALL xcloc_finalize()
         RETURN
      ENDIF
      ! Initialize the cross-correlogram calculator.
      CALL xcloc_fdxc_initialize(npts, nptsPad, nxcs, xcPairs, &
                                 verbose, prec, accuracy, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         CALL xcloc_finalize()
         RETURN
      ENDIF
      ! Get the length of the correlograms.
      CALL xcloc_fdxc_getCorrelogramLength(nptsInXCs, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,901)
         CALL xcloc_finalize()
         RETURN
      ENDIF
      ! Initialize the processing module.
      CALL xcloc_setXCFilter(nfcoeffs, ftype, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,902)
         CALL xcloc_finalize()
         RETURN
      ENDIF
      ! Initialize the diffraction stack migration module.
      CALL xcloc_dsmxc_initialize(ngrd, nxcs, nptsInXCs,    &
                                  dt, xcPairs, verbose, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,905)
         CALL xcloc_finalize()
         RETURN
      ENDIF
      CALL xcloc_fdxc_getNumberOfCorrelograms(nxcs_, ierr)
      CALL xcloc_fdxc_getCorrelogramLength(nptsInXCs_, ierr)
      CALL xcloc_dsmxc_getNumberOGridPointsInTable(ngrd_)
      ! Set some basic info
      precision_ = prec
      verbose_ = verbose
  895 FORMAT('xcloc_initialize: Failed to set signal migration type')
  900 FORMAT('xcloc_initialize: Failed to initialize fdxc')
  901 FORMAT('xcloc_initialize: Failed to get length of correlograms')
  902 FORMAT('xcloc_initialize: Failed to initialize spxc') 
  905 FORMAT('xcloc_initialize: Failed to initialize dsmxc')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Initializes the correlogram filter.
!>    @param[in] nfcoeffs  The number of filter coefficients.  If the processing is
!>                         to be done then this must be odd.
!>    @param[in] ftype     XCLOC_SPXC_DONOT_FILTER will not filter correlograms.
!>    @param[in] ftype     XCLOC_SPXC_ENVELOPE_FILTER will apply envelope to correlograms.
!>    @param[in] ftype     XCLOC_SPXC_RMS_FILTER will comptue RMS of correlograms.
!>    @param[out] ierr     0 indicates success.
!>    @ingroup xcloc_xcloc
      SUBROUTINE xcloc_setXCFilter(nfcoeffs, ftype, ierr) &
      BIND(C, NAME='xcloc_setXCFilter')
      INTEGER(C_INT), VALUE, INTENT(IN) :: nfcoeffs, ftype
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL xcloc_spxc_finalize()
      CALL xcloc_spxc_initialize(nfcoeffs, ftype, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
      ENDIF
  900 FORMAT('xcloc_setXCFilter: Failed to set filter for correlograms')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Maps a signal number to the table index.
!>    @param[in] is     Signal number.  This must be in the xcPairs table.
!>    @param[out] it    The table index corresponding to is.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup xcloc_xcloc 
      SUBROUTINE xcloc_signalToTableIndex(is, it, ierr) &
      BIND(C, NAME='xcloc_signalToTableIndex')
      INTEGER(C_INT), VALUE, INTENT(IN) :: is
      INTEGER(C_INT), INTENT(OUT) :: it, ierr
      CALL xcloc_dsmxc_signalToTableIndex(is, it, ierr) 
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900) is
         it = 0
      ENDIF
  900 FORMAT('xcloc_signalToTableIndex: Failed map signal', I4, ' to table')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the it'th signal corresponding to the is'th signal.  One should call
!>           xcloc_signalToTableIndex to get the appropriate table index for this signal.
!>    @param[in] tableNumber  Table index to set.  The table number for the is'th
!>                            signal can be determined with xcloc_signalToTableIndex.
!>    @param[in] ngrd         Number of grid points in the table.  This must equal ngrd_.
!>    @param[in] table        The travel times from the source to all points in the
!>                            grid in seconds.  This has dimension [ngrd].
!>    @param[out] ierr        0 indicates success.
!>    @ingroup xcloc_xcloc 
      SUBROUTINE xcloc_setTable64f(tableNumber, ngrd, table, ierr) &
      BIND(C, NAME='xcloc_setTable64f')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: tableNumber, ngrd
      REAL(C_DOUBLE), INTENT(IN) :: table(ngrd)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL xcloc_dsmxc_setTable64f(tableNumber, ngrd, table, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900) tableNumber
         ierr = 1
      ENDIF
      CALL xcloc_dsmxc_haveAllTables(lhaveAllTables_)
      IF (lhaveAllTables_ .AND. verbose_ > XCLOC_PRINT_WARNINGS) WRITE(OUTPUT_UNIT,905)
  900 FORMAT('xcloc_setTable64f: Failed to set table number', I4)
  905 FORMAT('xcloc_setTable64f: All tables set')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the it'th signal corresponding to the is'th signal.  One should call
!>           xcloc_signalToTableIndex to get the appropriate table index for this signal.
!>    @param[in] tableNumber  Table index to set.  The table number for the is'th
!>                            signal can be determined with xcloc_signalToTableIndex.
!>    @param[in] ngrd         Number of grid points in the table.  This must equal ngrd_.
!>    @param[in] table        The travel times from the source to all points in the
!>                            grid in seconds.  This has dimension [ngrd].
!>    @param[out] ierr        0 indicates success.
!>    @ingroup xcloc_xcloc 
      SUBROUTINE xcloc_setTable32f(tableNumber, ngrd, table, ierr) &
      BIND(C, NAME='xcloc_setTable32f')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: tableNumber, ngrd
      REAL(C_FLOAT), INTENT(IN) :: table(ngrd)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      CALL xcloc_dsmxc_setTable32f(tableNumber, ngrd, table, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900) tableNumber
         ierr = 1
      ENDIF
      CALL xcloc_dsmxc_haveAllTables(lhaveAllTables_)
      IF (lhaveAllTables_ .AND. verbose_ > XCLOC_PRINT_WARNINGS) WRITE(OUTPUT_UNIT,905)
  900 FORMAT('xcloc_setTable64f: Failed to set table number', I4)
  905 FORMAT('xcloc_setTable64f: All tables set')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @param[in] signalToMigrate  XCLOC_MIGRATE_PHASE_XCS indicates that phase
!>                                correlograms will be migrated.
!>    @param[in] signalToMigrate  XCLOC_MIGRATE_XCS indicates that cross correlograms
!>                                will be migrated.
!>    @param[out] ierr            0 indicates success.
!>    @ingroup xcloc_xcloc
      SUBROUTINE xcloc_setSignalToMigrate(signalToMigrate, ierr) &
      BIND(C, NAME='xcloc_setSignalToMigrate')
      INTEGER(C_INT), VALUE, INTENT(IN) :: signalToMigrate
      INTEGER(C_INT), INTENT(OUT) :: ierr
      LOGICAL lisValid
      ierr = 0
      signalMigrationType_ = XCLOC_MIGRATE_PHASE_XCS
      lisValid = xcloc_constants_isValidSignalToMigrate(signalToMigrate)
      IF (.NOT.lisValid) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      signalMigrationType_ = signalToMigrate
  900 FORMAT('xcloc_setSignalToMigrate: Invalid type of signal to migrate')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief After setting the travel time tables and signals this will:
!>           1. Compute the phase or cross-correlograms depnding on the signal migration
!>              policy.
!>           2. Either compute the envelope or RMS of the correlograms or leave the
!>              correlograms unchanged.
!>           3. Migrate the (processed) correlograms.
!>    @param[out] ierr   0 indicates successs.
!>    @ingroup xcloc_xcloc
      SUBROUTINE xcloc_compute(ierr) &
      BIND(C, NAME='xcloc_compute')
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xcs64f
      REAL, ALLOCATABLE, DIMENSION(:) :: xcs32f
      lhaveImage_ = .FALSE.
      lhaveXCs_ = .FALSE.
      IF (.NOT.linit_) THEN
         WRITE(ERROR_UNIT,850)
         ierr = 1 
         RETURN
      ENDIF
      IF (.NOT.lhaveSignals_) THEN
         WRITE(ERROR_UNIT,855)
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.lhaveAllTables_) THEN
         WRITE(ERROR_UNIT,860)
         ierr = 1
         RETURN
      ENDIF
      ! Compute the correlograms.
      IF (signalMigrationType_ == XCLOC_MIGRATE_PHASE_XCS) THEN
         CALL xcloc_fdxc_computePhaseCorrelograms(ierr)
      ELSE
         CALL xcloc_fdxc_computeCrossCorrelograms(ierr)
      ENDIF
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,865)
         ierr = 1
         RETURN
      ENDIF
      lhaveXCs_ = .TRUE.
      ! Get the correlograms, filter them, and set the filters xcs on the DSM module 
      IF (precision_ == XCLOC_SINGLE_PRECISION) THEN
         ALLOCATE(xcs32f(nxcs_*nptsInXCs_))
         CALL xcloc_fdxc_getCorrelograms32f(nptsInXCs_, nxcs_, xcs32f, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,870)
            ierr = 1
            RETURN
         ENDIF
         CALL xcloc_spxc_filterXCsInPlace32f(nptsInXCs_, nptsInXCs_, nxcs_, xcs32f, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,875)
            ierr = 1
            RETURN
         ENDIF
         CALL xcloc_dsmxc_setCorrelograms32f(nptsInXCs_, nptsInXCs_, nxcs_, xcs32f, ierr)
      ELSE
         ALLOCATE(xcs64f(nxcs_*nptsInXCs_))
         CALL xcloc_fdxc_getCorrelograms64f(nptsInXCs_, nxcs_, xcs64f, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,870)
            ierr = 1
            RETURN
         ENDIF
         CALL xcloc_spxc_filterXCsInPlace64f(nptsInXCs_, nptsInXCs_, nxcs_, xcs64f, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,875)
            ierr = 1
            RETURN
         ENDIF
         CALL xcloc_dsmxc_setCorrelograms64f(nptsInXCs_, nptsInXCs_, nxcs_, xcs64f, ierr)
      ENDIF
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,880)
         ierr = 1
         RETURN
      ENDIF
      ! Migrate the signals.
      CALL xcloc_dsmxc_compute(ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,885)
         ierr = 1
         RETURN
      ENDIF 
      lhaveImage_ = .TRUE.
  850 FORMAT('xcloc_compute: Module not initialized')
  855 FORMAT('xcloc_compute: Signals not yet set')
  860 FORMAT('xcloc_compute: Travel time tables not yet set')
  865 FORMAT('xcloc_compute: Failed to compute correlograms')
  870 FORMAT('xcloc_compute: Failed to get correlograms')
  875 FORMAT('xcloc_compute: Failed to process correlograms')
  880 FORMAT('xcloc_compute: Failed to get correlograms')
  885 FORMAT('xcloc_compute: Failed to compute image')
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
!>    @ingroup xcloc_xcloc
      SUBROUTINE xcloc_setSignals64f(ldx, npts, nsignals, x, ierr) &
      BIND(C, NAME='xcloc_setSignals64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldx, npts, nsignals
      REAL(C_DOUBLE), INTENT(IN) :: x(ldx*nsignals)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      lhaveSignals_ = .FALSE.
      IF (.NOT.linit_) THEN
         WRITE(ERROR_UNIT,850)
         ierr = 1
         RETURN
      ENDIF
      CALL xcloc_fdxc_setSignals64f(ldx, npts, nsignals, x, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      lhaveSignals_ = .TRUE.
  850 FORMAT('xcloc_setSignals64f: Module not initialized')
  900 FORMAT('xcloc_setSignals64f: Failed to set signals') 
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
!>    @ingroup xcloc_xcloc
      SUBROUTINE xcloc_setSignals32f(ldx, npts, nsignals, x, ierr) &
      BIND(C, NAME='xcloc_setSignals32f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldx, npts, nsignals
      REAL(C_DOUBLE), INTENT(IN) :: x(ldx*nsignals)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      lhaveSignals_ = .FALSE.
      IF (.NOT.linit_) THEN
         WRITE(ERROR_UNIT,850)
         ierr = 1 
         RETURN
      ENDIF
      CALL xcloc_fdxc_setSignals64f(ldx, npts, nsignals, x, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      lhaveSignals_ = .TRUE.
  850 FORMAT('xcloc_setSignals64f: Module not initialized')
  900 FORMAT('xcloc_setSignals32f: Failed to set signals')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases memory on the xcloc module.
!>    @ingroup xcloc_xcloc
      SUBROUTINE xcloc_finalize() &
      BIND(C, NAME='xcloc_finalize')
      CALL xcloc_fdxc_finalize()
      CALL xcloc_spxc_finalize()
      CALL xcloc_dsmxc_finalize()
      nxcs_ = 0
      nptsInXCs_ = 0
      ngrd_ = 0
      signalMigrationType_ = XCLOC_MIGRATE_PHASE_XCS
      precision_ = XCLOC_SINGLE_PRECISION
      verbose_ = XCLOC_PRINT_WARNINGS
      lhaveSignals_ = .FALSE.
      lhaveImage_ = .FALSE.
      lhaveXCs_ = .FALSE.
      linit_ = .FALSE.
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the precision of the module.
!>    @param[out] prec   Precision of module - e.g., XCLOC_SINGLE_PRECISION or 
!>                       XCLOC_DOUBLE_PRECISION.
!>    @ingroup xcloc_xcloc
      SUBROUTINE xcloc_getPrecision(prec) &
      BIND(C, NAME='xcloc_getPrecision')
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: prec
      prec = precision_
      RETURN
      END
END MODULE
