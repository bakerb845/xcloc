!> @defgroup xcloc_mpi xclocMPI
!> @breif The parallel xcloc libary.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.

!> @defgroup xcloc_xcloc_mpi Parallel Cross-Correlation Based Event Location
!> @brief Parallel cross-correlation event location utilities.
!> @ingroup xcloc_mpi
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_MPI
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      USE MPI_F08
      USE XCLOC_SPXC
      USE XCLOC_CONSTANTS
      !> MPI communicator.
#if defined(__INTEL_COMPILER)
      TYPE(MPI_Comm), PRIVATE, SAVE :: comm_ = MPI_COMM_WORLD
#else
      TYPE(MPI_Comm), PRIVATE, SAVE :: comm_
#endif
      !> @ingroup xcloc_mpi
      !> Root ID on communicator.
      INTEGER(C_INT), PRIVATE, SAVE :: root_ = 0
      !> Number of processes on communicator.
      INTEGER(C_INT), PRIVATE, SAVE :: nprocs_ = 0
      !> @ingroup xcloc_mpi
      INTEGER(C_INT), PRIVATE, SAVE :: ftype_ = XCLOC_SPXC_DONOT_FILTER
      !> @ingroup xcloc_mpi
      !> Number of filtering taps.
      INTEGER(C_INT), PRIVATE, SAVE :: nfcoeffs_ = 0
      !> Accuracy of MKL computations.
      INTEGER, PRIVATE, SAVE :: accuracy_ = XCLOC_HIGH_ACCURACY
      !> Signal to migrate.
      INTEGER, PRIVATE, SAVE :: xcTypeToMigrate_ = XCLOC_MIGRATE_PHASE_XCS
      !> Controls verbosity.
      INTEGER, PRIVATE, SAVE :: verbose_ = XCLOC_PRINT_WARNINGS
      !> The precision of the module.
      INTEGER, PRIVATE, SAVE :: precision_ = XCLOC_SINGLE_PRECISION
      !> My rank on the communicator.
      INTEGER, PRIVATE, SAVE :: myid_ = MPI_UNDEFINED
      !> Sampling period (seconds) of signals.
      REAL(C_DOUBLE), PRIVATE, SAVE :: dt_ = 0.d0
      !> Flag indicating that communicator has to be destroyed.
      LOGICAL, PRIVATE, SAVE :: lfreeComm_ = .FALSE.
      !> Flag indicating that the module is initialized.
      LOGICAL, PRIVATE, SAVE :: linit_ = .FALSE.


      PUBLIC :: xclocMPI_initialize
      PUBLIC :: xclocMPI_finalize
      CONTAINS
!========================================================================================!
!                                     Begin the Code                                     !
!========================================================================================!
!>    @brief Initializes the MPI-based xcloc module.
!>    @param[in] comm          The MPI communicator.  
!>    @param[in] root          The root process ID on the MPI communicator.
!>                             This will likely be 0.
!>    @param[in] dsmGroupSize  asdf
!>    @param[in] npts          Number of points in the input signals.  This is defined
!>                             on the root process.
!>    @param[in] nptsPad       A tuning parameter.  If 2*npts - 1 is a semi-prime
!>                             number then the DFT will be inefficient.  nptsPad 
!>                             will post-pad the input signals with nptsPad - npts
!>                             zeros to mitigate this.  This is defined on the root
!>                             process.
!>    @param[in] nxcs          The number of cross-correlograms.  
!>    @param[in] verbose       Controls the verbosity.  This is defined on the root
!>                             process.
!>    @param[in] prec          Defines the precision as float or double.  This is
!>                             defined on the root process.
!>    @param[in] accuracy      Controls the accuracy of the vectorized MKL calculations.
!>    @param[out] ierr         0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_initialize(comm, root,                     &
                                     dsmGroupSize,                   &
                                     npts, nptsPad, nxcs,            &
                                     s2m, dt, ngrd,                  &
                                     nfcoeffs, ftype,                &
                                     xcPairs,                        &
                                     verbose, prec, accuracy, ierr ) &
      BIND(C, NAME='xclocMPI_initialize')
      IMPLICIT NONE
      TYPE(MPI_Comm), VALUE, INTENT(IN) :: comm
      INTEGER(C_INT), VALUE, INTENT(IN) :: root, npts, dsmGroupSize,            &
                                           nptsPad, nxcs, ngrd, nfcoeffs,       &
                                           ftype, s2m, verbose, prec, accuracy
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt
      INTEGER(C_INT), INTENT(IN) :: xcPairs(2*nxcs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER mpierr
      ierr = 0
      root_ = root
      CALL MPI_COMM_RANK(comm, myid_,   mpierr) 
      CALL MPI_COMM_SIZE(comm, nprocs_, mpierr)
      IF (myid_ == root_) THEN
         root_ = root
         IF (npts < 1 .OR. nptsPad < npts .OR. nxcs < 1 .OR. &
             ngrd < 1 .OR. dt <= 0.d0) THEN
            IF (npts < 1) WRITE(ERROR_UNIT,905) npts 
            IF (nptsPad < npts) WRITE(ERROR_UNIT,907) nptsPad, npts
            IF (nxcs < 1) WRITE(ERROR_UNIT,908) nxcs
            IF (ngrd < 1) WRITE(ERROR_UNIT,909) ngrd
            IF (dt <= 0.d0) WRITE(ERROR_UNIT,910) dt
            ierr = 1 
         ENDIF
         IF (.NOT.xcloc_constants_isValidSignalToMigrate(s2m)) THEN
            WRITE(ERROR_UNIT,911)
            ierr = 1
         ENDIF
         IF (.NOT.xcloc_constants_isValidPrecision(prec)) THEN
            WRITE(ERROR_UNIT,912)
            ierr = 1
         ENDIF
         IF (.NOT.xcloc_constants_isValidAccuracy(accuracy)) THEN
            WRITE(ERROR_UNIT,913)
            ierr = 1
         ENDIF
         nfcoeffs_ = nfcoeffs
         IF (.NOT.xcloc_constants_isValidFilteringType(ftype)) THEN
            WRITE(ERROR_UNIT,914)
            ierr = 1
         ELSE
            IF (nfcoeffs < 1) THEN
               WRITE(ERROR_UNIT,915)
               ierr = 1
            ELSE
               IF (MOD(nfcoeffs, 2) /= 1) THEN
                  nfcoeffs_ = nfcoeffs_ + 1
                  WRITE(OUTPUT_UNIT,916) nfcoeffs_ 
               ENDIF
            ENDIF
         ENDIF 
         IF (dsmGroupSize > nprocs_ .OR. MOD(nprocs_, dsmGroupSize) /= 0) THEN
            WRITE(ERROR_UNIT,917) dsmGroupSize, nprocs_
            ierr = 1
         ENDIF
         IF (ierr /= 0) GOTO 500
         ftype_ = ftype
         xcTypeToMigrate_ = s2m
         dt_ = dt
         accuracy_ = accuracy
         precision_ = prec
         verbose_ = verbose
         IF (ierr /= 0) GOTO 500
      ENDIF
  500 CONTINUE
      CALL MPI_BCAST(ierr, 1, MPI_INTEGER, root_, comm, mpierr)
      IF (ierr /= 0) RETURN
      ! Copy the communicator
      CALL MPI_COMM_DUP_WITH_INFO(comm, MPI_INFO_NULL, comm_, mpierr)
      lfreeComm_ = .TRUE.
      CALL MPI_BCAST(ftype_,           1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_BCAST(nfcoeffs_,        1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_BCAST(accuracy_,        1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_BCAST(precision_,       1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_BCAST(xcTypeToMigrate_, 1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_BCAST(verbose_,         1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_BCAST(dt_, 1, MPI_DOUBLE_PRECISION, root_, comm_, mpierr)

  905 FORMAT("xclocMPI_initialize: npts=", I8, " must be positive")
! 906 FORMAT("xclocMPI_initialize: nsignals=", I8, " must be at least 2")
  907 FORMAT("xclocMPI_initialize: nptsPad=", I8, " must be greater than npts=", I8)
  908 FORMAT("xclocMPI_initialize: No correlation pairs=", I8)
  909 FORMAT("xclocMPI_initialize: No grid points=", I8)
  910 FORMAT("xclocMPI_initialize: Sampling period=", E12.4, " must be positive")
  911 FORMAT("xclocMPI_initialize: Invalid type of xc signal to migrate")
  912 FORMAT("xclocMPI_initialize: Invalid precision")
  913 FORMAT("xclocMPI_initialize: Invalid accuracy")
  914 FORMAT("xclocMPI_initialize: Filtering type=", I4, " is invalid")
  915 FORMAT("xclocMPI_initialize: Number of filter taps=", I4, "must be positive")
  916 FORMAT("xclocMPI_initialize: Number of taps should be odd; setting to=", I4)
  917 FORMAT("xclocMPI_initialize: dsmGroupSize=", I4, &
             " cannot exceed and must equally divide nprocs_=", I4)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases memory on the module.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_finalize() &
      BIND(C, NAME='xclocMPI_finalize')
      INTEGER mpierr

      accuracy_ = XCLOC_HIGH_ACCURACY
      precision_ = XCLOC_SINGLE_PRECISION 
      verbose_ = XCLOC_PRINT_ERRORS 
      ftype_ = XCLOC_SPXC_DONOT_FILTER
      nfcoeffs_ = 0
      dt_ = 0.d0
      IF (lfreeComm_) CALL MPI_COMM_FREE(comm_, mpierr)
      lfreeComm_ = .FALSE.
      linit_ = .FALSE.
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the cross-correlogram signal type to migrate.
!>    @param[in] root             Root process ID that is broadcasting xcTypeToMigrate.
!>    @param[in] xcTypeToMigrate  Cross-correlogram signal type to migrate defined on
!>                                the root process.
!>    @param[in] xcTypeToMigrate  XCLOC_MIGRATE_PHASE_XCS indicates that phase
!>                                correlograms will be migrated.
!>    @param[in] xcTypeToMigrate  XCLOC_MIGRATE_XCS indicates that cross correlograms
!>                                will be migrated.
!>    @param[out] ierr            0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_setXCTypeToMigrate(root, xcTypeToMigrate, ierr) &
      BIND(C, NAME='xclocMPI_setXCTypeToMigrate')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: root, xcTypeToMigrate
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER mpierr
      LOGICAL lisValid
      ierr = 0
      xcTypeToMigrate_ = XCLOC_MIGRATE_PHASE_XCS
      IF (.NOT.linit_) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (myid_ == root) THEN
         lisValid = xcloc_constants_isValidSignalToMigrate(xcTypeToMigrate)
         IF (.NOT.lisValid) THEN
            WRITE(ERROR_UNIT,905)
            ierr = 1
         ELSE
            xcTypeToMigrate_ = xcTypeToMigrate
         ENDIF
      ENDIF
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root, comm_, mpierr)
      IF (ierr /= 0) RETURN 
      CALL MPI_Bcast(xcTypeToMigrate_, 1, MPI_INTEGER, root, comm_, mpierr)
  900 FORMAT("xclocMPI_setXCTypeToMigrate: Module not initialized")
  905 FORMAT('xclocMPI_setXCTypeToMigrate: Invalid type of xc signal to migrate')
      RETURN
      END

END MODULE
