!> @defgroup xcloc_mpi xclocMPI
!> @brief The parallel xcloc libary.
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
      USE XCLOC_UTILS_MPI
      USE XCLOC_DSMXC_MPI
      USE XCLOC_SPXC
      USE XCLOC_CONSTANTS
#if defined(__INTEL_COMPILER)
      !> @ingroup xcloc_mpi
      !> Global MPI communicator
      TYPE(MPI_Comm), PRIVATE, SAVE :: xclocGlobalComm_ = MPI_COMM_NULL !MPI_COMM_WORLD
      !> @ingroup xcloc_mpi
      !> Handle for communication between the processes for each XC/DSM group.
      TYPE(MPI_Comm), PRIVATE, SAVE :: xcdsmInterComm_ = MPI_COMM_NULL
      !> @ingroup xcloc_mpi
      !> Handle for communication within any XC/DSM group. 
      TYPE(MPI_Comm), PRIVATE, SAVE :: xcdsmIntraComm_ = MPI_COMM_NULL
#else
      TYPE(MPI_Comm), PRIVATE, SAVE :: xclocGlobalComm_
      TYPE(MPI_Comm), PRIVATE, SAVE :: xcdsmInterComm_
      TYPE(MPI_Comm), PRIVATE, SAVE :: xcdsmIntraComm_
#endif
      !> @ingroup xcloc_mpi
      !> Maps from the idsmxc'th group to the start index of l2gSignal_.  This has
      !> dimension [nxcdsmGroups_].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: l2gSignalPtr_(:)
      !> @ingroup xcloc_mpi
      !> Maps from the idsmxc'th group's local signal to the global signal index.
      !> This has dimension [l2gSignalPtr_(nxcdsmGroups_+1)-1].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: l2gSignal_(:)
      !> @ingroup xcloc_mpi
      !> Maps from the is'th signal to the start index of signal2Group_.  This has
      !> dimension [nsignalsTotal_+1].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: signalToGroupPtr_(:) 
      !> Maps from the a signal to the DSM/XC group to which this signal belongs.
      !> This has dimension [signalToGroupPtr_(nsignalsTotal_+1)-1].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: signalToGroup_(:)
      !> @ingroup xcloc_mpi
      !> Root ID on communicator.
      INTEGER, PRIVATE, SAVE :: root_ = 0
      !> @ingroup xcloc_mpi
      !> Number of processes on global communicator.  This is equal to:
      !> nxcdsmGroups_ x xcGroupSize_.
      INTEGER, PRIVATE, SAVE :: nprocs_ = 0
      !> @ingroup xcloc_mpi
      !> Defines a process' rank on the xcdsm-intra communicator.
      INTEGER, PRIVATE, SAVE :: xcdsmIntraCommID_ = 0
      !> @ingroup xcloc_mpi
      !> Defines a process' rank on the xcdsm-inter communicator.
      INTEGER, PRIVATE, SAVE :: xcdsmInterCommID_ = 0
      !> @ingroup xcloc_mpi
      !> Total number of cross-correlograms.
      INTEGER, PRIVATE, SAVE :: nXCsTotal_ = 0
      !> @ingroup xcloc_mpi
      !> Number of XC/DSM groups.
      INTEGER, PRIVATE, SAVE :: nxcdsmGroups_ = 0
      !> @ingroup xcloc_mpi
      !> Number of processes in each cross-correloation group.
      INTEGER, PRIVATE, SAVE :: xcGroupSize_ = 0
      !> @ingroup xcloc_mpi
      !> Determines the filter to be applied to the correlograms prior to migrating.
      INTEGER, PRIVATE, SAVE :: ftype_ = XCLOC_SPXC_DONOT_FILTER
      !> @ingroup xcloc_mpi
      !> Number of filtering taps.
      INTEGER, PRIVATE, SAVE :: nfcoeffs_ = 0
      !> @ingroup xcloc_mpi
      !> Total number of travel time tables.
      INTEGER, PRIVATE, SAVE :: nTablesTotal_ = 0
      !> @ingroup xcloc_mpi
      !> Total number of signals.
      INTEGER, PRIVATE, SAVE :: nsignalsTotal_ = 0

      !> @ingroup xcloc_mpi
      !> Accuracy of MKL computations.
      INTEGER, PRIVATE, SAVE :: accuracy_ = XCLOC_HIGH_ACCURACY
      !> @ingroup xcloc_mpi
      !> Signal to migrate.
      INTEGER, PRIVATE, SAVE :: xcTypeToMigrate_ = XCLOC_MIGRATE_PHASE_XCS
      !> @ingroup xcloc_mpi
      !> Controls verbosity.
      INTEGER, PRIVATE, SAVE :: verbose_ = XCLOC_PRINT_WARNINGS
      !> @ingorup xcloc_mpi
      !> The precision of the module.
      INTEGER, PRIVATE, SAVE :: precision_ = XCLOC_SINGLE_PRECISION
      !> @ingroup xcloc_mpi
      !> My rank on the communicator.
      INTEGER, PRIVATE, SAVE :: myid_ = MPI_UNDEFINED
      !> @ingroup xcloc_mpi
      !> Sampling period (seconds) of signals.
      REAL(C_DOUBLE), PRIVATE, SAVE :: dt_ = 0.d0
      !> @ingroup xcloc_mpi
      !> Flag indicating that communicator has to be destroyed.
      LOGICAL, PRIVATE, SAVE :: lfreeComm_ = .FALSE.
      !> @ingroup xcloc_mpi
      !> Flag indicating that the module is initialized.
      LOGICAL, PRIVATE, SAVE :: linit_ = .FALSE.

integer, private, allocatable, save :: nDSMXCsPerGroup_(:)

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
!>    @param[in] nxcdsmGroups  Number of XC/DSM groups.
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
                                     nxcdsmGroups,                   &
                                     npts, nptsPad, nxcs,            &
                                     s2m, dt, ngrd,                  &
                                     nfcoeffs, ftype,                &
                                     xcPairs,                        &
                                     verbose, prec, accuracy, ierr ) &
      BIND(C, NAME='xclocMPI_initialize')
      IMPLICIT NONE
      TYPE(MPI_Comm), VALUE, INTENT(IN) :: comm
      INTEGER(C_INT), VALUE, INTENT(IN) :: root, npts, nxcdsmGroups,            &
                                           nptsPad, nxcs, ngrd, nfcoeffs,       &
                                           ftype, s2m, verbose, prec, accuracy
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt
      INTEGER(C_INT), DIMENSION(2*nxcs), INTENT(IN) :: xcPairs
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER, ALLOCATABLE :: myDSMXCs(:), myDSMXCPtr(:), iwork(:), jwork(:)
      INTEGER i1, i2, ig, indx, is, mpierr, nsignalPairs, nsignals, nsg, nunique, nxcsLocal
      ierr = 0
      root_ = root
      CALL MPI_Comm_rank(comm, myid_,   mpierr) 
      CALL MPI_Comm_size(comm, nprocs_, mpierr)
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
         IF (MINVAL(xcPairs(1:2*nxcs)) < 1) THEN
            WRITE(ERROR_UNIT,904)
            ierr = 1
         ENDIF
         nsignals = MAXVAL(xcPairs(1:2*nxcs))
         IF (nsignals < 2) THEN
            WRITE(ERROR_UNIT,906) nsignals
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
         IF (nxcdsmGroups > nprocs_ .OR. MOD(nprocs_, nxcdsmGroups) /= 0) THEN
            WRITE(ERROR_UNIT,917) nxcdsmGroups, nprocs_
            ierr = 1
         ENDIF
         IF (ierr /= 0) GOTO 500
         ! Copy
         nsignalsTotal_ = nsignals
         nxcdsmGroups_ = nxcdsmGroups
         xcGroupSize_ = nprocs_/nxcdsmGroups_
         nXCsTotal_ = nxcs
         ftype_ = ftype
         xcTypeToMigrate_ = s2m
         dt_ = dt
         accuracy_ = accuracy
         precision_ = prec
         verbose_ = verbose
         ! Divide up cross-correlations and DSM's
         ALLOCATE(myDSMXCs(nXCsTotal_))
         ALLOCATE(myDSMXCPtr(nxcdsmGroups_)) 
         CALL xcloc_utils_partitionTasks(nXCsTotal_, nxcdsmGroups_, &
                                         myDSMXCPtr, myDSMXCs, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,950)
            GOTO 500
         ENDIF
         ALLOCATE(nDSMXCsPerGroup_(nxcdsmGroups_))
         DO ig=1,nxcdsmGroups_
            nDSMXCsPerGroup_(ig) = myDSMXCPtr(ig+1) - myDSMXCPtr(ig)
            IF (verbose_ > XCLOC_PRINT_WARNINGS) THEN
               WRITE(OUTPUT_UNIT,800) ig-1, nDSMXCsPerGroup_(ig)
            ENDIF
         ENDDO
         ! Make a map of the unique signals for each group
         ALLOCATE(l2gSignalPtr_(nxcdsmGroups_+1)) 
         ALLOCATE(jwork(nxcdsmGroups_*nsignals)); jwork(:) = 0
         l2gSignalPtr_(1) = 1
         DO ig=1,nxcdsmGroups_
            ! Get the unique signals in this group
            i1 = 2*(myDSMXCPtr(ig) - 1) + 1
            i2 = 2*(myDSMXCPtr(ig+1) - 1)
            nsignalPairs = i2 - i1 + 1
            CALL xcloc_utils_unique32s(nsignalPairs, xcPairs(i1:i2), nunique, iwork, ierr)
            IF (ierr /= 0) GOTO 500
            l2gSignalPtr_(ig+1) = l2gSignalPtr_(ig) + nunique
            i1 = l2gSignalPtr_(ig)
            i2 = l2gSignalPtr_(ig+1) - 1
            jwork(i1:i2) = iwork(1:nunique)
            IF (ALLOCATED(iwork)) DEALLOCATE(iwork)
         ENDDO
         nsignalPairs = l2gSignalPtr_(nxcdsmGroups_+1)-1 
         ALLOCATE(l2gSignal_(nsignalPairs))
         l2gSignal_(1:nsignalPairs) = jwork(1:nsignalPairs)
         ! Make a map of the signals that belong to each group 
         ALLOCATE(signalToGroupPtr_(nsignalsTotal_+1))
         jwork(:) = 0
         signalToGroupPtr_(1) = 1
         nsg = 0
         DO is=1,nsignals
            DO ig=1,nxcdsmGroups_
               i1 = l2gSignalPtr_(ig)
               i2 = l2gSignalPtr_(ig+1) - 1
               ! Search for element  
               CALL xcloc_utils_bsearch32i(i2-i1+1, is, l2gSignal_(i1:i2), &
                                           indx, ierr, .FALSE.)
               IF (ierr == 0) THEN
                  nsg = nsg + 1
                  jwork(nsg) = is
               ENDIF
            ENDDO
            signalToGroupPtr_(is+1) = nsg + 1
         ENDDO
         ALLOCATE(signalToGroup_(nsg))
         signalToGroup_(1:nsg) = jwork(1:nsg)
         DEALLOCATE(jwork)

         IF (ierr /= 0) GOTO 500
      ENDIF
  500 CONTINUE
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root_, comm, mpierr)
      IF (ierr /= 0) GOTO 5555
      ! Copy the input communicator
      CALL MPI_Comm_dup_with_info(comm, MPI_INFO_NULL, xclocGlobalComm_, mpierr)
      IF (mpierr /= MPI_SUCCESS) WRITE(ERROR_UNIT, 960)
      ! Share the group sizes
print *, nprocs_
      CALL MPI_Bcast(nxcdsmGroups_, 1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(xcGroupSize_,  1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      ! Create a communicator on which to do the parallel cross-correlations and DSMs
      CALL xcloc_utilsMPI_splitComm(xclocGlobalComm_, nxcdsmGroups_,      &
                                    xcdsmIntraCommID_, xcdsmIntraComm_, &
                                    xcdsmInterCommID_, xcdsmInterComm_, &
                                    ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,961)
         GOTO 5555
      ENDIF 
      lfreeComm_ = .TRUE.
      ! Split the communicator
      CALL MPI_Bcast(ftype_,           1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(nfcoeffs_,        1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(nsignalsTotal_,   1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(nXCsTotal_,       1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(accuracy_,        1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(precision_,       1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(xcTypeToMigrate_, 1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(verbose_,         1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(dt_, 1, MPI_DOUBLE_PRECISION, root_, xclocGlobalComm_, mpierr)

      IF (ALLOCATED(myDSMXCs))   DEALLOCATE(myDSMXCs)
      IF (ALLOCATED(myDSMXCPtr)) DEALLOCATE(myDSMXCPtr)
      ! INfo
  800 FORMAT("xclocMPI_initialize: Group", I4, " will process ", I6, " correlograms")
      ! Errors/warnings
  904 FORMAT("xclocMPI_initialize: Minimum signal index must be positive")
  905 FORMAT("xclocMPI_initialize: npts=", I8, " must be positive")
  906 FORMAT("xclocMPI_initialize: nsignals=", I8, " must be at least 2")
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
  917 FORMAT("xclocMPI_initialize: Number of XC/DSM groups=", I4,  &
             " cannot exceed and must divide equally nprocs_=", I4)
  950 FORMAT("xclocMPI_initialize: Failed to partition work")
  960 FORMAT("xclocMPI_initialize: Failed to duplicate communicator")
  961 FORMAT("xclocMPI_initialize: Error splitting communicators")
 5555 IF (ierr /= 0) CALL xclocMPI_finalize()
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

      nxcdsmGroups_ = 0
      accuracy_ = XCLOC_HIGH_ACCURACY
      precision_ = XCLOC_SINGLE_PRECISION 
      verbose_ = XCLOC_PRINT_ERRORS 
      ftype_ = XCLOC_SPXC_DONOT_FILTER
      nfcoeffs_ = 0
      dt_ = 0.d0
      IF (lfreeComm_) THEN
         CALL MPI_Comm_free(xcdsmIntraComm_,  mpierr)
         CALL MPI_Comm_free(xcdsmInterComm_,  mpierr)
         CALL MPI_Comm_free(xclocGlobalComm_, mpierr)
      ENDIF
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
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN 
      CALL MPI_Bcast(xcTypeToMigrate_, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
  900 FORMAT("xclocMPI_setXCTypeToMigrate: Module not initialized")
  905 FORMAT('xclocMPI_setXCTypeToMigrate: Invalid type of xc signal to migrate')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the ID of the root process on the global xcloc communicator.
!>    @param[out] root    Root process ID.
!>    @param[out] ierr    0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_getGlobalRootProcessID(root, ierr) &
      BIND(C, NAME='xclocMPI_getGlobalRootProcessID')
      INTEGER(C_INT), INTENT(OUT) :: root, ierr
      ierr = 0
      root = 0
      IF (.NOT.linit_) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      root_ = root
  900 FORMAT("xclocMPI_getGlobalRootProcessID: Module not initialized")
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @param[in] tableNumberIn  Global table index.  This is defined on the root process.
!>    @param[in] ngrdIn         Number of grid points in table.  This is defined on the
!>                              the root process.
!>    @param[in] root           The root process ID on the xcloc global communicator.
!>    @param[in] table          The
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_setTable64f(tableNumberIn, ngrdIn, root, table, ierr) &
      BIND(C, NAME='xclocMPI_setTable64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: tableNumberIn, ngrdIn, root
      REAL(C_DOUBLE), DIMENSION(ngrdIn), INTENT(IN) :: table
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: work(:)
      INTEGER commRoot, ig, mpierr, ngrd, tableNumber
      TYPE(MPI_Status) :: stat
      LOGICAL lneedTable
      ierr = 0
      lneedTable = .FALSE.
      IF (.NOT.linit_) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (myid_ == root) THEN
         CALL xcloc_dsmxcMPI_getNumberOGridPointsInTable(ngrd)
         IF (ngrd /= ngrdIn) THEN
            WRITE(ERROR_UNIT,905) ngrdIn, ngrd          
            ierr = 1
         ENDIF
         IF (tableNumberIn < nTablesTotal_ .OR. tableNumberIn > ntablesTotal_) THEN
            WRITE(ERROR_UNIT,910) tableNumberIn, nTablesTotal_
            ierr = 1
         ENDIF
         commRoot = xcdsmIntraCommID_
         tableNumber = tableNumberIN
      ENDIF
      CALL MPI_Bcast(ierr,        1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN
      CALL MPI_Bcast(ngrd,        1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(tableNumber, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(commRoot,    1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
! TODO do i need this table?

      ! Distribute the table to the roots
      IF (xcdsmIntraCommID_ /= commRoot) THEN
         DO ig=0,nxcdsmGroups_-1
! TODO do i need to send this table?
            IF (myid_ == root) THEN
               CALL MPI_Send(table, ngrd, MPI_DOUBLE_PRECISION, ig, 0, &
                             xcdsmInterComm_, mpierr) 
            ELSE
               IF (lneedTable) THEN
                  ALLOCATE(work(ngrd))
                  CALL MPI_Recv(work, ngrd, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                                MPI_ANY_TAG, xcdsmInterComm_, stat, mpierr)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      IF (ALLOCATED(work)) DEALLOCATE(work)
  900 FORMAT("xclocMPI_setTable64f: Module not initialized")
  905 FORMAT("xclocMPI_setTable64f: ngrdIn=", I6, "should be ngrd=", I6)
  910 FORMAT("xclocMPI_setTable64f: table number=", I6, "should be in range[1,", I6, "]")
      RETURN
      END

END MODULE
