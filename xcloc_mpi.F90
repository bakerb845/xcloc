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
      USE XCLOC_FDXC_MPI
      USE XCLOC_DSMXC_MPI
      USE XCLOC_SPXC
      USE XCLOC_CONSTANTS
      IMPLICIT NONE
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
      !> RMA window for setting signals.
      TYPE(MPI_Win), PRIVATE, SAVE :: signalRMAWindow_
      !> @ingroup xcloc_mpi
DOUBLE PRECISION, POINTER, DIMENSION(:), PRIVATE, SAVE :: signalsRMA64f_
      !> @ingroup xcloc_mpi
REAL, POINTER, DIMENSION(:), PRIVATE, SAVE :: signalsRMA32f_ 
      !> @ingroup xcloc_mpi
TYPE(C_PTR), SAVE :: dataPtr_ = C_NULL_PTR
      !> @ingroup xcloc_mpi
      !> Maps from the XC/DSM's starting block of correlograms.  This has dimension
      !> [nxcdsmGroups_ + 1].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: myDSMXCPtr_(:)
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
      !> @ingroup xcloc_mpi
      !> Maps from the a signal to the DSM/XC group to which this signal belongs.
      !> This has dimension [signalToGroupPtr_(nsignalsTotal_+1)-1].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: signalToGroup_(:)
      !> @ingroup xcloc_mpi
      !> Tabulates the unique signals.  This has dimension [nSignalsTotal_].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: uniqueSignals_(:)
      !> @ingroup xcloc_mpi
      !> A local cross-correlation pairs table for a an XC/DSM group.
      !> This has dimension [2 x nxcsLocal_]. 
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: xcPairsLocal_(:)
      !> @ingroup xcloc_mpi
      !> Determines the ownernship ranges for each grid.  This has dimension
      !> [xcGroupSize_+1].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: grdPointsPtr_(:)
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
      !> Sampling period (seconds) of signals.
      REAL(C_DOUBLE), PRIVATE, SAVE :: dt_ = 0.d0

      !> @ingroup xcloc_mpi
      !> Number of filtering taps.
      INTEGER, PRIVATE, SAVE :: nfcoeffs_ = 0
      !> @ingroup xcloc_mpi
      !> Total number of grid points.
      INTEGER, PRIVATE, SAVE :: ngrdTotal_ = 0
      !> @ingroup xcloc_mpi
      !> Number of points in input signals.
      INTEGER, PRIVATE, SAVE :: npts_ = 0
      !> @ingroup xcloc_mpi
      !> Padded number of points in input signals.  The correlogram length will be
      !> 2*nptsPad_ - 1.  This can be used to avoid correlation lengths whose DFT is
      !> inefficient because 2*npts - 1 would be a large semi-prime number.
      INTEGER, PRIVATE, SAVE :: nptsPad_ = 0
      !> @ingroup xcloc_mpi
      !> Number of points in each correlogram.  This is equal to 2*nptsPad_ - 1.
      INTEGER, PRIVATE, SAVE :: nptsInXCs_ = 0
      !> @ingroup xcloc_mpi
      !> Number of signals pertinent to the xcdsm intra-communicator.
      INTEGER, PRIVATE, SAVE :: nsignalsLocal_ = 0
      !> @ingroup xcloc_mpi
      !> Total number of signals.  This should equal the total number of travel 
      !> time tables.
      INTEGER, PRIVATE, SAVE :: nsignalsTotal_ = 0
      !> @ingroup xcloc_mpi
      !> Total number of travel time tables.  This should equal the total number of
      !> signals.
      INTEGER, PRIVATE, SAVE :: nTablesTotal_ = 0
      !> @ingroup xcloc_mpi
      !> Number of correlograms that process is responsible for.
      INTEGER, PRIVATE, SAVE :: nxcsLocal_ = 0
      !> @ingroup xcloc_mpi
      !> Total number of cross-correlograms.
      INTEGER, PRIVATE, SAVE :: nxcsTotal_ = 0
      !> @ingroup xcloc_mpi
      !> Number of XC/DSM groups.
      INTEGER, PRIVATE, SAVE :: nxcdsmGroups_ = 0
      !> @ingroup xcloc_mpi
      !> Number of processes in each cross-correlation group.
      INTEGER, PRIVATE, SAVE :: xcGroupSize_ = 0

      !> @ingroup xcloc_mpi
      !> Determines the filter to be applied to the correlograms prior to migrating.
      INTEGER, PRIVATE, SAVE :: ftype_ = XCLOC_SPXC_DONOT_FILTER
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
      !> Flag indicating that all the signals are set.
      LOGICAL, PRIVATE, SAVE :: lhaveSignals_ = .FALSE.
      !> @ingroup xcloc_mpi
      !> @Flag indicating that cross-correlograms have been computed.
      LOGICAL, PRIVATE, SAVE :: lhaveXCs_ = .FALSE.
      !> @ingroup xcloc_mpi
      !> Flag indicating that the migration image is computed.
      LOGICAL, PRIVATE, SAVE :: lhaveImage_ = .FALSE.
      !> @ingroup xcloc_mpi
      !> Flag indicating that communicator has to be destroyed.
      LOGICAL, PRIVATE, SAVE :: lfreeComm_ = .FALSE.
      !> @ingroup xcloc_mpi
      !> Flag indicating that all the travel time tables have been set.
      LOGICAL, PRIVATE, SAVE :: lhaveAllTables_ = .FALSE.
      !> @ingroup xcloc_mpi
      !> Flag indiating that the RMA signal window has to be destroyed.
      LOGICAL, PRIVATE, SAVE :: lfreeWindow_ = .FALSE.
      !> @ingroup xcloc_mpi
      !> Flag indicating that the module is initialized.
      LOGICAL, PRIVATE, SAVE :: linit_ = .FALSE.

      !> @ingroup xcloc_mpi
      !> Make all the submodules use 0 as the root process ID.
      INTEGER, PRIVATE, PARAMETER :: subModuleRoot_ = 0

integer, private, allocatable, save :: nDSMXCsPerGroup_(:)

      PUBLIC :: xclocMPI_initialize
      PUBLIC :: xclocMPI_getCorrelogramLength
      PUBLIC :: xclocMPI_getNumberOfCorrelograms
      PUBLIC :: xclocMPI_setXCTypeToMigrate
      PUBLIC :: xclocMPI_setTable64f
      PUBLIC :: xclocMPI_setSignals64f
      PUBLIC :: xclocMPI_setSignals32f
      PUBLIC :: xclocMPI_finalize
      CONTAINS
!========================================================================================!
!                                     Begin the Code                                     !
!========================================================================================!
!>    @brief Initializes the MPI-based xcloc module.
!>    @param[in] fcomm         The MPI communicator.  
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
      SUBROUTINE xclocMPI_initialize(fcomm, root,                    &
                                     nxcdsmGroups,                   &
                                     npts, nptsPad, nxcs,            &
                                     s2m, dt, ngrd,                  &
                                     nfcoeffs, ftype,                &
                                     xcPairs,                        &
                                     verbose, prec, accuracy, ierr ) &
      BIND(C, NAME='xclocMPI_initialize')
      IMPLICIT NONE
      !TYPE(MPI_Comm), VALUE, INTENT(IN) :: fcomm
      INTEGER(C_INT64_T), VALUE, INTENT(IN) :: fcomm
      INTEGER(C_INT), VALUE, INTENT(IN) :: root, npts, nxcdsmGroups,            &
                                           nptsPad, nxcs, ngrd, nfcoeffs,       &
                                           ftype, s2m, verbose, prec, accuracy
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt
      INTEGER(C_INT), DIMENSION(2*nxcs), INTENT(IN) :: xcPairs
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT64_T) commTemp
      INTEGER, ALLOCATABLE :: myDSMXCs(:), iwork(:), jwork(:), xcPairsWork(:)
      INTEGER i1, i2, ig, ierrLoc, indx, is, ixc, mpierr, nalloc, nsignalPairs, &
              nsignals, nsg, nunique
      TYPE(MPI_Comm) comm
      INTEGER(KIND=MPI_ADDRESS_KIND) nallocMPI
      INTEGER arrayShape(1), sizeDouble
      comm%MPI_VAL = INT(fcomm, KIND(comm%MPI_VAL))
      ierr = 0
      CALL xclocMPI_finalize()
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
         !nsignals = MAXVAL(xcPairs(1:2*nxcs))
         CALL xcloc_utils_unique32s(2*nxcs, xcPairs,             &
                                   nsignals, uniqueSignals_, ierr)
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
         nTablesTotal_ = nsignals
         nxcdsmGroups_ = nxcdsmGroups
         xcGroupSize_ = nprocs_/nxcdsmGroups_
         ngrdTotal_ = ngrd
         npts_ = npts
         nptsPad_ = nptsPad
         nxcsTotal_ = nxcs
         ftype_ = ftype
         xcTypeToMigrate_ = s2m
         dt_ = dt
         accuracy_ = accuracy
         precision_ = prec
         verbose_ = verbose
         ! Divide up cross-correlations and DSM's
         ALLOCATE(myDSMXCs(nxcsTotal_))
         ALLOCATE(myDSMXCPtr_(nxcdsmGroups_+1)) 
         CALL xcloc_utils_partitionTasks(nxcsTotal_, nxcdsmGroups_, &
                                         myDSMXCPtr_, myDSMXCs, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,950)
            GOTO 500
         ENDIF
         ALLOCATE(nDSMXCsPerGroup_(nxcdsmGroups_))
         DO ig=1,nxcdsmGroups_
            nDSMXCsPerGroup_(ig) = myDSMXCPtr_(ig+1) - myDSMXCPtr_(ig)
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
            i1 = 2*(myDSMXCPtr_(ig) - 1) + 1
            i2 = 2*(myDSMXCPtr_(ig+1) - 1)
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
                  jwork(nsg) = ig - 1
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
      CALL MPI_Bcast(ngrdTotal_,       1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(npts_,            1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(nptsPad_,         1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(nsignalsTotal_,   1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(nTablesTotal_,    1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(nxcdsmGroups_,    1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(nxcsTotal_,       1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(accuracy_,        1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(precision_,       1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(xcTypeToMigrate_, 1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(verbose_,         1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(dt_, 1, MPI_DOUBLE_PRECISION, root_, xclocGlobalComm_, mpierr)
      ! Send vectors
      nalloc = 0
      ALLOCATE(xcPairsWork(2*nxcsTotal_)); xcPairsWork(:) = 0
      IF (myid_ == root_) xcPairsWork(1:2*nxcsTotal_) = xcPairs(1:2*nxcsTotal_)
      CALL MPI_Bcast(xcPairsWork, 2*nxcsTotal_, MPI_INTEGER, root_, &
                     xclocGlobalComm_, mpierr)
      IF (.NOT.ALLOCATED(myDSMXCPtr_)) ALLOCATE(myDSMXCPtr_(nxcdsmGroups_+1))
      CALL MPI_Bcast(myDSMXCPtr_, nxcdsmGroups_+1, MPI_INTEGER, root_, &
                     xclocGlobalComm_, mpierr) 
      IF (.NOT.ALLOCATED(l2gSignalPtr_)) ALLOCATE(l2gSignalPtr_(nxcdsmGroups_+1))
      CALL MPI_Bcast(l2gSignalPtr_, nxcdsmGroups_+1, MPI_INTEGER, root_, &
                     xclocGlobalComm_, mpierr)
      nalloc = l2gSignalPtr_(nxcdsmGroups_+1) - 1
      IF (.NOT.ALLOCATED(l2gSignal_)) ALLOCATE(l2gSignal_(nalloc))
      CALL MPI_Bcast(l2gSignal_, nalloc, MPI_INTEGER, root_, xclocglobalComm_, mpierr)
      IF (.NOT.ALLOCATED(signalToGroupPtr_)) ALLOCATE(signalToGroupPtr_(nSignalsTotal_+1))
      CALL MPI_Bcast(signalToGroupPtr_, nSignalsTotal_+1, MPI_INTEGER, root_, &
                     xclocGlobalComm_, mpierr)
      nalloc = signalToGroupPtr_(nSignalsTotal_+1) - 1
      IF (.NOT.ALLOCATED(signalToGroup_)) ALLOCATE(signalToGroup_(nalloc))
      CALL MPI_Bcast(signalToGroup_, nalloc, MPI_INTEGER, root_, &
                     xclocGlobalComm_, mpierr) 
      IF (.NOT.ALLOCATED(nDSMXCsPerGroup_)) ALLOCATE(nDSMXCsPerGroup_(nxcdsmGroups_))
      CALL MPI_Bcast(nDSMXCsPerGroup_, nxcdsmGroups_, MPI_INTEGER, root_, &
                     xclocGlobalComm_, mpierr)
      IF (.NOT.ALLOCATED(uniqueSignals_)) ALLOCATE(uniqueSignals_(nSignalsTotal_))
      CALL MPI_Bcast(uniqueSignals_, nSignalsTotal_, MPI_INTEGER, root_,  &
                     xclocGlobalComm_, mpierr)
      ! Compute some local information from pointers
      nsignalsLocal_ = l2gSignalPtr_(xcdsmInterCommID_ + 2) &
                     - l2gSignalPtr_(xcdsmInterCommID_ + 1)
      nxcsLocal_ = myDSMXCPtr_(xcdsmInterCommID_ + 2) - myDSMXCPtr_(xcdsmInterCommID_ + 1)
      ierrLoc = 0
      ! Make xcPairsWork local
      i1 = 2*(myDSMXCPtr_(xcdsmInterCommID_ + 1) - 1) + 1
      i2 = 2*(myDSMXCPtr_(xcdsmInterCommID_ + 2) - 1)
      CALL xcloc_utils_unique32s(2*nxcsLocal_, xcPairsWork(i1:i2), nunique, &
                                 iwork, ierr)
      ALLOCATE(xcPairsLocal_(MAX(1, 2*nxcsLocal_))); xcPairsLocal_(:) = 0
      DO ixc=1,2*nxcsLocal_
         CALL xcloc_utils_bsearch32i(nunique, xcPairsWork(i1-1+ixc), iwork, &
                                     xcPairsLocal_(ixc), ierr, .TRUE.)
      ENDDO
      IF (ALLOCATED(iwork)) DEALLOCATE(iwork)
      ! Initialize the cross-correlation engine
      ierr = 0
      DO ig=0,nxcdsmGroups_-1
         IF (ig == xcdsmInterCommID_ .AND. nxcsLocal_ > 0) THEN
            IF (xcdsmIntraCommID_ == subModuleRoot_ .AND. &
                verbose_ >= XCLOC_PRINT_INFO) THEN
               WRITE(OUTPUT_UNIT,810) ig
            ENDIF
            commTemp = INT(xcdsmIntraComm_%MPI_VAL, KIND(commTemp))
            CALL xcloc_fdxcMPI_initialize(commTemp, & !xcdsmIntraComm_.MPI_VAL, &
                                          subModuleRoot_,  &
                                          npts_, nptsPad_,                     &
                                          nxcsLocal_, xcPairsLocal_,           &
                                          verbose_, precision_, accuracy_, ierr)
            IF (ierr /= 0) THEN
               WRITE(ERROR_UNIT,990) myid_
               ierrLoc = 1
            ENDIF
            CALL xcloc_fdxcMPI_getCorrelogramLength(nptsInXCs_, ierr)
         ENDIF
         CALL MPI_Barrier(xclocGlobalComm_, mpierr)
      ENDDO
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) GOTO 5555 
      ! Initialize the signals processing engine 
      CALL xcloc_spxc_initialize(nfcoeffs_, ftype_, ierr) 
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,991) myid_
         ierrLoc = 1
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) GOTO 5555
      ! Initialize the DSM engine
      DO ig=0,nxcdsmGroups_-1
         IF (ig == xcdsmInterCommID_ .AND. nxcsLocal_ > 0) THEN
            IF (xcdsmIntraCommID_ == subModuleRoot_ .AND. &
                verbose_ >= XCLOC_PRINT_INFO) THEN
               WRITE(OUTPUT_UNIT,815) ig
            ENDIF
            commTemp = INT(xcdsmIntraComm_%MPI_VAL, KIND(commTemp))
            CALL xcloc_dsmxcMPI_initialize(commTemp, subModuleRoot_,            &
                                           ngrdTotal_, nxcsLocal_, nptsInXCs_,  &
                                           dt_, xcPairsLocal_, verbose_, ierr)
            IF (ierr /= 0) THEN
               WRITE(ERROR_UNIT,992) myid_
               ierrLoc = 1
            ENDIF
            CALL xcloc_dsmxcMPI_getGridOwnershipRange(grdPointsPtr_)
         ENDIF
         CALL MPI_Barrier(xclocGlobalComm_, mpierr)
      ENDDO
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) GOTO 5555 
      ! Create an RMA window to set the signals.  Double > float so this is plenty.
      CALL MPI_Sizeof(1.d0, sizeDouble, mpierr)
      arrayShape(1) = 1
      nallocMPI = 1
      IF (myid_ == root_) THEN
         arrayShape(1) = npts_*nSignalsTotal_
         nallocMPI = INT(arrayShape(1), KIND(nallocMPI)) &
                    *INT(sizeDouble,    KIND(nallocMPI))
      ENDIF
      CALL MPI_Alloc_mem(nallocMPI, MPI_INFO_NULL, dataPtr_, mpierr)
      CALL C_F_POINTER(dataPtr_, signalsRMA32f_, arrayShape)
      CALL C_F_POINTER(dataPtr_, signalsRMA64f_, arrayShape)
      CALL MPI_Win_create(signalsRMA64f_, nallocMPI, sizeDouble, MPI_INFO_NULL, &
                          xclocGlobalComm_, signalRMAWindow_, mpierr)
      lfreeWindow_ = .TRUE.
      linit_ = .TRUE.
      ! Finish and clean up
 5555 CONTINUE
      IF (ALLOCATED(myDSMXCs)) DEALLOCATE(myDSMXCs)
      IF (ALLOCATED(xcPairsWork)) DEALLOCATE(xcPairsWork)
      IF (ierr /= 0) THEN
         IF (myid_ == root_) WRITE(ERROR_UNIT,995)
         CALL xclocMPI_finalize()
      ENDIF
      ! Info
  800 FORMAT("xclocMPI_initialize: Group ", I0, " will process ", I0, " correlograms")
  810 FORMAT("xclocMPI_initialize: Initializing fdxcMPI for XC/DSM group ", I0)
  815 FORMAT("xclocMPI_initialize: Initializing dsmxcMPI for XC/DSM group ", I0)
      ! Errors/warnings
  904 FORMAT("xclocMPI_initialize: Minimum signal index must be positive")
  905 FORMAT("xclocMPI_initialize: npts=", I0, " must be positive")
  906 FORMAT("xclocMPI_initialize: nsignals=", I0, " must be at least 2")
  907 FORMAT("xclocMPI_initialize: nptsPad=", I0, " must be greater than npts=", I0)
  908 FORMAT("xclocMPI_initialize: No correlation pairs=", I0)
  909 FORMAT("xclocMPI_initialize: No grid points=", I0)
  910 FORMAT("xclocMPI_initialize: Sampling period=", E12.4, " must be positive")
  911 FORMAT("xclocMPI_initialize: Invalid type of xc signal to migrate")
  912 FORMAT("xclocMPI_initialize: Invalid precision")
  913 FORMAT("xclocMPI_initialize: Invalid accuracy")
  914 FORMAT("xclocMPI_initialize: Filtering type=", I0, " is invalid")
  915 FORMAT("xclocMPI_initialize: Number of filter taps=", I0, " must be positive")
  916 FORMAT("xclocMPI_initialize: Number of taps should be odd; setting to=", I0)
  917 FORMAT("xclocMPI_initialize: Number of XC/DSM groups=", I0,  &
             " cannot exceed and must divide equally nprocs_=", I0)
  950 FORMAT("xclocMPI_initialize: Failed to partition work")
  960 FORMAT("xclocMPI_initialize: Failed to duplicate communicator")
  961 FORMAT("xclocMPI_initialize: Error splitting communicators")
  990 FORMAT("xclocMPI_initialize: Failed to initialize fdxcMPI on process ", I0)
  991 FORMAT("xclocMPI_initialize: Failed to initialize spxc on process ", I0)
  992 FORMAT("xclocMPI_initialize: Failed to initialize dsmxcMPI on process ", I0)
  995 FORMAT("xclocMPI_initialize: Errors detected during initialization")
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

      CALL xcloc_fdxcMPI_finalize()
      CALL xcloc_spxc_finalize()
      CALL xcloc_dsmxcMPI_finalize()

      ngrdTotal_ = 0
      npts_ = 0
      nptsPad_ = 0
      nptsInXCs_ = 0
      nxcsLocal_ = 0
      nxcsTotal_ = 0
      nxcdsmGroups_ = 0
      xcGroupSize_ = 0
      ftype_ = XCLOC_SPXC_DONOT_FILTER
      accuracy_ = XCLOC_HIGH_ACCURACY
      precision_ = XCLOC_SINGLE_PRECISION 
      verbose_ = XCLOC_PRINT_ERRORS 
      nfcoeffs_ = 0
      nTablesTotal_ = 0
      nsignalsTotal_ = 0
      dt_ = 0.d0
      IF (lfreeComm_) THEN
         CALL MPI_Comm_free(xcdsmIntraComm_,  mpierr)
         CALL MPI_Comm_free(xcdsmInterComm_,  mpierr)
         CALL MPI_Comm_free(xclocGlobalComm_, mpierr)
      ENDIF
      IF (lfreeWindow_) CALL MPI_Win_free(signalRMAWindow_, mpierr)
      IF (ALLOCATED(myDSMXCPtr_))       DEALLOCATE(myDSMXCPtr_)
      IF (ALLOCATED(l2gSignalPtr_))     DEALLOCATE(l2gSignalPtr_)
      IF (ALLOCATED(l2gSignal_))        DEALLOCATE(l2gSignal_)
      IF (ALLOCATED(signalToGroupPtr_)) DEALLOCATE(signalToGroupPtr_)
      IF (ALLOCATED(signalToGroup_))    DEALLOCATE(signalToGroup_)
      IF (ALLOCATED(uniqueSignals_))    DEALLOCATE(uniqueSignals_)
      IF (ALLOCATED(xcPairsLocal_))     DEALLOCATE(xcPairsLocal_)
      IF (ALLOCATED(grdPointsPtr_))     DEALLOCATE(grdPointsPtr_)
!     IF (ASSOCIATED(signalsRMA32f_))   CALL MPI_Free_mem(signalsRMA32f_, mpierr)
      IF (ASSOCIATED(signalsRMA64f_))   CALL MPI_Free_mem(signalsRMA64f_, mpierr)
if (allocated(nDSMXCsPerGroup_)) deallocate(nDSMXCsPerGroup_)
      root_ = 0
      nprocs_ = 0
      xcdsmIntraCommID_ = 0
      xcdsmInterCommID_ = 0
      lfreeComm_ = .FALSE.
      lfreeWindow_ = .FALSE.
      lhaveAllTables_ = .FALSE.
      lhaveSignals_ = .FALSE.
      lhaveXCs_ = .FALSE.
      lhaveImage_ = .FALSE.
      linit_ = .FALSE.
      dataPtr_ = C_NULL_PTR
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
!>    @brief Returns the number of correlograms to be computed.
!>    @param[out] nxcs   The number of correlograms.
!>    @param[out] ierr   0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_getNumberOfCorrelograms(nxcs, ierr) &
      BIND(C, NAME='xclocMPI_getNumberOfCorrelograms')
      INTEGER(C_INT), INTENT(OUT) :: nxcs, ierr
      INTEGER ierrLoc, mpierr, nxcsLoc
      nxcs = 0
      ierr = 0
      ierrLoc = 0
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      CALL xcloc_fdxcMPI_getNumberOfCorrelograms(nxcsLoc, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierrLoc = 1
      ENDIF
      IF (xcdsmIntraCommID_ == root_) THEN
         CALL MPI_Reduce(nxcsLoc, nxcs, 1, MPI_INTEGER, MPI_SUM, root_, &
                         xcdsmInterComm_, mpierr)
         IF (myid_ == root_ .AND. nxcs /= nxcsTotal_) THEN
            WRITE(ERROR_UNIT,901)
            ierrLoc = 1
            nxcs = 0
         ENDIF
      ENDIF
      CALL MPI_Bcast(nxcs, 1, MPI_INTEGER, root_, xclocGlobalComm_, mpierr)
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
  900 FORMAT('xclocMPI_getNumberOfCorrelograms: Failed to get number of local xcs')
  901 FORMAT('xclocMPI_getNumberOfCorrelograms: Internal error')
      RETURN
      END
!>    @brief Returns the number of points in the time domain correlations.
!>    @param[out] nptsInXCs  Number of points in the correlograms.
!>    @param[out] ierr       0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_getCorrelogramLength(nptsInXCs, ierr) &
      BIND(C, NAME='xclocMPI_getCorrelogramLength')
      INTEGER(C_INT), INTENT(OUT) :: nptsInXCs, ierr
      INTEGER ierrLoc, mpierr
      nptsInXCs = 0
      ierr = 0
      ierrLoc = 0
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here 
      CALL xcloc_fdxcMPI_getCorrelogramLength(nptsInXCs, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900)
         ierrLoc = 1
      ENDIF
      IF (nptsInXCs /= nptsInXCs_) THEN
         WRITE(ERROR_UNIT,901)
         ierrLoc = 1 
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
  900 FORMAT('xclocMPI_getCorrelogramLength: Failed to get correlogram length')
  901 FORMAT('xclocMPI_getCorrelogramLength: Internal error')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gathers all the correlograms onto the root process on the global 
!>           xcloc communicator.
!>    @param[in] ldxcIn  The leading dimension of xcs.  This is defined on the root
!>                       process and must be at least nptsInXCs_.
!>    @param[in] nxcs    The number of correlograms.  This is defined on the root
!>                       process and should equal nxcsTotal_.
!>    @param[in] root    The root process ID on the xcloc global communicator.
!>                       This must be defined on all processes.
!>    @param[out] xcs    The [ldxcIn x nxcs] matrix of correlograms stored in column 
!>                       major format.  This is only accessed on the root process.
!>    @param[out] ierr   0 indicates success.
!>    @result 0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_getCorrelograms64f(ldxcIn, nxcs, root, xcs, ierr) &
      BIND(C, NAME='xclocMPI_getCorrelograms64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldxcIn, nxcs, root
      REAL(C_DOUBLE), DIMENSION(ldxcIn*nxcs), INTENT(OUT) :: xcs 
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: xcsLocal(:)
      INTEGER, DIMENSION(:), ALLOCATABLE :: displs, recvCount
      INTEGER ierrLoc, interCommTarget, intraCommTarget, irecv, ldxc, mpierr, &
              nxcsInGroup, sendCount
      ierr = 0 
      ierrLoc = 0
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      ! Did I compute the XCs yet?
      IF (.NOT.lhaveXCs_) THEN
         IF (myid_ == root_) WRITE(ERROR_UNIT,900)
         ierrLoc = 1
      ENDIF
      ! Verify the input arguments and set dest info 
      IF (myid_ == root) THEN
         IF (ldxcIn < nptsInXCs_ .AND. nxcs /= nxcsTotal_) THEN
            IF (ldxcIn < nptsInXCs_) WRITE(ERROR_UNIT,901) nptsInXCs_
            IF (nxcs /= nxcsTotal_) WRITE(ERROR_UNIT,902) nxcsTotal_
            ierrLoc = 1
         ENDIF
         ! Set the destination information and initialize the result to nothing 
         xcs(:) = 0.d0
         interCommTarget = xcdsmInterCommID_
         intraCommTarget = xcdsmIntraCommID_
         ldxc = ldxcIn
      ENDIF
      ! Leave now if there is an error
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN
      ! Send around target information
      CALL MPI_Bcast(interCommTarget, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(intraCommTarget, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(ldxc,            1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      ! Gather the correlograms into this group's `root' process
      IF (xcdsmIntraCommID_ == intraCommTarget) ALLOCATE(xcsLocal(MAX(1,ldxc*nxcsLocal_)))
      CALL xcloc_fdxcMPI_gatherCorrelograms64f(ldxc, nxcsLocal_, intraCommTarget, &
                                               xcsLocal, ierrLoc)
      IF (ierrLoc /= 0) THEN
         IF (xcdsmIntraCommID_ == intraCommTarget) WRITE(ERROR_UNIT,903) xcdsmInterCommID_
         ierrLoc = 1
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN
      ! Gather the correlograms onto the root. 
      IF (xcdsmIntraCommID_ == interCommTarget) THEN
         sendCount = ldxc*nxcsLocal_
         IF (xcdsmIntraCommID_ == intraCommTarget) THEN
            ALLOCATE(recvCount(nxcdsmGroups_)); recvCount(:) = 0
            ALLOCATE(displs(nxcdsmGroups_)); displs(:) = 0
            displs(1) = 0 ! This is C indexed
            DO irecv=1,nxcdsmGroups_
               nxcsInGroup = myDSMXCPtr_(irecv+1) - myDSMXCPtr_(irecv)
               recvCount(irecv) = nxcsInGroup*nptsInXCs_
               IF (irecv < nxcdsmGroups_) THEN
                  displs(irecv+1) = displs(irecv) + recvCount(irecv)
               ENDIF
            ENDDO
         ENDIF
         CALL MPI_Gatherv(xcsLocal, sendCount, MPI_DOUBLE_PRECISION,                    &
                          xcs, recvCount, displs,                                       &
                          MPI_DOUBLE_PRECISION, interCommTarget, xcdsmInterComm_, mpierr) 
         IF (ALLOCATED(displs))    DEALLOCATE(displs)
         IF (ALLOCATED(recvCount)) DEALLOCATE(recvCount) 
      ENDIF
      IF (ALLOCATED(xcsLocal)) DEALLOCATE(xcsLocal)
  900 FORMAT('xclocMPI_getCorrelograms64f: xcs not yet computed')
  901 FORMAT('xclocMPI_getCorrelograms64f: ldxcIn must be at least ', I0)
  902 FORMAT('xclocMPI_getCorrelograms64f: nxcs must be ', I0)
  903 FORMAT('xclocMPI_getCorrelograms64f: Failed to get xcs on group ', I0)
      RETURN
      END
!>    @brief Gathers all the correlograms onto the root process on the global 
!>           xcloc communicator.
!>    @param[in] ldxcIn  The leading dimension of xcs.  This is defined on the root
!>                       process and must be at least nptsInXCs_.
!>    @param[in] nxcs    The number of correlograms.  This is defined on the root
!>                       process and should equal nxcsTotal_.
!>    @param[in] root    The root process ID on the xcloc global communicator.
!>                       This must be defined on all processes.
!>    @param[out] xcs    The [ldxcIn x nxcs] matrix of correlograms stored in column 
!>                       major format.  This is only accessed on the root process.
!>    @param[out] ierr   0 indicates success.
!>    @result 0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_getCorrelograms32f(ldxcIn, nxcs, root, xcs, ierr) &
      BIND(C, NAME='xclocMPI_getCorrelograms32f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldxcIn, nxcs, root
      REAL(C_FLOAT), DIMENSION(ldxcIn*nxcs), INTENT(OUT) :: xcs
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL, ALLOCATABLE :: xcsLocal(:)
      INTEGER, DIMENSION(:), ALLOCATABLE :: displs, recvCount
      INTEGER ierrLoc, interCommTarget, intraCommTarget, irecv, ldxc, mpierr, &
              nxcsInGroup, sendCount
      ierr = 0
      ierrLoc = 0
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      ! Did I compute the XCs yet?
      IF (.NOT.lhaveXCs_) THEN
         IF (myid_ == root_) WRITE(ERROR_UNIT,900)
         ierrLoc = 1
      ENDIF
      ! Verify the input arguments and set dest info 
      IF (myid_ == root) THEN
         IF (ldxcIn < nptsInXCs_ .AND. nxcs /= nxcsTotal_) THEN
            IF (ldxcIn < nptsInXCs_) WRITE(ERROR_UNIT,901) nptsInXCs_
            IF (nxcs /= nxcsTotal_) WRITE(ERROR_UNIT,902) nxcsTotal_
            ierrLoc = 1
         ENDIF
         ! Set the destination information and initialize the result to nothing 
         xcs(:) = 0.0
         interCommTarget = xcdsmInterCommID_
         intraCommTarget = xcdsmIntraCommID_
         ldxc = ldxcIn
      ENDIF
      ! Leave now if there is an error
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN
      ! Send around target information
      CALL MPI_Bcast(interCommTarget, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(intraCommTarget, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(ldxc,            1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      ! Gather the correlograms into this group's `root' process
      IF (xcdsmIntraCommID_ == intraCommTarget) ALLOCATE(xcsLocal(MAX(1,ldxc*nxcsLocal_)))
      CALL xcloc_fdxcMPI_gatherCorrelograms32f(ldxc, nxcsLocal_, intraCommTarget, &
                                               xcsLocal, ierrLoc)
      IF (ierrLoc /= 0) THEN
         IF (xcdsmIntraCommID_ == intraCommTarget) WRITE(ERROR_UNIT,903) xcdsmInterCommID_
         ierrLoc = 1
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN
      ! Gather the correlograms onto the root. 
      IF (xcdsmIntraCommID_ == interCommTarget) THEN
         sendCount = ldxc*nxcsLocal_
         IF (xcdsmIntraCommID_ == intraCommTarget) THEN
            ALLOCATE(recvCount(nxcdsmGroups_)); recvCount(:) = 0
            ALLOCATE(displs(nxcdsmGroups_)); displs(:) = 0
            displs(1) = 0 ! This is C indexed
            DO irecv=1,nxcdsmGroups_
               nxcsInGroup = myDSMXCPtr_(irecv+1) - myDSMXCPtr_(irecv)
               recvCount(irecv) = nxcsInGroup*nptsInXCs_
               IF (irecv < nxcdsmGroups_) THEN
                  displs(irecv+1) = displs(irecv) + recvCount(irecv)
               ENDIF
            ENDDO
         ENDIF
         CALL MPI_Gatherv(xcsLocal, sendCount, MPI_REAL,                    &
                          xcs, recvCount, displs,                           &
                          MPI_REAL, interCommTarget, xcdsmInterComm_, mpierr)
         IF (ALLOCATED(displs))    DEALLOCATE(displs)
         IF (ALLOCATED(recvCount)) DEALLOCATE(recvCount)
      ENDIF
      IF (ALLOCATED(xcsLocal)) DEALLOCATE(xcsLocal)
  900 FORMAT('xclocMPI_getCorrelograms32f: xcs not yet computed')
  901 FORMAT('xclocMPI_getCorrelograms32f: ldxcIn must be at least ', I0)
  902 FORMAT('xclocMPI_getCorrelograms32f: nxcs must be ', I0)
  903 FORMAT('xclocMPI_getCorrelograms32f: Failed to get xcs on group ', I0)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the input signals to cross-correlate.
!>    @param[in] ldxIn       Leading dimension of x.  This must be set on the root
!>                           process.
!>    @param[in] nptsIn      Number of points in each signal.  This must match npts_.
!>                           This must be defined on the root process.
!>    @param[in] nsignalsIn  Number of signals.  This must match nsignalsTotal_.
!>                           This must be defined the root process.
!>    @param[in] root        Root process ID on global xcloc communicator.  This must be
!>                           set on all processes.
!>    @param[in] x           This is a matrix [ldxIn x nsignalsIn] of signals stored in
!>                           column major format.  This must be set on the root process.
!>    @param[out] ierr       0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_setSignals64f(ldxIn, nptsIn, nsignalsIn, root, x, ierr) &
      BIND(C, NAME='xclocMPI_setSignals64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldxIn, nptsIn, nsignalsIn, root
      REAL(C_DOUBLE), DIMENSION(1:ldxIn*nsignalsIn), INTENT(IN) :: x
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: xloc(:)
      INTEGER i1, i2, ierrLoc, isl, isg, j1, j2, ldx, mpierr, npts, nsignals, &
              rootGroup, source
      INTEGER(KIND=MPI_ADDRESS_KIND) :: offset
      ierr = 0
      lhaveSignals_ = .FALSE. 
      lhaveXCs_ = .FALSE.
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      IF (myid_ == root) THEN
         ldx = ldxIn
         npts = nptsIn
         nsignals = nsignalsIn
         IF (.NOT.linit_ .OR. ldx < npts .OR. &
             npts /= npts_ .OR. nsignals /= nSignalsTotal_) THEN
            IF (.NOT.linit_) WRITE(ERROR_UNIT,850)
            IF (ldx < npts) WRITE(ERROR_UNIT,900) ldx, npts
            IF (npts /= npts_) WRITE(ERROR_UNIT,901) npts_
            IF (nsignals /= nSignalsTotal_) WRITE(ERROR_UNIT,902) nSignalsTotal_
            ierr = 1
         ENDIF
         IF (root /= root_) THEN
            WRITE(ERROR_UNIT,903)
            ierr = 1
         ENDIF
         ! Have the root set the signals on the RMA window
         DO isg=1,nsignals
            i1 = (isg - 1)*ldx + 1
            i2 = i1 + npts - 1
            j1 = (isg - 1)*npts + 1
            j2 = j1 + npts - 1
            signalsRMA64f_(j1:j2) = x(i1:i2)
         ENDDO
      ENDIF
      ! Ditch early because of errors
      CALL MPI_Bcast(ierr,  1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN
      ! Figure out the DSM/XC group and the rank on that group that will send the data
      IF (myid_ == root) THEN
         rootGroup = xcdsmInterCommID_
         source = xcdsmIntraCommID_
      ENDIF
      CALL MPI_Bcast(ldx,       1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(npts,      1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(nsignals,  1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(rootGroup, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(source,    1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      ! Have the roots get the data
      CALL MPI_Win_fence(0, signalRMAWindow_, mpierr)
      IF (xcdsmIntraCommID_ == source) THEN
         ALLOCATE(xloc(nsignalsLocal_*npts)); xloc(:) = 0.d0
         DO isl=1,nsignalsLocal_
            i1 = (isl - 1)*npts + 1 
            i2 = i1 + npts - 1 
            isg = l2gSignal_(l2gSignalPtr_(xcdsmInterCommID_+1)-1+isl)
            offset = (isg - 1)*ldx  ! C numbering
            j1 = (isg - 1)*ldx + 1 
            j2 = j1 + npts - 1 
            IF (myid_ == root) THEN
               xloc(i1:i2) = x(j1:j2)
            ELSE
               CALL MPI_Get(xloc(i1), npts,                                     &
                            MPI_DOUBLE_PRECISION, root, offset,                 &
                            npts, MPI_DOUBLE_PRECISION, signalRMAWindow_, mpierr)
            ENDIF
         ENDDO
      ENDIF
      CALL MPI_Win_fence(0, signalRMAWindow_, mpierr)
      ! Set the signals
      ierrLoc = 0 
      IF (nsignalsLocal_ > 0) THEN
         CALL MPI_Barrier(xcdsmIntraComm_, mpierr)
         CALL xcloc_fdxcMPI_setSignals64f(npts, npts, nsignalsLocal_,     &
                                          source, xloc, ierrLoc)
         IF (ierrLoc /= 0) THEN
            IF (xcdsmIntraCommID_ == source) WRITE(ERROR_UNIT, 910) xcdsmInterCommId_ + 1 
            ierrLoc = 1 
         ENDIF
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) THEN
         IF (myid_ == root_) WRITE(ERROR_UNIT,915)
      ELSE
         lhaveSignals_ = .TRUE.
      ENDIF
      IF (ALLOCATED(xloc)) DEALLOCATE(xloc)
  850 FORMAT("xclocMPI_setSignals64f: Module not yet initialized")
  900 FORMAT('xclocMPI_setSignals64f: Error ldx=', I0, ' < ', 'npts=', I0)
  901 FORMAT('xclocMPI_setSignals64f: Error expecting npts=', I0)
  902 FORMAT('xclocMPI_setSignals64f: Error expecting nsignals=', I0)
  903 FORMAT('xclocMPI_setSignals64f: Can only send from root on global communicator')
  910 FORMAT('xclocMPI_setSignals64f: Error setting signals on group', I0)
  915 FORMAT('xclocMPI_setSignals64f: Errors encountered while setting signals')
      RETURN
      END
!>    @brief Sets the input signals to cross-correlate.
!>    @param[in] ldxIn       Leading dimension of x.  This must be set on the root
!>                           process.
!>    @param[in] nptsIn      Number of points in each signal.  This must match npts_.
!>                           This must be defined on the root process.
!>    @param[in] nsignalsIn  Number of signals.  This must match nsignalsTotal_.
!>                           This must be defined the root process.
!>    @param[in] root        Root process ID on global xcloc communicator.  This must be
!>                           set on all processes.
!>    @param[in] x           This is a matrix [ldxIn x nsignalsIn] of signals stored in
!>                           column major format.  This must be set on the root process.
!>    @param[out] ierr       0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_setSignals32f(ldxIn, nptsIn, nsignalsIn, root, x, ierr) &
      BIND(C, NAME='xclocMPI_setSignals32f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldxIn, nptsIn, nsignalsIn, root
      REAL(C_FLOAT), DIMENSION(1:ldxIn*nsignalsIn), INTENT(IN) :: x
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL, ALLOCATABLE :: xloc(:)
      INTEGER i1, i2, ierrLoc, isl, isg, j1, j2, ldx, mpierr, npts, nsignals, &
              rootGroup, source
      INTEGER(KIND=MPI_ADDRESS_KIND) :: offset
      ierr = 0
      lhaveSignals_ = .FALSE.
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      IF (myid_ == root) THEN
         ldx = ldxIn
         npts = nptsIn
         nsignals = nsignalsIn
         IF (.NOT.linit_ .OR. ldx < npts .OR. &
             npts /= npts_ .OR. nsignals /= nSignalsTotal_) THEN
            IF (.NOT.linit_) WRITE(ERROR_UNIT,850)
            IF (ldx < npts) WRITE(ERROR_UNIT,900) ldx, npts
            IF (npts /= npts_) WRITE(ERROR_UNIT,901) npts_
            IF (nsignals /= nSignalsTotal_) WRITE(ERROR_UNIT,902) nSignalsTotal_
            ierr = 1
         ENDIF
         IF (root /= root_) THEN
            WRITE(ERROR_UNIT,903)
            ierr = 1
         ENDIF
         ! Have the root set the signals on the RMA window
         DO isg=1,nsignals
            i1 = (isg - 1)*ldx + 1
            i2 = i1 + npts - 1
            j1 = (isg - 1)*npts + 1
            j2 = j1 + npts - 1
            signalsRMA32f_(j1:j2) = x(i1:i2)
         ENDDO
      ENDIF
      ! Ditch early because of errors
      CALL MPI_Bcast(ierr,  1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN
      ! Figure out the DSM/XC group and the rank on that group that will send the data
      IF (myid_ == root) THEN
         rootGroup = xcdsmInterCommID_
         source = xcdsmIntraCommID_
      ENDIF
      CALL MPI_Bcast(ldx,       1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(npts,      1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(nsignals,  1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(rootGroup, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(source,    1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      ! Have the roots get the data
      CALL MPI_Win_fence(0, signalRMAWindow_, mpierr)
      IF (xcdsmIntraCommID_ == source) THEN
         ALLOCATE(xloc(nsignalsLocal_*npts)); xloc(:) = 0.0
         DO isl=1,nsignalsLocal_
            i1 = (isl - 1)*npts + 1 
            i2 = i1 + npts - 1 
            isg = l2gSignal_(l2gSignalPtr_(xcdsmInterCommID_+1)-1+isl)
            offset = (isg - 1)*ldx  ! C numbering
            j1 = (isg - 1)*ldx + 1 
            j2 = j1 + npts - 1 
            IF (myid_ == root) THEN
               xloc(i1:i2) = x(j1:j2)
            ELSE
               CALL MPI_Get(xloc(i1), npts,                         &
                            MPI_REAL, root, offset,                 &
                            npts, MPI_REAL, signalRMAWindow_, mpierr)
            ENDIF
         ENDDO
      ENDIF
      CALL MPI_Win_fence(0, signalRMAWindow_, mpierr)
      ! Set the signals
      ierrLoc = 0 
      IF (nsignalsLocal_ > 0) THEN
         CALL MPI_Barrier(xcdsmIntraComm_, mpierr)
         CALL xcloc_fdxcMPI_setSignals32f(npts, npts, nsignalsLocal_,     &
                                          source, xloc, ierrLoc)
         IF (ierrLoc /= 0) THEN
            IF (xcdsmIntraCommID_ == source) WRITE(ERROR_UNIT, 910) xcdsmInterCommId_ + 1 
            ierrLoc = 1 
         ENDIF
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) THEN
         IF (myid_ == root_) WRITE(ERROR_UNIT,915)
      ELSE
         lhaveSignals_ = .TRUE.
      ENDIF
      IF (ALLOCATED(xloc)) DEALLOCATE(xloc)
  850 FORMAT("xclocMPI_setSignals32f: Module not yet initialized")
  900 FORMAT('xclocMPI_setSignals32f: Error ldx=', I0, ' < ', 'npts=', I0)
  901 FORMAT('xclocMPI_setSignals32f: Error expecting npts=', I0)
  902 FORMAT('xclocMPI_setSignals32f: Error expecting nsignals=', I0)
  903 FORMAT('xclocMPI_setSignals32f: Can only send from root on global communicator')
  910 FORMAT('xclocMPI_setSignals32f: Error setting signals on group', I0)
  915 FORMAT('xclocMPI_setSignals32f: Errors encountered while setting signals')
      RETURN
      END
!>    @brief Maps a signal number to the table index.
!>    @param[in] is     Signal number.  This must be in the xcPairs table.
!>    @param[in] root   ID of root process on global communicator.
!>    @param[out] it    Table index corresponding to the signal number.  This is
!>                      Fortran indexed.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup xcloc_mpi 
      SUBROUTINE xclocMPI_signalToTableIndex(is, root, it, ierr) &
      BIND(C, NAME='xclocMPI_signalToTableIndex')
      INTEGER(C_INT), VALUE, INTENT(IN) :: is, root
      INTEGER(C_INT), INTENT(OUT) :: it, ierr
      INTEGER mpierr
      it = 0
      ierr = 0
      IF (.NOT.linit_) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1 
         RETURN
      ENDIF
      IF (root == myid_) THEN
         IF (is < 1 .OR. is > nsignalsTotal_) THEN
            WRITE(ERROR_UNIT,905) is, nsignalsTotal_
            ierr = 1 
            GOTO 500
         ENDIF
         CALL xcloc_utils_bsearch32i(nSignalsTotal_, is, uniqueSignals_, &
                                     it, ierr, .FALSE.)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,910)
            ierr = 1
            it = 0
         ENDIF
      ENDIF
  500 CONTINUE
      CALL MPI_Bcast(it,   1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) it = 0
  900 FORMAT("xclocMPI_signalToTableIndex: Module not initialized")
  905 FORMAT("xclocMPI_signalToTableIndex: Source index = ", I0, &
            " must be in range [1,", I0, "]")
  910 FORMAT("xclocMPI_signalToTableIndex: Failed to find signal index ", I0)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the travel time table on the module.
!>    @param[in] tableNumberIn  Global table index.  This is defined on the root process.
!>    @param[in] ngrdIn         Number of grid points in table.  This is defined on the
!>                              the root process.
!>    @param[in] root           The root process ID on the xcloc global communicator.
!>    @param[in] table          The
!>    @param[out] ierr          0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_setTable64f(tableNumberIn, ngrdIn, root, table, ierr) &
      BIND(C, NAME='xclocMPI_setTable64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: tableNumberIn, ngrdIn, root
      REAL(C_DOUBLE), DIMENSION(ngrdIn), TARGET, INTENT(IN) :: table
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE, TARGET :: work(:)
      DOUBLE PRECISION, CONTIGUOUS, POINTER :: ttptr(:)
      TYPE(MPI_Status) :: stat
      INTEGER commRoot, i1, i2, ierrLoc, ihaveAll, ihaveAllLoc, ig, indx, mpierr, ngrd, &
              tableNumber, localTableNumber 
      LOGICAL lneedTable
      LOGICAL(C_BOOL) lhaveAll
      ierr = 0
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      lneedTable = .FALSE.
      localTableNumber = 0
      NULLIFY(ttptr)
      IF (.NOT.linit_) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (myid_ == root) THEN
         ngrd = ngrdIn
         CALL xcloc_dsmxcMPI_getNumberOGridPointsInTable(ngrd)
         IF (ngrd /= ngrdIn) THEN
            WRITE(ERROR_UNIT,905) ngrdIn, ngrd
            ierr = 1
         ENDIF
         IF (tableNumberIn < 1 .OR. tableNumberIn > ntablesTotal_) THEN
            WRITE(ERROR_UNIT,910) tableNumberIn, nTablesTotal_
            ierr = 1
         ENDIF
         commRoot = xcdsmIntraCommID_
         tableNumber = tableNumberIn
         ttptr => table(1:ngrd) ! Associate pointer 
      ENDIF
      CALL MPI_Bcast(ierr,        1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN
      CALL MPI_Bcast(ngrd,        1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(tableNumber, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(commRoot,    1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      ! Is this signal coming my way?
      lneedTable = .FALSE.
      i1 = signalToGroupPtr_(tableNumber)
      i2 = signalToGroupPtr_(tableNumber+1) - 1
      CALL xcloc_utils_bsearch32i(i2-i1+1, xcdsmInterCommID_, signalToGroup_(i1:i2), &
                                  indx, ierr, .FALSE.)
      IF (ierr == 0) lneedTable = .TRUE.
      ! Distribute the table to the roots (intra-communicator IDs must match)
      IF (xcdsmIntraCommID_ == commRoot) THEN
         ! Pass the table to all groups that need it
         DO ig=0,nxcdsmGroups_-1
            IF (myid_ == root) THEN
               ! Do I need to send this table?
               CALL xcloc_utils_bsearch32i(i2-i1+1, ig, signalToGroup_(i1:i2), &
                                           indx, ierr, .FALSE.)
               IF (ierr == 0 .AND. ig /= xcdsmInterCommID_) THEN
                  CALL MPI_Send(ttptr, ngrd, MPI_DOUBLE_PRECISION, ig, 0, &
                                xcdsmInterComm_, mpierr) 
               ENDIF
            ELSE
               IF (ig == xcdsmInterCommID_ .AND. lneedTable) THEN
                  ALLOCATE(work(ngrd)); work(:) = 0.d0
                  CALL MPI_Recv(work, ngrd, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                                MPI_ANY_TAG, xcdsmInterComm_, stat, mpierr)
                  ttptr => work(1:ngrd) ! Associate pointer
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      ! Roots distribute the table to members of its DSM/XC group
      ierrLoc = 0
      IF (lneedTable) THEN
         IF (xcdsmIntraCommID_ == commRoot) THEN
            i1 = l2gSignalPtr_(xcdsmInterCommID_+1)
            i2 = l2gSignalPtr_(xcdsmInterCommID_+2) - 1
            CALL xcloc_utils_bsearch32i(i2-i1+1, tableNumber, l2gSignal_(i1:i2), &
                                        localTableNumber, ierrLoc, .TRUE.)
            IF (ierrLoc /= 0) THEN
               ierrLoc = 1
               localTableNumber =-1
            ENDIF
         ENDIF
         CALL MPI_Bcast(localTableNumber, 1, MPI_INTEGER, commRoot, &
                        xcdsmIntraComm_, mpierr)
         IF (localTableNumber > 0) THEN
            CALL xcloc_dsmxcMPI_setTable64f(localTableNumber, ngrd, commRoot, &
                                            ttptr, ierrLoc)
            IF (ierrLoc /= 0) THEN
               WRITE(ERROR_UNIT,915) tableNumber, xcdsmInterCommID_
               ierrLoc = 1
            ENDIF
         ENDIF
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_SUM, xclocGlobalComm_, mpierr)
      IF (ASSOCIATED(ttptr)) NULLIFY(ttptr)
      IF (ALLOCATED(work)) DEALLOCATE(work)
      ! Do I have all the tables?
      CALL xcloc_dsmxcMPI_haveAllTables(lhaveAll)
      ihaveAllLoc = 0
      IF (lhaveAll) ihaveAllLoc = 1
      CALL MPI_Allreduce(ihaveAllLoc, ihaveAll, 1, MPI_INTEGER, MPI_MIN, &
                         xclocGlobalComm_, mpierr)
      IF (ihaveAll == 1) lhaveAllTables_ = .TRUE.
      IF (lhaveAllTables_ .AND. myid_ == root_ .AND. verbose_ > XCLOC_PRINT_WARNINGS) &
      WRITE(OUTPUT_UNIT,920)
  900 FORMAT("xclocMPI_setTable64f: Module not initialized")
  905 FORMAT("xclocMPI_setTable64f: ngrdIn = ", I0, " should be equal ngrd =", I0)
  910 FORMAT("xclocMPI_setTable64f: table number = ", I0, " should be in range[1,", &
             I0, "]")
  915 FORMAT("xclocMPI_setTable64f: Failed to set table = ", I0, " on group ", I0)
  920 FORMAT("xclocMPI_setTable64f: All tables set")
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE xclocMPI_getImage32f(ngrd, root, image, ierr) &
      BIND(C, NAME='xclocMPI_getImage32f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngrd, root
      REAL(C_FLOAT), DIMENSION(ngrd), INTENT(OUT) :: image
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL, POINTER, CONTIGUOUS, DIMENSION(:) :: localImagePtr
      REAL, ALLOCATABLE, DIMENSION(:) :: imageWork
      INTEGER, DIMENSION(:), ALLOCATABLE :: displs(:), recvcount(:)
      INTEGER ierrLoc, intraCommTarget, interCommTarget, irecv, mpierr, ngrdLocal, &
              sendCount
      
      ierr = 0
      ierrLoc = 0
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      IF (.NOT.lhaveImage_) ierrLoc = 1
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) THEN 
         IF (myid_ == root_) WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (myid_ == root) THEN
         IF (ngrd < ngrdTotal_) THEN
            WRITE(ERROR_UNIT,901) ngrd
            ierrLoc = 1
         ENDIF
         image(:) = 0.0
         interCommTarget = xcdsmInterCommID_
         intraCommTarget = xcdsmIntraCommID_
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN
      ! Get the pointer to the local image
      CALL MPI_Bcast(intraCommTarget, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      CALL MPI_Bcast(interCommTarget, 1, MPI_INTEGER, root, xclocGlobalComm_, mpierr)
      NULLIFY(localImagePtr)
      CALL xcloc_dsmxcMPI_getLocalImagePtr(localImagePtr, ierrLoc)
      ngrdLocal = SIZE(localImagePtr)
      IF (xcdsmInterCommID_ == interCommTarget) THEN
         ALLOCATE(imageWork(ngrdLocal)); imageWork(:) = 0.0
      ENDIF
      CALL MPI_Reduce(localImagePtr, imageWork, ngrdLocal, MPI_REAL, MPI_SUM, root_, &
                      xcdsmInterComm_, mpierr)
      NULLIFY(localImagePtr)
      ! Gather onto master
      IF (xcdsmInterCommID_ == interCommTarget) THEN
         sendCount = ngrdLocal
         IF (xcdsmIntraCommID_ == intraCommTarget) THEN
            ALLOCATE(displs(xcGroupSize_));    displs(:) = 0
            ALLOCATE(recvcount(xcGroupSize_)); recvcount(:) = 0
            displs(1) = 0 ! C indexed
            DO irecv=1,xcGroupSize_
               recvcount(irecv) = grdPointsPtr_(irecv+1) - grdPointsPtr_(irecv)
               IF (irecv < xcGroupSize_) &
               displs(irecv+1) = displs(irecv) + recvcount(irecv)
            ENDDO
         ENDIF
         CALL MPI_Gatherv(imageWork, sendCount, MPI_REAL,                   &
                          image, recvCount, displs,                         &
                          MPI_REAL, intraCommTarget, xcdsmIntraComm_, mpierr)
         IF (ALLOCATED(displs))    DEALLOCATE(displs)
         IF (ALLOCATED(recvcount)) DEALLOCATE(recvcount)
if (xcdsmIntraCommID_ == intraCommTarget) print *, maxloc(image), maxval(image)
      ENDIF

  900 FORMAT('xclocMPI_getImage32f: Image not yet computed')
  901 FORMAT('xclocMPI_getImage32f: ngrd must be at least ', I0)
      RETURN
      END

!>    @brief Gets the maximum value of the image.
!>    @param[out] maxIndex  The index of the maximum.  This is Fortran indexed.
!>    @param[out] maxValue  Maximum value corresponding to the maxIndex.
!>    @param[out] ierr      0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_getImageMax(maxIndex, maxValue, ierr) &
      BIND(C, NAME='xclocMPI_getImageMax')
      REAL(C_FLOAT), INTENT(OUT) :: maxValue
      INTEGER(C_INT), INTENT(OUT) :: maxIndex, ierr
      REAL, ALLOCATABLE, DIMENSION(:) :: imageWork, xmax, xmaxWork
      INTEGER, ALLOCATABLE, DIMENSION(:) :: imax, imaxWork
      REAL, POINTER, CONTIGUOUS, DIMENSION(:) :: localImagePtr
      INTEGER i, ierrLoc, mpierr, ngrdLocal
      ierr = 0
      maxIndex = 1
      maxValue = 0.0
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      ierrLoc = 0
      IF (.NOT.lhaveImage_) ierrLoc = 1
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) THEN
         IF (myid_ == root_) WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      ! Sum the images onto the root DSM group 
      NULLIFY(localImagePtr)
      CALL xcloc_dsmxcMPI_getLocalImagePtr(localImagePtr, ierrLoc)
      ngrdLocal = SIZE(localImagePtr)
 print *, maxval(localImagePtr)
      IF (xcdsmInterCommID_ == root_) THEN
         ALLOCATE(imageWork(ngrdLocal))
      ELSE
         ALLOCATE(imageWork(1))
      ENDIF
      CALL MPI_Reduce(localImagePtr, imageWork, ngrdLocal, MPI_REAL, MPI_SUM, root_, &
                      xcdsmInterComm_, mpierr) 
      NULLIFY(localImagePtr)
      ! Get the max
      IF (xcdsmInterCommID_ == root_) THEN
         ALLOCATE(imaxWork(xcGroupSize_)); imaxWork(:) = 0
         ALLOCATE(xmaxWork(xcGroupSize_)); xmaxWork(:) = 0.0
         ALLOCATE(imax(xcGroupSize_))
         ALLOCATE(xmax(xcGroupSize_))
         imaxWork(xcdsmIntraCommID_+1) = MAXLOC(imageWork, 1)
         xmaxWork(xcdsmIntraCommID_+1) = imageWork(imaxWork(xcdsmIntraCommID_+1))
print *, 'look', xmaxWork
         CALL MPI_Reduce(imaxWork, imax, xcGroupSize_, MPI_INTEGER, MPI_SUM, root_, &
                         xcdsmIntraComm_, mpierr)
         CALL MPI_Reduce(xmaxWork, xmax, xcGroupSize_, MPI_REAL,    MPI_SUM, root_, &
                         xcdsmIntraComm_, mpierr) 
         IF (xcdsmIntraCommID_ == root_) THEN
            i = MAXLOC(xmax, 1)
            maxIndex = grdPointsPtr_(i) + imax(i) - 1 
            maxValue = xmax(i)
print *, maxIndex, maxValue
         ENDIF
         IF (ALLOCATED(imaxWork)) DEALLOCATE(imaxWork)
         IF (ALLOCATED(xmaxWork)) DEALLOCATE(xmaxWork)
         IF (ALLOCATED(imax)) DEALLOCATE(imax)
         IF (ALLOCATED(xmax)) DEALLOCATE(xmax)
      ENDIF
      IF (ALLOCATED(imageWork)) DEALLOCATE(imageWork)
  900 FORMAT('xclocMPI_getImageMax: Image not yet computed')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief After setting the travel time tables and signals this will:
!>           1. Compute the phase or cross-correlograms depending on the signal migration
!>              policy.
!>           2. Either compute the envelope or RMS of the correlograms or leave the 
!>              correlograms unchanged.
!>           3. Migrate the (processed) correlograms.
!>    @param[out] ierr   0 indicates success.
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_compute(ierr) &
      BIND(C, NAME='xclocMPI_compute')
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xcs64f
      REAL, DIMENSION(:), ALLOCATABLE :: xcs32f
      DOUBLE PRECISION timer_end, timer_start
      INTEGER ierrLoc, mpierr
      ierr = 0
      lhaveImage_ = .FALSE.
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      ierrLoc = 0
      IF (.NOT.linit_) THEN
         WRITE(ERROR_UNIT,850) myid_
         ierrLoc = 1
      ENDIF
      IF (.NOT.lhaveSignals_) THEN
         WRITE(ERROR_UNIT,855) myid_
         ierrLoc = 1
      ENDIF
      IF (.NOT.lhaveAllTables_) THEN
         WRITE(ERROR_UNIT,860) myid_
         ierrLoc = 1
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) RETURN
      timer_start = MPI_Wtime()
      IF (xcTypeToMigrate_ == XCLOC_MIGRATE_PHASE_XCS) THEN
         CALL xcloc_fdxcMPI_computePhaseCorrelograms(ierrLoc)
      ELSE
         CALL xcloc_fdxcMPI_computeCrossCorrelograms(ierrLoc)
      ENDIF
      IF (ierrLoc /= 0) THEN
         WRITE(ERROR_UNIT,865) myid_
         ierrLoc = 1
         GOTO 500
      ENDIF
      ! This timing isn't the most accurate.  Technically I should block.
      timer_end = MPI_Wtime()
      IF (myid_ == root_ .AND. verbose_ > XCLOC_PRINT_WARNINGS) &
      WRITE(OUTPUT_UNIT,905) timer_end - timer_start
      timer_start = MPI_Wtime()
      IF (precision_ == XCLOC_SINGLE_PRECISION) THEN
         IF (xcdsmIntraCommID_ == root_) ALLOCATE(xcs32f(nxcsLocal_*nptsInXCs_))
         CALL xcloc_fdxcMPI_gatherCorrelograms32f(nptsInXCs_, nxcsLocal_, root_, &
                                                  xcs32f, ierrLoc)
         IF (ierrLoc /= 0) THEN 
            WRITE(ERROR_UNIT,870) myid_
            GOTO 500
         ENDIF
         IF (xcdsmIntraCommID_ == root_) THEN
            CALL xcloc_spxc_filterXCsInPlace32f(nptsInXCs_, nptsInXCs_, nxcsLocal_, &
                                                xcs32f, ierrLoc)
            IF (ierrLoc /= 0) THEN 
               WRITE(ERROR_UNIT,875) xcdsmInterCommID_
               ierrLoc = 1
            ENDIF
         ENDIF
         CALL MPI_Bcast(ierrLoc, 1, MPI_INTEGER, root_, xcdsmIntraComm_, mpierr)   
         IF (ierrLoc /= 0) GOTO 500
         ! Set the correlograms
         CALL xcloc_dsmxcMPI_setCorrelograms32f(nptsInXCs_, nptsInXCs_,           &
                                                nxcsLocal_, root_, xcs32f, ierrLoc)
         IF (ierrLoc /= 0) THEN
            IF (xcdsmIntraCommID_ == root_) WRITE(ERROR_UNIT,880) xcdsmInterCommID_
            ierrLoc = 1
            GOTO 500
         ENDIF
         IF (ALLOCATED(xcs32f)) DEALLOCATE(xcs32f)
      ELSE
         IF (xcdsmIntraCommID_ == root_) ALLOCATE(xcs64f(nxcsLocal_*nptsInXCs_))
         CALL xcloc_fdxcMPI_gatherCorrelograms64f(nptsInXCs_, nxcsLocal_, root_, &
                                                  xcs64f, ierrLoc)
         IF (ierrLoc /= 0) THEN
            WRITE(ERROR_UNIT,870) myid_
            GOTO 500
         ENDIF
         IF (xcdsmIntraCommID_ == root_) THEN
            CALL xcloc_spxc_filterXCsInPlace64f(nptsInXCs_, nptsInXCs_, nxcsLocal_, &
                                                xcs64f, ierrLoc)
            IF (ierrLoc /= 0) THEN
               WRITE(ERROR_UNIT,875) xcdsmInterCommID_
               ierrLoc = 1
            ENDIF
         ENDIF
         CALL MPI_Bcast(ierrLoc, 1, MPI_INTEGER, root_, xcdsmIntraComm_, mpierr)   
         IF (ierrLoc /= 0) GOTO 500
         ! Set the correlograms
         CALL xcloc_dsmxcMPI_setCorrelograms64f(nptsInXCs_, nptsInXCs_,           &
                                                nxcsLocal_, root_, xcs64f, ierrLoc)
         IF (ierrLoc /= 0) THEN 
            IF (xcdsmIntraCommID_ == root_) WRITE(ERROR_UNIT,880) xcdsmInterCommID_
            ierrLoc = 1
            GOTO 500
         ENDIF
         IF (ALLOCATED(xcs64f)) DEALLOCATE(xcs64f)
      ENDIF
      timer_end = MPI_Wtime()
      IF (verbose_ > XCLOC_PRINT_WARNINGS .AND. myid_ == root_) &
      WRITE(OUTPUT_UNIT,910) timer_end - timer_start
      ! Compute the DSM but leave the image distributed across all processes
      timer_start = MPI_Wtime()
      CALL xcloc_dsmxcMPI_compute(ierrLoc)
      IF (ierrLoc /= 0) THEN
         IF (xcdsmIntraCommID_ == root_) WRITE(ERROR_UNIT,885) xcdsmInterCommID_
         ierrLoc = 1
      ENDIF
      timer_end = MPI_Wtime()
      IF (verbose_ > XCLOC_PRINT_WARNINGS .AND. myid_ == root_) &
      WRITE(OUTPUT_UNIT,915) timer_end - timer_start
      ! Check for errors
  500 CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, xclocGlobalComm_, mpierr)
      IF (ierr /= 0) THEN
         IF (myid_ == root_) WRITE(ERROR_UNIT,950)
      ELSE
         lhaveXCs_ = .TRUE.
         lhaveImage_ = .TRUE.
      ENDIF
  850 FORMAT('xclocMPI_compute: Module not initialized on rank ', I0)
  855 FORMAT('xclocMPI_compute: Signals not yet set on rank ', I0)
  860 FORMAT('xclocMPI_compute: Travel time tables not yet set on rank ', I0)
  865 FORMAT('xclocMPI_compute: Failed to compute correlograms on rank ', I0)
  870 FORMAT('xclocMPI_compute: Failed to get correlograms on rank ', I0)
  875 FORMAT('xclocMPI_compute: Failed to process correlograms on group ', I0)
  880 FORMAT('xclocMPI_compute: Failed to set correlograms on group ', I0)
  885 FORMAT('xclocMPI_compute: Failed to compute DSM on group ', I0)
  905 FORMAT('xclocMPI_compute: Correlogram computation time=', F8.4, 's')
  910 FORMAT('xclocMPI_compute: Correlogram processing time=', F8.4, 's')
  915 FORMAT('xclocMPI_compute: Correlogram migration time=', F8.4, 's')
  950 FORMAT('xclocMPI_compute: Errors detected in processing')
      RETURN
      END

END MODULE
