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
      !> Maps from the a signal to the DSM/XC group to which this signal belongs.
      !> This has dimension [signalToGroupPtr_(nsignalsTotal_+1)-1].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: signalToGroup_(:)
 
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: xcPairsLocal_(:)
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
      !> Total number of signals.
      INTEGER, PRIVATE, SAVE :: nsignalsTotal_ = 0
      !> @ingroup xcloc_mpi
      !> Total number of travel time tables.
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
      !> Number of processes in each cross-correloation group.
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
      !> Flag indicating that communicator has to be destroyed.
      LOGICAL, PRIVATE, SAVE :: lfreeComm_ = .FALSE.
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
      INTEGER i1, i2, ig, ierrLoc, indx, is, ixc, mpierr, nalloc, nsignalPairs, nsignals, &
              nsg, nunique
      TYPE(MPI_Comm) comm
      INTEGER(KIND=MPI_ADDRESS_KIND) nallocMPI
      INTEGER arrayShape(1), sizeDouble
      comm%MPI_VAL = INT(fcomm, KIND(comm%MPI_VAL))
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
!print *, jwork(signalToGroupPtr_(is):signalToGroupPtr_(is+1)-1)
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
      IF (ALLOCATED(xcPairsLocal_))     DEALLOCATE(xcPairsLocal_)
!     IF (ASSOCIATED(signalsRMA32f_))   CALL MPI_Free_mem(signalsRMA32f_, mpierr)
      IF (ASSOCIATED(signalsRMA64f_))   CALL MPI_Free_mem(signalsRMA64f_, mpierr)
if (allocated(nDSMXCsPerGroup_)) deallocate(nDSMXCsPerGroup_)
      root_ = 0
      nprocs_ = 0
      xcdsmIntraCommID_ = 0
      xcdsmInterCommID_ = 0
      lfreeComm_ = .FALSE.
      lfreeWindow_ = .FALSE.
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
      IF (ierr /= 0 .AND. myid_ == root_) WRITE(ERROR_UNIT,915)
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
      IF (ierr /= 0 .AND. myid_ == root_) WRITE(ERROR_UNIT,915)
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
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the travel time table on the module.
!>    @param[in] tableNumberIn  Global table index.  This is defined on the root process.
!>    @param[in] ngrdIn         Number of grid points in table.  This is defined on the
!>                              the root process.
!>    @param[in] root           The root process ID on the xcloc global communicator.
!>    @param[in] table          The
!>    @ingroup xcloc_mpi
      SUBROUTINE xclocMPI_setTable64f(tableNumberIn, ngrdIn, root, table, ierr) &
      BIND(C, NAME='xclocMPI_setTable64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: tableNumberIn, ngrdIn, root
      REAL(C_DOUBLE), DIMENSION(ngrdIn), TARGET, INTENT(IN) :: table
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE, TARGET :: work(:)
      DOUBLE PRECISION, CONTIGUOUS, POINTER :: ttptr(:)
      TYPE(MPI_Status) :: stat
      INTEGER commRoot, i1, i2, ierrLoc, ig, indx, mpierr, ngrd, tableNumber, &
              localTableNumber 
      LOGICAL lneedTable
      ierr = 0
      lneedTable = .FALSE.
      localTableNumber = 0
      NULLIFY(ttptr)
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
               IF (lneedTable) THEN
                  ALLOCATE(work(ngrd)); work(:) = 0.d0
                  CALL MPI_Recv(ttptr, ngrd, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                                MPI_ANY_TAG, xcdsmInterComm_, stat, mpierr)
                  ttptr => work(1:ngrd) ! Associate pointer
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      ! Roots distribute the table to members of its DSM/XC group
      ierrLoc = 0
      IF (lneedTable) THEN
         CALL xcloc_dsmxcMPI_setTable64f(localTableNumber, ngrd, xcdsmIntraCommID_, &
                                         ttptr, ierrLoc)
         IF (ierrLoc /= 0) THEN
            ierrLoc = 1
         ENDIF
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_SUM, xclocGlobalComm_, mpierr)
      IF (ASSOCIATED(ttptr)) NULLIFY(ttptr)
      IF (ALLOCATED(work)) DEALLOCATE(work)
  900 FORMAT("xclocMPI_setTable64f: Module not initialized")
  905 FORMAT("xclocMPI_setTable64f: ngrdIn=", I6, "should be ngrd=", I6)
  910 FORMAT("xclocMPI_setTable64f: table number=", I6, "should be in range[1,", I6, "]")
      RETURN
      END

END MODULE
