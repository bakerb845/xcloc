!> @defgroup fdxcmpi Parallel Computation of Correlograms 
!> @ingroup xcloc
!> @brief Computes the cross-correlograms via the Fourier transform with MPI.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_FDXC_MPI
      USE ISO_C_BINDING
      USE ISO_FORTRAN_ENV
      USE MPI_F08
      USE XCLOC_CONSTANTS
      USE XCLOC_MEMORY
      USE XCLOC_FDXC
      USE XCLOC_UTILS
#ifdef _OPENMP
      USE OMP_LIB
#endif
      IMPLICIT NONE
      !> MPI communicator.
#if defined(__INTEL_COMPILER)
      TYPE(MPI_Comm), PRIVATE, SAVE :: comm_ = MPI_COMM_WORLD
#else
      TYPE(MPI_Comm), PRIVATE, SAVE :: comm_
#endif
!     !> RMA window for setting signals.
!     TYPE(MPI_Win), PRIVATE, SAVE :: signalRMAWindow_
!     !> Memory on the RMA signal window for holding the input signals.
!     DOUBLE PRECISION, POINTER, DIMENSION(:), PRIVATE, SAVE :: inputSignalsRMA64f_
!     !> Memory on the RMA signal window for holding the input signals.
!     REAL, POINTER, DIMENSION(:), PRIVATE, SAVE :: inputSignalsRMA32f_
!     !> Pointer to the data on the window.
!     TYPE(C_PTR), SAVE :: dataPtr_ = C_NULL_PTR
      !> Root process ID.
      INTEGER, PRIVATE, SAVE :: root_ = 0
      !> Number of processes on the communicator.
      INTEGER, PRIVATE, SAVE :: nprocs_ = 0
      !> My rank on the communicator.
      INTEGER, PRIVATE, SAVE :: myid_ = MPI_UNDEFINED
      !> A list of each processes list of unique signals global signal IDs.
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: uniqueGlobalSignalIDs_(:)
      !> This maps from the ip'th process to the start index of uniqueGlobalSignalIDs_.
      !> This has dimension [nprocs_ + 1].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: uniqueGlobalSignalIDPtr_(:)
      !> This is the number of cross-correlation the ip'th process must perform.
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: nXCsPerProcess_(:)
      !> This maps from the isLoc'th local signal number to the global signal number.
      !> This has dimension [nsignalsLocal_].
      INTEGER, ALLOCATABLE, PRIVATE, SAVE :: local2GlobalSignal_(:)
      !> Length of input signals.
      INTEGER, PRIVATE, SAVE :: npts_ = 0
      !> The leading dimension of the input signals RMA buffer.
      INTEGER, PRIVATE, SAVE :: lds_ = 0 
      !> The length of the signals to transform.  This can mitigate the pathologic
      !> case where the signal transform lengths are large (semi)prime numbers which
      !> make the DFT very expensive.
      INTEGER, PRIVATE, SAVE :: nptsPad_ = 0
      !> Length of the time domain cross-correloagrams.
      INTEGER, PRIVATE, SAVE :: nptsInXCs_ = 0
      !> Total number of input signals to correlate.
      INTEGER, PRIVATE, SAVE :: nSignalsTotal_ = 0
      !> The number of input signals specific to a process.
      INTEGER, PRIVATE, SAVE :: nSignalsLocal_ = 0 
      !> Total number of cross-correlations
      INTEGER, PRIVATE, SAVE :: nXCsTotal_ = 0
      !> Number of cross-correlations that I must perform.
      INTEGER, PRIVATE, SAVE :: nXCsLocal_ = 0
      !> Controls verbosity.
      INTEGER, PRIVATE, SAVE :: verbose_ = XCLOC_PRINT_WARNINGS
      !> Accuracy of the MKL computations.
      INTEGER, PRIVATE, SAVE :: accuracy_ = XCLOC_HIGH_ACCURACY
      !> Precision of module.
      INTEGER, PRIVATE, SAVE :: precision_ = XCLOC_SINGLE_PRECISION
      !> If true then this process will be computing cross-correlations. 
      LOGICAL, PRIVATE, SAVE :: ldoXC_ = .FALSE.
      !> If true then free the communicator.
      LOGICAL, PRIVATE, SAVE :: lfreeComm_ = .FALSE.
      !> If true then free the input signal RMA window
!     LOGICAL, PRIVATE, SAVE :: lfreeSignalRMA_ = .FALSE. 
      !> If true then free the cross-correlation RMA window
      INTEGER(KIND=8), PRIVATE, PARAMETER :: alignPage_ = 64
      !LOGICAL, PRIVATE, SAVE :: lfree
      PUBLIC :: xcloc_fdxcMPI_initialize
      CONTAINS
!========================================================================================!
!                                      Begin the Code                                    !
!========================================================================================!
!>    @brief Initializes the parallel frequency domain cross-correlation calculator.
!>    @param[in] comm      MPI communicator.  This must be defined on all processes.
!>    @param[in] root      ID of root process.  This will likely be 0 and must be 
!>                         defined on all processes.
!>    @param[in] npts      Number of points in each input signal.  This is defined on
!>                         the root process.
!>    @param[in] nptsPad   A tuning parameter to mitigate DFT lengths that could
!>                         potentially be large semi-prime numbers.  This is defined
!>                         on the root process.
!>    @param[in] nxcs      Number of cross-correlations.  This is defined on the
!>                         root process.
!>    @param[in] xcPairs   This is a [2 x nxcs] matrix in column major format where
!>                         the indices, (2*(ixc-1)+1, 2*(ixc-1)+2), map to the
!>                         (i,j)'th signal pair comprising a correlation.
!>                         This is defined on the root process.
!>    @param[in] verbose   Controls the verbosity of the module.  This is defined
!>                         on the root process.
!>    @param[in] prec      The precision of the module.  This is defined on the
!>                         root process.
!>    @param[in] accuracy  Controls the accuracy of the vector calculations in MKL.
!>                         This is defined on the root process.
!>    @param[out] ierr     0 indicates sucess.
!>    @ingroup fdxcmpi
      SUBROUTINE xcloc_fdxcMPI_initialize(comm, root,                     &
                                          npts, nptsPad,                  &
                                          nxcs, xcPairs,                  &
                                          verbose, prec, accuracy, ierr)  &
      BIND(C, NAME='xcloc_fdxcMPI_initialize')
      TYPE(MPI_Comm),  VALUE, INTENT(IN) :: comm
      !INTEGER(C_INT), VALUE, INTENT(IN) :: fcomm !TODO could be problematic w/ *finter.h
      INTEGER(C_INT), VALUE, INTENT(IN) :: root, npts, nptsPad, nxcs, &
                                           verbose, prec, accuracy
      INTEGER(C_INT), INTENT(IN) :: xcPairs(2*nxcs) !(*)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER, ALLOCATABLE, DIMENSION(:) :: l2gSignal, myXCs,  myXCPtr, signalList, &
                                            work, xcPairsWork, xcPairsLocal
      INTEGER i, i1, i2, ierrLocal, indx, ip, mpierr, nsignals, nsloc, nwork
      ! Release module in case someone called this twice.
      CALL xcloc_fdxcMPI_finalize()
      ! Get communicator size and my rank
      root_ = root
      CALL MPI_COMM_RANK(comm, myid_, mpierr)
      CALL MPI_COMM_SIZE(comm, nprocs_, mpierr)
      ! Some basic checks by root process
      ierr = 0
      IF (myid_ == root_) THEN
         root_ = root
         IF (npts < 1 .OR. nptsPad < npts .OR. nxcs < 1) THEN
            IF (npts < 1) WRITE(ERROR_UNIT,905) npts 
            IF (nptsPad < npts) WRITE(ERROR_UNIT,907) nptsPad, npts
            IF (nxcs < 1) WRITE(ERROR_UNIT,908) nxcs
            ierr = 1
         ENDIF
         nsignals = MAXVAL(xcPairs(1:2*nxcs))
         IF (nsignals < 2) THEN
            WRITE(ERROR_UNIT,906) nsignals
            ierr = 1
         ENDIF
         nsignals = MAXVAL(xcPairs(1:2*nxcs))
         IF (nsignals < 2) THEN
            WRITE(ERROR_UNIT,906) nsignals
            ierr = 1
         ENDIF
         IF (.NOT. xcloc_constants_isValidPrecision(prec)) ierr = 1 
         IF (.NOT. xcloc_constants_isValidAccuracy(accuracy)) ierr = 1
         IF (ierr /= 0) GOTO 500
         npts_     = npts
         nptsPad_  = nptsPad
         nSignalsTotal_ = nsignals
         verbose_ = verbose
         nptsInXCs_ = 2*nptsPad_ - 1   ! Length of the cross-correlations
         nXCsTotal_ = nxcs
         accuracy_ = accuracy
         precision_ = prec
         ! build up the workspaces
         IF (precision_ == XCLOC_SINGLE_PRECISION) THEN
            lds_ = xcloc_memory_padLength(alignPage_, SIZEOF(1.0),  nptsPad_)
         ELSE
            lds_ = xcloc_memory_padLength(alignPage_, SIZEOF(1.d0), nptsPad_)
         ENDIF
         ! Need to partition the work
         ALLOCATE(myXCs(nXCsTotal_))
         ALLOCATE(myXCPtr(nprocs_+1))
         CALL xcloc_utils_partitionTasks(nXCsTotal_, nprocs_, myXCPtr, myXCs, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,920)
            GOTO 500
         ENDIF
         ! Tabulate the number of XCs per process
         ALLOCATE(nXCsPerProcess_(nprocs_))
         DO ip=1,nprocs_
            nXCsPerProcess_(ip) = myXCPtr(ip+1) - myXCPtr(ip)
            IF (verbose_ > XCLOC_PRINT_WARNINGS) THEN
               WRITE(OUTPUT_UNIT,910) ip-1, nXCsPerProcess_(ip) !myXCPtr(ip+1) - myXCPtr(ip)
            ENDIF
         ENDDO
         ! Now that the work is partitioned figure out which processes get which signals.
         ALLOCATE(signalList(2*nXCsTotal_)); signalList(:) = 0
         ALLOCATE(uniqueGlobalSignalIDPtr_(nprocs_+1)); uniqueGlobalSignalIDPtr_(:) = 0
         ALLOCATE(work(nprocs_*nSignalsTotal_)); work(:) =-1 
         uniqueGlobalSignalIDPtr_(1) = 1
         DO ip=1,nprocs_
            i1 = myXCPtr(ip)
            i2 = myXCPtr(ip+1) - 1
            i1 = 2*(i1 - 1) + 1
            i2 = 2*i2
            nwork = i2 - i1 + 1
            nsloc = 0
            IF (nwork > 0) THEN
               signalList(1:nwork) = xcPairs(i1:i2)
               CALL xcloc_utils_unique32s(nwork, signalList, nsloc, &
                                          l2gSignal, ierr)
               indx = uniqueGlobalSignalIDPtr_(ip)
               work(indx:indx+nsloc-1) = l2gSignal(1:nsloc)
               IF (ALLOCATED(l2gSignal)) DEALLOCATE(l2gSignal)
            ENDIF
            uniqueGlobalSignalIDPtr_(ip+1) = uniqueGlobalSignalIDPtr_(ip) + nsloc
         ENDDO
         nwork = uniqueGlobalSignalIDPtr_(nprocs_+1) - 1
         ALLOCATE(uniqueGlobalSignalIDs_(nwork))
         uniqueGlobalSignalIDs_(:) =-1
         uniqueGlobalSignalIDs_(1:nwork) = work(1:nwork)
         IF (ALLOCATED(signalList)) DEALLOCATE(signalList)
         IF (ALLOCATED(work)) DEALLOCATE(work)
      ENDIF
  500 CONTINUE
      CALL MPI_BCAST(ierr, 1, MPI_INTEGER, root_, comm, mpierr)
      IF (ierr /= 0) RETURN
      ! Copy the communicator
      CALL MPI_COMM_DUP_WITH_INFO(comm, MPI_INFO_NULL, comm_, mpierr)
      lfreeComm_ = .TRUE.
      ! Set some basic information
      CALL MPI_Bcast(root_,          1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(npts_,          1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(nptsPad_,       1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(nSignalsTotal_, 1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(verbose_,       1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(accuracy_,      1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(precision_,     1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(nptsInXCs_,     1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(nXCsTotal_,     1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(lds_,           1, MPI_INTEGER, root, comm_, mpierr)
      IF (.NOT.ALLOCATED(myXCs))       ALLOCATE(myXCs(nXCsTotal_))
      IF (.NOT.ALLOCATED(myXCPtr))     ALLOCATE(myXCPtr(nprocs_+1))
      IF (.NOT.ALLOCATED(xcPairsWork)) ALLOCATE(xcPairsWork(2*nXCsTotal_))
      IF (.NOT.ALLOCATED(nXCsPerProcess_)) ALLOCATE(nXCsPerProcess_(nprocs_))
      IF (.NOT.ALLOCATED(uniqueGlobalSignalIDPtr_)) THEN
         ALLOCATE(uniqueGlobalSignalIDPtr_(nprocs_+1))
      ENDIF
      IF (myid_ == root_) xcPairsWork(1:2*nXCsTotal_) = xcPairs(1:2*nXCsTotal_)
      CALL MPI_Bcast(myXCs,       nXCsTotal_,   MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(myXCPtr,     nprocs_+1,    MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(xcPairsWork, 2*nXCsTotal_, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(nXCsPerProcess_, nprocs_,  MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(uniqueGlobalSignalIDPtr_, nprocs_+1, MPI_INTEGER, &
                     root_, comm_, mpierr)
      nwork = uniqueGlobalSignalIDPtr_(nprocs_+1) - 1
      IF (.NOT.ALLOCATED(uniqueGlobalSignalIDs_)) ALLOCATE(uniqueGlobalSignalIDs_(nwork))
      CALL MPI_Bcast(uniqueGlobalSignalIDs_, nwork, MPI_INTEGER, root_, comm_, mpierr)
      ! Figure out the number of cross-correlations that I need to perform
      nXCsLocal_ = myXCPtr(myid_+2) - myXCPtr(myid_+1)
      CALL MPI_ALLREDUCE(nXCsLocal_, nwork, 1, MPI_INTEGER, MPI_SUM, comm_, mpierr)
      IF (nwork /= nXCsTotal_) THEN
         IF (myid_ == root_) WRITE(ERROR_UNIT,911) nwork, nXCsTotal_
         ierr = 1
         RETURN
      ENDIF
      ! Figure out the unique signals I must concern myself with
      ALLOCATE(signalList(2*MAX(1, nXCsLocal_))); signalList(:) = 0
      ierrLocal = 0
      nSignalsLocal_ = 0
      IF (nXCsLocal_ > 0) THEN
         i1 = uniqueGlobalSignalIDPtr_(myid_+1)
         i2 = uniqueGlobalSignalIDPtr_(myid_+2) - 1
         nSignalsLocal_ = i2 - i1 + 1
         ALLOCATE(local2GlobalSignal_(nSignalsLocal_))
         local2GlobalSignal_(1:nSignalsLocal_) = uniqueGlobalSignalIDs_(i1:i2)
         i1 = myXCPtr(myid_+1)
         i2 = myXCPtr(myid_+2) - 1
         i1 = 2*(i1 - 1) + 1
         i2 = 2*i2
         signalList(1:2*nXCsLocal_) = xcPairsWork(i1:i2)
         IF (MINVAL(local2GlobalSignal_) < MINVAL(signalList) .OR. &
             MAXVAL(local2GlobalSignal_) > MAXVAL(signalList) .OR. ierr /= 0) THEN
             WRITE(ERROR_UNIT,912) myid_
             ierrLocal = 1
         ENDIF
         ! Make a local xcPairs 
         ALLOCATE(xcPairsLocal(2*nxcsLocal_)); xcPairsLocal(:) =-1
         DO i=i1,i2
            ! The index of global signal list that matches the signal index to
            ! is the local signal index.
            CALL xcloc_utils_bsearch32i(nSignalsLocal_, xcPairsWork(i),  &
                                        local2GlobalSignal_, indx, ierr)
            IF (ierr /= 0) THEN
               WRITE(ERROR_UNIT,913) myid_
               ierrLocal = 1
            ENDIF
            xcPairsLocal(i-i1+1) = indx
         ENDDO
         IF (MINVAL(xcPairsLocal) < 1 .OR. MINVAL(xcPairsLocal) > nSignalsLocal_) THEN
            WRITE(ERROR_UNIT,914)
            ierrLocal = 1
         ENDIF
         ! Can now initialize
         CALL xcloc_fdxc_initialize(npts_, nptsPad_,                     &
                                    nxcsLocal_, xcPairsLocal,            &
                                    verbose_, precision_, accuracy_, ierr) 
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,915) myid_
            ierrLocal = 1
         ENDIF
         ldoXC_ = .TRUE.
      ELSE
         ldoXC_ = .FALSE.
      ENDIF
      CALL MPI_ALLREDUCE(ierrLocal, ierr, 1, MPI_INTEGER, MPI_MAX, comm_, mpierr) 
      IF (ierr /= 0) RETURN
      ! Free workspace
      IF (ALLOCATED(myXCs))       DEALLOCATE(myXCs)
      IF (ALLOCATED(myXCPtr))     DEALLOCATE(myXCPtr)
      IF (ALLOCATED(signalList))  DEALLOCATE(signalList)
      IF (ALLOCATED(xcPairsWork)) DEALLOCATE(xcPairsWork)
      ! Format statements
  905 FORMAT("xcloc_fdxcMPI_initialize: npts=", I8, " must be positive")
  906 FORMAT("xcloc_fdxcMPI_initialize: nsignals=", I8, " must be at least 2")
  907 FORMAT("xcloc_fdxcMPI_initialize: nptsPad=", I8, " must be greater than npts=", I8)
  908 FORMAT("xcloc_fdxcMPI_initialize: No correlation pairs=", I8)
  910 FORMAT('xcloc_fdxcMPI_initialize: Process', I4, ' has ', I6, ' correlations')
  911 FORMAT('xcloc_fdxcMPI_initialize: Have ', I6, ' xcs but need ', I6, ' xcs')
  912 FORMAT('xcloc_fdxcMPI_initialize: Local to global signal map is wrong on rank', I4)
  913 FORMAT('xcloc_fdxcMPI_initialize: bsearch failure on rank:', I4, ' check l2g map')
  914 FORMAT('xcloc_fdxcMPI_initialize: failed to generate xcPairsLocal on rank', I4)
  915 FORMAT('xcloc_fdxcMPI_initialize: fdxc initialization failed on rank', I4)
  920 FORMAT('xcloc_fdxcMPI_initialize: Failed to partition XCs')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases memory on the module.
!>    @ingroup fdxcmpi
      SUBROUTINE xcloc_fdxcMPI_finalize()  &
      BIND(C, NAME='xcloc_fdxcMPI_finalize')
      INTEGER mpierr
      CALL xcloc_fdxc_finalize()
      IF (ALLOCATED(uniqueGlobalSignalIDPtr_)) DEALLOCATE(uniqueGlobalSignalIDPtr_)
      IF (ALLOCATED(uniqueGlobalSignalIDs_)) DEALLOCATE(uniqueGlobalSignalIDs_)
      IF (lfreeComm_) CALL MPI_COMM_FREE(comm_, mpierr)
      IF (ALLOCATED(nXCsPerProcess_))     DEALLOCATE(nXCsPerProcess_)
      IF (ALLOCATED(local2GlobalSignal_)) DEALLOCATE(local2GlobalSignal_)
      lfreeComm_ = .FALSE.
      ldoXC_ = .FALSE.
      nsignalsLocal_ = 0
      nsignalsTotal_ = 0
      npts_ = 0
      nptsPad_ = 0
      nptsInXCs_ = 0
      nXCsLocal_ = 0
      myid_ = MPI_UNDEFINED
      accuracy_ = XCLOC_HIGH_ACCURACY
      precision_ = XCLOC_SINGLE_PRECISION
      verbose_ = XCLOC_PRINT_WARNINGS
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the cross-correlograms.
!>    @param[out] ierr   0 indicates success.
!>    @ingroup fdxcmpi
      SUBROUTINE xcloc_fdxcMPI_computeCrossCorrelograms(ierr) &
      BIND(C, NAME='xcloc_fdxcMPI_computeCrossCorrelograms')
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER ierrLoc, mpierr 
      ierr = 0
      CALL xcloc_fdxc_computeCrossCorrelograms(ierrLoc)
      IF (ierrLoc /= 0) THEN
         WRITE(ERROR_UNIT,900) myid_
         ierr = 1 
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_SUM, comm_, mpierr)
  900 FORMAT('xcloc_fdxcMPI_computeCrossCorrelograms: Error computing xcs on rank', I4) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the phase correlograms.
!>    @param[out] ierr   0 indicates success.
!>    @ingroup fdxcmpi
      SUBROUTINE xcloc_fdxcMPI_computePhaseCorrelograms(ierr) &
      BIND(C, NAME='xcloc_fdxcMPI_computePhaseCorrelograms')
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER ierrLoc, mpierr
      ierr = 0
      CALL xcloc_fdxc_computePhaseCorrelograms(ierrLoc)
      IF (ierrLoc /= 0) THEN
         WRITE(ERROR_UNIT,900) myid_
         ierrLoc = 1
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_SUM, comm_, mpierr)
  900 FORMAT('xcloc_fdxcMPI_computePhaseCorrelograms: Error computing pxcs on rank', I4)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gathers the correlograms onto the root process on the module.
!>    @param[in] ldxcIn   The leading dimension of the correlograms.  This must be at
!>                        least nptsInXCs_.  This is defined on the root process.
!>    @param[in] nxcsIn   Total number of correlograms.  This must equal nXCsTotal_.
!>                        This is defined on the root process.
!>    @param[in] root     The root process ID on which to gather the correlograms.
!>    @param[out] xcs     All of the cross-correlograms.  This is an array of dimension
!>                        [ldxcIn x nxcsIn] with leading dimension ldxcIn.  This need
!>                        only be defined on the root process.
!>    @param[out] ierr    0 indicates succcess.
!>    @ingroup fdxcmpi
      SUBROUTINE xcloc_fdxcMPI_gatherCorrelograms64f(ldxcIn, nxcsIn, root, xcs, ierr) &
      BIND(C, NAME='xcloc_fdxcMPI_gatherCorrelograms64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldxcin, nxcsIn, root
      REAL(C_DOUBLE), INTENT(OUT) :: xcs(ldxcIn*nxcsIn)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: work
      INTEGER irecv, ldxc, mpierr, nxcs, nxcsLoc, sendCount
      INTEGER, ALLOCATABLE :: displs(:), recvCount(:)
      ierr = 0
      IF (myid_ == MPI_UNDEFINED) RETURN
      IF (root < 0 .OR. root > nprocs_ - 1) THEN
         IF (myid_ == root_) WRITE(ERROR_UNIT,895) root
         ierr = 1
      ENDIF
      IF (myid_ == root_) THEN
         IF (ldxcIn < nptsInXCs_) THEN
            WRITE(ERROR_UNIT,900) nptsInXCs_
            ierr = 1
         ENDIF
         IF (nxcsIn < nXCsTotal_) THEN 
            WRITE(ERROR_UNIT,905) nxcsIn
            ierr = 1
         ENDIF
         ldxc = ldxcIn
         nxcs = nxcsIn
      ENDIF
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root, comm_, mpierr)
      IF (ierr /= 0) RETURN
      CALL MPI_Bcast(ldxc, 1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(nxcs, 1, MPI_INTEGER, root, comm_, mpierr)
      ! Set space to gather result
      sendCount = nXCsLocal_*ldxc
      ALLOCATE(work(sendCount))
      ! Have master figure out who is sending what and how much
      IF (myid_ == root) THEN
         ALLOCATE(recvCount(nprocs_)); recvCount(:) = 0 
         ALLOCATE(displs(nprocs_)); displs(:) = 0 
         displs(1) = 0 ! This is C indexed
         DO irecv=1,nprocs_
            nxcsLoc = nXCsPerProcess_(irecv)
            recvCount(irecv) = nxcsLoc*ldxc
            IF (irecv < nprocs_) displs(irecv+1) = displs(irecv) + recvCount(irecv)
         ENDDO
      ENDIF
      ! Get the data
      CALL xcloc_fdxc_getCorrelograms64f(ldxc, nXCsLocal_, work, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,910) myid_
      ENDIF
      ! Gather the data
      CALL MPI_Gatherv(work, sendCount, MPI_DOUBLE_PRECISION,    &
                       xcs, recvCount, displs,                   &
                       MPI_DOUBLE_PRECISION, root, comm_, mpierr)
      IF (ALLOCATED(work))      DEALLOCATE(work)
      IF (ALLOCATED(recvCount)) DEALLOCATE(recvCount) 
      IF (ALLOCATED(displs))    DEALLOCATE(displs)
  895 FORMAT('xcloc_fdxcMPI_gatherCorrelograms64f: Invalid root', I6)
  900 FORMAT('xcloc_fdxcMPI_gatherCorrelograms64f: ldxcIn must be at least', I6)
  905 FORMAT('xcloc_fdxcMPI_gatherCorrelograms64f: nxcsIn must be at least', I6)
  910 FORMAT("xcloc_fdxcMPI_gatherCorrelograms64f: Failed to get data on rank", I4)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the signals to cross-correlate.
!>    @param[in] ldx       Leading dimension of x.  This cannot be less than npts.
!>                         This is defined on the root process.
!>    @param[in] npts      Number of points in the signals.  This must match npts_.
!>                         This is defined on the root process.
!>    @param[in] nsignals  Number of total signals to set.  This must match
!>                         nSignalsTotal_.  This is defined on the root process.
!>    @param[in] root      Root process on that is sending the data.
!>    @param[in] x         The signals to set.  This an [ldx x nsignals] matrix.
!>                         This is defined on the root process.
!>    @param[out] ierr     0 indicates success.
!>    @ingroup fdxcmpi
      SUBROUTINE xcloc_fdxcMPI_setSignals64f(ldx, npts, nsignals, root, x, ierr) &
      BIND(C, NAME='xcloc_fdxcMPI_setSignals64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldx, npts, nsignals, root
      REAL(C_DOUBLE), INTENT(IN) :: x(ldx*nsignals)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: work(:)
      TYPE(MPI_Request), ALLOCATABLE :: send_request(:), recv_request(:)
      INTEGER dest, i, i1, i2, ierrLoc, indx, ip, is, j1, j2, mpierr, nrecv, nsend, source
      LOGICAL flag
      TYPE(MPI_Status) stat
      ierr = 0
      CALL xcloc_fdxc_haveNoSignals() ! Indicate that I have no signals set.
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      IF (myid_ == root) THEN
         IF (ldx < npts .OR. npts /= npts_ .OR. nsignals /= nSignalsTotal_) THEN
            IF (ldx < npts) WRITE(ERROR_UNIT,900) ldx, npts
            IF (npts /= npts_) WRITE(ERROR_UNIT,901) npts_
            IF (nsignals /= nSignalsTotal_) WRITE(ERROR_UNIT,902) nSignalsTotal_
            ierr = 1
         ENDIF
         source = root
      ENDIF
      CALL MPI_Bcast(ierr,   1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(source, 1, MPI_INTEGER, root, comm_, mpierr)
      IF (ierr /= 0) RETURN
      ierrLoc = 0
      ! Have the root send the signals 
      IF (myid_ == root) THEN
         ! Have the root send the input data to the other processes
         nsend = uniqueGlobalSignalIDPtr_(nprocs_+1) &
               - uniqueGlobalSignalIDPtr_(2) + 1
         ALLOCATE(send_request(MAX(1, nsend))); send_request(:) = MPI_REQUEST_NULL 
         is = 1
         DO ip=2,nprocs_
            i1 = uniqueGlobalSignalIDPtr_(ip)
            i2 = uniqueGlobalSignalIDPtr_(ip+1) - 1
            DO i=i1,i2
               dest = uniqueGlobalSignalIDs_(i)
               j1 = (dest - 1)*ldx + 1
               call MPI_Isend(x(j1), npts_, MPI_DOUBLE_PRECISION, ip-1, root, &
                              comm_, send_request(is), mpierr)
               is = is + 1
            ENDDO
         ENDDO 
         IF (is /= nsend) THEN
            WRITE(ERROR_UNIT,910)
            ierrLoc = 1
         ENDIF
         ! Set what I have to set
         i1 = uniqueGlobalSignalIDPtr_(myid_+1)
         i2 = uniqueGlobalSignalIDPtr_(myid_+2) - 1
         DO i=i1,i2
            is = local2GlobalSignal_(i)
            j1 = (is - 1)*ldx + 1
            j2 = j1 + npts_ - 1
            CALL xcloc_fdxc_setSignal64fF(i, npts_, x(j1:j2), ierr)
            IF (ierr /= 0) THEN
               WRITE(ERROR_UNIT,915) i, myid_
               ierrLoc = 1
            ENDIF 
         ENDDO
         DEALLOCATE(send_request)
      ELSE
         ALLOCATE(work(npts_*nSignalsLocal_)); work(:) = 0.d0
         ALLOCATE(recv_request(MAX(1, nSignalsLocal_))); recv_request(:) = MPI_REQUEST_NULL
         ! Get the data from the root
         DO i=1,nsignalsLocal_
            indx = (i - 1)*npts_ + 1
            CALL MPI_Irecv(work(indx), npts_, MPI_DOUBLE_PRECISION, source, &
                           MPI_ANY_TAG, comm_, recv_request(i), mpierr)
         ENDDO
         ! As I receive things put them into the signal buffer
         nrecv = 0
         DO WHILE (nrecv < nSignalsLocal_)
            DO i=1,nSignalsLocal_
               IF (recv_request(i) == MPI_REQUEST_NULL) CYCLE
               CALL MPI_Test(recv_request(i), flag, stat, mpierr) 
               IF (flag) THEN
                  j1 = (i - 1)*npts_ + 1
                  CALL xcloc_fdxc_setSignal64fF(i, npts_, work(j1), ierr)
                  IF (ierr /= 0) THEN
                     WRITE(ERROR_UNIT,915) i, myid_
                     ierrLoc = 1
                  ENDIF
                  !print *, local2GlobalSignal_(i), work(j1)
                  nrecv = nrecv + 1
               ENDIF
            ENDDO
         ENDDO
         ! Free my workspace
         DEALLOCATE(work)
         DEALLOCATE(recv_request)
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_SUM, comm_, mpierr) 
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,920)
         RETURN
      ENDIF
  900 FORMAT('xcloc_fdxc_setSignals64f: Error ldx=', I6, '<', 'npts=', I6)
  901 FORMAT('xcloc_fdxc_setSignals64f: Error expecting npts=', I6)
  902 FORMAT('xcloc_fdxc_setSignals64f: Error expecting nsignals=', I6)
  910 FORMAT('xcloc_fdxc_setSignals64f: Some signals may not be sent')
  915 FORMAT('xcloc_fdxc_setSignals64f: Error setting signal:', I4, ' on process ', I4)
  920 FORMAT('xcloc_fdxc_setSignals64f: Errors detected')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the signals to cross-correlate.
!>    @param[in] ldx       Leading dimension of x.  This cannot be less than npts.
!>                         This is defined on the root process.
!>    @param[in] npts      Number of points in the signals.  This must match npts_.
!>                         This is defined on the root process.
!>    @param[in] nsignals  Number of total signals to set.  This must match
!>                         nSignalsTotal_.  This is defined on the root process.
!>    @param[in] root      Root process on that is sending the data.
!>    @param[in] x         The signals to set.  This an [ldx x nsignals] matrix.
!>                         This is defined on the root process.
!>    @param[out] ierr     0 indicates success.
!>    @ingroup fdxcmpi
      SUBROUTINE xcloc_fdxcMPI_setSignals32f(ldx, npts, nsignals, root, x, ierr) &
      BIND(C, NAME='xcloc_fdxcMPI_setSignals32f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldx, npts, nsignals, root
      REAL(C_FLOAT), INTENT(IN) :: x(ldx*nsignals)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL, ALLOCATABLE :: work(:)
      TYPE(MPI_Request), ALLOCATABLE :: send_request(:), recv_request(:)
      INTEGER dest, i, i1, i2, ierrLoc, indx, ip, is, j1, j2, mpierr, nrecv, nsend, source
      LOGICAL flag
      TYPE(MPI_Status) stat
      ierr = 0
      CALL xcloc_fdxc_haveNoSignals() ! Indicate that I have no signals set.
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      IF (myid_ == root) THEN
         IF (ldx < npts .OR. npts /= npts_ .OR. nsignals /= nSignalsTotal_) THEN
            IF (ldx < npts) WRITE(ERROR_UNIT,900) ldx, npts
            IF (npts /= npts_) WRITE(ERROR_UNIT,901) npts_
            IF (nsignals /= nSignalsTotal_) WRITE(ERROR_UNIT,902) nSignalsTotal_
            ierr = 1
         ENDIF
         source = root
      ENDIF
      CALL MPI_Bcast(ierr,   1, MPI_INTEGER, root, comm_, mpierr)
      CALL MPI_Bcast(source, 1, MPI_INTEGER, root, comm_, mpierr)
      IF (ierr /= 0) RETURN
      ierrLoc = 0
      ! Have the root send the signals 
      IF (myid_ == root) THEN
         ! Have the root send the input data to the other processes
         nsend = uniqueGlobalSignalIDPtr_(nprocs_+1) &
               - uniqueGlobalSignalIDPtr_(2) + 1
         ALLOCATE(send_request(MAX(1, nsend))); send_request(:) = MPI_REQUEST_NULL 
         is = 1
         DO ip=2,nprocs_
            i1 = uniqueGlobalSignalIDPtr_(ip)
            i2 = uniqueGlobalSignalIDPtr_(ip+1) - 1
            DO i=i1,i2
               dest = uniqueGlobalSignalIDs_(i)
               j1 = (dest - 1)*ldx + 1
               call MPI_Isend(x(j1), npts_, MPI_REAL, ip-1, root, &
                              comm_, send_request(is), mpierr)
               is = is + 1
            ENDDO
         ENDDO 
         IF (is /= nsend) THEN
            WRITE(ERROR_UNIT,910)
            ierrLoc = 1
         ENDIF
         ! Set what I have to set
         i1 = uniqueGlobalSignalIDPtr_(myid_+1)
         i2 = uniqueGlobalSignalIDPtr_(myid_+2) - 1
         DO i=i1,i2
            is = local2GlobalSignal_(i)
            j1 = (is - 1)*ldx + 1
            j2 = j1 + npts_ - 1
            CALL xcloc_fdxc_setSignal32fF(i, npts_, x(j1:j2), ierr)
            IF (ierr /= 0) THEN
               WRITE(ERROR_UNIT,915) i, myid_
               ierrLoc = 1
            ENDIF 
         ENDDO
         DEALLOCATE(send_request)
      ELSE
         ALLOCATE(work(npts_*nSignalsLocal_)); work(:) = 0.d0
         ALLOCATE(recv_request(MAX(1, nSignalsLocal_))); recv_request(:) = MPI_REQUEST_NULL
         ! Get the data from the root
         DO i=1,nsignalsLocal_
            indx = (i - 1)*npts_ + 1
            CALL MPI_Irecv(work(indx), npts_, MPI_REAL, source,       &
                           MPI_ANY_TAG, comm_, recv_request(i), mpierr)
         ENDDO
         ! As I receive things put them into the signal buffer
         nrecv = 0
         DO WHILE (nrecv < nSignalsLocal_)
            DO i=1,nSignalsLocal_
               IF (recv_request(i) == MPI_REQUEST_NULL) CYCLE
               CALL MPI_Test(recv_request(i), flag, stat, mpierr) 
               IF (flag) THEN
                  j1 = (i - 1)*npts_ + 1
                  CALL xcloc_fdxc_setSignal32fF(i, npts_, work(j1), ierr)
                  IF (ierr /= 0) THEN
                     WRITE(ERROR_UNIT,915) i, myid_
                     ierrLoc = 1
                  ENDIF
                  !print *, local2GlobalSignal_(i), work(j1)
                  nrecv = nrecv + 1
               ENDIF
            ENDDO
         ENDDO
         ! Free my workspace
         DEALLOCATE(work)
         DEALLOCATE(recv_request)
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_SUM, comm_, mpierr) 
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,920)
         RETURN
      ENDIF
  900 FORMAT('xcloc_fdxc_setSignals64f: Error ldx=', I6, '<', 'npts=', I6)
  901 FORMAT('xcloc_fdxc_setSignals64f: Error expecting npts=', I6)
  902 FORMAT('xcloc_fdxc_setSignals64f: Error expecting nsignals=', I6)
  910 FORMAT('xcloc_fdxc_setSignals64f: Some signals may not be sent')
  915 FORMAT('xcloc_fdxc_setSignals64f: Error setting signal:', I4, ' on process ', I4)
  920 FORMAT('xcloc_fdxc_setSignals64f: Errors detected')
      RETURN
      END

END MODULE
