!> @brief Computes the cross-correlograms via the Fourier transform with MPI.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_FDXC_MPI
      USE ISO_C_BINDING
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
      !> RMA window for setting signals.
      TYPE(MPI_Win), PRIVATE, SAVE :: signalRMAWindow_
      !> RMA window for gathering correlograms.
      TYPE(MPI_Win), PRIVATE, SAVE :: xcRMAWindow_
      !> Holds all the input signals.

      !> Memory on the RMA signal window for holding the input signals.
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE, SAVE :: inputSignalsRMA64f_(:)
      !> Memory on the RMA signal window for holding the input signals.
      REAL, ALLOCATABLE, PRIVATE, SAVE :: inputSignalsRMA32f_(:)
      !> Root (master) process ID.
      INTEGER, PRIVATE, SAVE :: root_ = 0
      !> Number of processes on the communicator.
      INTEGER, PRIVATE, SAVE :: nprocs_ = 0
      !> My rank on the process.
      INTEGER, PRIVATE, SAVE :: myid_ = MPI_UNDEFINED
      !> represents the (i,j)'th signal (globally indexed).
      !> This maps from the isLoc'th local signal number to the global signal number.
      !> This has dimension [nsignalsLocal_].
      INTEGER, ALLOCATABLE, PRIVATE, SAVE :: local2GlobalSignal_(:)
      !> This is a [2 x nXCsLocal_] matrix where (2
      INTEGER, ALLOCATABLE, PRIVATE, SAVE :: xcPairsLocal_localNumbering_(:)
      INTEGER, ALLOCATABLE, PRIVATE, SAVE :: xcPairsLocalNumbering_(:)
      !> Maps from
      INTEGER, ALLOCATABLE, PRIVATE, SAVE :: myXCs_(:)
      INTEGER, ALLOCATABLE, PRIVATE, SAVE :: myXCPtr_(:)
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
      PUBLIC :: xcloc_fdxcMPI_initialize
      CONTAINS
!========================================================================================!
!                                      Begin the Code                                    !
!========================================================================================!
!>    @brief Initializes the parallel frequency domain cross-correlation calculator.
!>    @param[in] comm      MPI communicator.
!>    @param[in] master    ID of master process.  This will likely be 0.
!>    @param[in] npts      Number of points in each input signal.  This is defined on
!>                         the master process.
!>    @param[in] nptsPad   A tuning parameter to mitigate DFT lengths that could
!>                         potentially be large semi-prime numbers.  This is defined
!>                         on the master process.
!>    @param[in] nxcs      Number of cross-correlations.  This is defined on the
!>                         master process.
!>    @param[in] xcPairs   This is a [2 x nxcs] matrix in column major format where
!>                         the indices, (2*(ixc-1)+1, 2*(ixc-1)+2), map to the
!>                         (i,j)'th signal pair comprising a correlation.
!>                         This is defined on the master process.
!>    @param[in] verbose   Controls the verbosity of the module.  This is defined
!>                         on the master process.
!>    @param[in] prec      The precision of hte module.  This is defined on the
!>                         master process.
!>    @param[in] accuracy  Controls the accuracy of the vector calculations in MKL.
!>                         This is defined on the master process.
!>    @param[out] ierr     0 indicates sucess.
!>
      SUBROUTINE xcloc_fdxcMPI_initialize(comm, master,                   &
                                          npts, nptsPad,                  &
                                          nxcs, xcPairs,                  &
                                          verbose, prec, accuracy, ierr)  &
      BIND(C, NAME='xcloc_fdxcMPI_initialize')
      TYPE(MPI_Comm), VALUE, INTENT(IN) :: comm
      INTEGER(C_INT), VALUE, INTENT(IN) :: master, npts, nptsPad, nxcs, &
                                           verbose, prec, accuracy
      INTEGER(C_INT), INTENT(IN) :: xcPairs(*)
      INTEGER ierr, indx, mpierr, myid, nsignals
      INTEGER(KIND=8) nwork
      ! Copy the communicator
      CALL MPI_COMM_DUP_WITH_INFO(comm, MPI_INFO_NULL, comm_)
      lfreeComm_ = .TRUE.
      ! Get information 
      CALL MPI_COMM_RANK(comm_, myid_, mpierr)
      CALL MPI_COMM_SIZE(comm_, nprocs_, mpierr)
      ! Some basic checks by master process
      ierr = 1
      IF (myid == master) THEN
         ierr = 0
         root_ = master
         IF (npts < 1 .OR. nsignals < 2 .OR. nptsPad < npts .OR. nxcs < 1) THEN
            IF (npts < 1) WRITE(*,905) npts 
            IF (nptsPad < npts) WRITE(*,907) nptsPad, npts
            IF (nxcs < 1) WRITE(*,908) nxcs
            ierr = 1
         ENDIF
         nsignals = MAXVAL(xcPairs(1:2*nxcs))
         IF (nsignals < 2) THEN
            WRITE(*,906) nsignals
            ierr = 1
         ENDIF
         IF (.NOT. xcloc_constants_isValidPrecision(prec)) ierr = 1 
         IF (.NOT. xcloc_constants_isValidAccuracy(accuracy)) ierr = 1
         IF (ierr == 0) THEN
            npts_     = npts
            nptsPad_  = nptsPad
            nSignalsTotal_ = nsignals
            verbose_ = verbose
            nptsInXCs_ = 2*nptsPad_ - 1   ! Length of the cross-correlations
            nXCsTotal_ = nxcs
            accuracy_ = accuracy
            precision_ = prec
            ! build up the workspaces
            !lds_ = xcloc_memory_padLength(4096, SIZEOF(1.d0), nptsPad_) 
            !nwork = lds_*
         ENDIF
      ENDIF
      CALL MPI_BCAST(ierr, 1, MPI_INTEGER, master, comm_, mpierr)
      IF (ierr /= 0) RETURN
      ! Set some basic information
      CALL MPI_BCAST(root_,          1, MPI_INTEGER, master, comm_, mpierr)
      CALL MPI_BCAST(npts_,          1, MPI_INTEGER, master, comm_, mpierr)
      CALL MPI_BCAST(nptsPad_,       1, MPI_INTEGER, master, comm_, mpierr)
      CALL MPI_BCAST(nSignalsTotal_, 1, MPI_INTEGER, master, comm_, mpierr)
      CALL MPI_BCAST(verbose_,       1, MPI_INTEGER, master, comm_, mpierr)
      CALL MPI_BCAST(accuracy_,      1, MPI_INTEGER, master, comm_, mpierr)
      CALL MPI_BCAST(nptsInXCs_,     1, MPI_INTEGER, master, comm_, mpierr)
      CALL MPI_BCAST(nXCsTotal_,     1, MPI_INTEGER, master, comm_, mpierr)
      ! Load balance
      ALLOCATE(myXCs_(nXCsTotal_))
      ALLOCATE(myXCPtr_(nprocs_+1))
      CALL xcloc_utils_partitionTasks(nXCsTotal_, nprocs_, myXCs_, myXCPtr_, ierr)
      ! Figure out the number of cross-correlations that I need to perform
      indx = myXCs_(myid+1)
      nXCsLocal_ = myXCPtr_(indx+1) - myXCPtr_(indx)
!     ALLOCATE(xcPairsLocal_(MAX(2*nLocalXCs_, 1))); xcPairsLoc_(:) = 0
      IF (nXCsLocal_ > 0) THEN
!        CALL xcloc_fdxc_initialize(npts_,      nptsPad_,    &
!                                   nXCsLocal_, xcPairsLocal_, &
!                                   verbose_, prec_, accuracy_, ierr) 
      ENDIF
      ! Set space
      !IF (precision_ == 
      IF (myid_ == root_) THEN

      ENDIF
      !CALL MPI_Win_create( )
      ! Format statements
  905 FORMAT("xcloc_fdxcMPI_initialize: npts=", I8, "must be positive")
  906 FORMAT("xcloc_fdxcMPI_initialize: nsignals=", I8, "must be at least 2")
  907 FORMAT("xcloc_fdxcMPI_initialize: nptsPad=", I8, "must be greater than npts=", I8)
  908 FORMAT("xcloc_fdxc_initialize: No correlation pairs=", I8)
      RETURN
      END

      SUBROUTINE xcloc_fdxcMPI_finalize() &
      BIND(C, NAME='xcloc_fdxcMPI_finalize')
      IF (lfreeComm_) CALL MPI_COMM_FREE(comm_)
      lfreeComm_ = .FALSE.
      myid_ = MPI_UNDEFINED
      accuracy_ = XCLOC_HIGH_ACCURACY
      precision_ = XCLOC_SINGLE_PRECISION
      verbose_ = XCLOC_PRINT_WARNINGS
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
!>    @param[in] x         The signals to set.  This an [ldx x nsignals] matrix.
!>                         This is defined on the root process.
!>    @param[out] ierr     0 indicates success.
      SUBROUTINE xcloc_fdxcMPI_setSignals64f(ldx, npts, nsignals, x, ierr) &
      BIND(C, NAME='xcloc_fdxcMPI_setSignals64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldx, npts, nsignals
      REAL(C_DOUBLE), INTENT(IN) :: x(*)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: xwork(:)
      INTEGER(KIND=MPI_ADDRESS_KIND) disp
      INTEGER i1, i2, ierrLoc, is, isLoc, j1, j2, mpierr, rank
      ierr = 0
      CALL xcloc_fdxc_haveNoSignals() ! Indicate that I have no signals set.
      IF (myid_ == MPI_UNDEFINED) RETURN ! I don't belong here
      CALL MPI_Comm_rank(comm_, rank, mpierr)
      IF (rank == root_) THEN
         IF (ldx < npts .OR. npts /= npts_ .OR. nsignals /= nSignalsTotal_) THEN
            IF (ldx < npts) WRITE(*,900) ldx, npts
            IF (npts /= npts_) WRITE(*,901) npts_
            IF (nsignals /= nSignalsTotal_) WRITE(*,902) nSignalsTotal_
            ierr = 1
         ENDIF
         ! Have the master copy the input data
         IF (ierr == 0) THEN
            DO is=1,nSignalsTotal_
               i1 = (is - 1)*ldx + 1
               i2 = i1 + npts - 1
               j1 = (is - 1)*lds_ + 1
               j2 = j1 + npts - 1
               inputSignalsRMA64f_(j1:j2) = x(i1:i2) 
            ENDDO
         ENDIF
      ENDIF
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root_, comm_, mpierr)
      IF (ierr /= 0) RETURN
      CALL MPI_Barrier(comm_, mpierr)
      ! Have the processes fetch the signals that they require
      CALL MPI_Win_fence(MPI_MODE_NOPRECEDE, signalRMAWindow_)
      ! Get the signals that I need and set them
      ALLOCATE(xwork(npts_))
      DO isLoc=1,nSignalsLocal_
         is = local2GlobalSignal_(isLoc) ! Get global signal number
         disp = (is - 1)*lds_ + 1
         IF (myid_ == root_ ) THEN
            j1 = (is - 1)*ldx + 1
            CALL xcloc_fdxc_setSignal64fF(isLoc, npts, x(j1), ierrLoc)
         ELSE
            CALL MPI_Get(xwork, npts_, MPI_DOUBLE_PRECISION, &
                         root_, disp,                        &
                         npts, MPI_DOUBLE_PRECISION,         &
                         signalRMAWindow_, mpierr)
            CALL xcloc_fdxc_setSignal64fF(isLoc, npts, xwork, ierrLoc)
         ENDIF
         IF (ierrLoc /= 0) THEN
            WRITE(*,915) myid_, is
            ierr = ierr + 1
         ENDIF
      ENDDO
      IF (ALLOCATED(xwork)) DEALLOCATE(xwork)
      ! Ensure all RMA operations complete before proceeding.
      CALL MPI_Win_fence(MPI_MODE_NOSTORE + MPI_MODE_NOPUT + MPI_MODE_NOSUCCEED, &
                         signalRMAWindow_, mpierr)
  900 FORMAT('xcloc_fdxc_setSignals64f: Error ldx=', I6, '<', 'npts=', I6)
  901 FORMAT('xcloc_fdxc_setSignals64f: Error expecting npts=', I6)
  902 FORMAT('xcloc_fdxc_setSignals64f: Error expecting nsignals=', I6)
  910 FORMAT('xcloc_fdxc_setSignals64f: Error setting signal index', I4)
  915 FORMAT('xcloc_fdxc_setSignals64f: Error setting signal:', I4, ' on process ', I4)
      RETURN
      END

END MODULE
