!> @brief Computes the cross-correlograms via the Fourier transform with MPI.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_FDXC_MPI
      USE ISO_C_BINDING
      USE MPI_F08
      !USE MPI
      USE XCLOC_CONSTANTS
      USE XCLOC_MEMORY
      USE XCLOC_FDXC
      USE XCLOC_UTILS
#ifdef _OPENMP
      USE OMP_LIB
#endif
      !> Maps from
      INTEGER, ALLOCATABLE, PRIVATE, SAVE :: myXCs_(:)
      INTEGER, ALLOCATABLE, PRIVATE, SAVE :: myXCPtr_(:)
      !> Length of input signals.
      INTEGER, PRIVATE, SAVE :: npts_ = 0
      !> The length of the signals to transform.  This can mitigate the pathologic
      !> case where the signal transform lengths are large (semi)prime numbers which
      !> make the DFT very expensive.
      INTEGER, PRIVATE, SAVE :: nptsPad_ = 0
      !> Length of the time domain cross-correloagrams.
      INTEGER, PRIVATE, SAVE :: nptsInXCs_ = 0
      !> Total number of cross-correlations
      INTEGER, PRIVATE, SAVE :: nTotalXCs_ = 0
      !> Controls verbosity.
      INTEGER, PRIVATE, SAVE :: verbose_ = XCLOC_PRINT_WARNINGS
      !> If true then this process will be computing cross-correlations. 
      LOGICAL, PRIVATE, SAVE :: ldoXC_ = .FALSE.
      
      PUBLIC :: xcloc_fdxcMPI_initialize
      CONTAINS
!========================================================================================!
!                                      Begin the Code                                    !
!========================================================================================!
!>    @brief Initializes the Foureir transform
!>    @param[in] comm      MPI communicator.
!>    @param[in] master    ID of master process.  This should probably be 0.
!>    @param[in] npts      Number of points in each input signal.
!>    @param[in] nptsPad   Number of  
      SUBROUTINE xcloc_fdxcMPI_initialize(comm, master,                   &
                                          npts, nptsPad,                  &
                                          nxcs, xcPairs,                  &
                                          verbose, prec, accuracy, ierr)  &
      BIND(C, NAME='xcloc_fdxcMPI_initialize')
      TYPE(MPI_Comm), VALUE, INTENT(IN) :: comm
      INTEGER(C_INT), VALUE, INTENT(IN) :: master, npts, nptsPad, nxcs, &
                                           verbose, prec, accuracy
      INTEGER(C_INT), INTENT(IN) :: xcPairs(*)
      INTEGER ierr, mpierr, myid, nsignals
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      CALL MPI_COMM_SIZE(comm, nprocs, mpierr)
      ! Some basic checks by master process
      ierr = 1
      IF (myid == master) THEN
         ierr = 0
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
            RETURN
         ENDIF
         IF (.NOT. xcloc_constants_isValidPrecision(prec)) ierr = 1 
         IF (.NOT. xcloc_constants_isValidAccuracy(accuracy)) ierr = 1
         IF (ierr == 0) THEN
            npts_     = npts
            nptsPad_  = nptsPad
            nTotalSignals_ = nsignals
            verbose_ = verbose
            nptsInXCs_ = 2*nptsPad_ - 1   ! Length of the cross-correlations
            nTotalXCs_ = nxcs
         ENDIF
      ENDIF
      CALL MPI_BCAST(ierr, 1, MPI_INTEGER, master, comm, mpierr)
      IF (ierr /= 0) RETURN
      ! Set some basic information
      CALL MPI_BCAST(npts_,          1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(nptsPad_,       1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(nTotalSignals_, 1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(verbose_,       1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(nptsInXCs_,     1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(nTotalXCs_,     1, MPI_INTEGER, master, comm, mpierr)
      ! Load balance
      ALLOCATE(myXCs_(nTotalXCs_))
      ALLOCATE(myXCPtr_(nprocs+1))
      CALL xcloc_utils_partitionTasks(nTotalXCs_, nprocs, myXCs_, myXCPtr_, ierr)
      ! Initialize

      ! Format statements
  905 FORMAT("xcloc_fdxcMPI_initialize: npts=", I8, "must be positive")
  906 FORMAT("xcloc_fdxcMPI_initialize: nsignals=", I8, "must be at least 2")
  907 FORMAT("xcloc_fdxcMPI_initialize: nptsPad=", I8, "must be greater than npts=", I8)
  908 FORMAT("xcloc_fdxc_initialize: No correlation pairs=", I8)
      RETURN
      END


END MODULE
