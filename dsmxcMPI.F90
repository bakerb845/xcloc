!> @brief Computes the diffraction stack migration images in parallel with MPI.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_DSMXC_MPI
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      USE MPI_F08
      USE XCLOC_CONSTANTS 
      USE XCLOC_MEMORY
      USE XCLOC_DSMXC
      USE XCLOC_UTILS
#if defined(__INTEL_COMPILER)
      TYPE(MPI_Comm), PRIVATE, SAVE :: comm_ = MPI_COMM_WORLD
#else
      TYPE(MPI_Comm), PRIVATE, SAVE :: comm_
#endif

      DOUBLE PRECISION, PRIVATE, SAVE :: dt_ = 0.d0
      !> Root process ID.
      INTEGER, PRIVATE, SAVE :: root_ = 0
      !> Number of processes on the communicator.
      INTEGER, PRIVATE, SAVE :: nprocs_ = 0
      !> My rank on the process.
      INTEGER, PRIVATE, SAVE :: myid_ = MPI_UNDEFINED
      !> Total number of grid points in domain.
      INTEGER, PRIVATE, SAVE :: ngrdTotal_ = 0
      !> Number of grid points specific to this process.
      INTEGER, PRIVATE, SAVE :: ngrdLocal_ = 0
      !> Number of grid points belonging to each process.
      INTEGER, PRIVATE, DIMENSION(:), ALLOCATABLE, SAVE :: nGridPtsPerProcess_
      !> If true then free the communicator.
      LOGICAL, PRIVATE, SAVE :: lfreeComm_ = .FALSE.

      PUBLIC :: xcloc_dsmxcMPI_initialize
      PUBLIC :: xcloc_dsmxcMPI_finalize

      CONTAINS
!========================================================================================!
!                                       Begin the Code                                   !
!========================================================================================!
!>    @brief Initializes the MPI-based module to compute the diffraction stack migration
!>           of correlograms. 
!>    @param[in] comm     MPI communicator.  This must be defined on all processes.
!>    @param[in] root     The root process ID on the communicator.  This likely will be 0.
!>    @param[in] ntables  
!>    @param[out] ierr    0 indicates success.
      SUBROUTINE xcloc_dsmxcMPI_initialize(comm, root,                         &
                                           ntables, ngrd, nxcPairs, nptsInXCs, &
                                           dt, xcPairs, ierr)                  &
      BIND(C, NAME='xcloc_dsmxcMPI_initialize')
      TYPE(MPI_Comm),  VALUE, INTENT(IN) :: comm
      !INTEGER(C_INT), VALUE, INTENT(IN) :: fcomm !TODO could be problematic w/ *finter.h
      INTEGER(C_INT), VALUE, INTENT(IN) :: root, ntables, ngrd, nxcPairs, nptsInXCs
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt
      INTEGER(C_INT), INTENT(IN) :: xcPairs(2*nxcPairs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER, ALLOCATABLE :: myGrid(:), myGridPtr(:)
      INTEGER ip, mpierr
      ! iRelease module in case someone called this twice.
      CALL xcloc_dsmxcMPI_finalize()
      ! Get communicator size and my rank
      root_ = root
      CALL MPI_COMM_RANK(comm, myid_, mpierr)
      CALL MPI_COMM_SIZE(comm, nprocs_, mpierr)
      ! Check the inputs
      ierr = 0
      IF (myid_ == root_) THEN
         ! Do my checks
         IF (ngrd < 1 .OR. nxcPairs < 1 .OR. dt <= 0.d0) THEN
            IF (ngrd < 1) WRITE(ERROR_UNIT,900)
            IF (nxcPairs < 1) WRITE(ERROR_UNIT,901)
            IF (dt <= 0.d0) WRITE(ERROR_UNIT,902)
            ierr = 1
         ENDIF
         IF (MINVAL(xcPairs) < 1 .OR. MAXVAL(xcPairs) > ntables) THEN
            IF (MINVAL(xcPairs) < 1) WRITE(ERROR_UNIT,903) 
            IF (MAXVAL(xcPairs) > ntables) WRITE(ERROR_UNIT,904) MAXVAL(xcPairs), ntables
            ierr = 1 
         ENDIF
         IF (nptsInXCs < 1 .OR. MOD(nptsInXCs, 2) /= 1) THEN
            IF (nptsInXCs < 1) WRITE(ERROR_UNIT,905)
            IF (MOD(nptsInXCs, 2) /= 1) WRITE(ERROR_UNIT,906) nptsInXCs
            ierr = 1 
         ENDIF
         IF (ierr /= 0) GOTO 500
         ! Copy some input variables 
         ngrdTotal_ = ngrd
         ntables_ = ntables
         dt_ = dt
         nxcPairs_ = nxcPairs
         ! Decompose the grid
         ALLOCATE(myGrid(ngrdTotal_))
         ALLOCATE(myGridPtr(nprocs_+1))
         CALL xcloc_utils_partitionTasks(ngrdTotal_, nprocs_, myGridPtr, myGrid, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,920)
            GOTO 500
         ENDIF
         ALLOCATE(nGridPtsPerProcess_(nprocs_))
         DO ip=1,nprocs_
            nGridPtsPerProcess_(ip) = myGridPtr(ip+1) - myGridPtr(ip)
         ENDDO
      ENDIF
  500 CONTINUE
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root_, comm, mpierr)
      IF (ierr /= 0) RETURN
      ! Copy the communicator
      CALL MPI_COMM_DUP_WITH_INFO(comm, MPI_INFO_NULL, comm_, mpierr)
      lfreeComm_ = .TRUE.
      ! Send some information to the other processes
      CALL MPI_Bcast(ngrdTotal_, 1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(ntables_,   1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(dt_, 1, MPI_DOUBLE_PRECISION, root_, comm_, mpierr)
      IF (.NOT.ALLOCATED(nGridPtsPerProcess_)) ALLOCATE(nGridPtsPerProcess_(nprocs_))
      CALL MPI_Bcast(nGridPtsPerProcess_, nprocs_, MPI_INTEGER, root_, comm_, mpierr)
      IF (ngrdLocal_ > 0) THEN
!        CALL xcloc_dsmxc_initialize(ntables_, ngrdLocal_, nxcPairs_, nptsInXCs_, &
!                                    dt_, xcPairsWork, ierr)
      ENDIF
      ! Free workspace
      IF (ALLOCATED(myGrid))    DEALLOCATE(myGrid)
      IF (ALLOCATED(myGridPtr)) DEALLOCATE(myGridPtr)
  900 FORMAT('xcloc_dsmxcMPI_initialize: No grid points')
  901 FORMAT('xcloc_dsmxcMPI_initialize: No correlation pairs')
  902 FORMAT('xcloc_dsmxcMPI_initialize: Sampling period must be positive')
  903 FORMAT('xcloc_dsmxcMPI_initialize: All table indices must be positive')
  904 FORMAT('xcloc_dsmxcMPI_initialize: Max table index=', I6, 'exceeds ntables=', I6)
  905 FORMAT('xcloc_dsmxcMPI_initialize: No points in xcs')
  906 FORMAT('xcloc_dsmxcMPI_initialize: Number of points in xcs is even', I6)
  920 FORMAT('xcloc_dsmxcMPI_initialize: Failed to partition grid')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases memory on the module.
!>
      SUBROUTINE xcloc_dsmxcMPI_finalize() &
      BIND(C, NAME='xcloc_dsmxcMPI_finalize')
      CALL xcloc_dsmxc_finalize()
      IF (ALLOCATED(nGridPtsPerProcess_)) DEALLOCATE(nGridPtsPerProcess_)
      IF (lfreeComm_) CALL MPI_COMM_FREE(comm_, mpierr)
      lfreeComm_ = .FALSE.
      ngrdTotal_ = 0
      ngrdLocal_ = 0
      ntables_ = 0
      myid_ = MPI_UNDEFINED
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gathers the diffraction stack migration of correlograms onto the root
!>           process.
!>    @param[in] nwork   Size of image.  This is defined on the root process and must
!>                       be at least ngrdTotal_.
!>    @param[in] root    The process ID onto which the image will be gathered.
!>    @param[out] image  The migrated image.  This is only pertinent to the root process
!>                       and on the root process this must have dimension [nwork].
!>    @param[out] ierr   0 indicates success.
!>
      SUBROUTINE xcloc_dsmxcMPI_gatherImage32f(nwork, root, image, ierr) &
      BIND(C, NAME='xcloc_dsmxcMPI_gatherImage32f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: nwork, root
      REAL(C_FLOAT), INTENT(OUT) :: image(*)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_FLOAT), CONTIGUOUS, DIMENSION(:), POINTER :: imagePtr
      INTEGER, ALLOCATABLE :: displs(:), recvCount(:)
      INTEGER ierrLoc, irecv, mpierr
      ierr = 0
      IF (myid_ == root) THEN
         IF (nwork < ngrdTotal_) THEN
            WRITE(ERROR_UNIT,900) nwork, ngrdTotal_ 
            ierr = 1
         ENDIF
      ENDIF
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root, comm_, mpierr)
      IF (ierr /= 0) RETURN
      ! Get a pointer to the image
      NULLIFY(imagePtr)
      CALL xcloc_dsmxc_getImagePtr(imagePtr, ierrLoc)
      IF (ierrLoc /= 0) THEN
          WRITE(ERROR_UNIT,905) myid_
          ierrLoc = 1
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, comm_, mpierr)
      IF (ierrLoc /= 0) RETURN
      IF (myid_ == root) THEN
         ALLOCATE(recvCount(nprocs_)); recvCount(:) = 0
         ALLOCATE(displs(nprocs_)); displs(:) = 0
         displs(1) = 0 ! This is C indexed
         DO irecv=1,nprocs_
            recvCount(irecv) = nGridPtsPerProcess_(irecv) 
            IF (irecv < nprocs_) displs(irecv+1) = displs(irecv) + recvCount(irecv) 
         ENDDO
      ENDIF
      ! Gather the image onto the root
      CALL MPI_Gatherv(image, ngrdLocal_, MPI_REAL, & 
                       image, recvCount, displs,    &
                       MPI_REAL, root, comm_, mpierr) 
      ! Dereference pointers and free workspace
      NULLIFY(imagePtr)
      IF (ALLOCATED(displs)) DEALLOCATE(displs)
  900 FORMAT('xcloc_dsmxcMPI_getImage32f: nwork=', I6, ' must be at least =', I6)
  905 FORMAT('xcloc_dsmxcMPI_getImage32f: Failed to get local image on rank', I4)
      RETURN
      END
 
END MODULE
