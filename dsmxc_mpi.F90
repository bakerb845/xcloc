!> @defgroup dsmxcMPI Parallel Diffraction Stack Migration of Correlograms
!> @ingroup xcloc_mpi
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
      IMPLICIT NONE
#if defined(__INTEL_COMPILER)
      !> @ingroup dsmxcMPI
      !> The MPI communicator for the module.
      TYPE(MPI_Comm), PRIVATE, SAVE :: comm_ = MPI_COMM_NULL !MPI_COMM_WORLD
#else
      TYPE(MPI_Comm), PRIVATE, SAVE :: comm_
#endif
      !> @ingroup dsmxcMPI
      !> The sampling period of the correlograms.
      DOUBLE PRECISION, PRIVATE, SAVE :: dt_ = 0.d0
      !> @ingroup dsmxcMPI
      !> Root process ID.
      INTEGER, PRIVATE, SAVE :: root_ = 0
      !> @ingroup dsmxcMPI
      !> Number of processes on the communicator.
      INTEGER, PRIVATE, SAVE :: nprocs_ = 0
      !> @ingroup dsmxcMPI
      !> My rank on the process.
      INTEGER, PRIVATE, SAVE :: myid_ = MPI_UNDEFINED
      !> @ingroup dsmxcMPI
      !> Total number of grid points in domain.
      INTEGER, PRIVATE, SAVE :: ngrdTotal_ = 0
      !> @ingroup dsmxcMPI
      !> Number of grid points specific to this process.
      INTEGER, PRIVATE, SAVE :: ngrdLocal_ = 0
      !> @ingroup dsmxcMPI
      !> Number of points in correlograms.
      INTEGER, PRIVATE, SAVE :: nptsInXCs_ = 0
      !> @ingroup dsmxcMPI
      !> The number of tables.
      INTEGER, PRIVATE, SAVE :: ntables_ = 0
      !> @ingroup dsmxcMPI
      !> The number of cross-correlatin pairs.
      INTEGER, PRIVATE, SAVE :: nxcPairs_ = 0
      !> @ingroup dsmxcMPI
      !> Verbosity on the module.
      INTEGER, PRIVATE, SAVE :: verbose_ = XCLOC_PRINT_WARNINGS
      !> @ingroup dsmxcMPI
      !> Number of grid points belonging to each process.
      INTEGER, PRIVATE, DIMENSION(:), ALLOCATABLE, SAVE :: nGridPtsPerProcess_
      !> @ingroup dsmxcMPI
      !> If true then all tables are set.
      LOGICAL, PRIVATE, SAVE :: lhaveAllTables_ = .FALSE.
      !> @ingroup dsmxcMPI
      !> If true then free the communicator.
      LOGICAL, PRIVATE, SAVE :: lfreeComm_ = .FALSE.

      PUBLIC :: xcloc_dsmxcMPI_initialize
      PUBLIC :: xcloc_dsmxcMPI_finalize
      PUBLIC :: xcloc_dsmxcMPI_setTable64f
      PUBLIC :: xcloc_dsmxcMPI_getNumberOGridPointsInTable
      PUBLIC :: xcloc_dsmxcMPI_getNumberOfTables
      PUBLIC :: xcloc_dsmxcMPI_haveAllTables

      CONTAINS
!========================================================================================!
!                                       Begin the Code                                   !
!========================================================================================!
!>    @brief Initializes the MPI-based module to compute the diffraction stack migration
!>           of correlograms. 
!>    @param[in] fcomm      MPI communicator.  This must be defined on all processes.
!>    @param[in] root       The root process ID on the communicator.  This likely will
!>                          be 0.
!>    @param[in] ngrd       Number of points in migration grid.  This is defined on the
!>                          root process.
!>    @param[in] nxcPairs   The nmber of cross-correlation pairs.  This is defined on the
!>                          root process.
!>    @param[in] nptsInXCs  The number of points in the correlograms.  This must be
!>                          positive and an odd number.  This is defined on the root
!>                          process.
!>    @param[in] dt         The sampling period of the correlograms in seconds.  This
!>                          must be defined on the root process.
!>    @param[in] xcPairs    This is a [2 x ncxPairs] matrix stored in column major format
!>                          that defines the signal indices comprising a correlation
!>                          pair.  This is defined on the root process. 
!>    @param[in] verbose    Controls the verbosity. 
!>    @param[out] ierr      0 indicates success.
!>    @ingroup dsmxcMPI
!>    @bug Need to verify that MPI_Comm is the same size as MPI_Fcomm.
      SUBROUTINE xcloc_dsmxcMPI_initialize(fcomm, root,                &
                                           ngrd, nxcPairs, nptsInXCs,  &
                                           dt, xcPairs, verbose, ierr) &
      BIND(C, NAME='xcloc_dsmxcMPI_initialize')
      !TYPE(MPI_Comm),  VALUE, INTENT(IN) :: comm
      INTEGER(C_INT64_T), VALUE, INTENT(IN) :: fcomm
      !INTEGER(C_INT), VALUE, INTENT(IN) :: fcomm !TODO could be problematic w/ *finter.h
      INTEGER(C_INT), VALUE, INTENT(IN) :: root, ngrd, nxcPairs, nptsInXCs, &
                                           verbose
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt
      INTEGER(C_INT), INTENT(IN) :: xcPairs(2*nxcPairs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER, ALLOCATABLE :: uniqueSignals(:), myGrid(:), myGridPtr(:), xcPairsWork(:)
      INTEGER ierrLoc, ip, mpierr
      TYPE(MPI_Comm) comm
      comm%MPI_VAL = INT(fcomm, KIND(comm%MPI_VAL))
      ! iRelease module in case someone called this twice.
      CALL xcloc_dsmxcMPI_finalize()
      ! Get communicator size and my rank
      root_ = root
      CALL MPI_Comm_rank(comm, myid_, mpierr)
      CALL MPI_Comm_size(comm, nprocs_, mpierr)
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
         IF (nptsInXCs < 1 .OR. MOD(nptsInXCs, 2) /= 1) THEN
            IF (nptsInXCs < 1) WRITE(ERROR_UNIT,905)
            IF (MOD(nptsInXCs, 2) /= 1) WRITE(ERROR_UNIT,906) nptsInXCs
            ierr = 1
         ENDIF
         IF (ierr /= 0) GOTO 500
         CALL xcloc_utils_unique32s(2*nxcPairs, xcPairs, ntables_, uniqueSignals, ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,907)
            GOTO 500
         ENDIF 
         IF (ALLOCATED(uniqueSignals)) DEALLOCATE(uniqueSignals)
         ! Copy some input variables 
         ngrdTotal_ = ngrd
         dt_ = dt
         nptsInXCs_ = nptsInXCs
         nxcPairs_ = nxcPairs
         verbose_ = verbose
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
            IF (verbose_ > XCLOC_PRINT_WARNINGS) THEN
               WRITE(OUTPUT_UNIT,910) ip-1, nGridPtsPerProcess_(ip)
            ENDIF
         ENDDO
         ALLOCATE(xcPairsWork(2*nxcPairs_))
         xcPairsWork(1:2*nxcPairs_) = xcPairs(1:2*nxcPairs_) 
      ENDIF
  500 CONTINUE
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root_, comm, mpierr)
      IF (ierr /= 0) RETURN
      ! Copy the communicator
      CALL MPI_Comm_dup_with_info(comm, MPI_INFO_NULL, comm_, mpierr)
      lfreeComm_ = .TRUE.
      ! Send some information to the other processes
      CALL MPI_Bcast(ngrdTotal_, 1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(nxcPairs_,  1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(nptsInXCs_, 1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(verbose_,   1, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(dt_, 1, MPI_DOUBLE_PRECISION, root_, comm_, mpierr)
      IF (.NOT.ALLOCATED(nGridPtsPerProcess_)) ALLOCATE(nGridPtsPerProcess_(nprocs_))
      IF (.NOT.ALLOCATED(xcPairsWork)) ALLOCATE(xcPairsWork(2*nxcPairs_))
      CALL MPI_Bcast(nGridPtsPerProcess_, nprocs_, MPI_INTEGER, root_, comm_, mpierr)
      CALL MPI_Bcast(xcPairsWork, 2*nxcPairs_, MPI_INTEGER, root_, comm_, mpierr)
      ngrdLocal_ = nGridPtsPerProcess_(myid_+1)
      ierrLoc = 0
      IF (ngrdLocal_ > 0) THEN
         CALL xcloc_dsmxc_initialize(ngrdLocal_, nxcPairs_, nptsInXCs_, &
                                     dt_, xcPairsWork, verbose_, ierrLoc)
         IF (ierrLoc /= 0) THEN
            WRITE(ERROR_UNIT,930) myid_
            ierrLoc = 1
         ENDIF
         CALL xcloc_dsmxc_getNumberOfTables(ntables_)
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, comm_, mpierr) 
      ! Free workspace
      IF (ALLOCATED(myGrid))      DEALLOCATE(myGrid)
      IF (ALLOCATED(myGridPtr))   DEALLOCATE(myGridPtr)
      IF (ALLOCATED(xcPairsWork)) DEALLOCATE(xcPairsWork)
  900 FORMAT('xcloc_dsmxcMPI_initialize: No grid points')
  901 FORMAT('xcloc_dsmxcMPI_initialize: No correlation pairs')
  902 FORMAT('xcloc_dsmxcMPI_initialize: Sampling period must be positive')
  905 FORMAT('xcloc_dsmxcMPI_initialize: No points in xcs')
  906 FORMAT('xcloc_dsmxcMPI_initialize: Number of points in xcs is even ', I0)
  907 FORMAT('xcloc_dsmxcMPI_initialize: Failed to make unique table list')
  910 FORMAT('xcloc_dsmxcMPI_initialize: Process ', I0, ' has ', I0, ' grid points')
  920 FORMAT('xcloc_dsmxcMPI_initialize: Failed to partition grid')
  930 FORMAT('xcloc_dsmxcMPI_initialize: Failed to initialize dsmxc on rank', I0)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Maps a signal number to the table index.
!>    @param[in] is     Signal number.  This must be in the xcPairs table.
!>    @param[out] it    The table index corresponding to is.
!>    @param[out] ierr  0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxcMPI_signalToTableIndex(is, it, ierr) &
      BIND(C, NAME='xcloc_dsmxcMPI_signalToTableIndex')
      INTEGER, INTENT(IN) :: is
      INTEGER, INTENT(OUT) :: it, ierr
      CALL xcloc_dsmxc_signalToTableIndex(is, it, ierr)
      IF (ierr /= 0) WRITE(ERROR_UNIT,900)
  900 FORMAT('xcloc_dsmxcMPI_signalToTableIndex: Failed to find signal number', I0, &
             ' in table')
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the travel time table on the module.
!>    @param[in] tableNumber  Index of the table to set.  This must be defined on the
!>                            root process. 
!>    @param[in] ngrd         Number of grid points in table.  This must match
!>                            ngrdTotal_.   This must be defined on the root process.
!>    @param[in] root         The root process ID on the DSM XC communicator.
!>    @param[in] table        Travel time (seconds) from the receiver to all points.
!>                            This is an array of dimension [ngrd] and must be defined
!>                            on the root process.
!>    @param[out] ierr        0 indicates success.
!>    @ingroup dsmxcMPI
      SUBROUTINE xcloc_dsmxcMPI_setTable64f(tableNumber, ngrd, root, table, ierr) &
      BIND(C, NAME='xcloc_dsmxcMPI_setTable64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: tableNumber, ngrd, root
      REAL(C_DOUBLE), DIMENSION(ngrd), INTENT(IN) :: table
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: twork(:)
      INTEGER, ALLOCATABLE :: displs(:), sendCounts(:)
      INTEGER i, ibool, ierrLoc, itemp, mpierr, recvCount, tnumber
      LOGICAL(C_BOOL) :: ltemp
      ierr = 0
      IF (myid_ == root) THEN
         tnumber = tableNumber
         IF (tnumber < 1 .OR. tnumber > ntables_) THEN
            WRITE(ERROR_UNIT,900) tnumber, ntables_
            ierr = 1
         ENDIF
         IF (ngrd /= ngrdTotal_) THEN
            WRITE(ERROR_UNIT,905) ngrd, ngrdTotal_
            ierr = 1 
         ENDIF
         IF (ierr /= 0) GOTO 500
         ALLOCATE(displs(nprocs_)); displs(:) = 0
         ALLOCATE(sendCounts(nprocs_)); sendCounts(:) = nGridPtsPerProcess_(:)
         DO i=1,nprocs_-1
            displs(i+1) = displs(i) + nGridPtsPerProcess_(i)
         ENDDO
      ENDIF
  500 CONTINUE
      CALL MPI_Bcast(ierr, 1, MPI_INTEGER, root, comm_, mpierr)
      IF (ierr /= 0) RETURN
      ! Distribute the table
      CALL MPI_Bcast(tnumber, 1, MPI_INTEGER, root, comm_, mpierr)
      recvCount = ngrdLocal_
      ALLOCATE(twork(ngrdLocal_))
      CALL MPI_Scatterv(table, sendCounts, displs,                &
                        MPI_DOUBLE_PRECISION, twork, recvcount,   &
                        MPI_DOUBLE_PRECISION, root, comm_, mpierr)
      CALL xcloc_dsmxc_setTable64f(tnumber, ngrdLocal_, twork, ierrLoc)
      IF (ALLOCATED(displs))     DEALLOCATE(displs)
      IF (ALLOCATED(sendCounts)) DEALLOCATE(sendCounts)
      IF (ALLOCATED(twork))      DEALLOCATE(twork)
      IF (ierrLoc /= 0) THEN
         WRITE(ERROR_UNIT,910) myid_
         ierrLoc = 1
      ENDIF
      CALL MPI_Allreduce(ierrLoc, ierr, 1, MPI_INTEGER, MPI_MAX, comm_, mpierr)
      IF (ierr /= 0) RETURN
      CALL xcloc_dsmxc_haveAllTables(ltemp) 
      itemp = 0
      IF (ltemp) itemp = 1
      CALL MPI_Allreduce(itemp, ibool, 1, MPI_INTEGER, MPI_MIN, comm_, mpierr)
      IF (ibool == 1) lhaveAllTables_ = .TRUE.
  900 FORMAT('xcloc_dsmxcMPI_setTable64f: tableNumber = ', I0, &
             ' must be in range [1,',I0,']')
  905 FORMAT('xcloc_dsmxcMPI_setTable64f: ngrd = ', I0, ' expecting ngrdTotal_ =', I0)
  910 FORMAT('xcloc_dsmxcMPI_setTable64f: Failed to set table on rank ', I0)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Convenience utility to determine if all the tables have been set.
!>    @param[out] lhaveAllTables  If true then all tables have been set.
!>    @ingroup dsmxcMPI
      SUBROUTINE xcloc_dsmxcMPI_haveAllTables(lhaveAllTables) &
      BIND(C, NAME='xcloc_dsmxcMPI_haveAllTables')
      LOGICAL(C_BOOL), INTENT(OUT) :: lhaveAllTables
      lhaveAllTables = lhaveAllTables_
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the total number of grid points in a table.
!>    @param[out] ngrd   Number of grid points in each table.
!>    @ingroup dsmxcMPI
      SUBROUTINE xcloc_dsmxcMPI_getNumberOGridPointsInTable(ngrd) &
      BIND(C, NAME='xcloc_dsmxcMPI_getNumberOGridPointsInTable')
      INTEGER(C_INT), INTENT(OUT) :: ngrd
      ngrd = ngrdTotal_
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the number of travel time tables on the module.
!>    @param[out] ntables  The number of travel time tables.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxcMPI_getNumberOfTables(ntables) &
      BIND(C, NAME='xcloc_dsmxcMPI_getNumberOfTables')
      INTEGER(C_INT), INTENT(OUT) :: ntables
      ntables = ntables_
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases memory on the module.
!>    @ingroup dsmxcMPI
      SUBROUTINE xcloc_dsmxcMPI_finalize() &
      BIND(C, NAME='xcloc_dsmxcMPI_finalize')
      INTEGER mpierr
      CALL xcloc_dsmxc_finalize()
      IF (ALLOCATED(nGridPtsPerProcess_)) DEALLOCATE(nGridPtsPerProcess_)
      IF (lfreeComm_) CALL MPI_Comm_free(comm_, mpierr)
      lfreeComm_ = .FALSE.
      lhaveAllTables_ = .FALSE.
      ngrdTotal_ = 0
      ngrdLocal_ = 0
      nptsInXCs_ = 0
      ntables_ = 0
      dt_ = 0.d0
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
!>    @ingroup dsmxcMPI
      SUBROUTINE xcloc_dsmxcMPI_gatherImage32f(nwork, root, image, ierr) &
      BIND(C, NAME='xcloc_dsmxcMPI_gatherImage32f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: nwork, root
      REAL(C_FLOAT), DIMENSION(:), INTENT(OUT) :: image
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
