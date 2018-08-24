!> @defgroup utilsMPI MPI Utilities
!> @ingroup xcloc_mpi
!> @brief Utilities for easing the use of MPI.
!> @copyright Ben Baker distributed under the MIT license. 
!> @date March 2015
MODULE XCLOC_UTILS_MPI
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
      USE MPI_F08
      IMPLICIT NONE

      PUBLIC :: xcloc_utilsMPI_splitComm
      CONTAINS
!=======================================================================================!
!                                       Begin the Code                                  !
!=======================================================================================!
!>    @brief Splits an MPI communicator with nprocs processes.  As an example the  
!>           module will split 12 processes into 4 rows and 3 columns.
!>           On the original communicator the processes have IDs: \n
!>              myid:  0  1  2  3  4  5  6  7  8  9 10 11 \n
!>           This is rearranged into the column-major square: 
!>              0  4   8
!>              1  5   9
!>              2  6  10
!>              3  7  11
!>           Each column represents an intra communicator so that the intracommunicator
!>           has ranks given by
!>              0  0  0
!>              1  1  1
!>              2  2  2 
!>              3  3  3
!>           Communication between the groups happens on intercommunicators whose
!>           ranks are given by
!>              0  1  2
!>              0  1  2
!>              0  1  2
!>              0  1  2 
!>           
!>           Note this is an efficacious strategy for nodes in that a common MPI default
!>           is to assign process IDs sequentially node by node, hence, if making 
!>           subgroups on a node this division strategy can decrease the likelihood that
!>           point to point will be performed on a node.
!>
!>    @param[in] mpiComm        This is the communicator to split.
!>    @param[in] ncols          Number of columns in which to split mpiComm.
!>
!>    @param[out] intraCommID   Process ID on intraComm.
!>    @param[out] intraComm     MPI communicator for communication along a column.  These
!>                              processes will likely physically be closer to one another.
!>    @param[out] interCommID   Process ID on interComm.
!>    @param[out] interComm     MPI communicator for communication along a row.  These
!>                              processes will likely be split across compute nodes.
!>    @ingroup utilsMPI
      SUBROUTINE xcloc_utilsMPI_splitComm(mpiComm, ncols,         &
                                          intraCommID, intraComm, &
                                          interCommID, interComm, &
                                          ierr)                   &
      BIND(C, NAME='xcloc_utilsMPI_splitComm')
      IMPLICIT NONE
      TYPE(MPI_Comm), INTENT(IN) :: mpiComm 
      INTEGER, INTENT(IN) :: ncols
      TYPE(MPI_Comm), INTENT(OUT) :: intraComm, interComm
      INTEGER, INTENT(OUT) :: ierr, intraCommID, interCommID
      INTEGER :: color, myid, mpierr, nprocs, nrows
      ! Set a null result
      ierr = 0
      intraCommID = MPI_UNDEFINED
      interCommID = MPI_UNDEFINED
      intraComm = MPI_COMM_NULL
      interComm = MPI_COMM_NULL
      IF (ncols < 1) THEN
         WRITE(ERROR_UNIT,850) ncols
         ierr = 1
         RETURN         
      ENDIF
      ! Check divisibility - TODO set processes out of square to MPI_UNDEFINED
      IF (MOD(nprocs, ncols) /= 0) WRITE(OUTPUT_UNIT,900)
      ! Determine size of communicator and process' ID
      CALL MPI_Comm_rank(mpiComm, myid, mpierr)
      CALL MPI_Comm_size(mpiComm, nprocs, mpierr)
      nrows = nprocs/ncols
      ! Edge case
      IF (nprocs == 1) THEN
         CALL MPI_Comm_dup_with_info(mpiComm, MPI_INFO_NULL, intraComm, mpierr)
         CALL MPI_Comm_rank(intraComm, intraCommID, mpierr)
         CALL MPI_Comm_dup_with_info(mpiComm, MPI_INFO_NULL, intercomm, mpierr) 
         CALL MPI_Comm_rank(mpiComm, interCommID, mpierr)
         RETURN
      ENDIF
      ! Split the intra communicator as described in preamble
      color = myid/nrows !MOD(myid, nrows) ! Processes with same color map to same group
      CALL MPI_Comm_split(mpiComm, color, myid, intraComm, mpierr) 
      CALL MPI_Comm_rank(intraComm, intraCommID, mpierr)
!print *, myid, intraCommID
      ! Figure out process' new ID on the split intra-communicator
  500 CONTINUE
      ! Now split the inter-communicator
      color = MOD(myid, nrows) !myid/nrows ! Processes with same color map to same group
      CALL MPI_Comm_split(mpiComm, color, myid, interComm, mpierr)
      CALL MPI_Comm_rank(interComm, interCommID, mpierr)
!print *, myid, interCommID
  850 FORMAT('xcloc_utilsMPI_splitComm: At least one column required')
  900 FORMAT('xcloc_utilsMPI_splitComm: Warning groups groups may not be evenly split')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!
!>    @brief Creates a communicator from the head processes.
!>    @param[in] mpiGlobalComm   Communicator on which all processes reside.
!>    @param[in] myGlobalID      Process ID on global communicator.
!>    @param[in] ngroups         Number of groups in which mpiGlobalComm has been split.
!>    @param[out] mpiHeadComm    Communicator for MPI groups.  The communicator
!>                               is formed from each subgroup's root process ID
!>                               which is likely 0.
!>    @ingroup utilsMPI
      SUBROUTINE xcloc_utilsMPI_head2Comm(mpiGlobalComm, myGlobalID, &
                                          ngroups, mpiHeadComm)      &
      BIND(C, NAME='xcloc_utilsMPI_head2Comm')
      IMPLICIT NONE
      TYPE(MPI_Comm), VALUE, INTENT(IN) :: mpiGlobalComm
      INTEGER, INTENT(IN) :: myGlobalID, ngroups
      TYPE(MPI_Comm), INTENT(OUT) :: mpiHeadComm
      INTEGER, ALLOCATABLE :: ranks(:)
      TYPE(MPI_Group) myHeadGroup, myHeadRow
      INTEGER nprocs, nrows, i, mpierr, rem
      INTEGER, PARAMETER :: master = 0
      ! Set a null result
      mpiHeadComm = MPI_COMM_NULL 
      ! Completely serial
      CALL MPI_Comm_size(mpiGlobalComm, nprocs, mpierr)
      IF (nprocs == 1) THEN
         IF (myGlobalID == master) THEN
            CALL MPI_Comm_dup_with_info(mpiGlobalComm, MPI_INFO_NULL, mpiHeadComm, mpierr)
         ENDIF
         RETURN
      ENDIF
      ! Slice through the group the group in the other direction
      rem = MOD(myGlobalID, ngroups)
      ALLOCATE(ranks(ngroups))
      nrows = nprocs/ngroups
      DO 1 i=1,ngroups
         ranks(i) = (i - 1)*nrows + rem
    1 CONTINUE
      CALL MPI_Comm_group(mpiGlobalComm, myHeadGroup, mpierr)
      CALL MPI_Group_incl(myHeadGroup, ngroups, ranks, myHeadRow, mpierr)
      CALL MPI_Comm_create(mpiGlobalComm, myHeadRow, mpiHeadComm, mpierr)
      !IF (my_local_id == master) THEN
      !   CALL MPI_COMM_RANK(mpiHeadComm,my_head_id,mpierr)
      !   print *, myGlobalID, my_local_id, my_head_id
      !ENDIF
      DEALLOCATE(ranks)
      RETURN
      END
END MODULE
