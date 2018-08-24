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
      PUBLIC :: xcloc_utilsMPI_head2Comm
      CONTAINS
!=======================================================================================!
!                                       Begin the Code                                  !
!=======================================================================================!
!>    @brief Splits an MPI communicator with nprocs processes.  As an example the  
!>           module will split 12 processes into 3 groups.  On the original communicator
!>           the processes have ID: \n
!>              myid:  0  1  2  3  4  5  6  7  8  9 10 11 \n
!>           The three groups have `colors' given by: \n
!>             color:  0  0  0  0  1  1  1  1  2  2  2  2 \n
!>           The new process IDs in each color are then: \n
!>             newCommID:  0  1  2  3  0  1  2  3  0  1  2  3 
!>
!>           Note this is an efficacious strategy for nodes in that a common MPI default
!>           is to assign process IDs sequentially node by node, hence, if making 
!>           subgroups on a node this division strategy can decrease the likelihood that
!>           point to point will be performed on a node.
!>
!>    @param[in] mpiComm        This is the communicator to split.
!>    @param[in] ngroups        Number of groups in which to split mpiComm.
!>
!>    @param[out] color         Group (or color) to which the process belongs.
!>    @param[out] newCommID     Process ID on mpiSplitComm.
!>    @param[out] mpiSplitComm  MPI communicator for newCommID.
!>    @ingroup utilsMPI
      SUBROUTINE xcloc_utilsMPI_splitComm(mpiComm, ngroups,               &
                                          color, newCommID, mpiSplitComm) &
      BIND(C, NAME='xcloc_utilsMPI_splitComm')
      IMPLICIT NONE
      TYPE(MPI_Comm), VALUE, INTENT(IN) :: mpiComm 
      INTEGER, INTENT(IN) :: ngroups
      TYPE(MPI_Comm), INTENT(OUT) :: mpiSplitComm
      INTEGER, INTENT(OUT) :: newCommID, color
      INTEGER myid, mpierr, nprocs, nrows
      ! Determine size of communicator and process' ID
      CALL MPI_Comm_rank(mpiComm, myid, mpierr)
      CALL MPI_Comm_size(mpiComm, nprocs, mpierr)
      ! Edge case
      IF (nprocs == 1 .OR. ngroups == 1) THEN
         color = 0
         newCommID = myid
         mpiSplitComm = mpiComm
         RETURN
      ENDIF
      ! Not a good idea to make uneven groups (load balancing)
      IF (MOD(nprocs, ngroups) /= 0) WRITE(OUTPUT_UNIT,900)
      ! Split as described in preamble
      nrows = nprocs/ngroups
      color = myid/nrows  !processes with same color are in same group
      ! Split the communicator
      CALL MPI_Comm_split(mpiComm, color, myid, mpiSplitComm,mpierr) 
      ! Figure out process' new ID on the split communicator
      CALL MPI_Comm_rank(mpiSplitComm, newCommID, mpierr)
  900 FORMAT('xcloc_utilsMPI_splitComm: Warning groups groups may not be evenly split')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!
!>    @brief Creates a communicator from the head processes.
!>    @param[in] mpiGlobalComm   Communicator on which all processes reside.
!>    @param[in] myGlobalID      Process ID On global communicator.
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
      INTEGER nprocs, nrows, i, mpierr
      INTEGER, PARAMETER :: master = 0
      ! Completely serial 
      CALL MPI_Comm_size(mpiGlobalComm, nprocs, mpierr)
      IF (nprocs == 1) THEN 
         IF (myGlobalID == master) mpiHeadComm = mpiGlobalComm
         RETURN
      ENDIF
      ! otherwise have to split it up
      ALLOCATE(ranks(ngroups))
      nrows = nprocs/ngroups
      DO 1 i=1,ngroups
         ranks(i) = (i - 1)*nrows
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
