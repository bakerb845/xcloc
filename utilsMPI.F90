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
      TYPE(MPI_Comm), INTENT(IN) :: mpiComm 
      INTEGER(C_INT), INTENT(IN) :: ncols
      TYPE(MPI_Comm), INTENT(OUT) :: intraComm, interComm
      INTEGER(C_INT), INTENT(OUT) :: ierr, intraCommID, interCommID
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
      ! Determine size of communicator and process' ID
      CALL MPI_Comm_rank(mpiComm, myid, mpierr)
      CALL MPI_Comm_size(mpiComm, nprocs, mpierr)
      ! Check divisibility - TODO set processes out of square to MPI_UNDEFINED
      IF (MOD(nprocs, ncols) /= 0) WRITE(OUTPUT_UNIT,900) nprocs, ncols
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
      ! Now split the inter-communicator
      color = MOD(myid, nrows) !myid/nrows ! Processes with same color map to same group
      CALL MPI_Comm_split(mpiComm, color, myid, interComm, mpierr)
      CALL MPI_Comm_rank(interComm, interCommID, mpierr)
!print *, myid, interCommID
  850 FORMAT('xcloc_utilsMPI_splitComm: At least one column required')
  900 FORMAT('xcloc_utilsMPI_splitComm: Warning groups groups may not be evenly split', &
             'nprocs=', I0, 'ncols=', I0)
      RETURN
      END
END MODULE
