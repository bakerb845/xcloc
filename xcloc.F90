!> @brief Utilities for driving xcloc.
!> @author Ben Baker 
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC
      USE MPI_F08
      USE ISO_C_BINDING
      TYPE(MPI_Comm), PRIVATE, SAVE :: globalComm_ = MPI_COMM_WORLD

      LOGICAL, PRIVATE, SAVE :: mpi_isInit_ = 0
integer nprocs_, myid_
      PUBLIC :: xcloc_initialize
      PUBLIC :: xcloc_finalize
      CONTAINS

      SUBROUTINE xcloc_initialize() &
      BIND(C, NAME='xcloc_initializeF')
      IMPLICIT NONE
      INTEGER mpierr, provided
      CALL MPI_Initialized(mpi_isInit_, mpierr)
      IF (.NOT. mpi_isInit_) THEN
         CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, mpierr)
      ENDIF
      CALL MPI_COMM_SIZE(globalComm_, nprocs_, mpierr)
      CALL MPI_COMM_RANK(globalComm_, myid_,   mpierr)
      END

      SUBROUTINE xcloc_finalize() &
      BIND(C, NAME='xcloc_finalizeF')
      INTEGER mpierr
      IF (.NOT. mpi_isInit_) CALL MPI_FINALIZE(mpierr)
      END

END MODULE
