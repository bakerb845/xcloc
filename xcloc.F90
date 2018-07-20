!> @defgroup xcloc xcloc
!> @brief Utilities for driving xcloc.
!> @author Ben Baker 
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC
      USE ISO_C_BINDING
#if defined(XCLOC_USE_MPI_F)
      USE MPI_F08
      !USE MPI_F08_TYPES, ONLY : MPI_COMM_WORLD
      TYPE(MPI_Comm), PRIVATE, SAVE :: globalComm_ != MPI_COMM_WORLD
#endif
      LOGICAL, PRIVATE, SAVE :: mpi_isInit_ = .FALSE.
integer nprocs_, myid_
      PUBLIC :: xcloc_initialize
      PUBLIC :: xcloc_finalize
      CONTAINS

      SUBROUTINE xcloc_initialize() &
      BIND(C, NAME='xcloc_initializeF')
      IMPLICIT NONE
      INTEGER mpierr, provided
#if defined(XCLOC_USE_MPI_F)
      globalComm_ = MPI_COMM_WORLD
      CALL MPI_Initialized(mpi_isInit_, mpierr)
      IF (.NOT. mpi_isInit_) THEN
         CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, mpierr)
      ENDIF
      CALL MPI_COMM_SIZE(globalComm_, nprocs_, mpierr)
      CALL MPI_COMM_RANK(globalComm_, myid_,   mpierr)
#endif
      END

      SUBROUTINE xcloc_finalize() &
      BIND(C, NAME='xcloc_finalizeF')
#if defined(XCLOC_USE_MPI_F)
      INTEGER mpierr
      IF (.NOT. mpi_isInit_) CALL MPI_FINALIZE(mpierr)
#endif
      END

END MODULE
