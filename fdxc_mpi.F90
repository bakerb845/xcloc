!> @brief Computes the cross-correlograms via the Fourier transform with MPI.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_FDXC_MPI
      USE ISO_C_BINDING
      USE MPI_F08
      USE XCLOC_CONSTANTS
      USE XCLOC_MEMORY
      USE XCLOC_FDXC
#ifdef _OPENMP
      USE OMP_LIB
#endif
      CONTAINS
      !----------------------------------------------------------------------------------!
      !                                  Begin the Code                                  !
      !----------------------------------------------------------------------------------!
      !> @brief Initializes the Foureir transform
      SUBROUTINE xcloc_fdxc_mpi_initialize(comm, npts, nsignals, nptsPad,   &
                                           verbose, prec, accuracy, ierr  ) &
      BIND(C, NAME='xcloc_fdxc_mpi_initialize')
      TYPE(MPI_Comm) comm
      INTEGER mpierr, rank
      CALL MPI_COMM_RANK(comm, rank, mpierr)
      CALL MPI_COMM_SIZE(comm, nprocs, mpierr)
      RETURN
      END

END MODULE
