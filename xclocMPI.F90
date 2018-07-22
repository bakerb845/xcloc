!> @defgroup xcloc_mpi xcloc
!> @breif The parallel xcloc libary.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.

!> @defgroup xcloc_xcloc_mpi Parallel Cross-Correlation Based Event Location
!> @brief Parallel cross-correlation event location utilities.
!> @ingroup xcloc_mpi
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_MPI
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      USE MPI_F08
      USE XCLOC_SPXC
      USE XCLOC_CONSTANTS

      CONTAINS
!========================================================================================!
!                                     Begin the Code                                     !
!========================================================================================!
      SUBROUTINE xclocMPI_initialize( ) &
      BIND(C, NAME='xclocMPI_initialize')
      RETURN
      END

END MODULE
