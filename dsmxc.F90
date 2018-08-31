!> @defgroup dsmxc Diffraction Stack Migration of Correlograms
!> @ingroup xcloc
!> @brief Diffraction stack migration of cross-correlograms.
!> @author Ben Baker
!> @copyright Ben Baker distributed the MIT license.
MODULE XCLOC_DSMXC
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
      USE XCLOC_IPPS
      USE XCLOC_MEMORY
      USE XCLOC_UTILS
#ifdef _OPENMP
      USE OMP_LIB
#endif
      IMPLICIT NONE
      !> Holds the cross-correlograms.  This is an array of dimension
      !> [dataOffset_ x nxcs_] stored in column major format. 
      REAL(C_FLOAT), PRIVATE, TARGET, ALLOCATABLE, SAVE :: xcs32f_(:)
      !> Holds the travel time tables.  This is an array of dimension [ldg_ x ntables_]
      !> stored in column major format.
      INTEGER(C_INT), PRIVATE, TARGET, ALLOCATABLE, SAVE :: ttimes_(:)
      !> Holds the migration image.  This is an array of dimension [ngrd_].
      REAL(C_FLOAT), PRIVATE, TARGET, ALLOCATABLE, SAVE :: image32f_(:)
      !> Maps from the is'th signal number to the table index.  This is an array
      !> of dimension [ntables_].
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: signal2Table_(:)
      !> This is defines the two signal numbers comprising a correlation pair.
      !> This is an array of dimension [2 x nxcPairs_] stored in column major format.
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: xcPairs_(:)
      !> Maps from the ixc'th correlation to the table pairs.  This in an array of
      !> dimension [2 x nxcPairs_] stored in column major format.
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: xcTablePairs_(:)
      !> Determines if all tables are set.
      LOGICAL, PRIVATE, ALLOCATABLE, SAVE :: lhaveTable_(:)
      !> Sampling period (seconds) of signals.
      DOUBLE PRECISION, PRIVATE, SAVE :: dt_ = 0.d0
      !> Leading dimension of the travel time table.
      INTEGER, PRIVATE, SAVE :: ldg_ = 0
      !> Number of grid points in each table.
      INTEGER, PRIVATE, SAVE :: ngrd_ = 0
      !> Number of points in correlograms.
      INTEGER, PRIVATE, SAVE :: nptsInXCs_ = 0
      !> Number of travel time tables.
      INTEGER, PRIVATE, SAVE :: ntables_ = 0
      !> Number of cross-correlation pairs.
      INTEGER, PRIVATE, SAVE :: nxcPairs_ = 0
      !> Leading dimension of correlograms.
      INTEGER, PRIVATE, SAVE :: dataOffset_ = 0
      !> Controls the verbosity.
      INTEGER, PRIVATE, SAVE :: verbose_ = XCLOC_PRINT_WARNINGS
      !> Flag indicating that all tables are set.
      LOGICAL, PRIVATE, SAVE :: lhaveAllTables_ = .FALSE.
      !> Flag indicating that the image has been computed.
      LOGICAL, PRIVATE, SAVE :: lhaveImage_ = .FALSE.
      !> A tuning parameter that can improve performance by emphasizing cache coherency.
      INTEGER, PRIVATE, SAVE :: blockSize_ = XCLOC_DEFAULT_BLOCKSIZE 
      !> The data alignment.
      INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: dataAlignment = 64
      !> The size of a float.
      INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: sizeof_float = 4
      !> The size of an integer.
      INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: sizeof_int = 4

      PUBLIC :: xcloc_dsmxc_initialize
      PUBLIC :: xcloc_dsmxc_compute
      PUBLIC :: xcloc_dsmxc_finalize
      PUBLIC :: xcloc_dsmxc_setBlockSize
      PUBLIC :: xcloc_dsmxc_setTable64f
      PUBLIC :: xcloc_dsmxc_setTable32f
      PUBLIC :: xcloc_dsmxc_getImage64f
      PUBLIC :: xcloc_dsmxc_getImage32f
      PUBLIC :: xcloc_dsmxc_getNumberOGridPointsInTable
      PUBLIC :: xcloc_dsmxc_getNumberOfTables
      PUBLIC :: xcloc_dsmxc_setCorrelograms64f
      PUBLIC :: xcloc_dsmxc_setCorrelograms32f
      PUBLIC :: xcloc_dsmxc_setCorrelogram64fF 
      PUBLIC :: xcloc_dsmxc_setCorrelogram32fF
      PUBLIC :: xcloc_dsmxc_signalToTableIndex
      PUBLIC :: xcloc_dsmxc_haveAllTables
      PUBLIC :: xcloc_dsmxc_getImageMax
      ! Public but for Fortran only
      PUBLIC :: xcloc_dsmxc_getImagePtr

      PRIVATE :: xcloc_dsmxc_checkTables
      CONTAINS
!----------------------------------------------------------------------------------------!
!                                   Begin the Code                                       !
!----------------------------------------------------------------------------------------!
!>    @brief Initializes the diffraction statck migration cross-correlation module.
!>    @param[in] ngrd       Number of grid points in each table.
!>    @param[in] nxcPairs   Number of cross-correlation pairs.
!>    @param[in] nptsInXCs  Number of points in cross-correlograms.
!>    @param[in] dt         Sampling period of cross-correlograms (seconds).
!>    @param[in] xcPairs    This is a map from the ixc'th correlation to the table
!>                          pairs.  It is a column major matrix with dimension 
!>                          [2 x nxcPairs].
!>    @param[in] verbose    Controls the verbosity.
!>    @param[out] ierr      0 indicates success.  
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_initialize(ngrd, nxcPairs, nptsInXCs,   &
                                        dt, xcPairs, verbose, ierr)  &
      BIND(C, NAME='xcloc_dsmxc_initialize')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngrd, nptsInXCs, nxcPairs, verbose
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt
      INTEGER(C_INT), INTENT(IN) :: xcPairs(2*nxcPairs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER ixc
      ! Check the variables
      ierr = 0
      CALL xcloc_dsmxc_finalize()
      IF (ngrd < 1) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
      ENDIF
      IF (nxcPairs < 1) THEN
         WRITE(ERROR_UNIT,901)
         ierr = 1
      ENDIF
      IF (dt <= 0.d0) THEN
         WRITE(ERROR_UNIT,902)
         ierr = 1
      ENDIF
      ! Create a map of unique tables
      CALL xcloc_utils_unique32s(2*nxcPairs, xcPairs, ntables_, signal2Table_, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,905)
         ierr = 1
      ENDIF
      IF (ntables_ < 2) THEN
         WRITE(ERROR_UNIT,906)
         ierr = 1
      ENDIF
      IF (nptsInXCs < 1 .OR. MOD(nptsInXCs, 2) /= 1) THEN
         IF (nptsInXCs < 1) WRITE(ERROR_UNIT,907)
         IF (MOD(nptsInXCs, 2) /= 1) WRITE(ERROR_UNIT,908) nptsInXCs
         ierr = 1
      ENDIF
      IF (ierr /= 0) THEN
         CALL xcloc_dsmxc_finalize()
         RETURN
      ENDIF
      ! Set the variables 
      lhaveAllTables_ = .FALSE.
      ngrd_ = ngrd
      nptsInXCs_ = nptsInXCs
      nxcPairs_ = nxcPairs
      verbose_ = verbose
      dt_ = dt
      dataOffset_ = xcloc_memory_padLength(dataAlignment, sizeof_float, nptsInXCs_) 
      ldg_ = xcloc_memory_padLength(dataAlignment, sizeof_int, ngrd_)
      ALLOCATE(lhaveTable_(ntables_)); lhaveTable_(:) = .FALSE.
      ALLOCATE(xcs32f_(dataOffset_*nptsInXCs_)); xcs32f_(:) = 0.0
      ALLOCATE(image32f_(ldg_)); image32f_(:) = 0.0
      ALLOCATE(ttimes_(ldg_*ntables_)); ttimes_(:) = 0
      ALLOCATE(xcPairs_(2*nxcPairs_)); xcPairs_(:) = 0
      ALLOCATE(xcTablePairs_(2*nxcPairs_)); xcTablePairs_(:) = 0
      xcPairs_(1:2*nxcPairs_) = xcPairs(1:2*nxcPairs_)
      DO ixc=1,2*nxcPairs_
         CALL xcloc_dsmxc_signalToTableIndex(xcPairs_(ixc), xcTablePairs_(ixc), ierr) 
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,910) ixc
         ENDIF
      ENDDO
  900 FORMAT('xcloc_dsmxc_initialize: No grid points')
  901 FORMAT('xcloc_dsmxc_initialize: No correlation pairs')
  902 FORMAT('xcloc_dsmxc_initialize: Sampling period must be positive')
  905 FORMAT('xcloc_dsmxc_initialize: Failed to make signal to table map')
  906 FORMAT('xcloc_dsmxc_initialize: There should be at least two unique tables') 
  907 FORMAT('xcloc_dsmxc_initialize: No points in xcs')
  908 FORMAT('xcloc_dsmxc_initialize: Number of points in xcs is even', I6)
  910 FORMAT('xcloc_dsmxc_initialize: Failed to from signal', I4, ' to table index')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases the memory on the module.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_finalize( ) &
      BIND(C, NAME='xcloc_dsmxc_finalize')
      IF (ALLOCATED(lhaveTable_))   DEALLOCATE(lhaveTable_)
      IF (ALLOCATED(xcs32f_))       DEALLOCATE(xcs32f_)
      IF (ALLOCATED(ttimes_))       DEALLOCATE(ttimes_)
      IF (ALLOCATED(image32f_))     DEALLOCATE(image32f_)
      IF (ALLOCATED(xcPairs_))      DEALLOCATE(xcPairs_)
      IF (ALLOCATED(xcTablePairs_)) DEALLOCATE(xcTablePairs_)
      IF (ALLOCATED(signal2Table_)) DEALLOCATE(signal2Table_)
      dt_  = 0.d0
      ldg_ = 0 
      ngrd_ = 0 
      nptsInXCs_ = 0 
      ntables_ = 0 
      nxcPairs_ = 0 
      dataOffset_ = 0 
      verbose_ = XCLOC_PRINT_WARNINGS
      lhaveAllTables_ = .FALSE.
      lhaveImage_ = .FALSE.
      blockSize_ = XCLOC_DEFAULT_BLOCKSIZE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Convenience utility to determine if all the tables have been set.
!>    @param[out] lhaveAllTables  If true then all tables have been set.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_haveAllTables(lhaveAllTables) &
      BIND(C, NAME='xcloc_dsmxc_haveAllTables')
      LOGICAL(C_BOOL), INTENT(OUT) :: lhaveAllTables
      lhaveAllTables = lhaveAllTables_
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the total number of grid points in a table.
!>    @param[out] ngrd   Number of grid points in each table.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_getNumberOGridPointsInTable(ngrd) &
      BIND(C, NAME='xcloc_dsmxc_getNumberOGridPointsInTable')
      INTEGER(C_INT), INTENT(OUT) :: ngrd
      ngrd = ngrd_
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Returns the number of travel time tables on the module.
!>    @param[out] ntables  The number of travel time tables.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_getNumberOfTables(ntables) &
      BIND(C, NAME='xcloc_dsmxc_getNumberOfTables')
      INTEGER(C_INT), INTENT(OUT) :: ntables
      ntables = ntables_
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the maximum value of the image.  
!>    @param[out] maxIndex  The index of the maximum.  This is Fortran indexed.
!>    @param[out] maxValue  Maximum value corresponding the maxIndex.
!>    @param[out] ierr      0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_getImageMax(maxIndex, maxValue, ierr) &
      BIND(C, NAME='xcloc_dsmxc_getImageMax')
      REAL(C_FLOAT), INTENT(OUT) :: maxValue
      INTEGER(C_INT), INTENT(OUT) :: maxIndex, ierr
      ierr = 0
      maxIndex = 1
      maxValue = 0.0
      IF (.NOT.lhaveImage_) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      ierr = ippsMaxIndx_32f(image32f_, ngrd_, maxValue, maxIndex)
      IF (ierr /= ippStsNoErr) THEN
         WRITE(ERROR_UNIT,905)
         ierr = 1
         RETURN
      ENDIF
      maxIndex = maxIndex + 1 ! C to Fortran
      !maxIndex = MAXLOC(image32f_(1:ngrd_), 1)
      !maxValue = image32f_(maxIndex)      
  900 FORMAT('xcloc_dsmxc_getImageMax: Image not yet computed')
  905 FORMAT('xcloc_dsmxc_getImageMax: Failed to get image max')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets a pointer to the migration image.
!>    @param[out] imagePtr  A pointer to the image. 
!>    @param[out] ierr      0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_getImagePtr(imagePtr, ierr)
      REAL(C_FLOAT), CONTIGUOUS, POINTER, DIMENSION(:), INTENT(OUT) :: imagePtr
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (.NOT.lhaveImage_) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1 
      ENDIF
      imagePtr(1:ngrd_) => image32f_(1:ngrd_)
  900 FORMAT('xcloc_dsmxc_getImagPtr: Image not yet computed')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gets the diffraction stack migration image of the correlograms from the
!>           module.
!>    @param[in] nwork   Size of image.
!>    @param[out] image  This contains the DSM image.  This is an array of dimension
!>                       [nwork] but only the first ngrd_ points are accessed. 
!>    @param[out] ierr   0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_getImage64f(nwork, image, ierr) &
      BIND(C, NAME='xcloc_dsmxc_getImage64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: nwork
      REAL(C_DOUBLE), INTENT(OUT) :: image(nwork)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (.NOT.lhaveImage_) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1 
         RETURN
      ENDIF
      IF (nwork < ngrd_) THEN
         WRITE(ERROR_UNIT,905) nwork, ngrd_
         ierr = 1
         RETURN
      ENDIF
      image(1:ngrd_) = DBLE(image32f_(1:ngrd_))
  900 FORMAT('xcloc_dsmxc_getImage64f: Image not yet computed')
  905 FORMAT('xcloc_dsmxc_getImage64f: nwork=', I6, ' must be at least ngrd_=', I6)
      RETURN
      END
!>    Gets the diffraction stack migration image of the correlograms from the module.
!>    @param[in] nwork   Size of image.
!>    @param[out] image  This contains the DSM image.  This is an array of dimension
!>                       [nwork] but only the first ngrd_ points are accessed. 
!>    @param[out] ierr   0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_getImage32f(nwork, image, ierr) &
      BIND(C, NAME='xcloc_dsmxc_getImage32f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: nwork
      REAL(C_FLOAT), INTENT(OUT) :: image(nwork)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (.NOT.lhaveImage_) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (nwork < ngrd_) THEN
         WRITE(ERROR_UNIT,905) nwork, ngrd_
         ierr = 1
         RETURN
      ENDIF
      image(1:ngrd_) = image32f_(1:ngrd_)
  900 FORMAT('xcloc_dsmxc_getImage32f: Image not yet computed')
  905 FORMAT('xcloc_dsmxc_getImage32f: nwork=', I6, ' must be at least ngrd_=', I6)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the correlograms to migrate on the module.
!>    @param[in] ldxc       Leading dimension of nxcs.  This must be greater than or equal
!>                          to nptsInXCs.
!>    @param[in] nptsInXCs  The number of points in each correlogram.  This must
!>                          equal nptsInXCs_.
!>    @param[in] nxcPairs   The number of correlogram pairs.  This must equal nxcPairs_.
!>    @param[in] xcs        The cross-correlograms to migrate.  This an array of dimension
!>                          [ldxc x nxcPairs] in column major format with leading
!>                          dimension ldxc.  The order of the correlograms is dictated
!>                          by xcPairs_ - i.e., for the ixc'th correlation the correlogram
!>                          is expected to represent the cross-correlation between the
!>                          signal index pairs given by xcPairs_(2*(ixc-1)+1) and
!>                          xcPairs_(2*ixc).
!>    @param[out] ierr      0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_setCorrelograms64f(ldxc, nptsInXCs, nxcPairs, xcs, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setCorrelograms64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldxc, nptsInXCs, nxcPairs
      REAL(C_DOUBLE), INTENT(IN) :: xcs(ldxc*nxcPairs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i, indx
      ierr = 0
      IF (nxcPairs /= nxcPairs_) THEN
         WRITE(ERROR_UNIT,900) nxcPairs, nxcPairs_
         ierr = 1
         RETURN
      ENDIF
      IF (nptsInXCs /= nptsInXCs_) THEN
         WRITE(ERROR_UNIT,901) nptsInXCs, nptsInXCs_
         ierr = 1
         RETURN
      ENDIF
      DO i=1,nxcPairs_
         indx = (i - 1)*ldxc + 1 
         CALL xcloc_dsmxc_setCorrelogram64fF(i, nptsInXCs, xcs(indx), ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,902) i
            RETURN
         ENDIF
      ENDDO
  900 FORMAT('xcloc_dsmxc_setCorrelograms64f: nxcPairs=', I4,'- expecting nxcPairs=',I4)
  901 FORMAT('xcloc_dsmxc_setCorrelograms64f: nPtsInXCs=', I6,'- expecting nptsInXCS=',I6)
  902 FORMAT('xcloc_dsmxc_setCorrelograms64f: Error setting xc number', I4) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the correlograms to migrate on the module.
!>    @param[in] ldxc       Leading dimension of nxcs.  This must be greater than or equal
!>                          to nptsInXCs.
!>    @param[in] nptsInXCs  The number of points in each correlogram.  This must
!>                          equal nptsInXCs_.
!>    @param[in] nxcPairs   The number of correlograms pairs.  This must equal nxcPairs_.
!>    @param[in] xcs        The cross-correlograms to migrate.  This an array of dimension
!>                          [ldxc x nxcPairs] in column major format with leading
!>                          dimension ldxc.  The order of the correlograms is dictated
!>                          by xcPairs_ - i.e., for the ixc'th correlation the correlogram
!>                          is expected to represent the cross-correlation between the
!>                          signal index pairs given by xcPairs_(2*(ixc-1)+1) and
!>                          xcPairs_(2*ixc).
!>    @param[out] ierr      0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_setCorrelograms32f(ldxc, nptsInXCs, nxcPairs, xcs, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setCorrelograms32f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldxc, nptsInXCs, nxcPairs
      REAL(C_FLOAT), INTENT(IN) :: xcs(ldxc*nxcPairs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i, indx
      ierr = 0
      IF (nxcPairs /= nxcPairs_) THEN
         WRITE(ERROR_UNIT,900) nxcPairs, nxcPairs_
         ierr = 1
         RETURN
      ENDIF
      IF (nptsInXCs /= nptsInXCs_) THEN
         WRITE(ERROR_UNIT,901) nptsInXCs, nptsInXCs_
         ierr = 1
         RETURN
      ENDIF
      DO i=1,nxcPairs_
         indx = (i - 1)*ldxc + 1
         CALL xcloc_dsmxc_setCorrelogram32fF(i, nptsInXCs, xcs(indx), ierr)
         IF (ierr /= 0) THEN
            WRITE(ERROR_UNIT,902) i
            RETURN
         ENDIF
      ENDDO
  900 FORMAT('xcloc_dsmxc_setCorrelograms32f: nxcPairs=', I4,'- expecting nxcPairs=',I4)
  901 FORMAT('xcloc_dsmxc_setCorrelograms32f: nPtsInXCs=', I6,'- expecting nptsInXCS=',I6)
  902 FORMAT('xcloc_dsmxc_setCorrelograms32f: Error setting xc number', I4)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the xcIndex'th correlogram.
!>    @param[in] xcIndex    Cross-correlation index.  This must be in the range
!>                          [1,nxcPairs_].
!>    @param[in] nptsInXCs  Number of points in cross-correlation.  This must equal
!>                          nptsInXCs_.
!>    @param[in] xc         Cross-correlogram to set.  This has dimension [nptsInXCs].
!>    @param[out] ierr      0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_setCorrelogram64fF(xcIndex, nptsInXCs, xc, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setCorrelogram64fF')
      INTEGER(C_INT), VALUE, INTENT(IN) :: xcIndex, nptsInXCs
      REAL(C_DOUBLE), INTENT(IN) :: xc(nptsInXCs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i1, i2
      ierr = 0
      IF (nptsInXCs /= nptsInXCs_ .OR. xcIndex < 1 .OR. xcIndex > nxcPairs_) THEN
         IF (nptsInXCs /= nptsInXCs_) WRITE(ERROR_UNIT,900) nptsInXCs, nptsInXCs_
         IF (xcIndex < 1 .OR. xcIndex > nxcPairs_) WRITE(ERROR_UNIT,901) nxcPairs_
         ierr = 1
         RETURN
      ENDIF
      i1 = (xcIndex - 1)*dataOffset_ + 1
      i2 = i1 + nptsInXCs_ - 1 
      xcs32f_(i1:i2) = SNGL(xc(1:nptsInXCs_))
  900 FORMAT('xcloc_dsmxc_setCorrelogram64fF: nptsInXCs=', I6,'- expecting nptsInXCS=',I6)
  901 FORMAT('xcloc_dsmxc_setCorrelogram64fF: xcIndex must be in range [1,',I4,']')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the xcIndex'th correlogram.
!>    @param[in] xcIndex    Cross-correlation index.  This must be in the range
!>                          [1,nxcPairs_].
!>    @param[in] nptsInXCs  Number of points in cross-correlation.  This must equal
!>                          nptsInXCs_.
!>    @param[in] xc         Cross-correlogram to set.  This has dimension [nptsInXCs].
!>    @param[out] ierr      0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_setCorrelogram32fF(xcIndex, nptsInXCs, xc, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setCorrelogram32fF')
      INTEGER(C_INT), VALUE, INTENT(IN) :: xcIndex, nptsInXCs
      REAL(C_FLOAT), INTENT(IN) :: xc(nptsInXCs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i1, i2
      ierr = 0
      IF (nptsInXCs /= nptsInXCs_ .OR. xcIndex < 1 .OR. xcIndex > nxcPairs_) THEN
         IF (nptsInXCs /= nptsInXCs_) WRITE(ERROR_UNIT,900) nptsInXCs, nptsInXCs_
         IF (xcIndex < 1 .OR. xcIndex > nxcPairs_) WRITE(ERROR_UNIT,901) nxcPairs_
         ierr = 1
         RETURN
      ENDIF 
      i1 = (xcIndex - 1)*dataOffset_ + 1
      i2 = i1 + nptsInXCs_ - 1
      xcs32f_(i1:i2) = xc(1:nptsInXCs_)
  900 FORMAT('xcloc_dsmxc_setCorrelogram32fF: nptsInXCs=', I6,'- expecting nptsInXCS=',I6)
  901 FORMAT('xcloc_dsmxc_setCorrelogram32fF: xcIndex must be in range [1,',I4,']')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the block size in the migration loop.
!>    @param[in] blockSize  This many grid points will be updated during the migration.
!>                          This must be a power of 2.
!>    @param[out] ierr      0 indicates success.  
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_setBlockSize(blockSize, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setBlockSize')
      INTEGER(C_INT), VALUE, INTENT(IN) :: blockSize
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i, i2
      ierr = 0
      i2 = 1
      blockSize_ = XCLOC_DEFAULT_BLOCKSIZE
      DO i=1,30
         IF (blockSize == i2) GOTO 500
         i2 = i2*2
      ENDDO
      ierr = 1
      WRITE(ERROR_UNIT,900) blockSize
  900 FORMAT('xcloc_dsmxc_setBlockSize: blockSize=', I6, 'is not a power of 2')
      RETURN
  500 CONTINUE 
      blockSize_ = blockSize 
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
      SUBROUTINE xcloc_dsmxc_signalToTableIndex(is, it, ierr) &
      BIND(C, NAME='xcloc_dsmxc_signalToTableIndex')
      INTEGER(C_INT), VALUE, INTENT(IN) :: is
      INTEGER(C_INT), INTENT(OUT) :: it, ierr
      LOGICAL, PARAMETER :: lshowError = .TRUE.
      CALL xcloc_utils_bsearch32i(ntables_, is, signal2Table_, it, ierr, lshowError)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900) is
         it = 0
      ENDIF
  900 FORMAT('xcloc_dsmxc_signalToTableIndex: Failed to find signal number', I4, &
             ' in table')
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets a travel time table.
!>    @param[in] tableNumber  Table index to set.  This must be in the range 
!>                            [1, ntables_].
!>    @param[in] ngrd         Number of grid points in the table.  This must equal ngrd_.
!>    @param[in] table        The travel times from the source to all points in the
!>                            grid in seconds.  This has dimension [ngrd].
!>    @param[out] ierr        0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_setTable64f(tableNumber, ngrd, table, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setTable64f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: tableNumber, ngrd
      REAL(C_DOUBLE), INTENT(IN) :: table(ngrd)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i, igrd
      ierr = 0 
      IF (tableNumber < 1 .OR. tableNumber > ntables_) THEN
         WRITE(ERROR_UNIT,900) tableNumber, ntables_ 
         ierr = 1 
      ENDIF
      IF (ngrd /= ngrd_) THEN
         WRITE(ERROR_UNIT,905) ngrd, ngrd_
         ierr = 1 
      ENDIF
      IF (ierr /= 0) RETURN
      igrd = (tableNumber - 1)*ldg_
      DO i=1,ngrd_
         ttimes_(igrd+i) = INT(table(i)/dt_ + 0.5d0)
      ENDDO
      lhaveTable_(tableNumber) = .TRUE.
      ! Check if all tables are set
      DO i=1,ntables_
         IF (.NOT.lhaveTable_(i)) GOTO 500
      ENDDO
      lhaveAllTables_ = .TRUE.
  500 CONTINUE
      ! Check the tables
      IF (lhaveAllTables_) THEN
         CALL xcloc_dsmxc_checkTables(ierr)
         IF (ierr /= 0) WRITE(OUTPUT_UNIT,910)
      ENDIF
  900 FORMAT('xcloc_dsmxc_setTable64f: tableNumber=', I4, ' must be in range [1,',I4,']')
  905 FORMAT('xcloc_dsmxc_setTable64f: ngrd=', I6, ' expecting ngrd_=', I6) 
  910 FORMAT('xcloc_dsmxc_setTable64f: All tables set')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets a travel time table.
!>    @param[in] tableNumber  Table index to set.  This must be in the range 
!>                            [1, ntables_].
!>    @param[in] ngrd         Number of grid points in the table.  This must equal ngrd_.
!>    @param[in] table        The travel times from the source to all points in the
!>                            grid in seconds.  This has dimension [ngrd].
!>    @param[out] ierr        0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_setTable32f(tableNumber, ngrd, table, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setTable32f')
      INTEGER(C_INT), VALUE, INTENT(IN) :: tableNumber, ngrd
      REAL(C_FLOAT), INTENT(IN) :: table(ngrd)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_FLOAT) dt4
      INTEGER i, igrd
      ierr = 0
      IF (tableNumber < 1 .OR. tableNumber > ntables_) THEN
         WRITE(ERROR_UNIT,900) tableNumber, ntables_ 
         ierr = 1
      ENDIF
      IF (ngrd /= ngrd_) THEN
         WRITE(ERROR_UNIT,905) ngrd, ngrd_
         ierr = 1
      ENDIF
      IF (ierr /= 0) RETURN
      igrd = (tableNumber - 1)*ldg_
      dt4 = SNGL(dt_)
      DO i=1,ngrd_
         ttimes_(igrd+i) = INT(table(i)/dt4 + 0.5)
      ENDDO
      ! Check if all tables are set
      DO i=1,ntables_
         IF (.NOT.lhaveTable_(i)) GOTO 500 
      ENDDO
      lhaveAllTables_ = .TRUE.
  500 CONTINUE
      ! Check the tables
      IF (lhaveAllTables_) THEN
         CALL xcloc_dsmxc_checkTables(ierr)
         IF (ierr /= 0) WRITE(OUTPUT_UNIT,910) 
      ENDIF
  900 FORMAT('xcloc_dsmxc_setTable32f: tableNumber=', I4, ' must be in range [1,',I4,']')
  905 FORMAT('xcloc_dsmxc_setTable32f: ngrd=', I6, ' expecting ngrd_=', I6)
  910 FORMAT('xcloc_dsmxc_setTable32f: All tables set')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the diffraction stack migration image of the cross-correlograms.
!>           Note, that the correlations are defined by
!>           \f$
!>             s_1 \star s_2 (t) = \int_{-\infty}^\infty s_1(\tau+t) s_2(\tau) \, d\tau
!>           \f$
!>           whose maximum is the number of samples to move signal 2 so that it is
!>           best aligned with signal 1.
!>
!>    @param[out] ierr 0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_compute(ierr) &
      BIND(C, NAME='xcloc_dsmxc_compute')
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i, igrd, igrd1, igrd2, indxXC, it1, it2, ixc, ixc1, ixc2, jgrd1, jgrd2, &
              kgrd1, kgrd2, lxc2, ngrdLoc
      INTEGER(C_INT), CONTIGUOUS, POINTER :: tt1(:), tt2(:)
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: imagePtr32f(:)
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: xcPtr32f(:)
      ! This would be problematic if all tables aren't set
      ierr = 0
      lhaveImage_ = .FALSE.
      IF (.NOT.lhaveAllTables_) WRITE(OUTPUT_UNIT,905)
  905 FORMAT('xcloc_dsmxc_compute: Only a subset of tables were set')
      ! Compute
      lxc2 = nptsInXCs_/2 
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(blockSize_, dataOffset_, image32f_, ldg_, ngrd_) &
      !$OMP SHARED(nptsInXCs_, nxcPairs_, ttimes_, xcTablePairs_, xcs32f_) &
      !$OMP PRIVATE(i, igrd, igrd1, ixc, ixc1, ixc2, igrd2, it1, it2, indxXC) &
      !$OMP PRIVATE(jgrd1, jgrd2, kgrd1, kgrd2, imagePtr32f, ngrdLoc)   &
      !$OMP PRIVATE(tt1, tt2, xcPtr32f) &
      !$OMP FIRSTPRIVATE(lxc2)
      DO igrd=1,ngrd_,blockSize_
         ngrdLoc = MIN(ngrd_ - igrd + 1, blockSize_)
         igrd1 = igrd
         igrd2 = igrd + ngrdLoc - 1
         imagePtr32f => image32f_(igrd1:igrd2)
         imagePtr32f(1:ngrdLoc) = 0.0 ! 0 out block prior to stack
         ! Loop on correlation pairs
         DO ixc=1,nxcPairs_
            it1 = xcTablePairs_(2*(ixc-1)+1)
            it2 = xcTablePairs_(2*(ixc-1)+2)
            jgrd1 = (it1 - 1)*ldg_ + igrd1
            jgrd2 = (it1 - 1)*ldg_ + igrd2
            kgrd1 = (it2 - 1)*ldg_ + igrd1 
            kgrd2 = (it2 - 1)*ldg_ + igrd2
            ixc1 = (ixc - 1)*dataOffset_ + 1
            ixc2 = (ixc - 1)*dataOffset_ + nptsInXCs_
            tt1 => ttimes_(jgrd1:jgrd2)
            tt2 => ttimes_(kgrd1:kgrd2)
            xcPtr32f => xcs32f_(ixc1:ixc2)
            !$OMP SIMD ALIGNED(imagePtr32f, tt1, tt2: 64)
            DO i=1,ngrdLoc
               indxXC = lxc2 + tt1(i) - tt2(i)
               imagePtr32f(i) = imagePtr32f(i) + xcPtr32f(indxXC)
            ENDDO
            NULLIFY(tt1)
            NULLIFY(tt2)
            NULLIFY(xcPtr32f)
         ENDDO
         NULLIFY(imagePtr32f)
      ENDDO
      lhaveImage_ = .TRUE.
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Verifies that the migration will not segfault wen accessing the correlograms.
!>    @param[out] ierr    0 indicates success.
!>    @ingroup dsmxc
      SUBROUTINE xcloc_dsmxc_checkTables(ierr)
      INTEGER, INTENT(OUT) :: ierr
      INTEGER igrd, igrd1, igrd2, indxXC, indxMaxXC, indxMinXC, &
              it1, it2, ixc, ixc1, ixc2, jgrd1, jgrd2, kgrd1, kgrd2, lxc2, ngrdLoc
      INTEGER(C_INT), CONTIGUOUS, POINTER :: tt1(:), tt2(:)
      ierr = 0
      IF (.NOT.lhaveAllTables_) THEN
         WRITE(ERROR_UNIT,905)
         ierr = 1
         RETURN
      ENDIF
      indxMaxXC = 0
      indxMinXC = HUGE(1)
      lxc2 = nptsInXCs_/2
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(blockSize_, dataOffset_, image32f_, ldg_, ngrd_) &
      !$OMP SHARED(nptsInXCs_, nxcPairs_, ttimes_, xcTablePairs_, xcs32f_) &
      !$OMP PRIVATE(igrd, igrd1, ixc, ixc1, ixc2, igrd2, it1, it2, indxXC) &
      !$OMP PRIVATE(jgrd1, jgrd2, kgrd1, kgrd2, ngrdLoc)   &
      !$OMP PRIVATE(tt1, tt2) &
      !$OMP FIRSTPRIVATE(lxc2) &
      !$OMP REDUCTION(max:indxMaxXC) REDUCTION(min:indxMinXC)
      DO igrd=1,ngrd_,blockSize_
         ngrdLoc = MIN(ngrd_ - igrd + 1, blockSize_)
         igrd1 = igrd
         igrd2 = igrd + ngrdLoc - 1
         DO ixc=1,nxcPairs_
            it1 = xcTablePairs_(2*(ixc-1)+1)
            it2 = xcTablePairs_(2*(ixc-1)+2)
            jgrd1 = (it1 - 1)*ldg_ + igrd1
            jgrd2 = (it1 - 1)*ldg_ + igrd2
            kgrd1 = (it2 - 1)*ldg_ + igrd1 
            kgrd2 = (it2 - 1)*ldg_ + igrd2
            ixc1 = (ixc - 1)*dataOffset_ + 1 
            ixc2 = (ixc - 1)*dataOffset_ + nptsInXCs_
            tt1(1:ngrdLoc) => ttimes_(jgrd1:jgrd2)
            tt2(1:ngrdLoc) => ttimes_(kgrd1:kgrd2)
            indxXC = lxc2 + MAXVAL(tt1(1:ngrdLoc) - tt2(1:ngrdLoc))
            indxMinXC = MIN(indxMinXC, lxc2 + MINVAL(tt1(1:ngrdLoc) - tt2(1:ngrdLoc)))
            indxMaxXC = MAX(indxMaxXC, lxc2 + MAXVAL(tt1(1:ngrdLoc) - tt2(1:ngrdLoc)))
            NULLIFY(tt1)
            NULLIFY(tt2)
         ENDDO
      ENDDO
      IF (indxMinXC < 1) THEN
         WRITE(ERROR_UNIT,910) indxMinXC
         ierr = 1
      ENDIF
      IF (indxMaxXC > nptsInXCs_) THEN
         WRITE(ERROR_UNIT,911) indxMaxXC, nptsInXCs_
         WRITE(ERROR_UNIT,912) nptsInXCs_
         ierr = 1
      ENDIF
      !print *, indxMinXC, indxMaxXC, nptsInXCs_
  905 FORMAT('xcloc_dsmxc_checkTables: Only a subset of tables were set')
  910 FORMAT('xcloc_dsmxc_checkTables: Min index=', I6, ' is less than 0')
  911 FORMAT('xcloc_dsmxc_checkTables: Max index=', I6, ' exceeds nptsInXCs_=', I6)
  912 FORMAT('xcloc_dsmxc_checktables: Segfault will occur - increase nptsInXCs to:', I6)
      RETURN
      END
END MODULE
