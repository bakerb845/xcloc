!> @brief Differential stack migration of cross-correlograms.
!> @author Ben Baker
!> @copyright Ben Baker distributed the MIT license.
MODULE XCLOC_DSMXC
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
      USE XCLOC_IPPS
      USE XCLOC_MEMORY
#ifdef _OPENMP
      USE OMP_LIB
#endif
      !> Holds the cross-correlograms.  This is an array of dimension
      !> [dataOffset_ x nxcs_] stored in column major format. 
      REAL(C_FLOAT), PRIVATE, TARGET, ALLOCATABLE, SAVE :: xcs32f_(:)
      !> Holds the travel time tables.  This is an array of dimension [ldg_ x ntables_]
      !> stored in column major format.
      INTEGER(C_INT), PRIVATE, TARGET, ALLOCATABLE, SAVE :: ttimes_(:)
      !> Holds the migration image.  This is an array of dimension [ngrd_].
      REAL(C_FLOAT), PRIVATE, TARGET, ALLOCATABLE, SAVE :: image32f_(:)
      !> Maps from the ixc'th correlation to the table pairs.  This in an array of
      !> dimension [2 x nxcPairs_] stored in column major format.
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: xcPairs_(:)
      !> Sampling period (seconds) of signals.
      REAL(C_DOUBLE), PRIVATE, SAVE :: dt_ = 0.d0
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
      !> A tuning parameter that can improve performance by emphasizing cache coherency.
      INTEGER, PRIVATE, SAVE :: blockSize_ = XCLOC_DEFAULT_BLOCKSIZE 
      !> The data alignment.
      INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: dataAlignment = 64
      !> The size of a float.
      INTEGER(C_SIZE_T), PRIVATE, PARAMETER :: sizeof_float = 4

      PUBLIC :: xcloc_dsmxc_initialize
      PUBLIC :: xcloc_dsmxc_finalize
      PUBLIC :: xcloc_dsmxc_setBlockSize
      PUBLIC :: xcloc_dsmxc_setTable64fF
      PUBLIC :: xcloc_dsmxc_setTable32fF
      PUBLIC :: xcloc_dsmxc_setCorrelograms32f
      PUBLIC :: xcloc_dsmxc_setCorrelogram32fF
      CONTAINS
!----------------------------------------------------------------------------------------!
!                                   Begin the Code                                       !
!----------------------------------------------------------------------------------------!
!>    @brief Initializes the diffraction statck migration cross-correlation module.
!>    @param[in] ntables    Number of travel time tables.
!>    @param[in] ngrd       Number of grid points in each table.
!>    @param[in] nxcPairs   Number of cross-correlation pairs.
!>    @param[in] nptsInXCs  Number of points in cross-correlograms.
!>    @param[in] dt         Sampling period of cross-correlograms (seconds).
!>    @param[in] xcPairs    This is a map from the ixc'th correlation to the table
!>                          pairs.  It is a column major matrix with dimension 
!>                          [2 x nxcPairs].
!>    @param[out] ierr      0 indicates success.  
!>
      SUBROUTINE xcloc_dsmxc_initialize(ntables, ngrd, nxcPairs, nptsInXCs, &
                                        dt, xcPairs, ierr)                  &
      BIND(C, NAME='xcloc_dsmxc_initialize')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ntables, ngrd, nptsInXCs, nxcPairs
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt
      INTEGER(C_INT), INTENT(IN) :: xcPairs(2*nxcPairs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ! Check the variables
      ierr = 0
      CALL xcloc_dsmxc_finalize()
      IF (ngrd < 1) THEN
         WRITE(*,900)
         ierr = 1
      ENDIF
      IF (nxcPairs < 1) THEN
         WRITE(*,901)
         ierr = 1
      ENDIF
      IF (dt <= 0.d0) THEN
         WRITE(*,902)
         ierr = 1
      ENDIF
      IF (MINVAL(xcPairs) < 1 .OR. MAXVAL(xcPairs) > ntables) THEN
         IF (MINVAL(xcPairs) < 1) WRITE(*,903) 
         IF (MAXVAL(xcPairs) > ntables) WRITE(*,904) MAXVAL(xcPairs), ntables 
         ierr = 1
      ENDIF
      IF (nptsInXCs < 1 .OR. MOD(nptsInXCs, 2) /= 1) THEN
         IF (nptsInXCs < 1) WRITE(*,905)
         IF (MOD(nptsInXCs, 2) /= 1) WRITE(*,906) nptsInXCs
         ierr = 1
      ENDIF
      IF (ierr /= 0) RETURN
      ! Set the variables 
      ngrd_ = ngrd
      nptsInXCs_ = nptsInXCs_
      ntables_ = ntables
      nxcPairs_ = nxcPairs
      dt_ = dt
      dataOffset_ = xcloc_memory_padLength(dataAlignment, sizeof_float, nptsInXCs_) 
      ALLOCATE(xcs32f_(dataOffset_*nptsInXCs_)); xcs32f_(:) = 0.0
      ALLOCATE(image32f_(ngrd_)); image32f_(:) = 0.0
  900 FORMAT('xcloc_dsmxc_initialize: No grid points')
  901 FORMAT('xcloc_dsmxc_initialize: No correlation pairs')
  902 FORMAT('xcloc_dsmxc_initialize: Sampling period must be positive')
  903 FORMAT('xcloc_dsmxc_initialize: All table indices must be positive')
  904 FORMAT('xcloc_dsmxc_initialize: Max table index=', I6, 'greater than ntables=', I6)
  905 FORMAT('xcloc_dsmxc_initialize: No points in xcs')
  906 FORMAT('xcloc_dsmxc_initialize: Number of points in xcs is even', I6)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Releases the memory on the module.
!>
      SUBROUTINE xcloc_dsmxc_finalize( ) &
      BIND(C, NAME='xcloc_dsmxc_finalize')
      IF (ALLOCATED(xcs32f_))   DEALLOCATE(xcs32f_)
      IF (ALLOCATED(ttimes_))   DEALLOCATE(ttimes_)
      IF (ALLOCATED(image32f_)) DEALLOCATE(image32f_)
      IF (ALLOCATED(xcPairs_))  DEALLOCATE(xcPairs_)
      dt_  = 0.d0
      ldg_ = 0 
      ngrd_ = 0 
      nptsInXCs_ = 0 
      ntables_ = 0 
      nxcPairs_ = 0 
      dataOffset_ = 0 
      blockSize_ = XCLOC_DEFAULT_BLOCKSIZE
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
!>    @param[in] nxcs       The number of correlograms.  This must equal nxcPairs_.
!>    @param[in] xcs        The cross-correlograms to migrate.  This an array of dimension
!>                          [ldxc x nxcPairs] in column major format with leading
!>                          dimension ldxc.  The order of the correlograms is dictated
!>                          by xcPairs_ - i.e., for the ixc'th correlation the correlogram
!>                          is expected to represent the cross-correlation between the
!>                          signal index pairs given by xcPairs_(2*(ixc-1)+1) and
!>                          xcPairs_(2*ixc).
!>    @param[out] ierr      0 indicates success.
!>
      SUBROUTINE xcloc_dsmxc_setCorrelograms32f(ldxc, nptsInXCs, nxcPairs, xcs, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setCorrelograms32f')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ldxc, nptsInXCs, nxcPairs
      REAL(C_FLOAT), INTENT(IN) :: xcs(ldxc*nxcPairs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i, indx
      ierr = 0
      IF (nxcPairs /= nxcPairs_) THEN
         WRITE(*,900) nxcPairs, nxcPairs_
         ierr = 1
         RETURN
      ENDIF
      IF (nptsInXCs /= nptsInXCs_) THEN
         WRITE(*,901) nptsInXCs, nptsInXCs_
         ierr = 1
         RETURN
      ENDIF
      DO i=1,nxcPairs_
         indx = (i - 1)*ldxc + 1
         CALL xcloc_dsmxc_setCorrelogram32fF(i, nptsInXCs, xcs(indx), ierr)
         IF (ierr /= 0) THEN
            WRITE(*,902) i
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
!>
      SUBROUTINE xcloc_dsmxc_setCorrelogram32fF(xcIndex, nptsInXCs, xc, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setCorrelogram32fF')
      INTEGER(C_INT), VALUE, INTENT(IN) :: xcIndex, nptsInXCs
      REAL(C_FLOAT), INTENT(IN) :: xc(nptsInXCs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i1, i2
      ierr = 0
      IF (nptsInXCs /= nptsInXCs_ .OR. xcIndex < 1 .OR. xcIndex > nxcPairs_) THEN
         IF (nptsInXCs /= nptsInXCs_) WRITE(*,900) nxcPairs, nxcPairs_
         IF (xcIndex < 1 .OR. xcIndex > nxcPairs_) WRITE(*,901) nxcPairs_
         ierr = 1
         RETURN
      ENDIF 
      i1 = (xcIndex - 1)*dataOffset_ + 1
      i2 = i1 + nptsInXCs_ - 1
      xcs32f_(i1:i2) = xc(1:nptsInXCs_)
  900 FORMAT('xcloc_dsmxc_setCorrelogram32fF: nPtsInXCs=', I6,'- expecting nptsInXCS=',I6)
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
!>
      SUBROUTINE xcloc_dsmxc_setBlockSize(blockSize, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setBlockSize')
      IMPLICIT NONE
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
      WRITE(*,900) blockSize
  900 FORMAT('xcloc_dsmxc_setBlockSize: blockSize=', I6, 'is not a power of 2')
      RETURN
  500 CONTINUE 
      blockSize_ = blockSize 
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
!>
      SUBROUTINE xcloc_dsmxc_setTable64fF(tableNumber, ngrd, table, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setTable64fF')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: tableNumber, ngrd
      REAL(C_DOUBLE), INTENT(IN) :: table(ngrd)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i, igrd
      ierr = 0 
      IF (tableNumber < 1 .OR. tableNumber > ntables_) THEN
         WRITE(*,900) tableNumber, ntables_ 
         ierr = 1 
      ENDIF
      IF (ngrd /= ngrd_) THEN
         WRITE(*,905) ngrd, ngrd_
         ierr = 1 
      ENDIF
      IF (ierr /= 0) RETURN
      igrd = (tableNumber - 1)*ldg_
      DO i=1,ngrd_
         ttimes_(igrd+i) = INT(table(i)/dt_ + 0.5d0)
      ENDDO
  900 FORMAT('xcloc_dsmxc_setTable64fF: tableNumber=', I4, ' must be in range [1,',I4,']')
  905 FORMAT('xcloc_dsmxc_setTable64fF: ngrd=', I6, ' expecting ngrd_=', I6) 
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
!>
      SUBROUTINE xcloc_dsmxc_setTable32fF(tableNumber, ngrd, table, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setTable32fF')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: tableNumber, ngrd
      REAL(C_FLOAT), INTENT(IN) :: table(ngrd)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      REAL(C_FLOAT) dt4
      INTEGER i, igrd
      ierr = 0
      IF (tableNumber < 1 .OR. tableNumber > ntables_) THEN
         WRITE(*,900) tableNumber, ntables_ 
         ierr = 1
      ENDIF
      IF (ngrd /= ngrd_) THEN
         WRITE(*,905) ngrd, ngrd_
         ierr = 1
      ENDIF
      IF (ierr /= 0) RETURN
      igrd = (tableNumber - 1)*ldg_
      dt4 = REAL(dt_)
      DO i=1,ngrd_
         ttimes_(igrd+i) = INT(table(i)/dt4 + 0.5)
      ENDDO
  900 FORMAT('xcloc_dsmxc_setTable32fF: tableNumber=', I4, ' must be in range [1,',I4,']')
  905 FORMAT('xcloc_dsmxc_setTable32fF: ngrd=', I6, ' expecting ngrd_=', I6)
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
!>
      SUBROUTINE xcloc_dsmxc_compute(ierr) &
      BIND(C, NAME='xcloc_dsmxc_compute')
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER i, igrd, igrd1, igrd2, indxXC, it1, it2, ixc, ixc1, ixc2, jgrd1, jgrd2, &
              kgrd1, kgrd2, lxc2, ngrdLoc
      INTEGER(C_INT), CONTIGUOUS, POINTER :: tt1(:), tt2(:)
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: imagePtr32f(:)
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: xcPtr32f(:)
      ierr = 0
      lxc2 = nptsInXCs_/2
      DO igrd=1,ngrd_,blockSize_
         ngrdLoc = MIN(ngrd_ - igrd + 1, blockSize_)
         igrd1 = igrd
         igrd2 = igrd + ngrdLoc - 1
         imagePtr32f => image32f_(igrd1:igrd2)
         DO ixc=1,nxcPairs_
            it1 = xcPairs_(2*(ixc-1)+1)
            it2 = xcPairs_(2*(ixc-1)+2)
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
               indxXC = lxc2 + tt1(i) - tt2(igrd)
               imagePtr32f(i) = imagePtr32f(i) + xcPtr32f(indxXC)
            ENDDO
            NULLIFY(tt1)
            NULLIFY(tt2)
            NULLIFY(xcPtr32f)
         ENDDO
         NULLIFY(imagePtr32f)
      ENDDO
      RETURN
      END

END MODULE
