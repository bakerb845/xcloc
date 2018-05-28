!> @brief Differential stack migration of cross-correlograms.
!> @author Ben Baker
!> @copyright Ben Baker distributed the MIT license.
MODULE XCLOC_DSMXC
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
#ifdef _OPENMP
      USE OMP_LIB
#endif
      !> Holds the travel time tables.  This is an array of dimension [ldg x ntables_]
      !> stored in column major format.
      REAL(C_INT), PRIVATE, TARGET, ALLOCATABLE, SAVE :: ttimes32f_(:)
      !> Holds the migration image.  This is an array of dimension [ngrd_].
      REAL(C_FLOAT), PRIVATE, TARGET, ALLOCATABLE, SAVE :: image32f_(:)
      !> Maps from the ixc'th correlation to the table pairs.  This in an array of
      !> dimension [2 x nxcPairs_] stored in column major format.
      INTEGER, PRIVATE, ALLOCATABLE, SAVE :: xcPairs_(:)
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
      !> A tuning parameter that can improve performance by emphasizing cache coherency.
      INTEGER, PRIVATE, SAVE :: blockSize_ = XCLOC_DEFAULT_BLOCKSIZE 

      PUBLIC :: xcloc_dsmxc_initialize
      PUBLIC :: xcloc_dsmxc_setBlockSize
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
      SUBROUTINE xcloc_dsmxc_initialize(ntables, ngrd, nxcPairs, nptsInXCs, &
                                        dt, xcPairs, ierr) &
      BIND(C, NAME='xcloc_dsmxc_initialize')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ntables, ngrd, nptsInXCs, nxcPairs
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt
      INTEGER(C_INT), INTENT(IN) :: xcPairs(2*nxcPairs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ! Check the variables
      ierr = 0
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

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the block size in the migration loop.
!>    @param[in] blockSize  This many grid points will be updated during the migration.
!>                          This must be a power of 2.
!>    @param[out] ierr      0 indicates success.  
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
      SUBROUTINE xcloc_dsmxc_setTable32f(ngrd, table, ierr) &
      BIND(C, NAME='xcloc_dsmxc_setTable32f') 
      INTEGER(C_INT), VALUE, INTENT(IN) :: ngrd
      REAL(C_FLOAT), INTENT(IN) :: table(ngrd)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ierr = 0
      IF (ngrd /= ngrd_) THEN
         WRITE(*,900) ngrd, ngrd_
         ierr = 1 
      ENDIF
  900 FORMAT('xcloc_dsmxc_setTable32f: ngrd=', I6, 'expecting ngrd_=', I6)
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
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER igrd, igrd1, igrd2, indxXC, it1, it2, ixc, jgrd1, jgrd2, &
              kgrd1, kgrd2, lxc2, ngrdLoc
      REAL(C_INT), CONTIGUOUS, POINTER :: tt1(:), tt2(:)
      REAL(C_FLOAT), CONTIGUOUS, POINTER :: imagePtr32f(:)
      ierr = 0
      lxc2 = nptsInXCs_/2
      DO igrd=1,ngrd,blockSize_
         ngrdLoc = MIN(ngrd - igrd + 1, blockSize_)
         igrd1 = igrd
         igrd2 = igrd + ngrdLoc - 1
         imagePtr32f => image32f_(igrd1:igrd2)
         DO ixc=1,nxcPairs_
            it1 = xcPairs_(2*(ixc-1)+1)
            it2 = xcPairs_(2*(ixc-1)+2) 
            jgrd1 = (it1 - 1)*ldg + igrd1
            jgrd2 = (it1 - 1)*ldg + igrd2
            kgrd1 = (it2 - 1)*ldg + igrd1 
            kgrd2 = (it2 - 1)*ldg + igrd2
            DO i=1,ngrdLoc
               indxXC = lxc2 + tt1(i) - tt2(igrd)
               !imagePtr32f(i) = imagePtr32f(i) + xcPtr32f(indxXC)
            ENDDO
            NULLIFY(tt1)
            NULLIFY(tt2)
         ENDDO
         NULLIFY(imagePtr32f)
      ENDDO
      RETURN
      END

END MODULE
