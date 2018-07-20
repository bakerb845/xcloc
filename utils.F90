!> @defgroup utils Utilities
!> @ingroup xcloc
!> @brief Generic utilities to simplify using the library.
!> @author Ben Baker
!> @copyright Ben Baker distributed under the MIT license.
MODULE XCLOC_UTILS
      USE ISO_FORTRAN_ENV
      USE ISO_C_BINDING
      USE XCLOC_CONSTANTS
      USE XCLOC_IPPS
      IMPLICIT NONE

      PUBLIC :: xcloc_utils_computeDefaultXCTable
      PUBLIC :: xcloc_utils_partitionTasks
      PUBLIC :: xcloc_utils_unique32s
      PUBLIC :: xcloc_utils_bsearch32i
      CONTAINS
!========================================================================================!
!                                     Begin the Code                                     !
!========================================================================================!
!>    @brief Computes the default cross-correlation table pairs.  The pairs can be viewed
!>           as a symmetric matrix.
!>
!>    @param[in] ldoAutoCorrs  If true then the auto-correlations are to be performed in
!>                             addition to the cross-correlations; i.e., the main diagonal
!>                             of the cross-correlation matrix is computed.
!>    @param[in] ldoAutoCorrs  If false then only the cross-correlations are computed.
!>    @param[in] nsignals      Number of signals.
!>    @param[in] nwork         The size of xcPairs which should be equal to 2*nxcs.
!>    @param[in] nwork         If nwork is negative then this is a space query and
!>                             xcPairs will not be accessed.
!>    @param[in] numbering     If numbering is XCLOC_FORTRAN_NUMBERING then the signal
!>                             indices begin at 1.
!>    @param[in] numbering     If numbering is XCLOC_C_NUMBERING then the signal indices
!>                             begin at 0.
!>
!>    @param[out] nxcs         Number of cross-correlations.
!>    @param[out] xcPairs      If nwork > 0 then this contains the pairs such that the
!>                             i'th signal is to be correlated with the j'th signal.
!>                             This has dimension [nwork] but only the first 2*nxcs 
!>                             indices will be set.
!>    @param[out] ierr         0 indicates success.
!>
      SUBROUTINE xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, &
                                                   nwork, numbering,       &
                                                   nxcs, xcPairs, ierr)    &
      BIND(C, NAME='xcloc_utils_computeDefaultXCTable')
      LOGICAL(C_BOOL), VALUE, INTENT(IN) :: ldoAutoCorrs
      INTEGER(C_INT), VALUE, INTENT(IN) :: nsignals, numbering, nwork
      INTEGER(C_INT), INTENT(OUT) :: nxcs, xcPairs(nwork), ierr
      INTEGER i, j, indx
      ierr = 0
      nxcs = 0
      IF (nsignals < 1) THEN
         WRITE(ERROR_UNIT,900)
         ierr = 1
         RETURN
      ENDIF
      IF (.NOT.ldoAutoCorrs .AND. nsignals < 2) THEN
         WRITE(ERROR_UNIT,901)
         ierr = 1
         RETURN
      ENDIF
      ! Only the superdiagonal
      IF (.NOT.ldoAutoCorrs) THEN
         nxcs = (nsignals*(nsignals - 1))/2
         IF (nwork < 0) RETURN ! space inquiry
         xcPairs(:) = 0
         IF (nwork < 2*nxcs) THEN
            WRITE(ERROR_UNIT,905) 2*nxcs
            ierr = 1
            RETURN
         ENDIF
         DO i=1,nsignals
            DO j=i+1,nsignals
               indx = nsignals*(i - 1) - (i*(i+1))/2 + j
               xcPairs(2*(indx-1)+1) = i
               xcPairs(2*(indx-1)+2) = j
            ENDDO
         ENDDO
      ! Upper triangle (including diagonal)
      ELSE
         nxcs = (nsignals*(nsignals + 1))/2
         IF (nwork < 0) RETURN ! space inquiry
         xcPairs(:) = 0
         IF (nwork < 2*nxcs) THEN
            WRITE(ERROR_UNIT,905) 2*nxcs
            ierr = 1
            RETURN
         ENDIF
         DO i=1,nsignals
            DO j=i,nsignals
               indx = nsignals*(i - 1) - ((i-1)*i)/2 + j
               xcPairs(2*(indx-1)+1) = i
               xcPairs(2*(indx-1)+2) = j
            ENDDO
         ENDDO
      ENDIF
      IF (numbering == XCLOC_C_NUMBERING) xcPairs(:) = xcPairs(:) - 1
  900 FORMAT('xcloc_utils_computeDefaultXCTable: No signals')
  901 FORMAT('xcloc_utils_computeDefaultXCTable: At least 2 signals required')
  905 FORMAT('xcloc_utils_computeDefaultXCTable: Error - nwork must be >=', I5)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Crudely attemps to evenly divide ntasks tasks between nprocs processes.
!>           This assumes that each task is approximately the same amount of work as
!>           any other task.
!>
!>    @param[in] ntasks    Number of tasks to divide.
!>    @param[in] nprocs    Number of processes to assign tasks to.
!>    @param[out] taskPtr  Given the p'th process ID this will return the chunk of
!>                         tasks in myTasks that the process is to complete - i.e.,
!>                         myTasks(taskPtr(p+1):taskPtr:p+2)-1).
!>    @param[out] myTasks  This is a map from the it'th task to the process ID.
!>    @param[out] ierr     0 indicates success.
!>
      SUBROUTINE xcloc_utils_partitionTasks(ntasks, nprocs,          &
                                            taskPtr, myTasks, ierr)  &
      BIND(C, NAME='xcloc_utils_partitionTasks')
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(IN) :: ntasks, nprocs
      INTEGER(C_INT), INTENT(OUT) :: myTasks(ntasks), taskPtr(nprocs+1), ierr
      INTEGER, ALLOCATABLE :: taskCtr(:)
      DOUBLE PRECISION dPart
      INTEGER i, ip, isum, low, high
      ierr = 0
      IF (ntasks < 1 .OR. nprocs < 1) THEN
         ierr = 1
      ENDIF
      ! Initialize
      myTasks(:) =-1 ! ID of my task
      ALLOCATE(taskCtr(nprocs)); taskCtr(:) = 0 ! Number of tasks process must do
      dPart = DBLE(ntasks)/DBLE(nprocs)
      ! Try to evenly divide the number of takss 
      DO ip=1,nprocs
         low  = MAX(1,      INT(DBLE(ip-1)*dPart + 0.5d0))
         high = MIN(ntasks, INT(DBLE(ip)*dPart   + 0.5d0))
         IF (ip == 1) low = 1 
         IF (ip == nprocs) high = ntasks + 1 
         DO i=1,ntasks
            IF (i >= low .AND. i < high) THEN
               myTasks(i) = ip - 1
               taskCtr(ip) = taskCtr(ip) + 1 
            ENDIF
         ENDDO
      ENDDO
      ! Verify
      DO i=1,ntasks
         IF (myTasks(i) < 0) THEN
            WRITE(ERROR_UNIT,900) i
            ierr = 1 
            RETURN
         ENDIF
      ENDDO
      ! Make a map from the ip'th process to the start index of myTasks
      isum = 0
      taskPtr(1) = 1
      DO ip=1,nprocs
         isum = isum + taskCtr(ip)
         taskPtr(ip+1) = isum + 1
      ENDDO
      IF (ALLOCATED(taskCtr)) DEALLOCATE(taskCtr)
      IF (isum /= ntasks) THEN
         WRITE(ERROR_UNIT,905) isum, ntasks 
         ierr = 1
      ENDIF
  900 FORMAT('xcloc_utils_partitionTasks: Failed to initialize element', I4) 
  905 FORMAT('xcloc_utils_partitionTasks: Counted', I4, ' ntasks; need ', I6, ' tasks')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Creates a list of the unique elements in a list.
!>    @param[in] n            Number of points in the input list.
!>    @param[in] list         List of numbers of which to find unique values.
!>    @param[out] nUnique     Number of unique elements in list.
!>    @param[out] listUnique  The unique elements of list in ascending order. 
!>                            This should have dimension of at  
!>    @param[out] ierr        0 indicates success.
      SUBROUTINE xcloc_utils_unique32s(n, list, nUnique, listUnique, ierr)
      INTEGER, VALUE, INTENT(IN ) :: n
      INTEGER, INTENT(IN) :: list(n)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: listUnique(:)
      INTEGER, INTENT(OUT) :: ierr, nUnique
      INTEGER, ALLOCATABLE :: work(:), workUnique(:)
      INTEGER i, np1
      ! Initialize and error check
      nUnique = 0
      IF (ALLOCATED(listUnique)) DEALLOCATE(listUnique)
      IF (n < 1) THEN
         WRITE(ERROR_UNIT,900) 
         ierr = 1
         RETURN
      ENDIF
      ! Copy input and sort
      np1 = n + 1
      ALLOCATE(work(np1))
      work(1:n) = list(1:n)
      work(np1) = MINVAL(list) - 1 ! Will avoid a match in sorted list
      ierr = ippsSortAscend_32s_I(work, n) ! Only sort first n elements
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,901) 
         ierr = 1
         DEALLOCATE(work)
         RETURN
      ENDIF
      ! Locate the unique elements
      ALLOCATE(workUnique(n))
      nUnique = 0
      DO i=1,n
         IF (work(i+1) /= work(i)) THEN
            nUnique = nUnique + 1
            workUnique(nUnique) = work(i)
         ENDIF
      ENDDO
      ! Copy the unique elements
      ALLOCATE(listUnique(nUnique))
      listUnique(1:nUnique) = workUnique(1:nUnique) 
      ! Clean up
      DEALLOCATE(work)
      DEALLOCATE(workUnique)
  900 FORMAT('xcloc_utils_unique32s: No points in input list') 
  901 FORMAT('xcloc_utils_unique32s: Failed to sort list')
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Searches an array sorted in increasing order for a key.
!>    @param[in] n       Number of elements in values.
!>    @param[in] key     Item to search for in values.
!>    @param[in] values  Array of items sorted in increasing order.
!>    @param[out] indx   Index in array values that matches key.
!>    @param[out] ierr   0 indicates success.
      SUBROUTINE xcloc_utils_bsearch32i(n, key, values, indx, ierr)
      INTEGER, VALUE, INTENT(IN) :: key, n
      INTEGER, INTENT(IN) :: values(n)
      INTEGER, INTENT(OUT) :: indx, ierr
      INTERFACE
         SUBROUTINE xcloc_sort_bsearch32i(key, values, n, indx, ierr) &
         BIND(C, NAME='xcloc_sort_bsearch32i')
         USE ISO_C_BINDING
         INTEGER(C_INT), VALUE, INTENT(IN) :: key, n
         INTEGER(C_INT), INTENT(IN) :: values(n)
         INTEGER(C_INT), INTENT(OUT) :: indx, ierr
         END SUBROUTINE
      END INTERFACE
      ierr = 0
      CALL xcloc_sort_bsearch32i(key, values, n, indx, ierr)
      IF (ierr /= 0) THEN
         WRITE(ERROR_UNIT,900) key
         ierr = 1
      ELSE
         indx = indx + 1 ! C to F indexing
      ENDIF 
  900 FORMAT('xcloc_utils_bsearch32i: Could not find key=', I6, ' in values')
      RETURN
      END
END MODULE
