#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <mpi.h>
#include "xcloc_xcfftMPI.h"
#include <ipps.h>

#define PADDING   DALES_MEM_PADDING 
#define ALIGNMENT DALES_MEM_ALIGNMENT

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif

static int cmp_int(const void *x, const void *y) 
{
    int xx = *(int *) x;
    int yy = *(int *) y;
    if (xx < yy) return -1; 
    if (xx > yy) return  1;  
    return 0;
}

int dales_xcfftMPI_computePartition(const bool verbose,
                                    const int n, const int nprocs,
                                    int *__restrict__ myOps,
                                    int *__restrict__ nOps);
int dales_xcfftMPI_getForwardTransformMap(const int nsignals,
                                          const int ntfSignals,
                                          const int ntfSignalsLoc,
                                          const int myid,
                                          const int *__restrict__ mytf,
                                          const int *__restrict__ mytfPtr,
                                          const int *__restrict__ myxc,
                                          const int *__restrict__ xcPairs,
                                          int *__restrict__ xcPairsLoc,
                                          int *ngetFwdSignalsLoc,
                                          int **targetOrigin,
                                          int **targetOffset);

/*!
 * @brief Initializes the parallel FFT structure.
 *
 * @param[in] nptsIn       Number of points in input signals.
 * @param[in] nptsPadIn    As a tuning parameter or necessity it can be helpful
 *                         to pad the signals.  The length of the cross
 *                         correlations will be 2*nptsPad - 1 where nptsPad
 *                         is greater than or equal to npts.
 * @param[in] nsignalsIn   Number of input signals.
 *
 * @param[out] xcfftMPI    On exit this contains requisite space and
 *                         FFT wisdom to forward and inverse transform the
 *                         signals when computing the correlations.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfftMPI_initialize(const int nptsIn, const int nptsPadIn,
                              const int nsignalsIn,
                              const MPI_Comm comm, const int master,
                              struct xcfftMPI_struct *xcfftMPI)
{
    int *mytf, *myxc, *xcPairs, *xcPairsLoc,
        myid, nprocs, npts, nptsPad, nsignals, nsignalsLoc,
        ntfSignals, ntfSignalsLoc;
    int i, ierr, ierrAll, sizeOfType;
    MPI_Aint localSize;
    const bool lwantDiag = false;
    const bool verbose = true;
    ierr = 0;
    memset(xcfftMPI, 0, sizeof(struct xcfftMPI_struct));
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &nprocs);
    // TODO - determine whether or not I want to split the communicator
    xcfftMPI->root = master;
    MPI_Comm_dup(comm, &xcfftMPI->comm);
    MPI_Comm_rank(xcfftMPI->comm, &xcfftMPI->rank);
    MPI_Comm_size(xcfftMPI->comm, &xcfftMPI->nprocs);
    // Have master verify the inputs
    if (myid == master)
    {
        npts = nptsIn;
        nptsPad = nptsPadIn;
        nsignals = nsignalsIn;
        if (npts < 1 || nptsPad < npts || nsignals < 2)
        {
            if (npts < 1)
            {
                fprintf(stderr, "%s: Signal length %d must be positive\n",
                        __func__, npts);
            }
            if (nptsPad < npts)
            {
                fprintf(stderr, "%s: Pad length %d < signals length %d\n",
                        __func__, nptsPad, npts);
            }
            if (nsignals < 2)
            {
                fprintf(stderr, "%s: At least 2 signals required \n", __func__);
            }
            ierr = 1;
        }
    }
    MPI_Bcast(&ierr,     1, MPI_INT, master, comm);
    if (ierr != 0){return -1;}
    MPI_Bcast(&npts,     1, MPI_INT, master, comm);
    MPI_Bcast(&nptsPad,  1, MPI_INT, master, comm);
    MPI_Bcast(&nsignals, 1, MPI_INT, master, comm);
    // Parallelization is on inverse transforms - so divvy up the signals.
    xcfftMPI->nsignals = nsignals;
    xcfftMPI->ntfSignals = ((nsignals*(nsignals-1)))/2; 
    if (nprocs > xcfftMPI->ntfSignals)
    {
        fprintf(stderr, "%s: nprocs=%d > ntfSignals=%d not yet checked\n",
                 __func__, nprocs, xcfftMPI->ntfSignals); 
    }
    // Make a global table of stuff to do
    xcPairs = (int *) calloc((size_t) (xcfftMPI->ntfSignals*2), sizeof(int));
    ierr = dales_xcfft_computeXCTable(lwantDiag,
                                      xcfftMPI->nsignals,
                                      xcfftMPI->ntfSignals,
                                      xcPairs);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing xc table on rank=%d\n",
                __func__, myid);
        return -1;
    }
    // Take ownership of portions of the inverse transforms
    myxc = (int *) calloc((size_t) xcfftMPI->ntfSignals, sizeof(int));
    mytf = (int *) calloc((size_t) xcfftMPI->nsignals, sizeof(int));
    xcfftMPI->proc2ftPtr = (int *) calloc((size_t) nprocs + 1, sizeof(int));
    xcfftMPI->proc2xcPtr = (int *) calloc((size_t) nprocs + 1, sizeof(int));
    if (myid == master)
    {
        // Have processes compute an approximately equal number of forward
        // transforms.
        ierr = dales_xcfftMPI_computePartition(verbose,
                                               xcfftMPI->nsignals,
                                               nprocs,
                                               mytf,
                                               xcfftMPI->proc2ftPtr);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error partitioning xcs\n", __func__);
            goto PARTITION_ERROR;
        }
        // Give approximately equal work in the inverse transforms to all
        // processes on the communicator
        ierr = dales_xcfftMPI_computePartition(verbose,
                                               xcfftMPI->ntfSignals,
                                               nprocs,
                                               myxc,
                                               xcfftMPI->proc2xcPtr);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error partitioning xcs\n", __func__);
            goto PARTITION_ERROR;
        }
        // Generate a map from the forward transforms to the distributed
        // cross-correlation pairs.

    }
PARTITION_ERROR:;
    MPI_Bcast(&ierr, 1, MPI_INT, master, comm);
    if (ierr != 0){return -1;}
    MPI_Bcast(myxc, xcfftMPI->ntfSignals,  MPI_INT, master, comm); 
    MPI_Bcast(mytf, xcfftMPI->nsignals,    MPI_INT, master, comm);
    MPI_Bcast(xcfftMPI->proc2ftPtr, nprocs+1, MPI_INT, master, comm);
    MPI_Bcast(xcfftMPI->proc2xcPtr, nprocs+1, MPI_INT, master, comm);
    // Set the space of the forward transforms
    nsignalsLoc   = xcfftMPI->proc2ftPtr[myid+1] - xcfftMPI->proc2ftPtr[myid];
    ntfSignalsLoc = xcfftMPI->proc2xcPtr[myid+1] - xcfftMPI->proc2xcPtr[myid];
    //ntfSignalsLoc = nxcs[myid+1] - nxcs[myid];
    MPI_Allreduce(&ntfSignalsLoc, &ntfSignals, 1, MPI_INT, MPI_SUM, comm);
    if (ntfSignals != xcfftMPI->ntfSignals)
    {
        fprintf(stderr, "%s: Lost count of tfSignals\n", __func__);
    }
    xcfftMPI->precision = XCLOC_SINGLE_PRECISION;
    xcfftMPI->sigFwd.lxc = nptsPad*2 - 1;    // Length of cross-correlations
    xcfftMPI->sigFwd.nsignals = nsignalsLoc; // Number of signals to transform
    xcfftMPI->sigFwd.npts = npts;            // Number of samples in input signal
    xcfftMPI->lxc = xcfftMPI->sigFwd.lxc;
    // Number of FFT samples
    xcfftMPI->sigFwd.ntfPts = xcfftMPI->sigFwd.lxc/2 + 1;
    xcfftMPI->sigFwd.ntfSignals = ntfSignalsLoc; // Number of inverse transforms
    xcfftMPI->nsignalsLoc = nsignalsLoc;  // Number of local signals
    xcfftMPI->ntfSignalsLoc = ntfSignalsLoc; // Number of local XC's
    //printf("%d %d %d\n", myid, nsignalsLoc, ntfSignalsLoc);
    ierr = dales_xcfft_setSpace(&xcfftMPI->sigFwd);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error allocating space on process %d\n",
                __func__, myid);
        ierr = 1;
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, comm); 
    if (ierrAll != 0){return -1;} 
    // Now that I know the space of my forward transforms I can allocate a 
    // buffer and place it on the memory window.  The memory window will be
    // a duplicate copy of my local transform accessible from all processes
    // on the communicator.
    localSize = (MPI_Aint) (xcfftMPI->sigFwd.ftOffset
                           *xcfftMPI->sigFwd.nsignals);
    MPI_Type_size(MPI_C_COMPLEX, &sizeOfType);
    xcfftMPI->localSizeOfFTs = (size_t) localSize;
    MPI_Alloc_mem(localSize*sizeOfType, MPI_INFO_NULL, &xcfftMPI->fts);
    MPI_Win_create(xcfftMPI->fts, localSize*sizeOfType, sizeOfType,
                   MPI_INFO_NULL, comm, &xcfftMPI->tfWindow);
    memset(xcfftMPI->fts, 0.0, (size_t) (localSize*sizeOfType));
    // Make a map from my xc transforms to their sources
    xcPairsLoc = (int *) calloc((size_t) (ntfSignalsLoc*2), sizeof(int));
    ierr = dales_xcfftMPI_getForwardTransformMap(nsignals,
                                                 ntfSignals,
                                                 ntfSignalsLoc,
                                                 myid,
                                                 mytf,
                                                 xcfftMPI->proc2ftPtr,
                                                 myxc,
                                                 xcPairs,
                                                 xcPairsLoc,
                                                 &xcfftMPI->ngetFwdSignals,
                                                 &xcfftMPI->getFwdTarget,
                                                 &xcfftMPI->getFwdOffset);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error getting forward transform map\n", __func__);
        ierr = 1;
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, comm);
    if (ierrAll != 0){return -1;}
    // Now set the space for the inverse transforms
    xcfftMPI->xcInv.ntfSignals = xcfftMPI->ntfSignalsLoc;
    xcfftMPI->xcInv.nsignals = xcfftMPI->ngetFwdSignals;
    xcfftMPI->xcInv.ntfPts = xcfftMPI->sigFwd.ntfPts;
    xcfftMPI->xcInv.npts = xcfftMPI->sigFwd.npts;
    xcfftMPI->xcInv.lxc = xcfftMPI->sigFwd.lxc;
    xcfftMPI->xcInv.xcPairs = xcPairsLoc;
    ierr = dales_xcfft_setSpace(&xcfftMPI->xcInv);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error allocating xcInv space on process %d\n",
                __func__, myid);
        ierr = 1;
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, comm);
    if (ierrAll != 0){return -1;} 
    // Create the descriptors
    if (myid == master)
    {
        fprintf(stdout, "%s: Making descriptors...\n", __func__);
    }
    ierr = dales_xcfft_makeDftiDescriptors(&xcfftMPI->sigFwd);
    ierr = xcloc_xcfft_makeFFTWDescriptors(&xcfftMPI->sigFwd);
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, comm);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error making forward descriptor on %d\n",
                __func__, myid);
        return -1;
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, comm);
    if (ierrAll != 0){return -1;}
    ierr = dales_xcfft_makeDftiDescriptors(&xcfftMPI->xcInv);
    ierr = xcloc_xcfft_makeFFTWDescriptors(&xcfftMPI->xcInv);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error making inverse descriptor on %d\n",
                __func__, myid);
        return -1;
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, comm);
    if (ierrAll != 0){return -1;}

    //printf("%d %d %d %d %d\n",  myid,
    //       xcfftMPI->xcInv.ntfSignals, xcfftMPI->xcInv.nsignals,
    //       xcfftMPI->xcInv.npts, xcfftMPI->xcInv.lxc);
    // Make a global table describing who owns the transform pairs and the pairs
    xcfftMPI->xcrPairs = (int *)
                         calloc((size_t) (3*xcfftMPI->ntfSignals), sizeof(int));
    for (i=0; i<xcfftMPI->ntfSignals; i++)
    {
        xcfftMPI->xcrPairs[3*i]   = xcPairs[2*i]; 
        xcfftMPI->xcrPairs[3*i+1] = xcPairs[2*i+1]; 
        xcfftMPI->xcrPairs[3*i+2] = myxc[i]; 
    }
    // Free workspace    
    free(xcPairs);
    free(myxc);
    free(mytf);
    //free(nxcs);
    //nsignals/nprocs;
    return 0;
}
//============================================================================//
/*!
 * @brief Gathers the forward transforms onto the cross-correlation Fourier
 *        transform structure.
 *
 * @param[in,out] xcfftMPI   On input contains the communication pattern
 *                           for the RMA window and, if applicable, the
 *                           transformed signals. \n
 *                           On exit, all the transformed signals for this
 *                           process have been obtained from the other
 *                           processes participating in the RMA window.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under Apache 2.
 *
 */
int xcloc_xcfftMPI_getForwardTransforms(struct xcfftMPI_struct *xcfftMPI)
{
    MPI_Aint targetDisp;
    Ipp64fc *pSrc64, *pDst64;
    Ipp32fc *pSrc32, *pDst32;
    int ftOffset, is, originCount, targetCount, targetRank;
    ftOffset = xcfftMPI->sigFwd.ftOffset;
    originCount = xcfftMPI->sigFwd.ntfPts;
    targetCount = xcfftMPI->sigFwd.ntfPts;
    //fts = (float complex *) calloc((xcfftMPI->ngetFwdSignals*ftOffset), sizeof(float complex));
    MPI_Win_fence(MPI_MODE_NOPRECEDE, xcfftMPI->tfWindow);
    for (is=0; is<xcfftMPI->ngetFwdSignals; is++)
    {
        targetRank = xcfftMPI->getFwdTarget[is];
        targetDisp = (MPI_Aint) (xcfftMPI->getFwdOffset[is]*ftOffset);
//printf("%d %d %d %d %d\n", is, ftOffset, originCount, targetRank, targetDisp);
        if (xcfftMPI->precision == XCLOC_SINGLE_PRECISION)
        {
            // Have the owner get its own data to lower contention
            if (targetRank == xcfftMPI->rank)
            {
                pSrc32 = (Ipp32fc *) &xcfftMPI->sigFwd.fts[targetDisp]; 
                pDst32 = (Ipp32fc *) &xcfftMPI->xcInv.fts[is*ftOffset];
                ippsCopy_32fc(pSrc32, pDst32, targetCount);
            }
            else
            {
                MPI_Get(&xcfftMPI->xcInv.fts[is*ftOffset], originCount,
                        MPI_C_COMPLEX, targetRank, targetDisp,
                        targetCount, MPI_C_COMPLEX, xcfftMPI->tfWindow);
            }
        }
        else 
        {
            if (targetRank == xcfftMPI->rank)
            {
                pSrc64 = (Ipp64fc *) &xcfftMPI->sigFwd.fts[targetDisp];
                pDst64 = (Ipp64fc *) &xcfftMPI->xcInv.fts[is*ftOffset];
                ippsCopy_64fc(pSrc64, pDst64, targetCount);
            }
            else
            {
                MPI_Get(&xcfftMPI->xcInv.fts[is*ftOffset], originCount,
                        MPI_C_DOUBLE_COMPLEX, targetRank, targetDisp,
                        targetCount, MPI_C_DOUBLE_COMPLEX, xcfftMPI->tfWindow);
            }
        }
    }
    MPI_Win_fence(MPI_MODE_NOSTORE + MPI_MODE_NOPUT + MPI_MODE_NOSUCCEED,
                  xcfftMPI->tfWindow);
    return 0;
}
//============================================================================//
/*!
 * @brief Creates the map so that this processes' requisite forward transforms
 *        may be plucked from the forward transform RMA memory window.
 *        Note, that this can be improved so that long continuous chunks
 *        of memory are copied from the origin process.
 *
 * @param[in] nsignals        Number of signals.
 * @param[in] ntfSignals      Number of transform signals.  This is a
 *                            triangular computed from nsignals.
 * @param[in] ntfSignalsLoc   Number of local transform signals - this will
 *                            be the number of cross-correlated signals for
 *                            this process.
 * @param[in] myid            Processes' rank on the communicator.
 * @param[in] mytf            Maps from the is'th global signal to the owner
 *                            of the input transform.  This is an array of
 *                            dimension [nsignals].
 * @param[in] mytfPtr         This maps from the signal'th signal to the
 *                            start index of input data transforms that
 *                            the ip'th process is reponsible for.  This
 *                            is an array of dimension [nprocs+1]. 
 * @param[in] myxc            This maps from the ixc'th global cross-correlation
 *                            to the process rank responsible for computing
 *                            that cross-correlation.  This is an array of
 *                            dimension [ntfSignals].
 * @param[in] xcPairs         This maps from the ixc'th global cross-correlation
 *                            to the global signal IDs that constitute the
 *                            cross-correlation pair.  This is an array of
 *                            dimension [2 x ntfSignals] with leading dimension
 *                            2.
 *
 * @param[out] xcPairsLoc     For the ixcLoc'th local cross-correlation this
 *                            returns the local signal IDs that constitute the
 *                            cross-correlation pair.  This is an array of
 *                            dimension [2 x ntfSignalsLoc] with leading 
 *                            dimension 2.
 * @param[out] nFwdSignalsLoc This is the number of forward transform signals
 *                            required for this process so that it may compute
 *                            all of its cross-correlations.
 * @param[out] xcPairsLoc     This maps from the local
 * @param[out] targetProcess  This maps from the isLoc'th local signal
 *                            to the process ID holding the Fourier transform
 *                            of the input signal.  On successful exit this will
 *                            be an array of dimension [nFwdSignalsLoc].
 * @param[out] targetOffset   This maps from the isLoc'th local signal
 *                            to the offset of the processID so that the 
 *                            isLoc'th forward transform can be obtained the
 *                            RMA memory window.  On successful exit this will
 *                            be an array of dimension [nFwdSignalsLoc].
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under Apache 2.
 *
 */
int dales_xcfftMPI_getForwardTransformMap(const int nsignals,
                                          const int ntfSignals,
                                          const int ntfSignalsLoc,
                                          const int myid,
                                          const int *__restrict__ mytf,
                                          const int *__restrict__ mytfPtr,
                                          const int *__restrict__ myxc,
                                          const int *__restrict__ xcPairs,
                                          int *__restrict__ xcPairsLoc,
                                          int *nFwdSignalsLoc,
                                          int **targetOrigin,
                                          int **targetOffset)
{
    int *work, *item, *origin, *offset,
        i, i0, i1, is, isLoc, ixc, jxc, locSignal0, locSignal1, nSignalsLoc;
    work = (int *) calloc((size_t) (nsignals + 1), sizeof(int));
    memset(work, INT_MAX, (size_t) (nsignals + 1)*sizeof(int));
    // Create a list of origins
    nSignalsLoc = 0;
    for (ixc=0; ixc<ntfSignals; ixc++)
    {
        if (myxc[ixc] == myid)
        {
            i0 = xcPairs[2*ixc];
            i1 = xcPairs[2*ixc+1];
            // Put this in the list
            for (i=0; i<nsignals; i++)
            {
                if (work[i] == i0){break;}
                if (work[i] ==-1)
                {
                    work[i] = i0;
                    nSignalsLoc = nSignalsLoc + 1;
                    break;
                }
            }
            for (i=0; i<nsignals; i++)
            {
                if (work[i] == i1){break;}
                if (work[i] ==-1)
                {
                    work[i] = i1;
                    nSignalsLoc = nSignalsLoc + 1;
                    break;
                }
            }
        }
    }
    // Sort the list
    qsort((void *) work, (size_t) nSignalsLoc, sizeof(int), cmp_int);
    // Map from the global signals to the local signals
    jxc = 0;
    for (ixc=0; ixc<ntfSignals; ixc++)
    {
        if (myxc[ixc] == myid)
        {
            i0 = xcPairs[2*ixc];
            i1 = xcPairs[2*ixc+1];
            item = (int *) bsearch((const void *) &i0, (const void *) work,
                                   nSignalsLoc, sizeof(int),
                                   cmp_int);
            if (item == NULL)
            {
                fprintf(stderr, "%s: Couldn't find %d in work\n", __func__, i0);
                return -1;
            }
            locSignal0 = item - work;
            item = (int *) bsearch((const void *) &i1, (const void *) work,
                                   nSignalsLoc, sizeof(int),
                                   cmp_int);
            if (item == NULL)
            {   
                fprintf(stderr, "%s: Couldn't find %d in work\n", __func__, i1);
                return -1;
            }
            locSignal1 = item - work;
            xcPairsLoc[2*jxc]   = locSignal0;
            xcPairsLoc[2*jxc+1] = locSignal1;
            jxc = jxc + 1;
        }
    }
    if (jxc != ntfSignalsLoc)
    {
        fprintf(stderr, "%s: Missed some local transforms\n", __func__);
        return -1;
    }
    // Compute the offset in the root process' memory space
    origin = (int *) calloc((size_t) nSignalsLoc, sizeof(int));
    offset = (int *) calloc((size_t) nSignalsLoc, sizeof(int));
    memset(origin, -1, (size_t) nSignalsLoc*sizeof(int));
    memset(origin, -1, (size_t) nSignalsLoc*sizeof(int));
    for (isLoc=0; isLoc<nSignalsLoc; isLoc++)
    {
        is = work[isLoc];
        origin[isLoc] = mytf[is];
        offset[isLoc] = is - mytfPtr[mytf[is]];
        //if (myid == 0){printf("%d %d %d\n", mytf[is], origin[isLoc], offset[isLoc]);}
    }
    *nFwdSignalsLoc = nSignalsLoc; 
    *targetOrigin = origin;
    *targetOffset = offset;
    free(work);
    return 0;
}
//============================================================================//
/*!
 * @brief Gathers all the forward transforms onto all processes.
 *
 * @param[in,out] xcfftMPI   On input this contains the local forward Fourier
 *                           transforms on xcfftMPI->xcfft. \n
 *                           On output all the Fourier transforms have been
 *                           gathered into xcfftMPI->fts.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *                          
 */
int dales_xcfftMPI_allGatherFFTs(struct xcfftMPI_struct *xcfftMPI)
{
    void *sendBuf, *recvBuf;
    int *displs, ftOffset, i, sendCount;
    ftOffset = xcfftMPI->sigFwd.ftOffset;
    sendBuf   = (void *) xcfftMPI->sigFwd.fts; 
    sendCount = xcfftMPI->sigFwd.nsignals*ftOffset;
    recvBuf   = (void *) xcfftMPI->fts; 
    displs = (int *) calloc((size_t) xcfftMPI->nprocs, sizeof(int));
    for (i=1; i<xcfftMPI->nprocs; i++)
    {
        displs[i] = (displs[i-1] + xcfftMPI->proc2ftPtr[i])*ftOffset;
    }
    MPI_Allgatherv(sendBuf, sendCount,  MPI_C_COMPLEX,
                   recvBuf, xcfftMPI->proc2ftPtr, displs,
                   MPI_C_COMPLEX, xcfftMPI->comm);
    free(displs);
    return 0; 
}
//============================================================================//
/*!
 * @brief Gathers the cross-correlations onto the root process.
 *
 * @param[in] root        Root process on which to collect signals.
 * @param[in] ntfSignals  Number of cross-correlations.  This must be defined 
 *                        on the root process and must match xcfftMPI.ntfSignals. 
 * @param[in] ldxc        Leading dimension x.  This must be >= lxc and defined
 *                        on the root process.
 * @param[in] lxc         Length of the cross-correlations.  This must equal
 *                        xcfftMPI.lxc and must be defined on the root process.
 * @param[in] recvType    Precision, MPI_FLOAT or MPI_DOUBLE, of x.
 * @param[in] xcfftMPI    Structure containing the distributed time-domain
 *                        cross-correlograms.
 *
 * @param[out] x          Contains the cross-correlograms.  This is an array
 *                        of dimension [ldxc x ntfSignals] and need only be
 *                        defined on the root process.  It's precision is
 *                        given by recvType.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 * 
 */
int xcloc_xcfftMPI_gatherXCs(const int root, const int ntfSignals,
                             const int ldxc, const int lxc,
                             const MPI_Datatype recvType,
                             const struct xcfftMPI_struct xcfftMPI,
                             void *__restrict__ x)
{
    void *xloc = NULL;
    int *displs = NULL;
    int *recvCounts = NULL;
    int i, ierr, ierrAll, ldxcUse, lxcUse, ntfSignalsLoc, sendCount;
    size_t nbytes;
    // Verify the root makes sense
    ierr = 0;
    if (recvType != MPI_DOUBLE && recvType != MPI_FLOAT)
    {
        fprintf(stderr, "%s: Can't gather this datatype\n", __func__);
        return -1;
    }
    if (root < 0 || root > xcfftMPI.nprocs - 1)
    {
        if (xcfftMPI.rank == 0)
        {
            fprintf(stderr, "%s: Error root=%d not on communicator\n",
                    __func__, root);
        }
        return -1;
    }
    if (root == xcfftMPI.rank)
    {
        if (ntfSignals != xcfftMPI.ntfSignals)
        {
            fprintf(stderr,
                    "%s: ntfSignals=%d != xcfftMPI.ntfSignals=%d on rank %d\n",
                    __func__, ntfSignals, xcfftMPI.ntfSignals, xcfftMPI.rank);
            ierr = 1;
        }
        if (lxc != xcfftMPI.lxc)
        {
            fprintf(stderr, "%s: lxc=%d != xcfftMPI.lxc=%d on rank %d\n",
                    __func__, lxc, xcfftMPI.lxc, xcfftMPI.rank);
            ierr = 1;
        }
        ldxcUse = ldxc;
        lxcUse = lxc;
    }   
    MPI_Bcast(&ierr, 1, MPI_INT, root, xcfftMPI.comm);
    // Make sure the targets know what they're sending w.r.t. sizes
    MPI_Bcast(&ldxcUse, 1, MPI_INT, root, xcfftMPI.comm);
    MPI_Bcast(&lxcUse,  1, MPI_INT, root, xcfftMPI.comm);
    // Figure out how much data I'm going to receive
    sendCount = xcfftMPI.ntfSignalsLoc*ldxcUse;
    if (recvType == MPI_DOUBLE)
    {
        nbytes = MAX((size_t) sendCount, 1)*sizeof(double);
        xloc = calloc(1, nbytes);
    }
    else
    {
        nbytes = MAX((size_t) sendCount, 1)*sizeof(float);
        xloc = calloc(1, nbytes);
    }
    // Have the root identify how much data it must receive from the
    // other processes
    if (root == xcfftMPI.rank)
    {
        recvCounts = (int *) calloc((size_t) xcfftMPI.nprocs, sizeof(int));
        displs = (int *) calloc((size_t) xcfftMPI.nprocs, sizeof(int));
        for (i=0; i<xcfftMPI.nprocs; i++) 
        {
            ntfSignalsLoc = xcfftMPI.proc2xcPtr[i+1] - xcfftMPI.proc2xcPtr[i];
            recvCounts[i] = ntfSignalsLoc*ldxcUse;
            displs[i] = xcfftMPI.proc2ftPtr[i]*ldxcUse;
        }
    }
    // Put the cross-correlations onto my local workspace
    ierr = 0;
    if (xcfftMPI.nsignalsLoc > 0)
    {
        if (recvType == MPI_DOUBLE)
        {
            ierr = xcloc_xcfft_getAllXCData(xcfftMPI.ntfSignalsLoc, ldxcUse,
                                            lxcUse, XCLOC_DOUBLE_PRECISION,
                                            xcfftMPI.xcInv, xloc);
        }
        else
        {
            ierr = xcloc_xcfft_getAllXCData(xcfftMPI.ntfSignalsLoc, ldxcUse,
                                            lxcUse, XCLOC_SINGLE_PRECISION,
                                            xcfftMPI.xcInv, xloc);
        }
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error setting data on process %d\n",
                     __func__, xcfftMPI.rank);
            ierr = 1;
        }
    }
    // Gather the data
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, xcfftMPI.comm);
    if (ierrAll == 0)
    {
        MPI_Gatherv(xloc, sendCount, recvType,
                    x, recvCounts, displs,
                    recvType, root, xcfftMPI.comm);
    }
    // Free workspace
    if (xloc != NULL){free(xloc);}
    if (recvCounts != NULL){free(recvCounts);}
    if (displs != NULL){free(displs);}
    return ierrAll;
}
//============================================================================//
/*!
 * @brief Scatters the input signals on the root process to the other processes
 *        on the xcfftMPI communicator.
 *
 * @param[in] root          ID of process with input signals.
 * @param[in] nsignals      Number of input signals.
 * @param[in] lds           Leading dimension of signal matrix.  This must
 *                          greater than or equal to npts.
 * @param[in] npts          Number of points in each signal.
 * @param[in] sendType      MPI data type defining data to scatter.
 *                          This can be MPI_FLOAT or MPI_DOUBLE.
 * @param[in] x             Input signals.  This is an array of dimension 
 *                          [lds x nsignals] with leading dimension lds
 *
 * @param[in,out] xcfftMPI  On input contains the initialized xcfftMPI 
 *                          structure. \n
 *                          On exit the input data has been scattered to 
 *                          the local xcfftMPI->sigFwd structures.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under Apache 2.
 *
 */
int xcloc_xcfftMPI_scatterData(
    const int root, const int nsignals,
    const int lds, const int npts, 
    const MPI_Datatype sendType,
    const void *__restrict__ x,
    struct xcfftMPI_struct *xcfftMPI)
{
    double *xloc64 = NULL;
    float *xloc32 = NULL;
    void *xloc = NULL;
    int *displs = NULL; 
    int *sendCounts = NULL;
    int i, ierr, ldsUse, nptsUse, nSignalsLoc, recvcount;
    size_t nbytes;
    // Verify the root makes sense
    ierr = 0;
    if (sendType != MPI_DOUBLE && sendType != MPI_FLOAT)
    {
        fprintf(stderr, "%s: Can't scatter this datatype\n", __func__);
        return -1;
    }
    if (root < 0 || root > xcfftMPI->nprocs - 1)
    {
        if (xcfftMPI->rank == 0)
        {
            fprintf(stderr, "%s: Error root=%d not on communicator\n",
                    __func__, root);
        }
        return -1;
    }
    if (root == xcfftMPI->rank)
    {
        if (nsignals != xcfftMPI->nsignals)
        {
            fprintf(stderr,
                    "%s: nsignals=%d != xcfftMPI->nsignals=%d on rank %d\n",
                    __func__, nsignals, xcfftMPI->nsignals, xcfftMPI->rank);
            ierr = 1;
        }
        ldsUse = lds;
        nptsUse = npts;
    }
    MPI_Bcast(&ierr, 1, MPI_INT, root, xcfftMPI->comm);
    if (ierr != 0){return -1;}
    // Make sure the targets know what's coming their way w.r.t. sizes
    MPI_Bcast(&ldsUse,  1, MPI_INT, root, xcfftMPI->comm);
    MPI_Bcast(&nptsUse, 1, MPI_INT, root, xcfftMPI->comm);
    // Figure out how much data I'm going to receive
    recvcount = xcfftMPI->nsignalsLoc*ldsUse;
    if (sendType == MPI_DOUBLE)
    {
        nbytes = MAX((size_t) recvcount, 1)*sizeof(double);
        xloc = calloc(1, nbytes);
    }
    else
    {
        nbytes = MAX((size_t) recvcount, 1)*sizeof(float);
        xloc = calloc(1, nbytes);
    }
    // Have the root identify how much data it must send to other processes
    if (root == xcfftMPI->rank)
    {
        sendCounts = (int *) calloc((size_t) xcfftMPI->nprocs, sizeof(int));
        displs = (int *) calloc((size_t) xcfftMPI->nprocs, sizeof(int));
        for (i=0; i<xcfftMPI->nprocs; i++)
        {
            nSignalsLoc = xcfftMPI->proc2ftPtr[i+1] - xcfftMPI->proc2ftPtr[i];
            sendCounts[i] = nSignalsLoc*ldsUse;
            displs[i] = xcfftMPI->proc2ftPtr[i]*ldsUse;
        }
    }
    // Send the requisite subset of data to each process
    MPI_Scatterv(x, sendCounts, displs, sendType,
                 xloc, recvcount, sendType, root, xcfftMPI->comm);
    // Now put the signals onto my local workspace
    if (xcfftMPI->nsignalsLoc > 0)
    {
        if (sendType == MPI_DOUBLE)
        {
            xloc64 = (double *) xloc; 
            ierr = dales_xcfft_setAllData64f(xcfftMPI->nsignalsLoc, ldsUse,
                                             nptsUse, xloc64,
                                             &xcfftMPI->sigFwd);
        }
        else
        {
            xloc32 = (float *) xloc;
            ierr = dales_xcfft_setAllData32f(xcfftMPI->nsignalsLoc, ldsUse,
                                             nptsUse, xloc32,
                                             &xcfftMPI->sigFwd);
        }
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error setting data on process %d\n",
                     __func__, xcfftMPI->rank);
        }
    }
    // Free workspace
    xloc64 = NULL;
    xloc32 = NULL;
    if (xloc != NULL){free(xloc);}
    if (sendCounts != NULL){free(sendCounts);}
    if (displs != NULL){free(displs);}
    return ierr;
}
//============================================================================//
/*!
 * @brief Crudely divides the work (load balances) n operations amount n
 *        processes.  There is no idea of weighting to reflect communication
 *        expense.
 *
 * @param[in] verbose  If true then print a summary of each process' work.
 * @param[in] n        Number of operations to split.
 * @param[in] nprocs   Number of processes.
 *
 * @param[out] myOps   Maps from the i'th operation to the process ID that
 *                     responsible for performing this operation.  This is
 *                     an array of dimension [n].
 * @param[out] opsPtr  Points from ip'th process to start index of myOps.
 *                     The ip'th process has opsPtr[ip+1] - opsPtr[ip] 
 *                     operations to perform.  This is an array of dimension
 *                     [nprocs + 1].
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 * 
 */
int dales_xcfftMPI_computePartition(const bool verbose,
                                    const int n, const int nprocs,
                                    int *__restrict__ myOps,
                                    int *__restrict__ opsPtr)
{
    double dPart;
    int high, i, ip, isum, low, *nOps;
    if (nprocs < 1 || n < 1 || myOps == NULL || opsPtr == NULL)
    {
        if (nprocs < 1){fprintf(stderr, "%s: No processes\n", __func__);}
        if (n < 1){fprintf(stderr, "%s: No operations\n", __func__);}
        if (myOps == NULL){fprintf(stderr, "%s: myOps is NULL\n", __func__);}
        if (opsPtr == NULL){fprintf(stderr, "%s: opsPtr is NULL\n", __func__);}
        return -1;
    }
    nOps = (int *) calloc((size_t) nprocs, sizeof(int));
    // Compute spacing of partitions
    dPart = (double) n/(double) nprocs; 
    // Initialize
    for (i=0; i<n; i++){myOps[i] =-1;}
    for (ip=0; ip<nprocs; ip++){nOps[ip] = 0;}
    // Loop on processes
    for (ip=0; ip<nprocs; ip++)
    {
        low  = (int) round((double) ip*dPart);
        high = (int) round((double) (ip + 1)*dPart);
        if (ip == 0){low = 0;} 
        if (ip == nprocs - 1){high = n;} 
        for (i=0; i<n; i++)
        {
            if (i >= low && i < high)
            {
               //printf("%d %d %d\n", i, ip, mySignal[i]);
               myOps[i] = ip;
               nOps[ip] = nOps[ip] + 1;
            }
        }
    }
    // Verify
    for (i=0; i<n; i++)
    {
        if (myOps[i] < 0)
        {
            fprintf(stderr, "%s: Failed to initialize element %d\n",
                    __func__, i);
            return -1;
        }
    }
    isum = 0;
    opsPtr[0] = 0;
    for (ip=0; ip<nprocs; ip++)
    {
        isum = isum + nOps[ip];
        opsPtr[ip+1] = isum;
    }
    if (isum != n) 
    {
        fprintf(stderr, "%s: Lost count in bins\n", __func__);
        return -1;
    }
    // Let user review how the partition went
    if (verbose)
    {
        fprintf(stdout, "%s: Partition summary:\n", __func__);
        for (ip=0; ip<nprocs; ip++)
        {
             fprintf(stdout, "Process %d has %d operations\n", ip, nOps[ip]);
        }
    }
    free(nOps);
    return 0;
}
//============================================================================//
/*!
 * @brief Deallocates the memory on the xcfftMPI structure, frees the 
 *        communicators, deletes the RMA windows, and sets all variables to
 *        0.
 *
 * @param[in,out] xcfftMPI   On input this contains and initialized xfftMPI
 *                           structure. \n
 *                           On exit all memory has been freed.
 *
 * @result 0 indicates success.
 * 
 * @author Ben Baker
 *
 * @copyright Ben Baker distirbuted under the MIT license.
 *
 */
int xcloc_xcfftMPI_finalize(struct xcfftMPI_struct *xcfftMPI)
{
    int ierr;
    // Release the MPI windows
    MPI_Win_free(&xcfftMPI->tfWindow);
    MPI_Free_mem(xcfftMPI->fts);
    // Finalize the transforms
    ierr = dales_xcfft_finalize(&xcfftMPI->sigFwd);
    if (ierr != 0)
    {   
        fprintf(stderr, "%s: Error finalizing xcInv on process %d\n",
                __func__, xcfftMPI->rank);
    }
    ierr = dales_xcfft_finalize(&xcfftMPI->xcInv);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error finalizing xcInv on process %d\n",
                __func__, xcfftMPI->rank);
    }
    // Release the memory on the pointers
    if (xcfftMPI->proc2ftPtr != NULL){free(xcfftMPI->proc2ftPtr);}
    if (xcfftMPI->proc2xcPtr != NULL){free(xcfftMPI->proc2xcPtr);}
    if (xcfftMPI->getFwdOffset != NULL){free(xcfftMPI->getFwdOffset);}
    if (xcfftMPI->getFwdTarget != NULL){free(xcfftMPI->getFwdTarget);}
    if (xcfftMPI->xcrPairs != NULL){free(xcfftMPI->xcrPairs);}
    // Free the communicator
    MPI_Comm_free(&xcfftMPI->comm);
    memset(xcfftMPI, 0, sizeof(struct xcfftMPI_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Transforms the input signals then distributes the requisite
 *        transforms to the other processes on the communicator.
 *
 * @param[in,out] xcfftMPI  On input contains the sigFwd transform information
 *                          and communication pattern for the fetching the
 *                          forward transforms from the other process on the
 *                          communicator. \n
 *                          On output the forward transforms have been computed
 *                          and the requisite transforms for this process has
 *                          gotten its portion of Fourier transforms so that it
 *                          may compute the cross-correlation.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfftMPI_fourierTransform(struct xcfftMPI_struct *xcfftMPI)
{
    void *ftsSrc __attribute__((aligned(ALIGNMENT)));
    void *ftsDst;
    int ierr, ierrLoc;
    ierrLoc = 0;
    // Compute the Fourier transforms
    if (xcfftMPI->sigFwd.nsignals > 0)
    {
        ierrLoc = dales_xcfft_forwardTransform(&xcfftMPI->sigFwd);
        //ierrLoc = xcloc_xcfft_forwardTransform(&xcfftMPI->sigFwd);
        if (ierrLoc != 0)
        {
            fprintf(stdout, "%s: Error computing FFTs on process %d\n",
                    __func__, xcfftMPI->rank);
        }
        else
        {
            // Copy to the transforms onto the memory window 
            ftsSrc = (void *) xcfftMPI->sigFwd.fts; 
            ftsDst = (void *) xcfftMPI->fts;
            memcpy(ftsDst, ftsSrc, xcfftMPI->localSizeOfFTs);
        }
    }
    MPI_Allreduce(&ierrLoc, &ierr, 1, MPI_INT, MPI_SUM, xcfftMPI->comm);
    if (ierr != 0){return -1;}
    // Fetch the requisite transforms
    ierr = xcloc_xcfftMPI_getForwardTransforms(xcfftMPI);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to get fwd transforms onto process %d\n",
                __func__, xcfftMPI->rank);
        return -1;
    }
    return 0;
} 
//============================================================================//
int xcloc_xcfftMPI_computePhaseCorrelation(struct xcfftMPI_struct *xcfftMPI)
{
    int ierr;
    const bool lphaseOnly = true;
    // Fourier transform and let the ranks get the data that they require
    // on the data
    ierr = xcloc_xcfftMPI_fourierTransform(xcfftMPI);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error FT'ing data on rank %d\n",
                __func__, xcfftMPI->rank);
        return -1; 
    }
    // Compute the phase correlations -> the data is now xcInv.data
    ierr = xcloc_xcfft_computeCorrelationsWithFFTData(lphaseOnly,
                                                     xcfftMPI->xcInv.ntfSignals,
                                                     xcfftMPI->xcInv.ntfPts,
                                                     xcfftMPI->xcInv.ftOffset,
                                                     xcfftMPI->xcInv.xcPairs,
                                                     xcfftMPI->xcInv.fts,
                                                     xcfftMPI->xcInv.xcfts);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing phase correlation on rank %d\n",
                __func__, xcfftMPI->rank);
    }
    // Inverse transform and finish with the time-domain cross-correlations 
    ierr = dales_xcfft_inverseTransformXC(&xcfftMPI->xcInv);
    //ierr = xcloc_xcfft_inverseTransformXC(&xcfftMPI->xcInv);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error inverse transforming on rank %d\n",
                __func__, xcfftMPI->rank);
    }
    return 0;
}
/*
int xcloc_xcfftMPI_apply(const bool lphaseOnly,
                         struct xcfftMPI_struct *xcfftMPI)
{
    int ierr;
    // Compute the Fourier transforms
    ierr = dales_xcfftMPI_fourierTransform(xcfftMPI);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing FTs on process %d\n",
                __func__, xcfftMPI->rank);
        return -1;
    }    
    // Apply the 
    //ierr = dales_xcfft_computeCorrelationsWithFFTData(
    return 0;
}
*/
