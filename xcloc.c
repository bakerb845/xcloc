#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "xcloc.h"
#include <ipps.h>
/*
#ifndef XCLOC_USE_MPI
# warning "This only works for MPI" 
#endif
*/
int xcloc_splitCommunicator(const int nprocs, const int nparts, 
                            int *__restrict__ myGroup);

/*!
 * @brief Checks the xcloc parameters.
 */
int xcloc_checkParameters(const struct xclocParms_struct xclocParms)
{
    int ierr, pMin, pMax;
    ierr = xcloc_xcfft_checkParameters(xclocParms.npts,
                                       xclocParms.nptsPad,
                                       xclocParms.nsignals);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: A problem detected with xcfft parameters\n",
                __func__);
        return -1;
    }
    if (xclocParms.dt <= 0.0)
    {
        fprintf(stderr, "%s: Error dt=%e must be positive\n",
                __func__, xclocParms.dt);
        ierr = 1;
    }
    if (xclocParms.chunkSize%DALES_MEM_ALIGNMENT != 0)
    {
        fprintf(stderr, "%s: chunkSize=%d not divisible by alignment=%d\n",
                __func__, xclocParms.chunkSize, DALES_MEM_ALIGNMENT);
        ierr = 1;
    }
    if (xclocParms.ngrd < 1)
    {
        fprintf(stderr, "%s: ngrd=%d must be positive\n", __func__,
                xclocParms.ngrd);
        ierr = 1;
    }
    if (xclocParms.signalGroup == NULL)
    {
        fprintf(stdout, "%s: Assuming all signals in group 0\n", __func__);
    }
    else
    {
        ippsMax_32s(xclocParms.signalGroup, xclocParms.nsignals, &pMax);
        ippsMin_32s(xclocParms.signalGroup, xclocParms.nsignals, &pMin);
        if (pMin != 0)
        {
            fprintf(stderr, "%s: pMin = %d != 0\n", __func__, pMin);
            ierr = 1;
        }
        fprintf(stdout, "%s: Number of signal groups: %d\n",
                __func__, pMax - pMin + 1);
    }
    return ierr;
}
//============================================================================//
int xcloc_initialize(const MPI_Comm comm,
                     const struct xclocParms_struct xclocParms,
                     struct xcloc_struct *xcloc)
{
    int color, ierr, ierrAll, ift, ig, is, key, mpiInit, pMin, pMax;
    enum xclocPrecision_enum precision;
    enum xclocAccuracy_enum accuracy;
    const int root = 0;
    ierr = 0;
    memset(xcloc, 0, sizeof(struct xcloc_struct));
xcloc->precision = XCLOC_SINGLE_PRECISION;
xcloc->accuracy = XCLOC_HIGH_ACCURACY;
    xcloc->root = root;
    xcloc->nSignalGroups = 1;
    MPI_Initialized(&mpiInit);
    if (!mpiInit)
    {
        //MPI_Init(NULL, NULL); // Would need to read environment variables
        fprintf(stderr, "%s: MPI not initialized\n", __func__);
        return -1;
    }
    MPI_Comm_dup(comm, &xcloc->globalComm); 
    MPI_Comm_rank(xcloc->globalComm, &xcloc->globalCommRank);
    MPI_Comm_size(xcloc->globalComm, &xcloc->globalCommSize);
    if (xcloc->globalCommRank == xcloc->root)
    {
        ierr = xcloc_checkParameters(xclocParms);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Invalid parameters\n", __func__);
            goto BCAST_ERROR;
        }
        xcloc->nfftProcs = xclocParms.nfftProcs;
        if (xcloc->globalCommSize%xcloc->nfftProcs != 0)
        {
            fprintf(stderr, "%s: globalCommSize=%d not divisible by fftProc=%d\n", 
                    __func__, xcloc->globalCommSize, xcloc->nfftProcs);
            ierr = 1;
            goto BCAST_ERROR;
        }
/*
        xcloc->nfftProcs = MIN(xclocParms.nfftProcs, xcloc->globalCommSize);
        xcloc->nfftProcs = MAX(1, xclocParms.nfftProcs);
        if (xcloc->nfftProcs != xclocParms.nfftProcs)
        {
            fprintf(stdout, "%s: Over-riding nfftProcs to %d\n", __func__,
                    xcloc->nfftProcs);
        }
*/
        xcloc->nmigrateProcs = xclocParms.ngridProcs;
        if (xcloc->nfftProcs*xcloc->nmigrateProcs != xcloc->globalCommSize)
        {
            fprintf(stderr,
                    "%s: nfftProcs*nmigrateProcs=%dx%d != globalCommSize=%d\n",
                    __func__, xcloc->nfftProcs, xcloc->nmigrateProcs, 
                    xcloc->globalCommSize);
            ierr = 1;
            goto BCAST_ERROR;
        }
        xcloc->dt = xclocParms.dt;
        xcloc->chunkSize = xclocParms.chunkSize; 
        xcloc->nmigrateGroups = xcloc->nfftProcs;
        xcloc->npts    = xclocParms.npts;
        xcloc->nptsPad = xclocParms.nptsPad;
        xcloc->ngrd    = xclocParms.ngrd;
        xcloc->lphaseXCs = xclocParms.lphaseXCs;
        xcloc->lxc = 2*xcloc->nptsPad - 1;
        if (xclocParms.envFIRLen > 0)
        {
            xcloc->lenvelope = true;
            xcloc->envFIRLen = xclocParms.envFIRLen;
        }
        ippsMax_32s(xclocParms.signalGroup, xclocParms.nsignals, &pMax);
        ippsMin_32s(xclocParms.signalGroup, xclocParms.nsignals, &pMin);
        xcloc->nSignalGroups = pMax - pMin + 1;
        xcloc->nsignals
          = (int *) calloc((size_t) xcloc->nSignalGroups, sizeof(int));
        for (ig=0; ig<xcloc->nSignalGroups; ig++)
        {
            for (is=0; is<xclocParms.nsignals; is++)
            {
                if (xclocParms.signalGroup[is] == ig)
                {
                    xcloc->nsignals[ig] = xcloc->nsignals[ig] + 1;
                }
            }
            if (xcloc->nsignals[ig] < 2)
            {
                fprintf(stderr,
                        "%s: Signal group %d has %d signals; 2 required\n",
                         __func__, ig, xcloc->nsignals[ig]);
                ierr = 1;
                goto BCAST_ERROR;
            }
            xcloc->nTotalSignals = xcloc->nTotalSignals + xcloc->nsignals[ig]; 
        }
    }
BCAST_ERROR:;
    MPI_Bcast(&ierr, 1, MPI_INT, xcloc->root, xcloc->globalComm);
    if (ierr != 0){return -1;}
    // Broadcast the parameters
    MPI_Bcast(&xcloc->dt, 1, MPI_DOUBLE, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->envFIRLen,      1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->chunkSize,      1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->nmigrateGroups, 1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->nmigrateProcs,  1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->nfftProcs,      1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->npts,           1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->nptsPad,        1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->lxc,            1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->ngrd,           1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->nTotalSignals,  1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->nSignalGroups,  1, MPI_INT, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->lphaseXCs, 1, MPI_C_BOOL, root, xcloc->globalComm);
    MPI_Bcast(&xcloc->lenvelope, 1, MPI_C_BOOL, root, xcloc->globalComm);
    if (xcloc->globalCommRank != xcloc->root)
    {
        xcloc->nsignals
          = (int *) calloc((size_t) xcloc->nSignalGroups, sizeof(int));
    }
    MPI_Bcast(xcloc->nsignals, xcloc->nSignalGroups, MPI_INT,
              xcloc->root, xcloc->globalComm);
    // Make the envelope structure
    if (xcloc->lenvelope)
    {
        ierr = xcloc_envelope_initialize(xcloc->envFIRLen,
                                         xcloc->precision,
                                         xcloc->accuracy,
                                         &xcloc->envelope);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error initializing envelope\n", __func__);
            ierr = 1;
        }
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, xcloc->globalComm);
        if (ierrAll != 0){return -1;}
    }
    // Make the FFT communicator
    xcloc->fftCommRank = MPI_UNDEFINED;
    color = MPI_UNDEFINED; 
    key = 0;
    for (ift=0; ift<xcloc->nfftProcs; ift++)
    {
        if (ift*xcloc->nmigrateProcs == xcloc->globalCommRank)
        {
            color = 0; // I'm in the group
            key = ift; // My FFT ID
            xcloc->ldoFFT = true;
            //printf("go with %d %d %d\n", color, key, xcloc->nfftProcs);
        }
    }
    ierr = MPI_Comm_split(xcloc->globalComm, color, key, &xcloc->fftComm);
    if (ierr != MPI_SUCCESS)
    {
        fprintf(stdout, "%s: Failed making fftComm\n", __func__);
        return -1;
    }
/*
    // Make the migration communicator
    color = MPI_UNDEFINED;
    key = 0;
    ir = 0;
    for (ift=0; ift<xcloc->nfftProcs; ift++)
    {
        for (ig=0; ig<xcloc->nmigrateProcs; ig++)
        {
            if (xcloc->globalCommRank == ir)
            {
                color = ift; // I'll get my XCs from this group
                key   = ig;  // I'll retain my rank
            } 
            ir = ir + 1;
        }
    }
    ierr = MPI_Comm_split(xcloc->globalComm, color, key, &xcloc->migrateComm);
    if (ierr != MPI_SUCCESS)
    {
        fprintf(stderr, "%s: Failed making migrateComm\n", __func__);
        return -1;
    }
    MPI_Comm_rank(xcloc->migrateComm, &xcloc->migrateCommRank);
    MPI_Comm_size(xcloc->migrateComm, &xcloc->migrateCommSize);
    // Make the global XC pairs
*/
    ierr = xcloc_makeXCPairs(xcloc);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to make XC pair on rank %d\n",
                __func__, xcloc->globalCommRank);
        return -1;
    }
    // Initialize the FFTs
    ierr = 0;
    if (xcloc->ldoFFT)
    {
        MPI_Comm_rank(xcloc->fftComm, &xcloc->fftCommRank);
        MPI_Comm_size(xcloc->fftComm, &xcloc->fftCommSize);
        if (xcloc->globalCommRank == 0)
        {
            fprintf(stdout, "%s: Initializing FFT...\n", __func__);
        }
        ierr = xcloc_xcfftMPI_initialize(xcloc->npts,
                                         xcloc->nptsPad,
                                         xcloc->nTotalSignals,
                                         xcloc->nTotalXCs,
                                         xcloc->fftComm,
                                         xcloc->root,
                                         xcloc->xcPairs,
                                         &xcloc->xcfftMPI);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to initialize fft structure\n",
                    __func__);
            ierr = 1;
        }
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, xcloc->globalComm);
    if (ierrAll != 0){return -1;} 
    // Initialize the migration communicator
    if (xcloc->globalCommRank == 0)
    {
        fprintf(stdout, "%s: Initializing migration...\n", __func__);
    }
    ierr = xcloc_migrateMPI_initialize(xcloc->globalComm,
                                       xcloc->nmigrateGroups,
                                       xcloc->nTotalSignals,
                                       xcloc->nTotalXCs,
                                       xcloc->chunkSize, 
                                       xcloc->ngrd,
                                       xcloc->lxc, xcloc->dt,
                                       xcloc->xcfftMPI.xcrPairs,
                                       &xcloc->migrateMPI);
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, xcloc->globalComm);
    if (ierrAll != 0){return -1;}
    // Set the temporary space for cross-correlations
    if (xcloc->ldoFFT)
    {
        xcloc->ntfSignalsLoc = xcloc->xcfftMPI.ntfSignalsLoc;
    }
    MPI_Bcast(&xcloc->ntfSignalsLoc, 1, MPI_INT, xcloc->root,
              xcloc->migrateMPI.intraComm);
    xcloc->leny = xcloc->lxc*xcloc->ntfSignalsLoc;
    if (xcloc->precision == XCLOC_SINGLE_PRECISION)
    {
        xcloc->y1 = aligned_alloc(XCLOC_MEM_ALIGNMENT, (size_t) xcloc->leny*sizeof(float)); 
        xcloc->y2 = aligned_alloc(XCLOC_MEM_ALIGNMENT, (size_t) xcloc->leny*sizeof(float));
    }
    else
    {
        xcloc->y1 = aligned_alloc(XCLOC_MEM_ALIGNMENT, (size_t) xcloc->leny*sizeof(double));
        xcloc->y2 = aligned_alloc(XCLOC_MEM_ALIGNMENT, (size_t) xcloc->leny*sizeof(double));
    }
    // All done
    xcloc->linit = true;
    return 0;
}
//============================================================================//
/*
int xcloc_setTravelTimeTableFromRoot(const int ngrd,
                                     const xcloc_struct *xcloc)
{

}
*/
//============================================================================//
/*!
 * @brief Computes the Fourier transforms, the envelopes, and the migration
 *        image.
 *
 * @coyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_apply(struct xcloc_struct *xcloc)
{
    float *yf1, *yf2;
    int ierr, ierrAll;
    if (!xcloc->linit)
    {
        fprintf(stderr, "%s: xcloc not initialized\n", __func__);
        return -1;
    }
    // Compute the correlograms
    if (xcloc->ldoFFT)
    {
        ierr = xcloc_xcfftMPI_computePhaseCorrelation(&xcloc->xcfftMPI);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing phase xcs on rank %d\n",
                    __func__, xcloc->globalCommRank);
            ierr = 1;
        }
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, xcloc->fftComm);
    }
    MPI_Bcast(&ierrAll, 1, MPI_INT, xcloc->root, xcloc->globalComm);
    if (ierrAll != 0){return -1;}
    // Distribute the correlograms to the process on the migration grid
    if (xcloc->precision == XCLOC_SINGLE_PRECISION)
    {
        if (xcloc->ldoFFT)
        {
            xcloc_xcfft_getAllXCData(xcloc->ntfSignalsLoc, 
                                     xcloc->lxc,
                                     xcloc->lxc,
                                     XCLOC_SINGLE_PRECISION,
                                     xcloc->xcfftMPI.xcInv,
                                     xcloc->y1);
        }
        MPI_Bcast(xcloc->y1, xcloc->leny, MPI_FLOAT,
                  xcloc->root, xcloc->migrateMPI.intraComm); 
    }
    else
    {
        if (xcloc->ldoFFT)
        {
            xcloc_xcfft_getAllXCData(xcloc->ntfSignalsLoc,
                                     xcloc->lxc,
                                     xcloc->lxc,
                                     XCLOC_DOUBLE_PRECISION,
                                     xcloc->xcfftMPI.xcInv,
                                     xcloc->y1);
        }
        MPI_Bcast(xcloc->y1, xcloc->leny, MPI_DOUBLE,
                  xcloc->root, xcloc->migrateMPI.intraComm);
    }
    // Apply the envelope
    if (xcloc->lenvelope)
    {
        float *work = aligned_alloc(XCLOC_MEM_ALIGNMENT, xcloc->leny*sizeof(float));
        memset(work, 0, xcloc->leny*sizeof(float));
/*
        double di = round((double) xcloc->ntfSignalsLoc/(double) xcloc->xcloc->nProcsPerGroup);
        for (int i=0; i<xcloc->migrateMPI.nProcsPerGroup; i++)
        { 

        }
*/
        MPI_Allreduce(work, xcloc->y2, xcloc->leny, MPI_FLOAT, MPI_SUM, 
                      xcloc->migrateMPI.intraComm);
        free(work); 
    }
    else
    {
        memcpy(xcloc->y2, xcloc->y1, xcloc->leny*sizeof(float));
    }
    xcloc_migrate_setCrossCorrelations(xcloc->lxc, xcloc->lxc,
                                       xcloc->ntfSignalsLoc,
                                       xcloc->y2,
                                       xcloc->precision,
                                       &xcloc->migrateMPI.migrate);
/*
int xcloc_migrate_setCrossCorrelation(const int lxc, const int ixc, 
                                      const void *__restrict__ xcs, 
                                      const enum xclocPrecision_enum precision,
                                      struct migrate_struct *migrate)
*/

    // Compute the migration image
    xcloc_migrateMPI_computeMigrationImage(&xcloc->migrateMPI);
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the it'th travel time table.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_setTableFromRoot(const int itIn, const int ngrdIn,
                           const enum xclocPrecision_enum precisionIn,
                           const void *__restrict__ ttimes,
                           struct xcloc_struct *xcloc)
{
    int ierr, ierrAll, it, ngrd, precIn;
    enum xclocPrecision_enum precision;
    ierr = 0;
    if (!xcloc->linit)
    {
        fprintf(stderr, "%s: xcloc not initialized\n", __func__);
        return -1;
    }
    if (xcloc->globalCommRank == xcloc->root)
    {
        it = itIn;
        ngrd = ngrdIn;
        precIn = (int) precisionIn;
    }
    MPI_Bcast(&it,     1, MPI_INT, xcloc->root, xcloc->globalComm);
    MPI_Bcast(&ngrd,   1, MPI_INT, xcloc->root, xcloc->globalComm);
    MPI_Bcast(&precIn, 1, MPI_INT, xcloc->root, xcloc->globalComm);
    precision = (enum xclocPrecision_enum) precIn;
    ierr = xcloc_migrateMPI_setTableFromRoot(it, ngrd, precision, ttimes, 
                                             &xcloc->migrateMPI);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error setting table\n", __func__);
        ierr = 1;
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, xcloc->globalComm); 
    return ierrAll;
} 
//============================================================================//
int xcloc_scatterDataFromRoot(const int nsignals,
                              const int lds, const int npts,
                              const MPI_Datatype sendType,
                              const void *__restrict__ x,
                              struct xcloc_struct *xcloc)
{
    int ierr, ierrAll;
    if (!xcloc->linit)
    {
        fprintf(stderr, "%s: xcloc not initialized\n", __func__);
        return -1;
    }
    ierr = 0;
    if (xcloc->ldoFFT)
    {
        ierr =  xcloc_xcfftMPI_scatterData(xcloc->root, nsignals,
                                           lds, npts, sendType, x,
                                           &xcloc->xcfftMPI);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error scattering data on globalRank=%d\n",
                    __func__, xcloc->globalCommRank);
        }
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, xcloc->globalComm);
    return ierrAll;
}
//============================================================================//
/*!
 * @brief Defines the cross-correlation pairs for the different signal groups.
 *
 * @param[in,out] xcloc    On input contains the global communicator.  On, the
 *                         root process the number of total signals, the number
 *                         of signal groups, the signal groups, and the number
 *                         of signals in each group is defined.
 *                         On exit, all processes have the XC pairs and the
 *                         total number of cross-correlations.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_makeXCPairs(struct xcloc_struct *xcloc)
{
    const bool lwantDiag = false;
    int *xcWork, i, ierr, indx, is, nxcLoc, nsignals,
        nsignalsLoc, nwork, offset;
    if (xcloc->globalCommRank == xcloc->root)
    {
        nsignals = xcloc->nTotalSignals;
        nwork = (nsignals*(nsignals-1))/2;
        xcWork = (int *) calloc((size_t) (2*nwork), sizeof(int));
        indx = 0;
        offset = 0;
        xcloc->nTotalXCs = 0;
        for (is=0; is<xcloc->nSignalGroups; is++)
        {
            nsignalsLoc = xcloc->nsignals[is];
            nxcLoc = (nsignalsLoc*(nsignalsLoc-1))/2;
            ierr = xcloc_xcfft_computeXCTable(lwantDiag, nsignalsLoc,
                                              nxcLoc, &xcWork[indx]);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error creating %d'th xc table\n",
                        __func__, is);
                goto BCAST_ERROR;
            }
            // Shift this table
            for (i=0; i<2*nxcLoc; i++)
            {
                xcWork[indx+i] = xcWork[indx+i] + offset;
            }
            // Update the pointer indices
            offset = offset + xcloc->nsignals[is];
            indx = indx + 2*nxcLoc;
            xcloc->nTotalXCs = xcloc->nTotalXCs + nxcLoc;
        }
        if (xcloc->nTotalXCs < 1)
        {
            fprintf(stderr, "%s: No XCs!\n", __func__);
            ierr = 1;
            goto BCAST_ERROR;
        }
        xcloc->xcPairs
           = (int *) calloc((size_t) xcloc->nTotalXCs*2, sizeof(int));
        ippsCopy_32s(xcWork, xcloc->xcPairs, 2*xcloc->nTotalXCs);
        free(xcWork);
        /*
        for (int i=0; i<xcloc->nTotalXCs; i++)
        {
            printf("%d %d\n", xcloc->xcPairs[2*i], xcloc->xcPairs[2*i+1]);
        }
        */
    }
BCAST_ERROR:;
    MPI_Bcast(&ierr, 1, MPI_INT, xcloc->root, xcloc->globalComm);
    if (ierr != 0){return -1;}
    MPI_Bcast(&xcloc->nTotalXCs, 1, MPI_INT, xcloc->root, xcloc->globalComm);
    if (xcloc->globalCommRank != xcloc->root)
    {
        xcloc->xcPairs
           = (int *) calloc((size_t) xcloc->nTotalXCs*2, sizeof(int));
    }
    MPI_Bcast(xcloc->xcPairs, 2*xcloc->nTotalXCs, MPI_INT,
              xcloc->root, xcloc->globalComm);
    return 0;
}
//============================================================================//
/*!
 * @brief Releases the memory on the xcloc structure.
 *
 * @param[in,out] xcloc   On input this is an initialized xcloc structure. \n
 *                        On exit all memory has been freed, all communicators
 *                        destroyed, and all variables set to 0.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_finalize(struct xcloc_struct *xcloc)
{
    if (!xcloc->linit){return 0;} // Never initialized
    if (xcloc->signalGroup != NULL){free(xcloc->signalGroup);}
    if (xcloc->xcPairs     != NULL){free(xcloc->xcPairs);}
    if (xcloc->xcPairsLoc  != NULL){free(xcloc->xcPairsLoc);}
    if (xcloc->y1          != NULL){free(xcloc->y1);}
    if (xcloc->y2          != NULL){free(xcloc->y2);}
    if (xcloc->ldoFFT)
    {
        xcloc_xcfftMPI_finalize(&xcloc->xcfftMPI);
        MPI_Comm_free(&xcloc->fftComm);
    }
    xcloc_migrateMPI_finalize(&xcloc->migrateMPI); 
    //MPI_Comm_free(&xcloc->migrateComm);
    MPI_Comm_free(&xcloc->globalComm);
    memset(xcloc, 0, sizeof(struct xcloc_struct));
    return 0;
}
