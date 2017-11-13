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
    int ierr, ierrAll, ig, is, groupMax, groupMin, mpiInit, pMin, pMax;
    ierr = 0;
    memset(xcloc, 0, sizeof(struct xcloc_struct));
    xcloc->root = 0;
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
        xcloc->nfftProcs = MIN(xclocParms.nfftProcs, xcloc->globalCommSize);
        xcloc->nfftProcs = MAX(1, xclocParms.nfftProcs);
        if (xcloc->nfftProcs != xclocParms.nfftProcs)
        {
            fprintf(stdout, "%s: Over-riding nfftProcs to %d\n", __func__,
                    xcloc->nfftProcs);
        }
        xcloc->npts    = xclocParms.npts;
        xcloc->nptsPad = xclocParms.nptsPad;
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
            if (xcloc->nsignals[ig] == 0)
            {
                fprintf(stderr, "%s: Signal group %d is empty\n", __func__, ig);
                ierr = 1;
                goto BCAST_ERROR;
            }
            xcloc->nTotalSignals = xcloc->nTotalSignals + xcloc->nsignals[ig]; 
        }
/*
        if (xclocParms.dt <= 0.0)
        {
            fprintf(stderr, "%s: dt=%e must be positive\n",
                    __func__, xcloc->dt) 
*/
    }
BCAST_ERROR:;
    MPI_Bcast(&ierr, 1, MPI_INT, xcloc->root, xcloc->globalComm);
    if (ierr != 0){return -1;}
    // Broadcast the parameters
    MPI_Bcast(&xcloc->nfftProcs, 1, MPI_INT, xcloc->root, xcloc->globalComm);
    // Make the global XC pairs
    ierr = xcloc_makeXCPairs(xcloc);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to make XC pair on rank %d\n",
                __func__, xcloc->globalCommRank);
    }

    // Initialize the FFTs
/*
    ierrAll = 0;
    if (xcloc->ldoFFT)
    {
        ierrAll = 0;
        for (ig=0; ig<xcloc->nSignalGroups; ig++)
        {
            ierr = xcloc_xcfftMPI_initialize(xcloc->npts[ig],
                                             xcloc->nptsPad[ig],
                                             xcloc->nsignals[ig],
                                             xcloc->fftComm,
                                             xcloc->root,
                                             NULL,
                                             &xcloc->xcfftMPI[ig]);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error initialzing xcfftMPI group %d\n",
                        __func__, ig+1);
                ierrAll = ierrAll + 1;
            } 
        }
        ierr = ierrAll;
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, xcloc->fftComm); 
    }
    MPI_Allreduce(&ierrAll, &ierr, 1, MPI_INT, MPI_SUM, xcloc->globalComm);
    if (ierr != 0){return -1;}
*/
 
    return 0;
}
//============================================================================//
/*!
 * @brief Defines the cross-correlation pairs for the different signal groups.
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
int xcloc_makeCommunicators(struct xcloc_struct *xcloc)
{
    MPI_Comm globalComm;
    int color, i, ierr;
    globalComm = xcloc->globalComm;
    // Split the global communicator into signal communicators 
    xcloc->signalGroup
        = (int *) calloc((size_t) xcloc->globalCommSize, sizeof(int));
    ierr = xcloc_splitCommunicator(xcloc->globalCommSize,
                                   xcloc->nSignalGroups,
                                   xcloc->signalGroup);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error spliting global communicator\n", __func__);
        return -1;
    }
    // For each signal communicator an FFT group is required
    color = xcloc->signalGroup[xcloc->globalCommRank];
    MPI_Comm_split(xcloc->globalComm, color, xcloc->globalCommRank,
                   &xcloc->signalComm);
    MPI_Comm_size(xcloc->signalComm, &xcloc->signalCommSize);
    MPI_Comm_size(xcloc->signalComm, &xcloc->signalCommRank);
    // The FFT groups are the masters of the local migration tables
    return 0;
}
//============================================================================//
int xcloc_splitCommunicator(const int nprocs, const int nparts, 
                            int *__restrict__ myGroup)
{
    double dPart;
    int i, ipart, low, high;
    if (nprocs < 1 || nparts < 1 || myGroup == NULL)
    {
        if (nprocs < 1)
        {
            fprintf(stderr, "%s: No processes %d\n", __func__, nprocs);
        }
        if (nparts < 1)
        {
            fprintf(stderr, "%s: No partitions %d\n", __func__, nparts);
        }
        if (myGroup == NULL)
        {
            fprintf(stderr, "%s: myGroup is NULL\n", __func__);
        }        
        return -1;
    }
    // Everyone is in the same group 
    if (nparts == 1)
    {
        for (i=0; i< nprocs; i++){myGroup[i] = 0;}
        return 0;
    }
    // Initialize
    for (i=0; i<nprocs; i++){myGroup[i] =-1;}
    // Computing spacing of partitions
    dPart = (double) nprocs/(double) nparts;
    for (ipart=0; ipart<nparts; ipart++)
    {
        low  = (int) round((double) ipart*dPart);
        high = (int) round((double) (ipart + 1)*dPart);
        if (ipart == 0){low = 0;}
        if (ipart == nparts - 1){high = nprocs;}
        for (i=low; i<high; i++)
        {
            myGroup[i] = ipart;
        }
    }
    // Verify
    for (i=0; i<nparts; i++)
    {
        if (myGroup[i] < 0)
        {
            fprintf(stderr, "%s: Failed to initialize process %d\n",
                    __func__, i);
            return -1;
        }
    }
    return 0;
}
//============================================================================//
int xcloc_finalize(struct xcloc_struct *xcloc)
{
    if (xcloc->signalGroup != NULL){free(xcloc->signalGroup);}
    if (xcloc->xcPairs != NULL){free(xcloc->xcPairs);}
    MPI_Comm_free(&xcloc->globalComm);
    memset(xcloc, 0, sizeof(struct xcloc_struct));
    return 0;
}
