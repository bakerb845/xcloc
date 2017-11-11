#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "xcloc.h"
/*
#ifndef XCLOC_USE_MPI
# warning "This only works for MPI" 
#endif
*/

struct xclocParms_struct
{
    double dt;        /*!< Sampling period of input signals. */
    int *signalGroup; /*!< It is possible to assign each signal to a group. 
                           e.g., P and S.  Consequently,  only signals in the
                           same group will be cross-correlated.  This is an
                           array of dimension [nsignals] where the is'th signal
                           returns the signal group number. */
    int nfftProcs;  /*!< Number of processes in that will be involved in
                         computing the cross-correlations.  Note, that these
                         processes will be the master ranks in the migration.
                         Therefore, this can and likely should be less than
                         the total number of processes available on the global
                         communicator. */
    int nsignals;   /*!< Number of signals. */
    int npts;       /*!< Number of points in input signals. */
    int nptsPad;    /*!< This is a tuning parameter.  The cross-correlations
                         will have dimension [2*npts+1] by default, thus, if
                         2*npts+1 is a large prime number the FFT pefrormance
                         will be greatly diminished.  In this case it may be
                         advantageous to pad the transforms. */
};

struct xcloc_struct
{
    MPI_Comm globalComm; /*!< Global communicator. */
    MPI_Comm fftComm;    /*!< FFT communicator. */
    int globalCommSize;  /*!< Size of global communicator. */
    int globalCommRank;  /*!< Rank on global communicator. */
    int nfftProcs;       /*!< Number of processes in FFT. */
    int root;            /*!< Rank of root process.  Will be 0. */
};

/*!
 * @brief Checks the xcloc parameters.
 */
int xcloc_checkParameters(const struct xclocParms_struct xclocParms)
{
    int ierr;
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
    return ierr;
}

int xcloc_initialize(const MPI_Comm comm,
                     const struct xclocParms_struct xclocParms,
                     struct xcloc_struct *xcloc)
{
    int ierr, groupMax, groupMin, mpiInit;
    ierr = 0;
    memset(xcloc, 0, sizeof(struct xcloc_struct));
    xcloc->root = 0;
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
    return 0;
}

int xcloc_finalize(struct xcloc_struct *xcloc)
{
    memset(xcloc, 0, sizeof(struct xcloc_struct));
    return 0;
}
