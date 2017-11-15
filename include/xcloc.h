#ifndef XCLOC_H__
#define XCLOC_H__ 1
#include "xcloc_config.h"
#include "xcloc_enum.h"
#include "xcloc_xcfft.h"
#include "xcloc_migrate.h"
#include "xcloc_hdf5.h"
#include "xcloc_xdmf.h"
#include "xcloc_envelope.h"
#include "xcloc_rmsFilter.h"
#ifdef XCLOC_USE_MPI
#include "xcloc_xcfftMPI.h"
//#include "xcloc_migrateMPI.h"
#endif

struct xclocParms_struct
{
    double dt;        /*!< Sampling period of input signals. */
    int *signalGroup; /*!< It is possible to assign each signal to a group. 
                           e.g., P and S.  Consequently,  only signals in the
                           same group will be cross-correlated.  This is an
                           array of dimension [nsignals] where the is'th signal
                           returns the signal group number. */
    int ngridProcs; /*!< Number of processes in each migration grid. */
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
    struct xcfftMPI_struct xcfftMPI; /*!< Array of FFT structures.  This has
                                           dimension [nSignalGroups]. */
    int *nsignals;       /*!< Number of signals in each group.  This is an array
                              of dimension [nSignalGroups]. */
    int *signalGroup;
    int *xcPairs;
    MPI_Comm globalComm; /*!< Global communicator. */
    MPI_Comm signalComm; /*!< Signal communicator. */
    MPI_Comm signalCommSize;
    MPI_Comm signalCommRank;
    MPI_Comm fftComm;    /*!< FFT communicator. */
    MPI_Comm migrateComm;/*!< Migration communicator. */
    int globalCommSize;  /*!< Size of global communicator. */
    int globalCommRank;  /*!< Rank on global communicator. */
    //int nmigrateProcs;   /*!< Number of processes in a migration. */
    int ngridProcs;      /*!< Number of processes in the migration
                              group. */
    int nfftProcs;       /*!< Number of processes in FFT. */
    int nTotalSignals;   /*!< Cumulative number of signals. */
    int npts;           /*!< Number of points in signals. */
    int nptsPad;        /*!< Number of points to pad signals prior to
                             cross-correlating.  Note, that the correlation
                             length will be 2*nptsPad-1. */
    int nTotalXCs;       /*!< Cumulative number of cross-correlations. */
    int nSignalGroups;   /*!< Number of signal groups.  For example if processing
                              P and S waves then this would be 2. */
    int root;            /*!< Rank of root process.  Will be 0. */
    bool ldoFFT;         /*!< If true then I will compute the FFTs. */
    bool linit;          /*!< Flag indicating whether or not the structure is
                              initialized. */
};

#ifdef __cplusplus
extern "C"
{
#endif

int xcloc_initialize(const MPI_Comm comm,
                     const struct xclocParms_struct xclocParms,
                     struct xcloc_struct *xcloc);
int xcloc_makeXCPairs(struct xcloc_struct *xcloc);
int xcloc_finalize(struct xcloc_struct *xcloc);

#ifdef __cplusplus
}
#endif
#endif
