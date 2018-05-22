#ifndef XCLOC_H__
#define XCLOC_H__ 1
#include "xcloc_config.h"
#include "xcloc_finter.h"
#include "xcloc_enum.h"
#include "xcloc_xcfft.h"
#include "xcloc_migrate.h"
#include "xcloc_hdf5.h"
#include "xcloc_xdmf.h"
#include "xcloc_envelope.h"
#include "xcloc_rmsFilter.h"
#ifdef XCLOC_USE_MPI
#include "xcloc_migrateMPI.h"
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
    int envFIRLen;  /*!< Length of the FIR envelope function */
    int nsignals;   /*!< Number of signals. */
    int npts;       /*!< Number of points in input signals. */
    int nptsPad;    /*!< This is a tuning parameter.  The cross-correlations
                         will have dimension [2*npts+1] by default, thus, if
                         2*npts+1 is a large prime number the FFT pefrormance
                         will be greatly diminished.  In this case it may be
                         advantageous to pad the transforms. */
    int chunkSize;  /*!< Chunksize in migration grid. */
    int ngrd;       /*!< Number of grid points in model. */
    bool lphaseXCs; /*!< If true then the compute phase correlations. \n
                          Otherwise, compute the cross-correlations. */
};

struct xcloc_struct
{
    struct xcfftMPI_struct xcfftMPI; /*!< Array of FFT structures.  This has
                                           dimension [nSignalGroups]. */
    struct migrateMPI_struct migrateMPI;
    struct migrate_struct migrate;   /*!< Migration structure. */
    struct xclocEnvelope_struct envelope; /*!< Envelope structure. */
    int *nsignals;       /*!< Number of signals in each group.  This is an array
                              of dimension [nSignalGroups]. */
    int *signalGroup;
    int *xcPairs;
    int *xcPairsLoc;
    void *y1;
    void *y2;
    MPI_Comm globalComm; /*!< Global communicator. */
    //MPI_Comm signalComm; /*!< Signal communicator. */
    MPI_Comm fftComm;    /*!< FFT communicator. */
    MPI_Comm migrateComm;/*!< Migration communicator. */
    double dt;           /*!< Sampling period. */
    int chunkSize;       /*!< Chunksize in migration grid. */ 
    int globalCommSize;  /*!< Size of global communicator. */
    int globalCommRank;  /*!< Rank on global communicator. */
    int fftCommSize;
    int fftCommRank;
    int migrateCommSize;
    int migrateCommRank;
    //int nmigrateProcs;   /*!< Number of processes in a migration. */
    int nmigrateGroups;  /*!< Number of migration groups. */
    int nmigrateProcs;   /*!< Number of processes in the migration
                              group. */
    int nfftProcs;       /*!< Number of processes in FFT. */
    int ntablesLoc;      /*!< Number of tables process will end up owning. */
    int ntfSignalsLoc;   /*!< Number of transforms process will end up owning. */
    int nTotalSignals;   /*!< Cumulative number of signals. */
    int npts;            /*!< Number of points in signals. */
    int nptsPad;         /*!< Number of points to pad signals prior to
                              cross-correlating.  Note, that the correlation
                              length will be 2*nptsPad-1. */
    int ldxc;            /*!< Leading dimension of cross-correlations. */
    int lxc;             /*!< Length of cross-correlations = 2*nptsPad-1. */
    int nTotalXCs;       /*!< Cumulative number of cross-correlations. */
    int nSignalGroups;   /*!< Number of signal groups.  For example if processing
                              P and S waves then this would be 2. */
    int root;            /*!< Rank of root process.  Will be 0. */
    int ngrd;            /*!< Total number of grid points in model. */
    int leny;
    int envFIRLen;
    enum xclocPrecision_enum precision;
    enum xclocAccuracy_enum accuracy;
    bool ldoFFT;         /*!< If true then I will compute the FFTs. */
    bool lphaseXCs;      /*!< If true then compute the phase correlations. \n
                              Otherwise, compute the cross-correlations. */
    bool lenvelope;      /*!< If true compute the envelope of the cross
                              correlations. */
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
int xcloc_apply(struct xcloc_struct *xcloc);
int xcloc_gatherMigrationImage(
    const int ngrd,
    const struct xcloc_struct xcloc,
    float *image);
int xcloc_setTableFromRoot(const int itIn, const int ngrdIn,
                           const enum xclocPrecision_enum precisionIn,
                           const void *__restrict__ ttimes,
                           struct xcloc_struct *xcloc);
int xcloc_scatterDataFromRoot(const int nsignals,
                              const int lds, const int npts,
                              const MPI_Datatype sendType, 
                              const void *__restrict__ x,
                              struct xcloc_struct *xcloc);

#ifdef __cplusplus
}
#endif
#endif
