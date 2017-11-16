#ifndef XCLOC_MIGRATE_MPI_H__
#define XCLOC_MIGRATE_MPI_H__ 1
#include <stdbool.h>
#include "xcloc_config.h"
#include "xcloc_enum.h"
#include "xcloc_migrate.h"
#ifdef XCLOC_USE_MPI
#include <mpi.h>

struct migrateMPI_struct
{
    struct migrate_struct migrate;
    MPI_Comm globalComm; /*!< Global communicator. */
    MPI_Comm interComm;  /*!< Communications between the ip'th processes
                              in each migration group. */
    MPI_Comm intraComm;  /*!< Communications within the ig'th migration group. */
    float *imageLoc;     /*!< Local image.  Only pertinent on the root
                              migration group where it is an array of
                              dimension [ngrdLoc]. */
    int *xcPairsLoc;     /*!< Local XC pairs.   This i san array of dimension
                              [2*nxcLoc]. */
    int *xcrPairs;       /*!< Global XC pairs and rank that will migrate the XC.
                              This is an array of dimension [3*nxc]. */
    int *xcPtr;          /*!< Maps from the ig'th cross-correlation group
                              to thstart of xcrPairs. This is an array of
                              dimension [nMigrateGroups]. */
    int *proc2GridPtr;   /*!< Maps from the process to the start index in the
                              global migration grid.  This is an array of 
                              dimension [nprocsPerGroup+1]. */
    bool *lsendTT2Grp;   /*!< If true then send the travel time table 
                              for the is'th signal to the ig'th group.
                              This is an array of dimension
                              [nMigrateGroups x nsignals] with leading 
                              dimension nMigrateGroups. */
    double dt;           /*!< Sampling period (seconds). */
    int globalCommRank;  /*!< Rank on global communicator. */
    int globalCommSize;  /*!< Size of global communicator. */
    int nxc;             /*!< Number of cross-correlations. */
    int lxc;             /*!< Length of cross-correlations. */
    int nsignals;        /*!< Number of signals. */
    int nsignalsLoc;     /*!< Number of local signals. */
    int nxcLoc;          /*!< Number of local cross-correlations.  This is
                              the number of cross-correlations that this
                              process will migrate. */
    int nProcsPerGroup; /*!< Number of processes per migration group.  The
                             globalCommSize = nProcsPerGroup x nMigrateGroups.*/
    int nMigrateGroups; /*!< Number of migraiton groups.  The 
                             globalCommSize = nProcsPerGroup x nMigrateGroups.*/
    int myMigrateGroup; /*!< My migration group. */
    int root;           /*!< Root process.  This will be 0. */
    int ngrd;           /*!< Number of grid points in total grid. */
    int ngrdLoc;        /*!< Number of grid poitns in local grid. */
    int chunkSize;      /*!< Chunk size in migration.  This is a tuning
                             parameter. */
    bool linit;         /*!< Flag indicating the MPI XC correlation structure
                             has been initialized. */
};

#ifdef __cplusplus
extern "C"
{
#endif

/*! Initializes the migrateMPI structure. */
int xcloc_migrateMPI_initialize(const MPI_Comm globalComm,
                                const int nMigrateGroups,
                                const int nsignals, const int nxc,
                                const int chunkSize, const int ngrd,
                                const int lxc, const double dt,
                                const int *__restrict__ xcrPairs,
                                struct migrateMPI_struct *migrateMPI);
/*! Frees memory and deletes comm on the migrateMPI structure. */
int xcloc_migrateMPI_finalize(struct migrateMPI_struct *migrateMPI);
/*! Sets the travel time tables. */
int xcloc_migrateMPI_setTableFromRoot(const int itIn, const int ngrdIn,
                                      const enum xclocPrecision_enum precision,
                                      const void *__restrict__ ttimes,
                                      struct migrateMPI_struct *migrateMPI);

/*! Private function */
int xcloc_migrateMPI_createLocalXCPairs(struct migrateMPI_struct *migrateMPI);

#ifdef __cplusplus
}
#endif
#endif /* XCLOC_USE_MPI */
#endif /* XCLOC_MIGRATE_MPI_H__ */
