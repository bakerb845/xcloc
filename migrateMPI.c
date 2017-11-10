#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "xcloc_migrate.h"

struct migrateMPI_struct
{
    struct migrate_struct migrate;
    MPI_Comm globalComm; /*!< Global communicator. */
    MPI_Comm imageInterComm;  /*!< Communication between master and slaves in
                                   table. */
    MPI_Comm imageIntraComm;  /*!< Communication between master processes in
                                   in tables. */
 
    int *g2lTable;  /*!< Maps from the global number to the groups that require
                         it. */
    int rank;
    int nprocs;
    int root;
    int ntables;   /*!< Number of tables in grid. */
};

//============================================================================//
/*!
 * @param[in] nxc       Number of cross-correlograms. 
 * @param[in] xcPairs   Maps from the ixc'th cross-correlation pair to the
 *                      global table indices.
 *
 */
int xcloc_migrateMPI_initialize(const int ntablesIn,
                                const int ngrdIn,
                                const double dtIn,
                                const int nSubGridsIn,
                                const int chunkSizeIn,
                                const MPI_Comm comm, const int master,
                                struct migrateMPI_struct *migrate)
{
    double dt;
    int ierr, chunkSize, ngrd, nSubGrids, ntables;
    memset(migrate, 0, sizeof(struct migrateMPI_struct));
    MPI_Comm_dup(comm, &migrate->globalComm);
    MPI_Comm_rank(migrate->globalComm, &migrate->rank);
    MPI_Comm_size(migrate->globalComm, &migrate->nprocs);
    migrate->root = master;
    ierr = 0;
    if (migrate->rank == master)
    {
        dt = dtIn;
        ntables = ntablesIn;
        ngrd = ngrdIn;
        nSubGrids = nSubGridsIn;
        chunkSize = chunkSizeIn;
        if (ngrd < 1 || ntables < 1 || dt <= 0.0)
        {
            if (ngrd < 1)
            {
                fprintf(stderr, "%s: number of grid pts=%d must be positive\n",
                        __func__, ngrd);
            }
            if (ntables < 2)
            {
                fprintf(stderr, "%s: ntables=%d must be at least 2\n",
                        __func__, ntables);
            }
            if (chunkSize%DALES_MEM_ALIGNMENT != 0)
            {
                fprintf(stderr, "%s: chunkSize=%d not divisible by %d\n",
                        __func__, chunkSize, DALES_MEM_ALIGNMENT);
            }

            if (dt <= 0.0)
            {
                fprintf(stderr, "%s: dt=%f must be positive\n", __func__, dt);
            }
            ierr = 1; 
            goto ERROR1;
        }
        // Partition the domain for use by vector scatter and gathers

    }
ERROR1:;
    MPI_Bcast(&ierr, 1, MPI_INT, migrate->root, migrate->globalComm);
    if (ierr != 0){return -1;}
    // Broadcast variables
    MPI_Bcast(&chunkSize, 1, MPI_INT, migrate->root, migrate->globalComm);
    MPI_Bcast(&ngrd,      1, MPI_INT, migrate->root, migrate->globalComm);
    MPI_Bcast(&ntables,   1, MPI_INT, migrate->root, migrate->globalComm);
    MPI_Bcast(&dt,      1, MPI_DOUBLE, migrate->root, migrate->globalComm);
    // Partition the grid into 
    return 0;
}
//============================================================================//
/*
int xcloc_migrateMPI_scatterTable(const int it, const int ngrdAll,
                                  const void *__restrict__ ttimesGrd,
                                  const MPI_Datatype dataType,
                                  struct migrateMPI_struct *migrate)
{
    double *ttimesGrd64;
    float *ttimesGrd32;
    // Look up the table ID

    if (dataType == MPI_DOUBLE)
    {
        ttimesGrd64 = (double *) ttimesGrd;
    }
    else if (dataType == MPI_FLOAT)
    {
        ttimesGrd32 = (float *) ttimesGrd;
    }
    else
    {
        fprintf(stderr, "%s: Invalid data type\n", __func__);
        return -1;
    }
    return 0;
}
//============================================================================//
int xcloc_migrateMPI_computeDSMImage( )
{
    return 0;
}
//============================================================================//
int xcloc_migrateMPI_gatherImage(const int ngrdAll,
                                 float *__restrict__ image,
                                 const struct migrateMPI_struct migrate)
{
    return 0;
}
*/
