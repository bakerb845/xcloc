#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include "xcloc_migrateMPI.h"
#include <ipps.h>

static int cmp_int(const void *x, const void *y)
{
    int xx = *(int *) x;
    int yy = *(int *) y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}
//============================================================================//
/*!
 * @brief Initialize the migrateMPI structure.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrateMPI_initialize(const MPI_Comm globalComm,
                                const int nMigrateGroups,
                                const int nsignals, const int nxc,
                                const int chunkSize, const int ngrd, 
                                const int lxc, const double dt,
                                const int *__restrict__ xcrPairs,
                                struct migrateMPI_struct *migrate)
{
    double dPart;
    int color, ierr, ierrAll, ig, ip, ir, is1, is2, ixc, group, group0, key,
        ngrdAll, ngrdLoc;
    const int root = 0;
    ierr = 0;
    memset(migrate, 0, sizeof(struct migrateMPI_struct));
    MPI_Comm_dup(globalComm, &migrate->globalComm);
    MPI_Comm_rank(migrate->globalComm, &migrate->globalCommRank);
    MPI_Comm_size(migrate->globalComm, &migrate->globalCommSize);
    migrate->root = root;
    if (migrate->globalCommRank == migrate->root)
    {
        if (nMigrateGroups < 1 || nxc < 1 || nsignals < 1 ||
            ngrd < 1 || chunkSize%XCLOC_MEM_ALIGNMENT != 0 ||
            lxc < 1 || dt <= 0.0 ||
            nMigrateGroups > ngrd || xcrPairs == NULL)
        {
printf("input error\n");
            if (nxc < 1)
            {
                fprintf(stderr, "%s: No correlations\n", __func__);
            }
            ierr = 1;
            goto BCAST_ERROR; 
        }
        if (migrate->globalCommSize%nMigrateGroups != 0)
        {
            ierr = 1;
            goto BCAST_ERROR;
        }
        migrate->dt = dt; 
        migrate->lxc = lxc;
        migrate->nProcsPerGroup = migrate->globalCommSize/nMigrateGroups;
        migrate->nMigrateGroups = nMigrateGroups;
        migrate->ngrd = ngrd;
        migrate->chunkSize = chunkSize;
        migrate->nsignals = nsignals;
        migrate->nxc = nxc;
        // Figure out how to scatter/gather the tables

        // Ensure the cross-correlation pairs make sense
        ig = 1;
        group0 = xcrPairs[2];
        migrate->xcPtr
            = (int *) calloc((size_t) (migrate->nMigrateGroups+1), sizeof(int));
        migrate->xcrPairs 
            = (int *) calloc((size_t) (3*migrate->nxc), sizeof(int));
        ippsCopy_32s(xcrPairs, migrate->xcrPairs, 3*migrate->nxc);
        for (ixc=0; ixc<nxc; ixc++)
        {
            is1   = xcrPairs[3*ixc];
            is2   = xcrPairs[3*ixc+1];
            group = xcrPairs[3*ixc+2]; 
            if (group0 != group)
            {
                migrate->xcPtr[ig] = ixc;
                ig = ig + 1;
            }
            if (is1 < 0 || is1 >= nsignals)
            {
                fprintf(stderr, "%s: signal_1=%d must be in range [0,%d]\n",
                        __func__, is1, nsignals-1);
                ierr = 1;
                goto BCAST_ERROR;
            }
            if (is2 < 0 || is2 >= nsignals)
            {
                fprintf(stderr, "%s: signal_2=%d must be in range [0,%d]\n",
                        __func__, is2, nsignals-1);
                ierr = 1;
                goto BCAST_ERROR;
            }
            if (group < 0 || group >= nMigrateGroups)
            {
                fprintf(stderr, "%s: group=%d must be in range [0,%d]\n",
                        __func__, group, nMigrateGroups-1);
                ierr = 1;
                goto BCAST_ERROR;
            }
            if (group < group0)
            {
                fprintf(stderr, "%s: xcrPairs must be sorted on groups\n",
                        __func__);
                ierr = 1;
                goto BCAST_ERROR;
            }
            group0 = group;
        }
        if (ig != nMigrateGroups)
        {
            fprintf(stderr, "%s: xcrPairs has %d groups; expecting %d groups\n",
                    __func__, ig, nMigrateGroups); 
            ierr = 1;
            goto BCAST_ERROR;
        }
        migrate->xcPtr[nMigrateGroups] = migrate->nxc;
        // Split the grid
        migrate->proc2GridPtr
            = (int *) calloc((size_t) (migrate->nProcsPerGroup+1), sizeof(int));
        dPart = (double) migrate->ngrd/(double) migrate->nProcsPerGroup;
        for (ip=0; ip<migrate->nProcsPerGroup; ip++)
        {
            migrate->proc2GridPtr[ip] = (int) round((double) ip*dPart);
            if (ip == 0){migrate->proc2GridPtr[ip] = 0;}
            migrate->proc2GridPtr[ip] 
                = MIN(migrate->proc2GridPtr[ip], migrate->ngrd);
        }
        migrate->proc2GridPtr[migrate->nProcsPerGroup] = migrate->ngrd;
        ngrdAll = 0;
        for (ip=0; ip<migrate->nProcsPerGroup; ip++)
        {
            ngrdLoc = migrate->proc2GridPtr[ip+1] - migrate->proc2GridPtr[ip]; 
            ngrdAll = ngrdAll + ngrdLoc;
            fprintf(stdout, "%s: Process %d on grid will have %d grid points\n",
                    __func__, ip, ngrdLoc);
        }
        if (ngrdAll != ngrd)
        {
            fprintf(stderr, "%s: ngrdAll=%d != ngrd=%d\n",
                    __func__, ngrdAll, ngrd);
            ierr = 1;
            goto BCAST_ERROR;
        }
    }
BCAST_ERROR:;
    MPI_Bcast(&ierr, 1, MPI_INT, root, migrate->globalComm);
    if (ierr != 0){return -1;}
    MPI_Bcast(&migrate->dt, 1, MPI_DOUBLE, root, migrate->globalComm);
    MPI_Bcast(&migrate->lxc,            1, MPI_INT, root, migrate->globalComm);
    MPI_Bcast(&migrate->nxc,            1, MPI_INT, root, migrate->globalComm);
    MPI_Bcast(&migrate->nsignals,       1, MPI_INT, root, migrate->globalComm);
    MPI_Bcast(&migrate->nMigrateGroups, 1, MPI_INT, root, migrate->globalComm);
    MPI_Bcast(&migrate->nProcsPerGroup, 1, MPI_INT, root, migrate->globalComm);
    MPI_Bcast(&migrate->ngrd,           1, MPI_INT, root, migrate->globalComm);
    MPI_Bcast(&migrate->chunkSize,      1, MPI_INT, root, migrate->globalComm);
    if (migrate->globalCommRank != migrate->root)
    {
        migrate->xcPtr
            = (int *) calloc((size_t) (migrate->nMigrateGroups+1), sizeof(int));
        migrate->proc2GridPtr
            = (int *) calloc((size_t) (migrate->nProcsPerGroup+1), sizeof(int));
        migrate->xcrPairs
            = (int *) calloc((size_t) (3*migrate->nxc), sizeof(int));
    }
    MPI_Bcast(migrate->xcPtr,        migrate->nMigrateGroups+1, MPI_INT,
              root, migrate->globalComm);
    MPI_Bcast(migrate->proc2GridPtr, migrate->nProcsPerGroup+1, MPI_INT,
              root, migrate->globalComm);
    MPI_Bcast(migrate->xcrPairs, 3*migrate->nxc,                MPI_INT,
              root, migrate->globalComm);
    // Make it possible for all processes in a migration group to communicate
    color = MPI_UNDEFINED;
    key = 0;
    ir = 0;
    for (ig=0; ig<migrate->nMigrateGroups; ig++)
    {
        for (ip=0; ip<migrate->nProcsPerGroup; ip++)
        {
            if (migrate->globalCommRank == ir)
            {
                color = ig; // My migration group
                key   = ip; // My rank (position) in the group
            }
            ir = ir + 1;
        }
    }
    ierr = MPI_Comm_split(migrate->globalComm, color, key, &migrate->intraComm);
    if (ierr != MPI_SUCCESS)
    {
        fprintf(stdout, "%s: Failed making intraComm\n", __func__);
        return -1; 
    }
    migrate->myMigrateGroup = color;
    migrate->ngrdLoc = migrate->proc2GridPtr[key+1]
                     - migrate->proc2GridPtr[key];
    // Make it possible for processes of the same rank each migration group
    // to communicate.
    color = MPI_UNDEFINED;
    key = 0;
    ir = 0;
    for (ig=0; ig<migrate->nMigrateGroups; ig++)
    {
        for (ip=0; ip<migrate->nProcsPerGroup; ip++)
        {
            if (migrate->globalCommRank == ir) 
            {
                color = ip; // My migration group
                key   = ig; // My rank (position) in the interComm
            }
            ir = ir + 1;
        }
    }
    ierr = MPI_Comm_split(migrate->globalComm, color, key, &migrate->interComm);
    if (ierr != MPI_SUCCESS)
    {
        fprintf(stdout, "%s: Failed making interComm\n", __func__);
        return -1;
    }
    // Figure out my local table
    ierr = xcloc_migrateMPI_createLocalXCPairs(migrate);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to make xc pairs\n", __func__);
        return -1;
    }
    // Initialize the MPI grids
    if (migrate->globalCommRank == migrate->root)
    {
        fprintf(stdout, "%s: Initializing local migrate structs\n", __func__);
    }
    ierr = xcloc_migrate_initialize(migrate->nsignalsLoc, // ntables
                                    migrate->ngrdLoc,     // ngrd
                                    migrate->nxcLoc,
                                    migrate->lxc,
                                    migrate->chunkSize,
                                    migrate->dt,
                                    migrate->xcPairsLoc,
                                    &migrate->migrate);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed initializing migrate struct on rank %d\n",
                __func__, migrate->globalCommRank);
        ierr = 1;
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, migrate->globalComm);
    if (ierr != 0){return -1;} 
    // Set space for image on the 0'th group
    if (migrate->myMigrateGroup == 0)
    {
        migrate->imageLoc
           = (float *) calloc((size_t) migrate->ngrdLoc, sizeof(float));
    }
    migrate->linit = true;
    return 0;
}
//============================================================================//
/*!
 * @brief Releases memory on the migrateMPI and destroyed communicators.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrateMPI_finalize(struct migrateMPI_struct *migrateMPI)
{
    if (!migrateMPI->linit){return 0;}
    if (migrateMPI->imageLoc     != NULL){free(migrateMPI->imageLoc);}
    if (migrateMPI->xcPairsLoc   != NULL){free(migrateMPI->xcPairsLoc);}
    if (migrateMPI->xcrPairs     != NULL){free(migrateMPI->xcrPairs);}
    if (migrateMPI->xcPtr        != NULL){free(migrateMPI->xcPtr);}
    if (migrateMPI->proc2GridPtr != NULL){free(migrateMPI->proc2GridPtr);}
    if (migrateMPI->lsendTT2Grp  != NULL){free(migrateMPI->lsendTT2Grp);}
    xcloc_migrate_finalize(&migrateMPI->migrate);
    MPI_Comm_free(&migrateMPI->interComm);
    MPI_Comm_free(&migrateMPI->intraComm);
    MPI_Comm_free(&migrateMPI->globalComm);
    memset(migrateMPI, 0, sizeof(struct migrateMPI_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Creates the local XC pairs on the migrateMPI structure and
 *        determines which travel time tables will be sent to which
 *        groups.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrateMPI_createLocalXCPairs(struct migrateMPI_struct *migrateMPI)
{
    int *item, *work, i, i0, i1, ierr, ierrAll, ig, indx, ixc, ixc1, ixc2,
        jxc, nsend, nwork, nxcLoc, nSignalsLoc;
    ierr = 0;
    if (migrateMPI->nxc < 1 || migrateMPI->xcrPairs == NULL)
    {
        if (migrateMPI->nxc < 1)
        {
            fprintf(stderr, "%s: Error no xcs\n", __func__);
        }
        if (migrateMPI->xcrPairs == NULL)
        {
            fprintf(stderr, "%s: xcrPairs is NULL\n", __func__);
        }
        return -1;
    }
    work = (int *) calloc((size_t) (2*migrateMPI->nxc + 1), sizeof(int));  
    memset(work, -1, (size_t) (2*migrateMPI->nxc + 1)*sizeof(int));
    // Create a unique list of origins
    nSignalsLoc = 0; 
    nxcLoc = 0;
    for (ixc=0; ixc<migrateMPI->nxc; ixc++)
    {
        if (migrateMPI->xcrPairs[3*ixc+2] == migrateMPI->myMigrateGroup)
        {
            i0 = migrateMPI->xcrPairs[3*ixc];
            i1 = migrateMPI->xcrPairs[3*ixc+1];
            // Put this in the list
            for (i=0; i<nSignalsLoc+1; i++) 
            {
                if (work[i] == i0){break;}
                if (work[i] ==-1)
                {
                    work[nSignalsLoc] = i0;
                    nSignalsLoc = nSignalsLoc + 1;
                    break;
                }
            }
            for (i=0; i<nSignalsLoc+1; i++) 
            {
                if (work[i] == i1){break;}
                if (work[i] ==-1)
                {
                    work[nSignalsLoc] = i1;
                    nSignalsLoc = nSignalsLoc + 1; 
                    break;
                }
            }
            nxcLoc = nxcLoc + 1;
        }
    }
    if (nxcLoc < 1)
    {
        fprintf(stderr, "%s: No cross-correlations\n", __func__);
        ierr = 1;
    }
    if (nSignalsLoc < 1)
    {
        fprintf(stderr, "%s: No local signals\n", __func__);
        ierr = 1;
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, migrateMPI->globalComm);
    if (ierrAll != 0){return -1;}
    // Sort the list
    migrateMPI->nxcLoc = nxcLoc;
    migrateMPI->nsignalsLoc = nSignalsLoc;
    qsort((void *) work, (size_t) nSignalsLoc, sizeof(int), cmp_int);
    // Map the global signals to the local signals
    migrateMPI->xcPairsLoc = (int *) calloc((size_t) (2*nxcLoc), sizeof(int));
    jxc = 0;
    for (ixc=0; ixc<migrateMPI->nxc; ixc++)
    {
        if (migrateMPI->xcrPairs[3*ixc+2] == migrateMPI->myMigrateGroup)
        {
            i0 = migrateMPI->xcrPairs[3*ixc];
            i1 = migrateMPI->xcrPairs[3*ixc+1];
            item = (int *) bsearch((const void *) &i0, (const void *) work,
                                   nSignalsLoc, sizeof(int),
                                   cmp_int);
            if (item == NULL)
            {
                fprintf(stderr, "%s: Couldn't find %d in work\n", __func__, i0);
                return -1;
            }
            migrateMPI->xcPairsLoc[2*jxc+0] = item - work;
            item = (int *) bsearch((const void *) &i1, (const void *) work,
                                   nSignalsLoc, sizeof(int),
                                   cmp_int);
            if (item == NULL)
            {
                fprintf(stderr, "%s: Couldn't find %d in work\n", __func__, i1);
                return -1;
            }
            migrateMPI->xcPairsLoc[2*jxc+1] = item - work;
            jxc = jxc + 1;
        }
    }
    // Have the master figure out if the it'th table will be sent to the ig'th
    // group
    nwork = migrateMPI->nMigrateGroups*migrateMPI->nsignals;
    migrateMPI->lsendTT2Grp = (bool *) calloc((size_t) nwork, sizeof(bool));
    if (migrateMPI->globalCommRank == migrateMPI->root)
    {
        // Look through the migration groups
        for (ig=0; ig<migrateMPI->nMigrateGroups; ig++)
        {
            nsend = 0;
            ixc1 = migrateMPI->xcPtr[ig]; 
            ixc2 = migrateMPI->xcPtr[ig+1];
            // And tag all the matching signals 
            for (ixc=ixc1; ixc<ixc2; ixc++)
            {
                i0 = migrateMPI->xcrPairs[3*ixc+0];
                i1 = migrateMPI->xcrPairs[3*ixc+1];
                indx = i0*migrateMPI->nMigrateGroups + ig;
                if (!migrateMPI->lsendTT2Grp[indx])
                {
                    migrateMPI->lsendTT2Grp[indx] = true; 
                    nsend = nsend + 1;
                }
                indx = i1*migrateMPI->nMigrateGroups + ig;
                if (!migrateMPI->lsendTT2Grp[indx])
                {
                    migrateMPI->lsendTT2Grp[indx] = true;
                    nsend = nsend + 1;
                }
            }
            fprintf(stdout, "%s: Will send %d tables to group %d\n",
                    __func__, nsend, ig); 
        }
    }
    MPI_Bcast(migrateMPI->lsendTT2Grp, nwork, MPI_C_BOOL,
              migrateMPI->root, migrateMPI->globalComm);
    return 0;
}

//============================================================================//
/*!
 * @brief Sets the travel-time table on the appropriate migration groups.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrateMPI_setTableFromRoot(const int itIn, const int ngrdIn,
                                      const enum xclocPrecision_enum precisionIn,
                                      const void *__restrict__ ttimes,
                                      struct migrateMPI_struct *migrateMPI)
{
    void *work = NULL;
    void *workLoc = NULL;
    MPI_Datatype dataType;
    MPI_Status status;
    enum xclocPrecision_enum precision;
    int *displs, *sendCounts, ierr, ierrAll, ig, indx, ip, it, itLoc,
        ixc, ixc1, ixc2, gridID, groupID, ngrd, prec, recvCount;
    ierr = 0;
    if (!migrateMPI->linit)
    {
        fprintf(stderr, "%s: migrate not initialized\n", __func__);
        return -1;
    }
    MPI_Comm_rank(migrateMPI->intraComm, &gridID);
    MPI_Comm_rank(migrateMPI->interComm, &groupID);
    if (migrateMPI->globalCommRank == migrateMPI->root)
    {
        ngrd = ngrdIn;
        it   = itIn;
        if (ngrd != migrateMPI->ngrd)
        {
            fprintf(stderr, "%s: ngrd=%d != migrateMPI->ngrd=%d\n", __func__,
                    ngrd, migrateMPI->ngrd);
            ierr = 1;
        }
        if (it < 0 || it >= migrateMPI->nsignals)
        {
            fprintf(stderr, "%s: it=%d must be in range [0,%d]\n", __func__,
                    it, migrateMPI->nsignals);
            ierr = 1;
        }
        prec = (int) precisionIn;
    }
    MPI_Bcast(&ierr, 1, MPI_INT, migrateMPI->root, migrateMPI->globalComm);
    if (ierr != 0){return -1;}
    MPI_Bcast(&it,   1, MPI_INT, migrateMPI->root, migrateMPI->globalComm);
    MPI_Bcast(&ngrd, 1, MPI_INT, migrateMPI->root, migrateMPI->globalComm);
    MPI_Bcast(&prec, 1, MPI_INT, migrateMPI->root, migrateMPI->globalComm);
    precision = (enum xclocPrecision_enum) prec;
    if (precision == XCLOC_SINGLE_PRECISION)
    {
        dataType = MPI_FLOAT;
        if (gridID == migrateMPI->root)
        {
             work = calloc((size_t) migrateMPI->ngrd, sizeof(float)); 
        }
        workLoc = calloc((size_t) migrateMPI->ngrdLoc, sizeof(float));
    }
    else if (dataType == XCLOC_DOUBLE_PRECISION)
    {
        dataType = MPI_DOUBLE; 
        if (gridID == migrateMPI->root)
        {
            work = calloc((size_t) migrateMPI->ngrd, sizeof(double)); 
        }
        workLoc = calloc((size_t) migrateMPI->ngrdLoc, sizeof(double));
    }
    else
    {
        fprintf(stderr, "%s: Invalid precision %d\n", __func__, (int)precision);
    }
    // Loop on the migration groups
    ierr = 0;
    for (ig=0; ig<migrateMPI->nMigrateGroups; ig++)
    {
        indx = it*migrateMPI->nMigrateGroups + ig;
        if (!migrateMPI->lsendTT2Grp[indx]){continue;}
        // Have the root send the table to the appropriate migration table root
        if (migrateMPI->globalCommRank == migrateMPI->root)
        {
            // Send it to the ig'th group
            if (ig > 0)
            {
                MPI_Send(ttimes, migrateMPI->ngrd, dataType, ig,
                         migrateMPI->root, migrateMPI->interComm);
            }
            // Keep this guy for myself
            else
            {
                if (precision == XCLOC_SINGLE_PRECISION)
                {
                    memcpy(work, ttimes, migrateMPI->ngrd*sizeof(float));
                }
                else
                {
                    memcpy(work, ttimes, migrateMPI->ngrd*sizeof(double));
                }
            }
        }
        else
        {
            if (ig == migrateMPI->myMigrateGroup)
            {
                if (gridID == migrateMPI->root)
                {
                    MPI_Recv(work, migrateMPI->ngrd, dataType, migrateMPI->root,
                             MPI_ANY_TAG, migrateMPI->interComm, &status);
                } 
            }
        }
        // Have the root scatter it then put it into the table 
        if (ig == migrateMPI->myMigrateGroup)
        {
            ixc1 = migrateMPI->xcPtr[ig];
            ixc2 = migrateMPI->xcPtr[ig+1]; 
            itLoc =-1;
            for (ixc=ixc1; ixc<ixc2; ixc++)
            {
                if (it == migrateMPI->xcrPairs[3*ixc+0])
                {
                    itLoc = migrateMPI->xcPairsLoc[2*(ixc-ixc1)];
                    break;
                }
                if (it == migrateMPI->xcrPairs[3*ixc+1])
                {
                    itLoc = migrateMPI->xcPairsLoc[2*(ixc-ixc1)+1];
                    break;
                }
            }
            // Have the master figure out the scatter offsets
            displs = NULL;
            sendCounts = NULL;
            if (gridID == migrateMPI->root)
            {
                displs = (int *) calloc((size_t) migrateMPI->nProcsPerGroup,
                                        sizeof(int));
                sendCounts = (int *) calloc((size_t) migrateMPI->nProcsPerGroup,
                                            sizeof(int));
                for (ip=0; ip<migrateMPI->nProcsPerGroup; ip++)
                {
                    displs[ip]     = migrateMPI->proc2GridPtr[ip];
                    sendCounts[ip] = migrateMPI->proc2GridPtr[ip+1]
                                   - migrateMPI->proc2GridPtr[ip]; 
                }
            }
            recvCount = migrateMPI->ngrdLoc;
            MPI_Scatterv(work, sendCounts, displs, dataType,
                         workLoc, recvCount, dataType, 
                         migrateMPI->root, migrateMPI->intraComm);
            if (displs != NULL){free(displs);}
            if (sendCounts != NULL){free(sendCounts);}
            // Set the table
            ierr = 0;
            if (itLoc >-1)
            {
                ierr = xcloc_migrate_setTable(itLoc, migrateMPI->ngrdLoc,
                                              precision, workLoc,
                                              &migrateMPI->migrate);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Error setting traveltime table\n",
                            __func__); 
                    ierr = 1;
                }
            }
            else
            {
                fprintf(stderr, "%s: Failed to find table %d\n", __func__, it);
                ierr = 1;
            }
        }
        // Block
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT,
                       MPI_SUM, migrateMPI->globalComm);
        if (ierr != 0){break;}
    }
    if (work != NULL){free(work);}
    if (workLoc != NULL){free(workLoc);}
    return ierr;
}
//============================================================================//
/*!
 * @brief Computes the migration image of the cross-correlograms.  Note, that 
 *        this will leave the image distributed on the 0'th migration group. 
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrateMPI_computeMigrationImage(struct migrateMPI_struct *migrateMPI)
{
    int ierr, ierrAll;
    if (!migrateMPI->linit)
    {
        fprintf(stderr, "%s: migrateMPI not initialized\n", __func__);
        return -1;
    }
    ierr = xcloc_migrate_computeMigrationImage(&migrateMPI->migrate);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing migration on globalRank %d\n",
                __func__, migrateMPI->globalCommRank);
        ierr = 1;
    }
    // Reduce the results into the root group
    MPI_Reduce(migrateMPI->migrate.migrate, migrateMPI->imageLoc,
               migrateMPI->migrate.ngrd, MPI_FLOAT,
               MPI_SUM, migrateMPI->root, migrateMPI->interComm);
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, migrateMPI->globalComm);
    return ierrAll;
}
