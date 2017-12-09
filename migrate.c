#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "xcloc_migrate.h"
#include <ipps.h>

int xcloc_migrate_checkParameters(const int ntables, const int ngrd,
                                  const int nxc, const int lxc,
                                  const int chunkSize, const double dt,
                                  const int xcPairs[])
{
    int ierr, maxIDx, minIDx;
    ierr = 0;
    if (ngrd < 1)
    {
        fprintf(stderr, "%s: ngrd=%d must be positive\n", __func__, ngrd);
        ierr = 1;
    }
    if (ntables < 2)
    {
        fprintf(stderr, "%s: ntables=%d must at least 2\n", __func__, ntables);
        ierr = 1;
    }
    if (dt <= 0.0)
    {
        fprintf(stderr, "%s: dt=%f must be positive\n", __func__, dt);
        ierr = 1;
    }
    if (nxc < 1)
    {
        fprintf(stderr, "%s: No cross-correlations\n", __func__);
        ierr = 1;
    }
    if (lxc < 1)
    {
        fprintf(stderr, "%s: xc length=%d must be positive\n",
                __func__, lxc);
        ierr = 1;
    }
    if (lxc%2 != 1)
    {
        fprintf(stdout, "%s: xc length=%d is even - watch for weird behavior\n",
                __func__, lxc);
    }
    if (xcPairs == NULL)
    {
        fprintf(stderr, "%s: xcPairs is NULL\n", __func__);
        ierr = 1;
    }
    if (chunkSize%XCLOC_MEM_ALIGNMENT != 0)
    {
        fprintf(stderr, "%s: chunkSize=%d not divisible by alignment=%d\n",
                __func__, chunkSize, XCLOC_MEM_ALIGNMENT);
        ierr = 1;
    }
    // Get the min and max of xcPairs and make sure they are in range
    ippsMax_32s(xcPairs, 2*nxc, &maxIDx);
    ippsMin_32s(xcPairs, 2*nxc, &minIDx);
    if (minIDx < 0)
    {
        fprintf(stderr, "%s: min index in xcPairs=%d must be positive\n",
                __func__, minIDx);
        ierr = 1;
    }
    if (maxIDx + 1 > ntables)
    {
        fprintf(stderr, "%s: xcPairs requires at least %d tables\n",
                __func__, maxIDx + 1);
        ierr = 1;
    }
    return ierr;
}
/*!
 * @brief Initializes the gridded migration structure.
 *
 * @param[in] ntables    Number of travel time tables.
 * @param[in] ngrd       Number of grid points.
 * @param[in] nxc        Number of cross-correlations.
 * @param[in] lxc        Length of cross-correlations.
 * @param[in] dt         Sampling period (s) of the signals to migrate.
 * @param[in] chunkSize  This is a tuning parameter and must be evenly
 *                       divisible by the memory alignment.  It is
 *                       recommended to use powers of 2 that are
 *                       of size x*pageWidth/sizeof(float) where x is a 
 *                       small positive integer and page width is the size
 *                       of a page on the target architecture.
 *
 * @param[out] migrate   Migration structure with space allocated for the 
 *                       migration image and the travel-time tables.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under MIT.
 *
 */
int xcloc_migrate_initialize(const int ntables, const int ngrd,
                             const int nxc, const int lxc,
                             const int chunkSize,
                             const double dt, 
                             const int *__restrict__ xcPairs,
                             struct migrate_struct *migrate)
{
    size_t nbytes;
    int ierr;
    memset(migrate, 0, sizeof(struct migrate_struct));
    ierr = xcloc_migrate_checkParameters(ntables, ngrd, nxc, lxc,
                                         chunkSize, dt, 
                                         xcPairs);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Invalid inputs\n", __func__);
        return -1;
    }
    // Copy variables
    xcloc_migrate_setChunkSize(chunkSize, migrate);
    migrate->ntables = ntables;
    migrate->dt = dt;
    migrate->nxc = nxc;
    migrate->lxc = lxc;
    migrate->ngrd = ngrd;
    migrate->ldxc = XCLOC_PAD_FLOAT32(XCLOC_MEM_ALIGNMENT, migrate->lxc);
    migrate->mgrd = XCLOC_PAD_FLOAT32(XCLOC_MEM_ALIGNMENT, migrate->ngrd);
    // Allocate space and NULL it out
    nbytes = (size_t) (migrate->nxc*2*sizeof(int));
    migrate->xcPairs = (int *) aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
    if (migrate->xcPairs == NULL)
    {
        fprintf(stderr, "%s: Failed to set space for xcPairs\n", __func__);
        return -1;
    }
    ippsCopy_32s(xcPairs, migrate->xcPairs, 2*migrate->nxc);

    nbytes = (size_t) (migrate->nxc*migrate->ldxc)*sizeof(float);
    migrate->xcs = (float *) aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
    if (migrate->xcs == NULL)
    {
        fprintf(stderr, "%s: Failed to set space for xcs\n", __func__);
        return -1;
    }
    memset(migrate->xcs, 0, nbytes);

    nbytes = (size_t) migrate->mgrd*sizeof(float);
    migrate->migrate = (float *) aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
    if (migrate->migrate == NULL)
    {
        fprintf(stderr, "%s: Failed to set space for migrate\n", __func__);
        return -1;
    }

    nbytes = (size_t) (migrate->mgrd*ntables)*sizeof(int);
    migrate->ttimesInt = (int *) aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
    if (migrate->ttimesInt == NULL)
    {
        fprintf(stderr, "%s: Failed to set space for tables\n", __func__);
        return -1;
    }
    memset(migrate->ttimesInt, 0, nbytes);
    return 0;
}
//============================================================================//
int xcloc_migrate_setChunkSize(const int chunkSize,
                               struct migrate_struct *migrate)
{
    migrate->chunkSize = 2048;
    if (chunkSize%XCLOC_MEM_ALIGNMENT != 0)
    {
        fprintf(stderr, "%s: chunkSize=%d not divisible by alignment=%d\n",
                __func__, chunkSize, XCLOC_MEM_ALIGNMENT);
        return -1;
    }
    migrate->chunkSize = chunkSize;
    return 0;
}
//============================================================================//
/*!
 * @brief Deallocates memory on the migration structure and sets all variables
 *        to 0.
 *
 * @param[in,out] migrate   On input this is the initialized migration
 *                          structure. \n
 *                          On exit all variables have been set to zero and
 *                          memory freed.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrate_finalize(struct migrate_struct *migrate)
{
    if (migrate->migrate != NULL){free(migrate->migrate);}
    if (migrate->ttimesInt != NULL){free(migrate->ttimesInt);}
    if (migrate->xcPairs != NULL){free(migrate->xcPairs);}
    if (migrate->xcs != NULL){free(migrate->xcs);}
    memset(migrate, 0, sizeof(struct migrate_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the ixc'th cross-correlation on the migration structure.
 *
 * @param[in] lxc          Length of the cross-correlation.  This must be
 *                         equal to migrate->lxc. 
 * @param[in] ixc          Index of cross-correlation to add to the migration
 *                         structure.  This must be in the inclusive range 
 *                         [0, migrate->nxc-1]. 
 * @param[in] precision    Precision of the cross-correlations.  This can be
 *                         single or double precision but internally it will
 *                         be converted to single precision.
 * @param[in,out] migrate  On input contains the initialized the migration
 *                         structure. \n
 *                         On exit the ixc'th cross-correlation has been 
 *                         added to the structure.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrate_setCrossCorrelation(const int lxc, const int ixc,
                                      const void *__restrict__ xcs, 
                                      const enum xclocPrecision_enum precision,
                                      struct migrate_struct *migrate)
{
    Ipp64f *pSrc64;
    Ipp32f *pSrc32;
    IppStatus status;
    int indx;
    if (lxc != migrate->lxc || xcs == NULL || ixc < 0 || ixc >= migrate->nxc)
    {   
        if (lxc != migrate->lxc)
        {
            fprintf(stderr, "%s: lxc=%d != migrate->lxc=%d\n",
                    __func__, lxc, migrate->lxc); 
        }
        if (xcs == NULL){fprintf(stderr, "%s: xcs is NULL\n", __func__);}
        if (ixc < 0 || ixc >= migrate->nxc)
        {
            fprintf(stderr, "%s: ixc=%d must be in range [0,%d]\n",
                    __func__, ixc, migrate->nxc-1);
        }
        return -1; 
    }
    indx = ixc*migrate->ldxc;
    if (precision == XCLOC_DOUBLE_PRECISION)
    {
        pSrc64 = (Ipp64f *) xcs;
        status = ippsConvert_64f32f(pSrc64, &migrate->xcs[indx], lxc);
        if (status != ippStsNoErr)
        {
            fprintf(stderr, "%s: Error converting ixc=%d\n", __func__, ixc);
            return -1;
        }
    }
    else if (precision == XCLOC_SINGLE_PRECISION)
    {
        pSrc32 = (Ipp32f *) xcs;
        status = ippsCopy_32f(pSrc32, &migrate->xcs[indx], lxc);
        if (status != ippStsNoErr)
        {
            fprintf(stderr, "%s: Error copying ixc=%d\n",__func__, ixc);
            return -1;
        }
    }
    else
    {
        fprintf(stderr, "%s: Invalid precision: %d\n",
                 __func__, (int) precision); 
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets all the cross-correlograms on the migration structure.
 *
 * @param[in] ldxc         Leading dimension of xcs.
 * @param[in] lxc          Number of points in each cross-correlation.  This
 *                         must equal migrate->lxc.
 * @param[in] nxc          Number of cross-correlations.  This should equal but
 *                         cannot exceed migrate->nxc.
 * @param[in] xcs          Cross-correlations to place on the migration
 *                         structure.  This is an array of dimension
 *                         [ldxc x nxc] with leading dimension ldxc.
 * @param[in] precision    Precision of cross-correlations.  This can be
 *                         single or double precision but will be converted
 *                         to single precision internally.
 * @param[in,out] migrate  On input this contains the initialized migration
 *                         structure. \n
 *                         On exit contains the cross-correlations.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 * 
 */
int xcloc_migrate_setCrossCorrelations(const int ldxc, const int lxc,
                                       const int nxc, 
                                       const void *__restrict__ xcs,
                                       const enum xclocPrecision_enum precision,
                                       struct migrate_struct *migrate)
{
    size_t offset;
    int ierr, ixc;
    ierr = 0;
    if (ldxc < lxc || nxc > migrate->nxc)
    {
        if (ldxc < lxc)
        {
            fprintf(stderr, "%s: ldxc=%d < lxc=%d\n", __func__, ldxc, lxc);
        }
        if (nxc > migrate->nxc)
        {
            fprintf(stderr, "%s: nxc=%d > migrate->nxc=%d\n",
                    __func__, nxc, migrate->nxc);
        } 
        if (xcs == NULL){fprintf(stderr, "%s: xcs is NULL\n", __func__);}
        return -1;
    }
    if (nxc != migrate->nxc)
    {
        fprintf(stdout, "%s: Warning some cross-correlations will not be set\n",
                __func__);
        ippsZero_32f(migrate->xcs, migrate->nxc*migrate->ldxc);
    }
    for (ixc=0; ixc<nxc; ixc++)
    {
        offset = (size_t) (ixc*ldxc)*sizeof(float);
        if (precision == XCLOC_DOUBLE_PRECISION)
        {
            offset = ixc*ldxc*sizeof(double);
        }
        ierr = xcloc_migrate_setCrossCorrelation(lxc, ixc, &xcs[offset],
                                                 precision, migrate);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to set %d'th correlation\n",
                     __func__, ixc);
            return -1;
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets a the travel time table on the grid.  Note that the table will
 *        be divided by the sampling period and converted to an integer.
 *
 * @param[in] it           Travel time table index.  This must be in range
 *                         [0, migrate->ntables - 1].
 * @param[in] ngrd         Number of grid points.  This should equal
 *                         migrate->ngrd.
 * @parma[in] precision    Precision of the table.
 * @param[in] ttimes       Travel times (second) to put on the structure.
 *                         This is an array of dimension [ngrd] whose type,
 *                         float or double, is defined by precision.
 *
 * @param[in,out] migrate  On input holds the number of tables and the
 *                         space to collect the tables. \n
 *                         On exit holds the it'th travel-time table.
 *
 * @result 0 indicates success.
 * 
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrate_setTable(const int it, const int ngrd,
                           const enum xclocPrecision_enum precision,
                           const void *__restrict__ ttimes,
                           struct migrate_struct *migrate)
{
    const double *ttimes64;
    const float *ttimes32;
    double dt;
    int i, indx;
    if (ngrd > migrate->mgrd || ttimes == NULL)
    {   
        if (ngrd > migrate->mgrd)
        {
            fprintf(stderr, "%s: Insufficient space for table %d (%d > %d)\n",
                    __func__, it, ngrd, migrate->mgrd);
        }
        if (ttimes == NULL){fprintf(stderr, "%s: ttimes is NULL\n", __func__);}
        return -1; 
    }   
    if (ngrd != migrate->ngrd)
    {   
        fprintf(stdout, "%s: ngrd=%d != migrate->ngrd=%d; expect problems\n",
                __func__, ngrd, migrate->ngrd);
    }
    indx = it*migrate->mgrd;
    dt = migrate->dt;
    if (precision == XCLOC_SINGLE_PRECISION)
    {
        ttimes32 = (const float *) ttimes;
        for (i=0; i<ngrd; i++)
        {
            migrate->ttimesInt[indx+i] = (int) (round((double) ttimes32[i]/dt));
        } 
    }
    else if (precision == XCLOC_DOUBLE_PRECISION)
    {
        ttimes64 = (const double *) ttimes;
        for (i=0; i<ngrd; i++)
        {
            migrate->ttimesInt[indx+i] = (int) (round(ttimes64[i]/dt));
        }
    }
    else
    {
        fprintf(stderr, "%s: Invalid precision %d\n", __func__, (int)precision);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets a double precision travel time table on the grid.  Note that
 *        the double precision grid will be converted to single preicison.
 *
 * @param[in] it           Travel time table index.  This must be in range
 *                         [0, migrate->ntables - 1].
 * @param[in] ngrd         Number of grid points.  This should equal
 *                         migrate->ngrd.
 * @param[in] ttimes       Travel times (second) to put on the structure.
 *                         This is an array of dimension [ngrd].
 * @param[in,out] migrate  On input holds the number of tables and the
 *                         space to collect the tables. \n
 *                         On exit holds the it'th travel-time table.
 *
 * @result 0 indicates success.
 * 
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrate_setTable64f(const int it, const int ngrd,
                              const double *__restrict__ ttimes, 
                              struct migrate_struct *migrate)
{
    double dt;
    int i, indx;
    if (ngrd > migrate->mgrd || ttimes == NULL)
    {
        if (ngrd > migrate->mgrd)
        {
            fprintf(stderr, "%s: Insufficient space for table %d (%d > %d)\n",
                    __func__, it, ngrd, migrate->mgrd);
        }
        if (ttimes == NULL){fprintf(stderr, "%s: ttimes is NULL\n", __func__);}
        return -1;
    }
    if (ngrd != migrate->ngrd)
    {
        fprintf(stdout, "%s: ngrd=%d != migrate->ngrd=%d; expect problems\n",
                __func__, ngrd, migrate->ngrd);
    }
    indx = it*migrate->mgrd;
    dt = migrate->dt;
    for (i=0; i<ngrd; i++)
    {
        migrate->ttimesInt[indx+i] = (int) (round(ttimes[i]/dt));
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets a single precision travel time table on the grid.
 *
 * @param[in] it           Travel time table index.  This must be in range
 *                         [0, migrate->ntables - 1].
 * @param[in] ngrd         Number of grid points.  This should equal
 *                         migrate->ngrd.
 * @param[in] ttimes       Travel times (second) to put on the structure.
 *                         This is an array of dimension [ngrd].
 * @param[in,out] migrate  On input holds the number of tables and the
 *                         space to collect the tables. \n
 *                         On exit holds the it'th travel-time table.
 *
 * @result 0 indicates success.
 * 
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrate_setTable32f(const int it, const int ngrd,
                              const float *__restrict__ ttimes, 
                              struct migrate_struct *migrate)
{
    double dt;
    int i, indx;
    if (ngrd > migrate->mgrd || ttimes == NULL)
    {   
        if (ngrd > migrate->mgrd)
        {
            fprintf(stderr, "%s: Insufficient space for table %d (%d > %d)\n",
                    __func__, it, ngrd, migrate->mgrd);
        }
        if (ttimes == NULL){fprintf(stderr, "%s: ttimes is NULL\n", __func__);}
        return -1; 
    }   
    if (ngrd != migrate->ngrd)
    {   
        fprintf(stdout, "%s: ngrd=%d != migrate->ngrd=%d; expect problems\n",
                __func__, ngrd, migrate->ngrd);
    }
    indx = it*migrate->mgrd;
    dt = migrate->dt;
    for (i=0; i<ngrd; i++)
    {
        migrate->ttimesInt[indx+i] = (int) (round((double) ttimes[i]/dt));
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Zeros out the migration image.  This should be called prior to
 *        every diffraction stack migration image generation phase.
 *
 * @param[in,out] migrate  On input this is the initialized migration
 *                         structure. \n
 *                         On exit the migration image has been set to 0.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrate_setImageToZero(struct migrate_struct *migrate)
{
    if (migrate->migrate == NULL || migrate->ngrd < 1)
    {
        fprintf(stderr, "%s: Migration structure not initialized\n",
                __func__);
        return -1;
    }
    ippsZero_32f(migrate->migrate, migrate->ngrd);
    return 0; 
}
//============================================================================//
/*!
 * @brief Computes the diffraction stack migration image of the
 *        cross-correlograms. 
 *
 * @param[in,out] migrate  On input holds the full initialized migration
 *                         structure with travel time tables set. \n
 *                         On exit holds the diffraction stack image 
 *                         of the cross-correlations.
 *  
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrate_computeMigrationImage(struct migrate_struct *migrate)
{
    float *image    __attribute__((aligned(XCLOC_MEM_ALIGNMENT)));
    float *imagePtr __attribute__((aligned(XCLOC_MEM_ALIGNMENT)));
    float *xc       __attribute__((aligned(XCLOC_MEM_ALIGNMENT)));
    float *xcs      __attribute__((aligned(XCLOC_MEM_ALIGNMENT)));
    int *xcPairs    __attribute__((aligned(XCLOC_MEM_ALIGNMENT)));
    int *ttimesInt, *ttimesIntPtr1, *ttimesIntPtr2, indxXC, lxc2; 
    int chunkSize, ib, igrd, it1, it2, ixc, ldxc, lxc, mgrd, ngrd, ngrdLoc, nxc;
    if (migrate->nxc == 0){return 0;} // Nothing to do
    if (migrate->chunkSize%XCLOC_MEM_ALIGNMENT != 0 || migrate->ngrd < 1 ||
        migrate->dt <= 0.0 || migrate->migrate == NULL ||
        migrate->xcPairs == NULL || migrate->ttimesInt == NULL)
    {
        if (migrate->ngrd < 1)
        {
            fprintf(stderr, "%s: No grid points in migration structure\n",
                    __func__);
        }
        if (migrate->chunkSize%XCLOC_MEM_ALIGNMENT != 0)
        {
            fprintf(stderr, "%s: chunkSize=%d not divisible by alignment=%d\n",
                    __func__, migrate->chunkSize, XCLOC_MEM_ALIGNMENT);
        }
        if (migrate->dt <= 0.0)
        {
            fprintf(stderr, "%s: Sampling period %e must be positive\n",
                    __func__, migrate->dt);
        }
        if (migrate->migrate == NULL)
        {
            fprintf(stderr, "%s: migration image is NULL\n", __func__);
        }
        if (migrate->xcPairs == NULL)
        {
            fprintf(stderr, "%s: xcPairs is NULL\n", __func__);
        }
        if (migrate->ttimesInt == NULL)
        {
            fprintf(stderr, "%s: ttimesInt is NULL\n", __func__);
        }
        return -1; 
    } 
    // Set some constants and get pointers
    mgrd = migrate->mgrd;
    chunkSize = migrate->chunkSize;
    lxc = migrate->lxc;
    ldxc = migrate->ldxc;
    xcs = migrate->xcs;
    xcPairs = migrate->xcPairs;
    nxc = migrate->nxc;
    lxc2 = lxc/2; // cross-correlation is lxc = 2*npts - 1; i.e. odd
    ngrd = migrate->ngrd;
    image = migrate->migrate;
    ttimesInt = migrate->ttimesInt;
    // The problem with a naive implementation is that the the inner most 
    // loop stacking loop is too efficient which makes the program memory
    // bound.  Thus, a cache-blocking strategy is adopted so that the
    // read/write memory location (imagePtr) is updated slowest.  This is
    // improves performance because imagePtr requires a read and write which
    // can be 2x slower than a read.  A further optimization may express
    // itself in the natural ordering of the XC pairs.  Recall, that the 
    // XC matrix is created by moving along the columns of an upper
    // triangular matrix.  Thus, the first travel-time table pointer also
    // tends to change slowly.  Consequently, the fastest memory update 
    // is related to the second travel-time pointer.  It is important that
    // the chunksize be at least a page width so that the threads do not
    // try to access the same page containing imagePtr (this can yield false
    // sharing).
    #pragma omp parallel for default(none) \
     firstprivate(lxc2, mgrd, ngrd, nxc) \
     private(ib, igrd, imagePtr, indxXC, it1, it2, ixc, \
             ngrdLoc, ttimesIntPtr1, ttimesIntPtr2, xc)  \
     shared(chunkSize, image, ldxc, ttimesInt, xcPairs, xcs)
    for (ib=0; ib<ngrd; ib=ib+chunkSize)
    {
        ngrdLoc = MIN(ngrd - ib, chunkSize);
        imagePtr = (float *) &image[ib];
        ippsZero_32f(imagePtr, ngrdLoc); // 0 out the image in this chunk 
        // Loop on cross-correlation pairs
        for (ixc=0; ixc<nxc; ixc++)
        {
            it1 = xcPairs[2*ixc];
            it2 = xcPairs[2*ixc+1];
            xc = (float *) &xcs[ldxc*ixc];
            ttimesIntPtr1 = (int *) &ttimesInt[it1*mgrd+ib];
            ttimesIntPtr2 = (int *) &ttimesInt[it2*mgrd+ib];
            // Update image
            #pragma omp simd \
             aligned(ttimesIntPtr1,ttimesIntPtr2,imagePtr: XCLOC_MEM_ALIGNMENT)
            for (igrd=0; igrd<ngrdLoc; igrd++)
            {
                // Compute the differential times.  Recall that the correlations
                // imply the time to shift trace 2 to trace 1.  To understand
                // this imagine  an arrival at trace 1 and trace 2 in a
                // homogeneous medium where trace 1 is closer to the source.  
                // The cross-correlation in this instance will yield achieve its
                // maximum in the negative or acausal part - because trace 2
                // must be shifted a negative number of points (i.e., slid back
                // in time) as to align with trace 1.  Because the
                // back-projection works through travel times the
                // cross-correlation must be viewed in differential times.  In
                // the scenario described, trace 2 has a larger travel time than
                // trace 1 and the maximum of the cross-correlation is in the
                // acausal part of the signal.  Thus, to access the acausal part
                // of the  cross-correlation travelTime1 - travelTime2 must be
                // computed where travelTime1 < travelTime2.  
                //
                // Compute the index in the cross-correlation.  Note, that this
                // approximates the following code: 
                //  findxXC = roundf(dti*(ttimesPtr1[igrd] - ttimesPtr2[igrd]));
                //  indxXC = lxc2 + (int) findxXC;
                // by computing
                //  findxXCnew = roundf(dti*ttimesPtr1[igrd])
                //             - roundf(dti*ttimesPtr2[igrd])
                //  indxXCnew = (int) findxXCnew 
                // This is approximate because the roundf function does not
                // distribute for pathological cases (e.g., consider
                // roundf(1.0 - 0.5) = 1 and roundf(1.0) - roundf(0.5) = 0.)
                indxXC = lxc2 + (ttimesIntPtr1[igrd] - ttimesIntPtr2[igrd]);
                // Stack the cross-correlation value
                imagePtr[igrd] = imagePtr[igrd] + xc[indxXC];
            }
        } // Loop on cross-correlation pairs
    } // Loop on cross-correlation blocks
    xcs = NULL;
    xcPairs = NULL;
    ttimesInt = NULL;
    return 0;
} 
//============================================================================//
/*!
 * @brief Convenience function to compute the stack of individual 
 *        diffraction migration images. 
 *
 * @param[in] nxc           Number of cross-correlations.
 * @param[in] xcPairs       Maps from the ixc'th cross-correlation pair the
 *                          transform pair table IDs.  This is an array of
 *                          dimension [2 x nxc] with leading dimension 2.
 * @param[in] ldxc          Leading dimension of cross-correlations.
 * @param[in] lxc           Length of the cross-correlations.  This should be
 *                          an odd number.
 * @param[in] dt            Sampling period (seconds) of the signals.
 * @param[in] xcs           Cross-correlograms to migrate.  This is an array
 *                          of dimension [ldxc x nxc] with leading dimension
 *                          ldxc.
 *
 * @param[in,out] migrate   On input holds the full initialized migration
 *                          structure with travel time tables set. \n
 *                          On exit holds the diffraction stack image 
 *                          of the cross-correlations.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */ 
int xcloc_migrate_computeXCDSMImage(struct migrate_struct *migrate)
{
    int ierr, ierrAll, it1, it2, ixc, nxc;
    // Null out the migration image
    nxc = migrate->nxc;
    ierr = xcloc_migrate_setImageToZero(migrate);
    ierrAll = 0;
    // Stack the images
    for (ixc=0; ixc<nxc; ixc++)
    {
        it1 = migrate->xcPairs[2*ixc];
        it2 = migrate->xcPairs[2*ixc+1];
        ierr = xcloc_migrate_updateXCDSMImage(it1, it2, ixc, migrate);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error updating xcpair %d\n", __func__, ixc);
            ierrAll = ierrAll + 1;
        }
    }
    if (ierrAll != 0)
    {
        fprintf(stderr, "%s: %d errors encountered\n", __func__, ierrAll);
    }
    return ierrAll;
}
//============================================================================//
/*!
 * @brief Updates the diffraction stack migration image of the 
 *        cross-correlogram signal xc.  In this instance the 
 *        cross-correlation corresponds to:
 *        \f$
 *          s_1 \star s_2 (\tau) = \int_\infty^\infty s_1(t) s_2(\tau + t) dt
 *        \f$.
 *        This indicates that the maximum of the cross-correlation
 *        is the shift that number of samples to move signal 2 so that it
 *        best aligns with signal 1.
 *
 * @param[in] it1          Index of first travel time table.  This must be
 *                         in the range [0, migrate->ntables - 1].
 * @param[in] it2          Index of second travel time table.  This must be
 *                         in the range [0, migrate->ntables - 1].
 *
 * @param[in,out] migrate  On input contains the initialized migration image. \n
 *                         On exit the migration image has been updated with
 *                         with the diffraction stack of the input signal.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *        
 */
int xcloc_migrate_updateXCDSMImage(const int it1, const int it2,
                                   const int ixc,
                                   struct migrate_struct *migrate)
{
    int *ttimesIntPtr1 __attribute__((aligned(XCLOC_MEM_ALIGNMENT)));
    int *ttimesIntPtr2 __attribute__((aligned(XCLOC_MEM_ALIGNMENT)));
    float *migratePtr __attribute__((aligned(XCLOC_MEM_ALIGNMENT)));
    float *xc __attribute__((aligned(XCLOC_MEM_ALIGNMENT)));
    int igrd, indxXC, lxc, lxc2, ngrd;
    if ((it1 < 0 || it1 > migrate->ntables - 1) ||
        (it2 < 0 || it2 > migrate->ntables - 1))
    {
        if (it1 < 0 || it1 > migrate->ntables - 1)
        {
            fprintf(stderr, "%s: it1=%d must be in range=[0,%d]\n", 
                    __func__, it1, migrate->ntables-1);
        }
        if (it2 < 0 || it1 > migrate->ntables - 1)
        {
            fprintf(stderr, "%s: it2=%d must be in range=[0,%d]\n", 
                    __func__, it2, migrate->ntables-1);
        }
        return -1;
    }
    lxc = migrate->lxc;
    xc = (float *) &migrate->xcs[migrate->ldxc*ixc];
    lxc2 = lxc/2; // cross-correlation is lxc = 2*npts - 1; i.e. odd
    ngrd = migrate->ngrd;
    migratePtr = (float *) &migrate->migrate[0];
    ngrd = migrate->ngrd;
    ttimesIntPtr1 = (int *) &migrate->ttimesInt[it1*migrate->mgrd];
    ttimesIntPtr2 = (int *) &migrate->ttimesInt[it2*migrate->mgrd];
/*
    //dti = (float) (1.0/(double) dt);
    lxc2 = lxc/2; // cross-correlation is lxc = 2*npts - 1; i.e. odd
    ngrd = migrate->ngrd;
    ttimesIntPtr1 = (int *) &migrate->ttimesInt[it1*migrate->mgrd];
    ttimesIntPtr2 = (int *) &migrate->ttimesInt[it2*migrate->mgrd];
    migratePtr = (float *) &migrate->migrate[0];
*/
    #pragma omp parallel for simd default(none) \
            shared(lxc, migratePtr, ngrd, ttimesIntPtr1, ttimesIntPtr2) \
            firstprivate(lxc2, xc) \
            private(igrd, indxXC) \
            aligned(ttimesIntPtr1,ttimesIntPtr2,migratePtr: XCLOC_MEM_ALIGNMENT)
    for (igrd=0; igrd<ngrd; igrd++)
    {
        // Compute the differential times.  Recall that the correlations imply
        // the time to shift trace 2 to trace 1.  To understand this imagine 
        // an arrival at trace 1 and trace 2 in a homogeneous medium where
        // trace 1 is closer to the source.  The cross-correlation in this
        // instance will yield achieve its maximum in the negative or
        // acausal part - because trace 2 must be shifted a negative number
        // of points (i.e., slid back in time) as to align with trace 1.
        // Because the back-projection works through travel times the
        // cross-correlation must be viewed in differential times. 
        // In the scenario described, trace 2 has a larger travel time than
        // trace 1 and the maximum of the cross-correlation is in the acausal
        // part of the signal.  Thus, to access the acausal part of the 
        // cross-correlation travelTime1 - travelTime2 must be computed where
        // travelTime1 < travelTime2.  
        // Original code:
        //diffTime = ttimesPtr1[igrd] - ttimesPtr2[igrd];
        //indxXC = lxc2 + (int) (roundf(diffTime*dti));
        // Compute the index in the cross-correlation.  Works well on AVX2.
        //findxXC = roundf(dti*(ttimesPtr1[igrd] - ttimesPtr2[igrd]));
        //indxXC = lxc2 + (int) findxXC;
        indxXC = lxc2 + (ttimesIntPtr1[igrd] - ttimesIntPtr2[igrd]);
        //if (indxXC < 0 || indxXC > lxc){printf("problem\n");}
        // Stack the cross-correlation value
        migratePtr[igrd] = migratePtr[igrd] + xc[indxXC];
    }
    migratePtr = NULL;
    ttimesIntPtr1 = NULL;
    ttimesIntPtr2 = NULL;
    return 0;
}
//============================================================================//
/*!
 * @brief Gets a pointer to the migration image on the mgirate structure.
 *
 * @param[in] ngrd      Number of grid points in image.
 * @param[in] migrate   Structure with migration image.
 *
 * @param[out] ierr     0 indicates success.
 *
 * @result Pointer to the migration image.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
const float *xcloc_migrate_getImagePointer(const int ngrd,
                                           const struct migrate_struct migrate,
                                           int *ierr)
{
    const float *ptr = NULL;
    *ierr = 0;
    if (ngrd != migrate.ngrd)
    {
        fprintf(stderr, "%s: Error ngrd=%d != migrate.ngrd=%d\n",
                __func__,  ngrd, migrate.ngrd); 
        return ptr;
    }
    ptr = migrate.migrate;
    return ptr;
}
//============================================================================//
/*!
 * @brief Gets a copy of the migration image.
 *
 * @param[in] ngrd      Number of grid points in image.
 * @param[in] migrate   Structure with migration image.
 *
 * @param[out] image    The migration image.  This is an array of dimension
 *                      [ngrd].
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrate_getImage32f(const int ngrd,
                              const struct migrate_struct migrate,
                              float *__restrict__ image)
{
    const float *ptr;
    int ierr;
    if (ngrd < 1 || image == NULL)
    {
        if (ngrd < 1){fprintf(stderr, "%s: No grid points\n", __func__);}
        if (image == NULL){fprintf(stderr, "%s: image is NULL\n", __func__);}
        return 0;
    }
    ptr = xcloc_migrate_getImagePointer(ngrd, migrate, &ierr);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to get image pointer\n", __func__);
        return -1;
    }
    ippsCopy_32f(ptr, image, ngrd);
    ptr = NULL;
    return 0;
}
//============================================================================//
/*!
 * @brief Gets a copy of the migration image.
 *
 * @param[in] ngrd      Number of grid points in image.
 * @param[in] migrate   Structure with migration image.
 *
 * @param[out] image    The migration image.  This is an array of dimension
 *                      [ngrd].
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrate_getImage64f(const int ngrd,
                              const struct migrate_struct migrate,
                              double *__restrict__ image)
{
    const float *ptr;
    int ierr;
    if (ngrd < 1 || image == NULL)
    {
        if (ngrd < 1){fprintf(stderr, "%s: No grid points\n", __func__);}
        if (image == NULL){fprintf(stderr, "%s: image is NULL\n", __func__);}
        return 0;
    }
    ptr = xcloc_migrate_getImagePointer(ngrd, migrate, &ierr);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to get image pointer\n", __func__);
        return -1; 
    }
    ippsConvert_32f64f(ptr, image, ngrd);
    ptr = NULL;
    return 0;
}
//============================================================================//
/*!
 * @brief This is a utility routine which will compute the maximum differential
 *        time.  This function will help a user determine the maximum amount
 *        of samples in their cross-correlations to avoid a segfault.
 *
 * @param[in] migrate    The initialized migration structure with the travel
 *                       time tables and cross-correlation pairs.
 *
 * @param[out] absMaxDT  Absolute value of the max differential travel-time
 *                       in seconds.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_migrate_computeMaxDifferentialTime(
    const struct migrate_struct migrate, double *absMaxDT)
{
    int *xcPairs    __attribute__((aligned(XCLOC_MEM_ALIGNMENT)));
    int *ttimesInt, *ttimesIntPtr1, *ttimesIntPtr2, indxXC;
    int chunkSize, ib, igrd, it1, it2, ixc, mgrd, ngrd, ngrdLoc, nxc;
    *absMaxDT = 0.0;
    if (migrate.nxc == 0){return 0;} // Nothing to do
    if (migrate.chunkSize%XCLOC_MEM_ALIGNMENT != 0 || migrate.ngrd < 1 ||
        migrate.dt <= 0.0 || migrate.xcPairs == NULL ||
        migrate.ttimesInt == NULL)
    {
        if (migrate.ngrd < 1)
        {
            fprintf(stderr, "%s: No grid points in migration structure\n",
                    __func__);
        }
        if (migrate.chunkSize%XCLOC_MEM_ALIGNMENT != 0)
        {
            fprintf(stderr, "%s: chunkSize=%d not divisible by alignment=%d\n",
                    __func__, migrate.chunkSize, XCLOC_MEM_ALIGNMENT);
        }
        if (migrate.dt <= 0.0)
        {
            fprintf(stderr, "%s: Sampling period %e must be positive\n",
                    __func__, migrate.dt);
        }
        if (migrate.xcPairs == NULL)
        {
            fprintf(stderr, "%s: xcPairs is NULL\n", __func__);
        }
        if (migrate.ttimesInt == NULL)
        {
            fprintf(stderr, "%s: ttimesInt is NULL\n", __func__);
        }
        return -1; 
    }
    // Set some constants and get pointers
    indxXC = 0;
    mgrd = migrate.mgrd;
    chunkSize = migrate.chunkSize;
    xcPairs = migrate.xcPairs;
    nxc = migrate.nxc;
    ngrd = migrate.ngrd;
    ttimesInt = migrate.ttimesInt;
    // Loop on grid chunks
    #pragma omp parallel for default(none) \
     firstprivate(mgrd, ngrd, nxc) \
     private(ib, igrd, it1, it2, ixc, \
             ngrdLoc, ttimesIntPtr1, ttimesIntPtr2)  \
     shared(chunkSize, ttimesInt, xcPairs) \
     reduction(max: indxXC)
    for (ib=0; ib<ngrd; ib=ib+chunkSize)
    {
        ngrdLoc = MIN(ngrd - ib, chunkSize);
        // Loop on cross-correlation pairs
        for (ixc=0; ixc<nxc; ixc++)
        {
            it1 = xcPairs[2*ixc];
            it2 = xcPairs[2*ixc+1];
            ttimesIntPtr1 = (int *) &ttimesInt[it1*mgrd+ib];
            ttimesIntPtr2 = (int *) &ttimesInt[it2*mgrd+ib];
            // Update max differential time 
            #pragma omp simd reduction(max:indxXC) \
             aligned(ttimesIntPtr1,ttimesIntPtr2: XCLOC_MEM_ALIGNMENT)
            for (igrd=0; igrd<ngrdLoc; igrd++)
            {
                indxXC = MAX(indxXC,
                             abs(ttimesIntPtr1[igrd] - ttimesIntPtr2[igrd]));
            }
        } // Loop on cross-correlation pairs
    } // Loop on cross-correlation blocks
    xcPairs = NULL;
    ttimesInt = NULL;
    *absMaxDT = (double) indxXC*migrate.dt;
    return 0;
}
