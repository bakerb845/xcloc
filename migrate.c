#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "xcloc_migrate.h"
#include <ipps.h>

#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

/*!
 * @brief Initializes the gridded migration structure.
 *
 * @param[in] ntables   Number of travel time tables.
 * @param[in] ngrd      Number of grid points.
 * @param[in] dt        Sampling period (s) of the signals to migrate.
 * @param[in] nxc       Number of cross-correlations.
 * @param[in] lxc       Length of cross-correlations.
 *
 * @param[out] migrate  Migration structure with space allocated for the 
 *                      migration image and the travel-time tables.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under MIT.
 *
 */
int xcloc_migrate_initialize(const int ntables, const int ngrd,
                             const int nxc, const int lxc,
                             const double dt, 
                             const int *__restrict__ xcPairs,
                             struct migrate_struct *migrate)
{
    size_t nbytes;
    int minIDx, maxIDx;
    memset(migrate, 0, sizeof(struct migrate_struct));
    if (ngrd < 1 || ntables < 2 || dt <= 0.0 || nxc < 1 ||
        lxc < 1 || xcPairs == NULL)
    {   
        if (ngrd < 1)
        {
            fprintf(stderr, "%s: ngrd=%d must be positive\n", __func__, ngrd);
        }
        if (ntables < 2)
        {
            fprintf(stderr, "%s: ntables=%d must at least 2\n",
                    __func__, ntables);
        }
        if (dt <= 0.0)
        {
            fprintf(stderr, "%s: dt=%f must be positive\n", __func__, dt);
        }
        if (nxc < 1)
        {
            fprintf(stderr, "%s: No cross-correlations\n", __func__);
        }
        if (lxc < 1)
        {
            fprintf(stderr, "%s: xc length=%d must be positive\n",
                    __func__, lxc);
        }
        if (xcPairs == NULL)
        {
            fprintf(stderr, "%s: xcPairs is NULL\n", __func__);
        }
        return -1;
    }
    // Get the min and max of xcPairs and make sure they are in range
    ippsMax_32s(xcPairs, 2*nxc, &maxIDx);
    ippsMin_32s(xcPairs, 2*nxc, &minIDx);
    if (minIDx < 0)
    {
        fprintf(stderr, "%s: min index in xcPairs=%d must be positive\n",
                __func__, minIDx);
        return -1;
    }
    if (maxIDx + 1 > ntables)
    {
        fprintf(stderr, "%s: xcPairs requires at least %d tables\n",
                __func__, maxIDx + 1);
        return -1;
    }
    // Copy variables
    migrate->ntables = ntables;
    migrate->dt = dt;
    migrate->nxc = nxc;
    migrate->lxc = lxc;
    migrate->ngrd = ngrd;
    migrate->ldxc = DALES_PAD_FLOAT32(DALES_MEM_ALIGNMENT, migrate->lxc);
    migrate->mgrd = DALES_PAD_FLOAT32(DALES_MEM_ALIGNMENT, migrate->ngrd);
    // Allocate space and NULL it out
    nbytes = (size_t) (migrate->nxc*2*sizeof(int));
    migrate->xcPairs = (int *) aligned_alloc(DALES_MEM_ALIGNMENT, nbytes);
    if (migrate->xcPairs == NULL)
    {
        fprintf(stderr, "%s: Failed to set space for xcPairs\n", __func__);
        return -1;
    }
    ippsCopy_32s(xcPairs, migrate->xcPairs, 2*migrate->nxc);

    nbytes = (size_t) (migrate->nxc*migrate->ldxc)*sizeof(float);
    migrate->xcs = (float *) aligned_alloc(DALES_MEM_ALIGNMENT, nbytes);
    if (migrate->xcs == NULL)
    {
        fprintf(stderr, "%s: Failed to set space for xcs\n", __func__);
        return -1;
    }
    memset(migrate->xcs, 0, nbytes);

    nbytes = (size_t) migrate->mgrd*sizeof(float);
    migrate->migrate = (float *) aligned_alloc(DALES_MEM_ALIGNMENT, nbytes);
    if (migrate->migrate == NULL)
    {
        fprintf(stderr, "%s: Failed to set space for migrate\n", __func__);
        return -1;
    }

    nbytes = (size_t) (migrate->mgrd*ntables)*sizeof(int);
    migrate->ttimesInt = (int *) aligned_alloc(DALES_MEM_ALIGNMENT, nbytes);
    if (migrate->ttimesInt == NULL)
    {
        fprintf(stderr, "%s: Failed to set space for tables\n", __func__);
        return -1;
    }
    memset(migrate->ttimesInt, 0, nbytes);
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
    if (precision == XCLOC_SINGLE_PRECISION)
    {
        pSrc64 = (Ipp64f *) xcs;
        ippsConvert_64f32f(pSrc64, &migrate->xcs[indx], lxc);
    }
    else if (precision == XCLOC_DOUBLE_PRECISION)
    {
        pSrc32 = (Ipp32f *) xcs;
        ippsCopy_32f(pSrc32, &migrate->xcs[indx], lxc);
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
        ierr = xcloc_migrate_setCrossCorrelation(lxc, ixc, &xcs[ixc*ldxc], 
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
    double dti;
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
    //ippsConvert_64f32f(ttimes, &migrate->ttimes[indx], ngrd);
    dti = 1.0/migrate->dt;
    for (i=0; i<ngrd; i++)
    {
        migrate->ttimesInt[indx+i] = (int) (round(ttimes[i]*dti));
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
    double dti;
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
    //ippsCopy_32f(ttimes, &migrate->ttimes[indx], ngrd);
    dti = 1.0/migrate->dt;
    for (i=0; i<ngrd; i++)
    {
        migrate->ttimesInt[indx+i] = (int) (round((double) ttimes[i]*dti));
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
 * @param[in] chunkSize     This is a tuning parameter and must be evenly
 *                          divisible by the memory alignment.  It is
 *                          recommended to use powers of 2 that are
 *                          of size x*pageWidth/sizeof(float) where x is a 
 *                          small positive integer and page width is the size
 *                          of a page on the target architecture. 
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
int xcloc_migrate_computeXCDSImage(const int chunkSize,
                                   const int nxcIn,
                                   const int *__restrict__ xcPairsDead,
                                   const int ldxc, const int lxc,
                                   const float dt,
                                   const float *__restrict__ xcs,
                                   struct migrate_struct *migrate)
{
    float *image    __attribute__((aligned(DALES_MEM_ALIGNMENT)));
    float *imagePtr __attribute__((aligned(DALES_MEM_ALIGNMENT)));
    float *xc       __attribute__((aligned(DALES_MEM_ALIGNMENT)));
    int *xcPairs    __attribute__((aligned(DALES_MEM_ALIGNMENT)));
    int *ttimesInt, *ttimesIntPtr1, *ttimesIntPtr2, indxXC, lxc2; 
    int ib, igrd, it1, it2, ixc, ngrd, ngrdLoc, nxc;
    if (chunkSize%DALES_MEM_ALIGNMENT != 0 || lxc < 1 ||
        dt <= 0.0f || xcs == NULL || migrate->ngrd <1 ||
        migrate->migrate == NULL)
    {
        if (chunkSize%DALES_MEM_ALIGNMENT != 0)
        {
            fprintf(stderr, "%s: chunkSize=%d not divisible by alignment=%d\n",
                    __func__, chunkSize, DALES_MEM_ALIGNMENT);
        }
        if (lxc < 1)
        {
            fprintf(stderr, "%s: lxc=%d must be positive\n", __func__, lxc);
        }
        if (dt <= 0.0f)
        {
            fprintf(stderr, "%s: dtf=%f must be positive\n", __func__, dt);
        }
        if (xcs == NULL)
        {
            fprintf(stderr, "%s: xcs is NULL\n", __func__);
        }
        if (migrate->ngrd < 1)
        {
            fprintf(stderr, "%s: No grid points in migration structure\n",
                    __func__);
        }
        if (migrate->migrate == NULL)
        {
            fprintf(stderr, "%s: migration image is NULL\n", __func__);
        }
        return -1; 
    } 
    // Set some constants and get pointers
    //dti = (float) (1.0/(double) dt);
    xcPairs = migrate->xcPairs;
    nxc = migrate->nxc;
    lxc2 = lxc/2; // cross-correlation is lxc = 2*npts - 1; i.e. odd
    ngrd = migrate->ngrd;
    image = migrate->migrate;
    //ttimes = migrate->ttimes;
    ttimesInt = migrate->ttimesInt;
    // The problem with a naive implementation is that the loop is the inner
    // most loop stacking loop is too efficient which makes the program memory.
    // Thus, we we chunk so that the read/write memory location (imagePtr)
    // is updated slowest (because a write is 2x slower than a read) then we
    // rely on the natural ordering of the XC pairs so that the row of the
    // XC pair matrix updates the memory address corresponding to ttimesPtr1
    // intermediately and the columns of ttimesPtr2 update memory addresses
    // fastest.  It is important that the chunkSize be at least a page width
    // so that the threads do not try to access the same page - particularly
    // when updating migratePtr (this is referred to as false sharing).
    #pragma omp parallel for default(none) \
     firstprivate(lxc2, ngrd, nxc) \
     private(ib, igrd, imagePtr, indxXC, it1, it2, ixc, \
             ngrdLoc, ttimesIntPtr1, ttimesIntPtr2, xc)  \
     shared(chunkSize, image, ttimesInt, xcPairs, migrate, xcs)
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
            ttimesIntPtr1 = (int *) &ttimesInt[it1*migrate->mgrd+ib];
            ttimesIntPtr2 = (int *) &ttimesInt[it2*migrate->mgrd+ib];
            // Update image
            #pragma omp simd \
             aligned(ttimesIntPtr1,ttimesIntPtr2,imagePtr: DALES_MEM_ALIGNMENT)
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
int xcloc_migrate_computeXCDSMImage(const int nxc,
                                    const int *__restrict__ xcPairs,
                                    const int ldxc, const int lxc,
                                    const float dt,
                                    const float *__restrict__ xcs,
                                    struct migrate_struct *migrate)
{
    const float *xc;
    int ierr, ierrAll, it1, it2, ixc;
    // Null out the migration image
    ierr = xcloc_migrate_setImageToZero(migrate);
    ierrAll = 0;
    // Stack the images
    for (ixc=0; ixc<nxc; ixc++)
    {
        it1 = xcPairs[2*ixc];
        it2 = xcPairs[2*ixc+1];
        xc = (float *) &xcs[ldxc*ixc];
        ierr = xcloc_migrate_updateXCDSMImage(it1, it2, lxc, dt, xc,
                                              migrate);
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
 * @param[in] lxc          Length of the cross-correlation.  This should be
 *                         an odd number.
 * @param[in] xc           Cross-correlogram of it1 correlated with it2.
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
                                   const int lxc, const float dt,
                                   const float *__restrict__ xc,
                                   struct migrate_struct *migrate)
{
    int *ttimesIntPtr1 __attribute__((aligned(DALES_MEM_ALIGNMENT)));
    int *ttimesIntPtr2 __attribute__((aligned(DALES_MEM_ALIGNMENT)));
    float *migratePtr __attribute__((aligned(DALES_MEM_ALIGNMENT)));
    int igrd, indxXC, lxc2, ngrd;
    if ((it1 < 0 || it1 > migrate->ntables - 1) ||
        (it2 < 0 || it2 > migrate->ntables - 1) ||
         lxc < 1 || dt <= 0.0f || xc == NULL)
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
        if (lxc < 1)
        {
            fprintf(stderr, "%s: lxc=%d must be positive\n", __func__, lxc);
        } 
        if (dt <= 0.0f)
        {
            fprintf(stderr, "%s: dtf=%f must be positive\n", __func__, dt);
        }
        if (xc == NULL)
        {
            fprintf(stderr, "%s: xc is NULL\n", __func__);
        }
        return -1;
    }
    //dti = (float) (1.0/(double) dt);
    lxc2 = lxc/2; // cross-correlation is lxc = 2*npts - 1; i.e. odd
    ngrd = migrate->ngrd;
    ttimesIntPtr1 = (int *) &migrate->ttimesInt[it1*migrate->mgrd];
    ttimesIntPtr2 = (int *) &migrate->ttimesInt[it2*migrate->mgrd];
    migratePtr = (float *) &migrate->migrate[0];
    #pragma omp parallel for simd default(none) \
            shared(lxc, migratePtr, ngrd, ttimesIntPtr1, ttimesIntPtr2) \
            firstprivate(lxc2, xc) \
            private(igrd, indxXC) \
            aligned(ttimesIntPtr1,ttimesIntPtr2,migratePtr: DALES_MEM_ALIGNMENT)
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
//                             Functions to deprecate                         //
//============================================================================//
int xcloc_migrate_migrateSignalsOnGrid(
    const int ld1, const int ld2,
    const int n1, const int n2, const int n3, 
    const int npts, const float dt, 
    const float *__restrict__ signal,
    const float *__restrict__ ttimes,
    float *__restrict__ migrate) 
{
    float *ttimesPtr, *migratePtr, dti;
    int i, indxModel, indx, j, k;
    dti = 1.0f/dt;
    // Loop on grid - e.g., z, y, x
    #pragma omp parallel default(none) \
            shared(migrate, ttimes, npts) \
            firstprivate(dti, signal) \
            private(ttimesPtr, migratePtr, i, j, k, indxModel, indx)
    {
    #pragma omp for collapse(2)
    for (k=0; k<n3; k=k++)
    {
        for (j=0; j<n2; j=j++)
        {
            indxModel = k*ld1*ld2 + j*ld2;
            ttimesPtr  = (float *) &ttimes[indxModel];
            migratePtr = (float *) &migrate[indxModel];
            #pragma omp simd aligned(ttimesPtr, migratePtr: 64)
            for (i=0; i<n1; i=i++)
            {
                 // Compute the index in the signal
                 indx = (int) roundf(ttimesPtr[i]*dti);
                 // Stack the cross-correlation value
                 migratePtr[i] = migratePtr[i] + signal[indx];
            } // Loop on dimension 1
        } // Loop on dimension 2
    } // Loop on dimension 3
    } // End parallel
    return 0;
}
