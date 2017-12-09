#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mkl.h>
#include <ipps.h>
#include "xcloc_rmsFilter.h"

/*!
 * @brief Initializes the RMS filter structure.
 *
 * @param[in] winLen     Number of samples in RMS window.  This should be
 *                       an odd number greater than 1.
 * @param[in] npts       Number of points in signal to filter.
 * @param[in] precision  Precision of filter to create - float or double.
 *
 * @param[out] rms       RMS FIR filter.
 *
 * @copyright Ben Baker distributed under MIT.
 *
 */
int xcloc_rmsFilter_initialize(const int winLen, const int npts,
                               const enum xclocPrecision_enum precision,
                               struct xcfftRMSFilter_struct *rms)
{
    IppsFIRSpec_64f *pSpec64 = NULL;
    IppsFIRSpec_32f *pSpec32 = NULL;
    Ipp64f *pTaps64 = NULL;
    Ipp32f *pTaps32 = NULL;
    Ipp8u *pBuf;
    double avg64;
    float avg32;
    int bufferSize, specSize, tapsLen;
    IppStatus status;
    memset(rms, 0, sizeof(struct xcfftRMSFilter_struct));
    if (winLen < 2)
    {
        fprintf(stdout, "%s: Will not low-pass filter signal\n", __func__);
        return 0;
    }
    tapsLen = winLen;
    // Make filter have an odd length for my sanity - phase shift is tapsLen/2
    if (winLen%2 == 0){tapsLen = tapsLen + 1;}
    if (tapsLen > npts)
    {
        fprintf(stderr, "%s: Error winLen=%d > ntps=%d\n",
                __func__, tapsLen, npts);
        return -1;
    }
    if (winLen < 1)
    {
        fprintf(stdout, "%s: Will not low-pass filter signal\n", __func__);
        return 0;
    }
    // Create the taps and normalize the FIR summation by 1/N
    if (precision == XCLOC_SINGLE_PRECISION)
    {
        pTaps32 = ippsMalloc_32f(tapsLen);
        avg32 = (float) (1.0/(double) tapsLen);
        ippsSet_32f(avg32, pTaps32, tapsLen);
        // Get the space
        status = ippsFIRSRGetSize(tapsLen, ipp32f, &specSize, &bufferSize);
        if (status != ippStsNoErr)
        {
            fprintf(stderr, "%s: Error getting state size\n", __func__);
            return -1;
        }
        // Initialize the filter
        pSpec32 = (IppsFIRSpec_32f *) ippsMalloc_8u(specSize);
        pBuf = ippsMalloc_8u(bufferSize);
        status = ippsFIRSRInit_32f(pTaps32, tapsLen, ippAlgAuto, pSpec32);
        if (status != ippStsNoErr)
        {
            fprintf(stderr, "%s: Error initializing FIR filter\n", __func__);
            //return -1;
        }
        // Set it on xcfft
        rms->pSpec = (void *) pSpec32;
        rms->pBuf = (void *) pBuf;
        rms->pTaps = (void *) pTaps32;
        rms->tapsLen = tapsLen;
        rms->bufferSize = bufferSize;
        rms->specSize = specSize;
    }
    else if (precision == XCLOC_DOUBLE_PRECISION)
    {
        pTaps64 = ippsMalloc_64f(tapsLen);
        avg64 = 1.0/(double) tapsLen;
        ippsSet_64f(avg64, pTaps64, tapsLen);
        // Get the space
        status = ippsFIRSRGetSize(tapsLen, ipp64f, &specSize, &bufferSize);
        if (status != ippStsNoErr)
        {
            fprintf(stderr, "%s: Error getting state size\n", __func__);
            return -1;
        }
        // Initialize the filter
        pSpec64 = (IppsFIRSpec_64f *) ippsMalloc_8u(specSize);
        pBuf = ippsMalloc_8u(bufferSize);
        status = ippsFIRSRInit_64f(pTaps64, tapsLen, ippAlgAuto, pSpec64);
        if (status != ippStsNoErr)
        {
            fprintf(stderr, "%s: Error initializing FIR filter\n", __func__);
            return -1;
        }
        // Set it on xcfft 
        rms->pSpec = (void *) pSpec64;
        rms->pBuf = (void *) pBuf;
        rms->pTaps = (void *) pTaps64;
        rms->tapsLen = tapsLen;
        rms->bufferSize = bufferSize;
        rms->specSize = specSize;
    }
    else
    {
        fprintf(stderr, "%s: Invalid precision: %d\n",
                __func__, (int) precision);
        return -1;
    }
    rms->npts = npts;
    rms->precision = precision;
    if (winLen > 1){rms->lfilter = true;}
    return 0;
}
//============================================================================//
/*!
 * @brief Releases memory on the RMS filter and sets variables to zero.
 *
 * @param[in,out] rms    On input this is the initialized RMS structure. \n
 *                       On exit memory has been freed and variables set to 0.
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_rmsFilter_finalize(struct xcfftRMSFilter_struct *rms)
{
    IppsFIRSpec_64f *pSpec64 = NULL;
    IppsFIRSpec_32f *pSpec32 = NULL;
    Ipp64f *pTaps64 = NULL;
    Ipp32f *pTaps32 = NULL;
    Ipp8u *pBuf;
    pBuf = (Ipp8u *) rms->pBuf;
    if (pBuf != NULL){ippsFree(pBuf);}
    if (rms->precision == XCLOC_SINGLE_PRECISION)
    {
        pTaps32 = (Ipp32f *) rms->pTaps;
        pSpec32 = (IppsFIRSpec_32f *) rms->pSpec;
        if (pTaps32 != NULL){ippsFree(pTaps32);}
        if (pSpec32 != NULL){ippsFree(pSpec32);}
    }
    else
    {
        pTaps64 = (Ipp64f *) rms->pTaps;
        pSpec64 = (IppsFIRSpec_64f *) rms->pSpec;
        if (pTaps64 != NULL){ippsFree(pTaps64);}
        if (pSpec64 != NULL){ippsFree(pSpec64);}
    }
    memset(rms, 0, sizeof(struct xcfftRMSFilter_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Applies RMS filter to signals.
 *
 * @param[in] nsignals    Number of signals to filter.
 * @param[in] lds         Leading dimension of signals.  This must be greater
 *                        than or equal to npts.
 * @param[in] npts        Number of points in signals.
 * @param[in] precision   Precision of signals.
 * @param[in,out] rms     RMS filter structure.
 * @param[in] x           Signals to filter.  This is an array of dimension
 *                        [lds x nsignals] with leading dimension lds whose
 *                        type is float or double as indicated by precision.
 * @param[out] xfilt      Filtered signals.  This is an array of dimension
 *                        [lds x nsignals] with leading dimension lds whose
 *                        type is float or double as indicated by precision.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 * @bug rms is destroyed.
 *
 */
int xcloc_rmsFilter_apply(const int nsignals,
                          const int lds, const int npts, 
                          const enum xclocPrecision_enum precision,
                          struct xcfftRMSFilter_struct *rms,
                          const void *__restrict__ x,
                          void *__restrict__ xfilt)
{
    IppsFIRSpec_64f *pSpec64 = NULL;
    IppsFIRSpec_32f *pSpec32 = NULL;
    Ipp8u *pBuf;
    double *xcFiltered64, *xcSqr64, *xwork64, pMaxAbs64;
    float *xcFiltered32, *xcSqr32, *xwork32, pMaxAbs32;
    size_t indx, is, nbytes;
    int filterLen, winLen2;
    bool lfilter;
    lfilter = rms->lfilter; 
    if (nsignals == 0 || npts == 0){return 0;} // Nothing to do
    if (x == NULL || xfilt == NULL)
    {
        if (x == NULL){fprintf(stderr, "%s: x is NULL\n", __func__);}
        if (xfilt == NULL){fprintf(stderr, "%s: xfilt is NULL\n", __func__);}
        return -1;
    }
    if (lds < npts)
    {
        fprintf(stderr, "%s: lds = %d < npts = %d\n", __func__, lds, npts);
        return -1; 
    }
    // Straight copy
    if (!lfilter)
    {
        if (precision == XCLOC_SINGLE_PRECISION)
        {
            nbytes = (size_t) (lds*nsignals)*sizeof(float); 
        }
        else if (precision == XCLOC_DOUBLE_PRECISION)
        {
            nbytes = (size_t) (lds*nsignals)*sizeof(double);
        }
        else
        {
            fprintf(stderr, "%s: Invalid precision: %d\n",
                     __func__, (int) precision);
            return -1;
        }
        memcpy(xfilt, x, nbytes); 
        return 0;
    }
    if (rms->precision != precision)
    {
        fprintf(stderr, "%s: Precision mismatch\n", __func__);
        return -1;
    }
    if (npts != rms->npts)
    {
        fprintf(stderr, "%s: npts=%d != rms->npts=%d\n",
                 __func__, npts, rms->npts);
        return -1;
    }
    winLen2 = rms->tapsLen/2;
    filterLen = npts + winLen2; // Length of filtered signal
    if (rms->precision == XCLOC_SINGLE_PRECISION)
    {
        nbytes = (size_t) filterLen*sizeof(float);
/*
        #pragma omp parallel default(none) \
         private(indx, is, filterLen, pMaxAbs32, pSpec32, xcFiltered32) \
         private(xcSqr32, xwork32) \
         shared(lds, lfilter, nbytes, npts, nsignals, winLen2, x, xfilt) \
         firstprivate(rms)
*/
        {
        pBuf = (Ipp8u *) rms->pBuf;
        pSpec32 = (IppsFIRSpec_32f *) rms->pSpec;
        xcSqr32 = (float *) aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
        xcFiltered32 = (float *) aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
        xwork32 = (float *) aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
        memset(xcSqr32, 0, nbytes);  // Pre-pad signal 
        // Loop on the transform signals and compute RMS
        //#pragma omp for
        for (is=0; is<nsignals; is++)
        {
            // Index of cross-correlation pair 
            indx = (size_t) (is*lds)*sizeof(float);
            // Compute RMS over sliding window
            if (lfilter)
            {
                // Deal with potential precision problem
                ippsMaxAbs_32f(&x[indx], npts, &pMaxAbs32);
                //ippsMean_32f(&x[indx], npts, &pMaxAbs, ippAlgHintFast);
                if (pMaxAbs32 == 0.0f){continue;} // Dead trace
                pMaxAbs32 = (float) (1.0/(double) pMaxAbs32); // Precision issue
                ippsMulC_32f(&x[indx], pMaxAbs32, xwork32, npts);
                // Square all elements of cross-correlation
                vsSqr(npts, xwork32, xcSqr32);
                // Low pass filter
                ippsFIRSR_32f(xcSqr32, xcFiltered32, filterLen, pSpec32,
                              NULL, NULL, pBuf);
                // Take sqrt understanding xcFiltered is delayed by winLen2
                // samples Note, the importance of winLen2 here.  The trace has
                // been delayed grpdelay(b) = nb/2 samples where nb is odd.  
                // Consequently, the phase distortion is shifted by ignoring
                // the first winLen2 elements.
                vsSqrt(npts, &xcFiltered32[winLen2], xwork32);
                // Undo the scaling
                pMaxAbs32 = (float) (1.0/(double) pMaxAbs32); // Precision issue
                ippsMulC_32f(xwork32, pMaxAbs32, &xfilt[indx], npts);
           }
           // Window is length 1 so square then square root means absolute value
           else
           {
                vsAbs(npts, &x[indx], &xfilt[indx]);
           }
       } // Loop on signals
       free(xcSqr32);
       free(xcFiltered32);
       } // End parallel for
    }
    else
    {
        pSpec64 = (IppsFIRSpec_64f *) rms->pSpec;
        nbytes = (size_t) filterLen*sizeof(float);
        #pragma omp parallel default(none) \
         private(indx, is, filterLen, pMaxAbs64, xcFiltered64) \
         private(xcSqr64, xwork64) \
         shared(lds, lfilter, nbytes, npts, nsignals, winLen2, x, xfilt) \
         firstprivate(pSpec64, pBuf)
        {
        xcSqr64 = (double *) aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
        xcFiltered64 = (double *) aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
        xwork64 = (double *) aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
        memset(xcSqr64, 0, nbytes);  // Pre-pad signal 
        // Loop on the transform signals and compute RMS
        #pragma omp for
        for (is=0; is<nsignals; is++)
        {
            // Index of cross-correlation pair 
            indx = (size_t) (is*lds)*sizeof(double);
            // Compute RMS over sliding window
            if (lfilter)
            {
                // Deal with potential precision problem
                ippsMaxAbs_64f(&x[indx], npts, &pMaxAbs64);
                if (pMaxAbs64 == 0.0){continue;} // Dead trace
                ippsDivC_64f(&x[indx], pMaxAbs64, xwork64, npts);
                // Square all elements of cross-correlation
                vdSqr(npts, xwork64, xcSqr64);
                // Low pass filter
                ippsFIRSR_64f(xcSqr64, xcFiltered64, filterLen, pSpec64,
                              NULL, NULL, pBuf);
                // Take sqrt understanding xcFiltered is delayed by winLen2
                // samples Note, the importance of winLen2 here.  The trace has
                // been delayed grpdelay(b) = nb/2 samples where nb is odd.  
                // Consequently, the phase distortion is shifted by ignoring
                // the first winLen2 elements.
                vdSqrt(npts, &xcFiltered64[winLen2], xwork64);
                // Undo the scaling
                ippsMulC_64f(xwork64, pMaxAbs64, &xfilt[indx], npts);
           }
           // Window is length 1 so square then square root means absolute value
           else
           {
                vsAbs(npts, &x[indx], &xfilt[indx]);
           }
       } // Loop on signals
       free(xcSqr64);
       free(xcFiltered64);
       }
    }
    return 0;
} 
