#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ipps.h>
#include <mkl.h>

int ippsFIRSRGetSize_finter64f(int tapsLen, int *pSpecSize, int *pBufSize)
{
    IppStatus status;
    status = ippsFIRSRGetSize(tapsLen, ipp64f, pSpecSize, pBufSize);
    if (status != ippStsNoErr)
    {
        fprintf(stderr, "%s: Error getting state size\n", __func__);
        return -1;
    }
    return 0;  
}

int ippsFIRSRGetSize_finter32f(int tapsLen, int *pSpecSize, int *pBufSize)
{
    IppStatus status;
    status = ippsFIRSRGetSize(tapsLen, ipp32f, pSpecSize, pBufSize);
    if (status != ippStsNoErr)
    {
        fprintf(stderr, "%s: Error getting state size\n", __func__);
        return -1;
    }
    return 0;  
}

/*!
 * @brief Applies the RMS filter to the signals.
 * @param[in] nsignals   Number of signals to filter.
 * @param[in] lds        Leading dimension of signals.  This must be >= npts.
 * @param[in] npts       Number of points in each signal.
 * @param[in] tapsLen    Number of filter coefficients.  This will be odd.
 * @param[in] taps       The FIR filter coefficients for computing the runinng
 *                       average.  This has dimension [tapsLen].
 * @param[in] x          The matrix of signals to filter.  This is a column
 *                       major matrix of dimension [lds x nsignals].
 * @param[out] xfilt     The RMS filtered signals.  This is a column major
 *                       major matrix of dimension [lds x nsignals]. 
 * @result 0 indicates success.
 * @copyright Ben Baker distributed under the MIT license.
 */
int xcloc_firFilter_rmsFilter64f(const int nsignals,
                                 const int lds,
                                 const int npts,
                                 const int tapsLen,
                                 const double taps[],
                                 const double x[],
                                 double xfilt[])
{
    int filterLen, is, pBufSize, pSpecSize, winLen2;
    IppStatus status;
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
    // Initialize
    winLen2 = tapsLen/2;
    filterLen = npts + winLen2; // Length of filtered signal
    status = ippsFIRSRGetSize(tapsLen, ipp64f, &pSpecSize, &pBufSize);
    if (status != ippStsNoErr)
    {
        fprintf(stderr, "%s: Error getting state size\n", __func__);
        return -1; 
    }
    #pragma omp parallel default(none) \
     shared(nsignals, pBufSize, pSpecSize, x, xfilt, taps), \
     private(is) \
     firstprivate(filterLen, tapsLen, winLen2)
    {
    // Allocate space and initialize the filter
    Ipp64f *pTaps = ippsMalloc_64f(tapsLen);
    ippsCopy_64f(taps, pTaps, tapsLen);
    Ipp8u *pBuf = ippsMalloc_8u(pBufSize);
    IppsFIRSpec_64f *pSpec = (IppsFIRSpec_64f *) ippsMalloc_8u(pSpecSize);
    ippsFIRSRInit_64f(pTaps, tapsLen, ippAlgAuto, pSpec);
    Ipp64f *xcFiltered = ippsMalloc_64f(filterLen);
    Ipp64f *xwork = ippsMalloc_64f(filterLen);
    // Loop on the transformed signals and compute RMS filter
    #pragma omp for
    for (is=0; is<nsignals; is++)
    {
        int indx = is*lds;
        // Deal with potential precision problem
        double pMaxAbs;
        ippsMaxAbs_64f(&x[indx], npts, &pMaxAbs);
        if (pMaxAbs == 0.0){continue;} // Dead trace
        ippsDivC_64f(&x[indx], pMaxAbs, xwork, npts);
        // Square all elements of cross-correlation
        ippsSqr_64f_I(xwork, npts);
        // Low pass filter
        ippsFIRSR_64f(xwork, xcFiltered, filterLen, pSpec,
                      NULL, NULL, pBuf);
        // Take sqrt understanding xcFiltered is delayed by winLen2
        // samples Note, the importance of winLen2 here.  The trace has
        // been delayed grpdelay(b) = nb/2 samples where nb is odd.  
        // Consequently, the phase distortion is shifted by ignoring
        // the first winLen2 elements. 
        ippsSqrt_64f_I(&xwork[winLen2], npts);
        // Undo the scaling
        ippsMulC_64f(xwork, pMaxAbs, &xfilt[indx], npts);
    }
    // Release memory    
    ippsFree(pBuf);
    ippsFree(pSpec);
    ippsFree(pTaps);
    ippsFree(xwork);
    ippsFree(xcFiltered);
    ippsFree(xwork);
    } // End parallel
    return 0;
}
/*!
 * @brief Applies the RMS filter to the signals.
 * @param[in] nsignals   Number of signals to filter.
 * @param[in] lds        Leading dimension of signals.  This must be >= npts.
 * @param[in] npts       Number of points in each signal.
 * @param[in] tapsLen    Number of filter coefficients.  This will be odd.
 * @param[in] taps       The FIR filter coefficients for computing the runinng
 *                       average.  This has dimension [tapsLen].
 * @param[in,out] x      The matrix of signals to filter.  This is a column
 *                       major matrix of dimension [lds x nsignals].
 * @param[in,out] xfilt  On exit these are the RMS filtered signal.  Again, this
 *                       is a column major matrix of dimension [lds x nsignals].
 * @result 0 indicates success.
 * @ingroup xcloc_spxc
 * @copyright Ben Baker distributed under the MIT license.
 */
int xcloc_firFilter_rmsFilter32f(const int nsignals,
                                 const int lds,
                                 const int npts,
                                 const int tapsLen,
                                 const float taps[],
                                 float x[])
{
    int filterLen, is, pBufSize, pSpecSize, winLen2;
    IppStatus status;
    if (nsignals == 0 || npts == 0){return 0;} // Nothing to do
    if (x == NULL)
    {
        if (x == NULL){fprintf(stderr, "%s: x is NULL\n", __func__);}
        return -1;
    }
    if (lds < npts)
    {
        fprintf(stderr, "%s: lds = %d < npts = %d\n", __func__, lds, npts);
        return -1; 
    }
    // Initialize
    winLen2 = tapsLen/2;
    filterLen = npts + winLen2; // Length of filtered signal
    status = ippsFIRSRGetSize(tapsLen, ipp32f, &pSpecSize, &pBufSize);
    if (status != ippStsNoErr)
    {
        fprintf(stderr, "%s: Error getting state size\n", __func__);
        return -1; 
    }
    #pragma omp parallel default(none) \
     shared(nsignals, pBufSize, pSpecSize, x, taps), \
     private(is) \
     firstprivate(filterLen, tapsLen, winLen2)
    {
    // Allocate space and initialize the filter
    Ipp32f *pTaps = ippsMalloc_32f(tapsLen);
    ippsCopy_32f(taps, pTaps, tapsLen);
    Ipp8u *pBuf = ippsMalloc_8u(pBufSize);
    IppsFIRSpec_32f *pSpec = (IppsFIRSpec_32f *) ippsMalloc_8u(pSpecSize);
    ippsFIRSRInit_32f(pTaps, tapsLen, ippAlgAuto, pSpec);
    Ipp32f *xcFiltered = ippsMalloc_32f(filterLen);
    Ipp32f *xwork = ippsMalloc_32f(filterLen);
    // Loop on the transformed signals and compute RMS filter
    #pragma omp for
    for (is=0; is<nsignals; is++)
    {
        int indx = is*lds;
        // Deal with potential precision problem
        float pMaxAbs;
        ippsMaxAbs_32f(&x[indx], npts, &pMaxAbs);
        if (pMaxAbs == 0.0f){continue;} // Dead trace
        pMaxAbs = (float) (1.0/(double) pMaxAbs); // Precision issue
        ippsMulC_32f(&x[indx], pMaxAbs, xwork, npts);
        // Square all elements of cross-correlation
        ippsSqr_32f_I(xwork, npts);
        // Low pass filter
        ippsFIRSR_32f(xwork, xcFiltered, filterLen, pSpec,
                      NULL, NULL, pBuf);
        // Take sqrt understanding xcFiltered is delayed by winLen2
        // samples Note, the importance of winLen2 here.  The trace has
        // been delayed grpdelay(b) = nb/2 samples where nb is odd.  
        // Consequently, the phase distortion is shifted by ignoring
        // the first winLen2 elements. 
        ippsSqrt_32f_I(&xwork[winLen2], npts);
        // Undo the scaling
        pMaxAbs = (float) (1.0/(double) pMaxAbs); // Precision issue
        ippsMulC_32f(xwork, pMaxAbs, &x[indx], npts);
    }
    // Release memory    
    ippsFree(pBuf);
    ippsFree(pSpec);
    ippsFree(pTaps);
    ippsFree(xwork);
    ippsFree(xcFiltered);
    ippsFree(xwork);
    } // End parallel
    return 0;
}
