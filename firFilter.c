#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ipps.h>

#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

/*!
 * @brief Computes the envelope of a set of signals.  It is assumed that the
 *        Hilbert transformer is a Type III FIR filter.
 * @param[in] lds        Leading dimension of signals.  This must be at least npts.
 * @param[in] npts       Number of points in signals.
 * @param[in] nsignals   Number of signals.
 * @param[in] nReCoeffs  Number of real FIR coefficients.
 * @param[in] reCoeffs   Real coefficients.  This is unity at nReCoeffs/2
 *                       and 0 elsewhere.  This has dimension [nReCoeffs].
 * @param[in] nImCoeffs  Number of imaginary FIR coefficients.
 * @param[in] imCoeffs   Imaginary coefficients.  Every other point is 0.
 *                       This has dimension [nImCoeffs].
 * @param[in,out] x      On input this is an [lds x nsignals] array containing
 *                       the signals to filter.
 * @param[in,out] x      On exit the envelope of the signals have been computed.
 * @result 0 indicates success.
 * @author Ben Baker
 * @ingroup xcloc_spxc
 * @copyright Ben Baker distributed under the MIT license.
 * @bug Currently IPP doesn't have a double flavor for the FIR filter so 
 *      the full Hilbert transformer must be specified.
 */
int xcloc_firFilter_envelope64f(const int lds,
                                const int npts,
                                const int nsignals,
                                const int nReCoeffs,
                                const double reCoeffs[],
                                const int nImCoeffs,
                                const double imCoeffs[],
                                double x[])
{
    if (nsignals == 0 || npts == 0){return 0;} // Nothing to do
    if (x == NULL || reCoeffs == NULL || imCoeffs == NULL)
    {
        if (x == NULL){fprintf(stderr, "%s: x is NULL\n", __func__);}
        if (reCoeffs == NULL)
        {
             fprintf(stderr, "%s: reCoeffs is NULL\n", __func__);
        }
        if (imCoeffs == NULL)
        {
             fprintf(stderr, "%s: imCoeffs is NULL\n", __func__);
        }
        return -1;
    }
    if (nReCoeffs < 0)
    {
        fprintf(stderr, "%s: This doesn't makse sense\n", __func__);
        return -1;
    }
    if (lds < npts)
    {
        fprintf(stderr, "%s: lds = %d < npts = %d\n", __func__, lds, npts);
        return -1;
    }
    int winLen, winLen2;
    winLen = nImCoeffs;
    winLen2 = winLen/2;
    // Begin the parallel computation 
    #pragma omp parallel default(none) \
     shared(imCoeffs, lds, nImCoeffs, npts, nsignals, x) \
     firstprivate(winLen, winLen2)
    {
    int filterLen = npts + winLen2;
    Ipp64f *xmean = NULL;
    Ipp64f *yfiltIm = NULL;
    Ipp8u *pBufIm = NULL;
    IppsFIRSpec_64f *pSpecIm = NULL;
    int bufSizeIm, is, specSizeIm;
    // Initialize the FIR filter
    ippsFIRSRGetSize(nImCoeffs, ipp64f, &specSizeIm, &bufSizeIm);
    pSpecIm = (IppsFIRSpec_64f *) ippsMalloc_8u(specSizeIm);
    pBufIm = ippsMalloc_8u(bufSizeIm);
    ippsFIRSRInit_64f(imCoeffs, nImCoeffs, ippAlgAuto, pSpecIm);
    // Set the workspace
    xmean   = ippsMalloc_64f(filterLen);
    yfiltIm = ippsMalloc_64f(filterLen);
    ippsZero_64f(xmean,   filterLen);
    ippsZero_64f(yfiltIm, filterLen);
    // Loop on signals
    #pragma omp for
    for (is=0; is<nsignals; is++)
    {
        // Demean the signal
        double pMean;
        int indx = is*lds;
        ippsMean_64f(&x[indx], npts, &pMean);
        ippsSubC_64f(&x[indx], pMean, xmean, npts);
        // The real part of the FIR filter is just a delay of winLen/2 samples
        //ippsCopy_64f(xmean32, &yfiltRe[winLen2], npts);
        // Imaginary filtering
        ippsFIRSR_64f(xmean, yfiltIm, filterLen, pSpecIm, NULL, NULL, pBufIm);
        Ipp64f* ptrRe = (Ipp64f *) &xmean[0];
        Ipp64f* ptrIm = (Ipp64f *) &yfiltIm[winLen2];
        // Compute the absolute value of the Hilbert transform.  Note, 
        // winLen2 accounts for the phase shift
        ippsMagnitude_64f(ptrRe, ptrIm, xmean, npts);
        // Reincorporate the mean into the signal
        ippsAddC_64f(xmean, pMean, &x[indx], npts);
    }
    if (xmean != NULL){ippsFree(xmean);}
    if (yfiltIm != NULL){ippsFree(yfiltIm);}
    if (pSpecIm != NULL){ippsFree(pSpecIm);}
    if (pBufIm  != NULL){ippsFree(pBufIm);}
    } // End parallel
    return 0;
}

/*!
 * @brief Computes the envelope of a set of signals.  It is assumed that the
 *        Hilbert transformer is a Type III FIR filter.
 * @param[in] lds          Leading dimension of signals.  This must be at least
 *                         npts.
 * @param[in] npts         Number of points in signals.
 * @param[in] nsignals     Number of signals.
 * @param[in] nnzReCoeffs  Number of non-zero real coefficients.  This
 *                         should be 1.
 * @param[in] nzReCoeffs   Indices of sparse filter that are non-zeros.  This
 *                         is an array of dimension [nnzReCoeffs].
 * @param[in] reCoeffs     Real coefficients of sparse FIR filter.  This
 *                         is an array of dimension [nnzReCoeffs].
 * @param[in] nnzImCoeffs  Number of non-zero imaginary coefficients.
 * @param[in] nzImCoeffs   Indices of sparse filter that are non-zeros.  This
 *                         is an array of dimension [nnzImCoeffs].
 * @param[in] imCoeffs     Imaginary coefficients of sparse FIR filter.  This
 *                         is an array of dimension [nnzImCoeffs].
 * @param[in,out] x        On input this is an [lds x nsignals] array containing
 *                         the signals to filter.
 * @param[in,out] x        On exit the envelope of the signals have been
 *                         computed.
 * @result 0 indicates success.
 * @author Ben Baker
 * @ingroup xcloc_spxc
 * @copyright Ben Baker distributed under the MIT license.
 */
int xcloc_firFilter_envelope32f(const int lds,
                                const int npts,
                                const int nsignals,
                                const int nnzReCoeffs,
                                const int nzReCoeffs[],
                                const float reCoeffs[],
                                const int nnzImCoeffs,
                                const int nzImCoeffs[],
                                const float imCoeffs[],
                                float x[])
{
    int ierr = 0;
    int i;
    if (nsignals == 0 || npts == 0){return 0;} // Nothing to do
    if (x == NULL || nzReCoeffs == NULL || nzImCoeffs == NULL ||
        reCoeffs == NULL || imCoeffs == NULL)
    {
        if (x == NULL){fprintf(stderr, "%s: x is NULL\n", __func__);}
        if (nzReCoeffs == NULL)
        {
             fprintf(stderr, "%s: nzReCoeffs is NULL\n", __func__);
        }
        if (nzImCoeffs == NULL)
        {
             fprintf(stderr, "%s: nzImCoeffs is NULL\n", __func__);
        }
        if (reCoeffs == NULL)
        {
             fprintf(stderr, "%s: reCoeffs is NULL\n", __func__);
        }
        if (imCoeffs == NULL)
        {
             fprintf(stderr, "%s: imCoeffs is NULL\n", __func__);
        }
        return -1;
    }
    if (lds < npts)
    {
        fprintf(stderr, "%s: lds = %d < npts = %d\n", __func__, lds, npts);
        return -1; 
    }
    if (nnzReCoeffs != 1)
    {
        fprintf(stderr, "%s: More general case not yet considered\n", __func__);
        return -1;
    }
    if (nnzImCoeffs < 1)
    {
        fprintf(stderr, "%s: No imaginary filter taps!\n", __func__);
        return -1;
    }
    int winLen, winLen2;
    winLen = 0;
    #pragma omp simd reduction(min: winLen)
    for (i=0; i<nnzImCoeffs; i++){winLen = MIN(winLen, nzImCoeffs[i]);}
    //ippsMax_32s(nzImCoeffs, nnzImCoeffs, &winLen);
    winLen2 = winLen/2 + 1; // Filter is symmetric so remove first half
    // Begin the parallel computation
    #pragma omp parallel default(none) \
     shared(imCoeffs, lds, npts, nnzImCoeffs, nzImCoeffs, nsignals, stderr, x) \
     firstprivate(winLen, winLen2) \
     reduction(max:ierr)
    {
    IppStatus status;
    int filterLen = npts + winLen2;
    int bSizeIm, is;
    int orderIm = nzImCoeffs[nnzImCoeffs-1];
    Ipp8u *pBufIm = NULL;
    Ipp32f *xmean = NULL;
    Ipp32f *yfiltIm = NULL;
    IppsFIRSparseState_32f *ppStateIm = NULL;
    status = ippsFIRSparseGetStateSize_32f(nnzImCoeffs, orderIm, &bSizeIm);
    if (status != ippStsNoErr)
    {   
        fprintf(stderr, "%s: Error getting state size: %d\n",
                "xcloc_firFilter_envelope32f",status);
        ierr = 1;
        goto ERROR;
    }
    pBufIm = ippsMalloc_8u(bSizeIm);
    //ppStateIm = (IppsFIRSparseState_32f *) ippsMalloc_8u(bSizeIm);
    status = ippsFIRSparseInit_32f(&ppStateIm,
                                   imCoeffs,
                                   nzImCoeffs,
                                   nnzImCoeffs,
                                   NULL, pBufIm);
    if (status != ippStsNoErr)
    {
        fprintf(stderr, "%s: Error initializing sparse filter: %d\n",
                "xcloc_firFilter_envelope32f", status);
        ierr = 1;
        goto ERROR;
    }
    xmean   = ippsMalloc_32f(filterLen);
    yfiltIm = ippsMalloc_32f(filterLen); 
    ippsZero_32f(xmean,   filterLen);
    ippsZero_32f(yfiltIm, filterLen);
    // Loop on signals
    #pragma omp for
    for (is=0; is<nsignals; is++)
    {
        // Demean the signal
        float pMean;
        int indx = is*lds;
        ippsMean_32f(&x[indx], npts, &pMean, ippAlgHintFast);//Accurate);
        ippsSubC_32f(&x[indx], pMean, xmean, npts);
        // The real part of the FIR filter is just a delay of winLen/2 samples
        //ippsCopy_32f(xmean32, &yfiltRe[winLen2], npts);
        // Imaginary filtering
        ippsFIRSparse_32f(xmean, yfiltIm, filterLen, ppStateIm);
        Ipp32f* ptrRe = (Ipp32f *) &xmean[0];
        Ipp32f* ptrIm = (Ipp32f *) &yfiltIm[winLen2];
        // Compute the absolute value of the Hilbert transform.  Note, 
        // winLen2 accounts for the phase shift
        ippsMagnitude_32f(ptrRe, ptrIm, xmean, npts);
        // Reincorporate the mean into the signal
        ippsAddC_32f(xmean, pMean, &x[indx], npts); 
    }
ERROR:;
    if (pBufIm != NULL)
    {
        ippsFree(pBufIm);
        pBufIm = NULL;
    }
    if (xmean != NULL){ippsFree(xmean);}
    if (yfiltIm != NULL){ippsFree(yfiltIm);}
    //if (ppStateIm != NULL){ippsFree(ppStateIm);}
    //omp_barrier();
    }

    return ierr;
}
/*!
 * @brief Applies the RMS filter to the signals.
 * @param[in] lds        Leading dimension of signals.  This must be >= npts.
 * @param[in] npts       Number of points in each signal.
 * @param[in] nsignals   Number of signals to filter.
 * @param[in] tapsLen    Number of filter coefficients.  This will be odd.
 * @param[in] taps       The FIR filter coefficients for computing the runinng
 *                       average.  This has dimension [tapsLen].
 * @param[in,out] x      The matrix of signals to filter.  This is a column
 *                       major matrix of dimension [lds x nsignals].
 * @param[in,out] x      The RMS filtered signals.  This is a column major
 *                       major matrix of dimension [lds x nsignals]. 
 * @result 0 indicates success.
 * @ingroup xcloc_spxc
 * @copyright Ben Baker distributed under the MIT license.
 */
int xcloc_firFilter_rmsFilter64f(const int lds,
                                 const int npts,
                                 const int nsignals,
                                 const int tapsLen,
                                 const double taps[],
                                 double x[])
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
    status = ippsFIRSRGetSize(tapsLen, ipp64f, &pSpecSize, &pBufSize);
    if (status != ippStsNoErr)
    {
        fprintf(stderr, "%s: Error getting state size\n", __func__);
        return -1; 
    }
    #pragma omp parallel default(none) \
     shared(lds, npts, nsignals, pBufSize, pSpecSize, x, taps), \
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
        ippsMulC_64f(xwork, pMaxAbs, &x[indx], npts);
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
 * @param[in] lds        Leading dimension of signals.  This must be >= npts.
 * @param[in] npts       Number of points in each signal.
 * @param[in] nsignals   Number of signals to filter.
 * @param[in] tapsLen    Number of filter coefficients.  This will be odd.
 * @param[in] taps       The FIR filter coefficients for computing the runinng
 *                       average.  This has dimension [tapsLen].
 * @param[in,out] x      The matrix of signals to filter.  This is a column
 *                       major matrix of dimension [lds x nsignals].
 * @param[in,out] x      On exit these are the RMS filtered signal.  Again, this
 *                       is a column major matrix of dimension [lds x nsignals].
 * @result 0 indicates success.
 * @ingroup xcloc_spxc
 * @copyright Ben Baker distributed under the MIT license.
 */
int xcloc_firFilter_rmsFilter32f(const int lds,
                                 const int npts,
                                 const int nsignals,
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
     shared(lds, nsignals, npts, pBufSize, pSpecSize, x, taps), \
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
