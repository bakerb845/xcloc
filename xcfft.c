#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "xcloc_xcfft.h"
#include <float.h>
#include <mkl.h>
#include <ipps.h>
#include <fftw/fftw3.h>
#include <fftw/fftw3_mkl.h>
#include <mkl_cblas.h>

#define PADDING   DALES_MEM_PADDING 
#define ALIGNMENT DALES_MEM_ALIGNMENT 

#define CHECK_FFT(status) \
if (status && !DftiErrorClass(status, DFTI_NO_ERROR)) \
{ \
  fprintf(stderr, "%s: Error: %s on line %d\n", __func__, \
           DftiErrorMessage(status), __LINE__); \
  return -1; \
};

static int computePadding32f(const int n);
static int computePadding64f(const int n);

int xcloc_xcfft_checkParameters(const int npts,
                                const int nptsPad,
                                const int nsignals)
{
    if (npts < 1 || nsignals < 2 || nptsPad < npts)
    {
        if (npts < 1) 
        {
            fprintf(stderr, "%s: Signal length %d must be positive\n",
                     __func__, npts);
        }
        if (nptsPad < npts)
        {
            fprintf(stderr, "%s: Pad length %d < signal length %d\n", 
                    __func__, nptsPad, npts);
        }
        if (nsignals < 2) 
        {
            fprintf(stderr, "%s: At least 2 signals required\n", __func__);
        }
        return -1;
    }
    return 0;
}
/*!
 * @brief Initializes the Fourier transformer for the nsignals each of length
 *        npts.  For more details on implementation:
 *        https://software.intel.com/en-us/mkl-developer-reference-c-configuring-and-computing-an-fft-in-c/c
 *
 * @param[in] npts      Number of points in input signals.
 * @param[in] nptsPad   As a tuning parameter or necessity it can be helpful
 *                      to pad the signals.  The length of the 
 *                      cross-correlations will be 2*nptsPad - 1 where nptsPad
 *                      is greater than or equal to npts.
 * @param[in] nsignals  Number of input signals.
 * @param[out] xcfft    On exit this contains the requisite space and
 *                      FFT wisdom to forward and inverse transform the
 *                      signals when computing the correlations.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_initialize(const int npts, const int nptsPad,
                           const int nsignals, struct xcfft_struct *xcfft)
{
    int ierr, lxc;
    memset(xcfft, 0, sizeof(struct xcfft_struct));
    if (npts < 1 || nsignals < 2 || nptsPad < npts)
    {
        if (npts < 1)
        {
            fprintf(stderr, "%s: Signal length %d must be positive\n",
                     __func__, npts);
        }
        if (nptsPad < npts)
        {
            fprintf(stderr, "%s: Pad length %d < signal length %d\n", 
                    __func__, nptsPad, npts);
        }
        if (nsignals < 2)
        {
            fprintf(stderr, "%s: At least 2 signals required\n", __func__);
        }
        return -1;
    }
    // Set some variables on structure
    lxc = nptsPad*2 - 1;              // Length of cross-correlation
    xcfft->nsignals = nsignals;       // Number of input signals.
    xcfft->npts = npts;               // Number of samples in input signals.
    xcfft->lxc = lxc;                 // Length of cross-correlation.
    xcfft->ntfPts = xcfft->lxc/2 + 1; // Number of FFT samples.
    xcfft->ntfSignals = ((nsignals*(nsignals-1)))/2; // Triangular number
    // Set the space
    ierr = xcloc_xcfft_setSpace(xcfft); 
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error allocating space\n", __func__);
        return -1;
    }
    // Create the descriptors
    ierr = xcloc_xcfft_makeFFTWDescriptors(xcfft);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error making FFTw descriptors\n", __func__);
        return -1;
    }
/*
    ierr = dales_xcfft_makeDftiDescriptors(xcfft);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error making descriptors\n", __func__);
        return -1;
    }
*/
    // Compute the cross correlation table
    ierr = dales_xcfft_createXCTable(xcfft);
    if (ierr != 0)
    {   
        fprintf(stderr, "%s: Error creating XC table\n", __func__);
        return -1; 
    }   
    return 0;
}
//============================================================================//
/*!
 * @brief Initializes memory on the xcfft structure.
 *
 * @param[in,out] xcfft  On input contains the length of the cross-correlations,
 *                       the number of transform points, the number of input
 *                       signals, and the number of transform signals. \n
 *                       On exit memory has been appropriately allocated
 *                       to all arrays.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */ 
int xcloc_xcfft_setSpace(struct xcfft_struct *xcfft)
{
    size_t nbytes;
    int nsignals, ntfSignals;
    // Check inputs
    if (xcfft->lxc < 1 || xcfft->ntfPts < 1)
    {
        if (xcfft->lxc < 1)
        {
            fprintf(stderr, "%s: Invalid xc length\n", __func__);
        }
        if (xcfft->ntfPts < 1)
        {
            fprintf(stderr, "%s: Invalid number of transform pts\n", __func__);
        }
        return -1;    
    }
    // Pad out the space for leading dimensions
    xcfft->dataOffset = xcfft->lxc    + computePadding32f(xcfft->lxc);
    xcfft->ftOffset   = xcfft->ntfPts + computePadding64f(xcfft->ntfPts);
    // MPI communicators that are too large may not have any signals
    nsignals = MAX(1, xcfft->nsignals);
    ntfSignals = MAX(1, xcfft->ntfSignals);
    // Set space for input signals
    if (xcfft->precision == XCLOC_SINGLE_PRECISION)
    {
        nbytes = (size_t) (nsignals*xcfft->dataOffset)*sizeof(float);
    }
    else
    {
        nbytes = (size_t) (nsignals*xcfft->dataOffset)*sizeof(double);
    }      
    xcfft->x = aligned_alloc(ALIGNMENT, nbytes);
    memset(xcfft->x, 0, nbytes); // Pre-pad it once and for all
    // Set space for Fourier transform of input signals
    if (xcfft->precision == XCLOC_SINGLE_PRECISION)
    {
        nbytes = (size_t) (nsignals*xcfft->ftOffset)*sizeof(float complex);
    }
    else
    {
        nbytes = (size_t) (nsignals*xcfft->ftOffset)*sizeof(double complex);
    }
    xcfft->fts = aligned_alloc(ALIGNMENT, nbytes);
    memset(xcfft->fts, 0, nbytes); 
    // Set space for all the output cross-correlated signals
    if (xcfft->precision == XCLOC_SINGLE_PRECISION)
    {
        nbytes = (size_t) (ntfSignals*xcfft->dataOffset)*sizeof(float);
    }
    else
    {
        nbytes = (size_t) (ntfSignals*xcfft->dataOffset)*sizeof(double);
    }
    xcfft->y = aligned_alloc(ALIGNMENT, nbytes);
    memset(xcfft->y, 0, nbytes);
    // Set space for all the filtered cross-correlated signals
    xcfft->yfilt = aligned_alloc(ALIGNMENT, nbytes);
    memset(xcfft->yfilt, 0, nbytes);

    // Set space for Fourier transform of correlated signals
    if (xcfft->precision == XCLOC_SINGLE_PRECISION)
    {
        nbytes = (size_t) (ntfSignals*xcfft->ftOffset)*sizeof(float complex);
    }
    else
    {
        nbytes = (size_t) (ntfSignals*xcfft->ftOffset)*sizeof(double complex);
    }
    xcfft->xcfts = aligned_alloc(ALIGNMENT, nbytes);
    memset(xcfft->xcfts, 0, nbytes);
    return 0;
}
//============================================================================//
/*!
 * @brief This creates the averaging portion of the RMS filter.
 *
 */ 
int dales_xcfft_createRMSFilter(const int winLen, struct xcfft_struct *xcfft)
{
    IppsFIRSpec_32f *pSpec = NULL;
    Ipp32f *pTaps;
    Ipp8u *pBuf;
    IppStatus status;
    float avg;
    int bufferSize, specSize, tapsLen;
    memset(&xcfft->fir, 0, sizeof(struct xcfftFIR_struct));
    if (winLen < 2)
    {
        fprintf(stdout, "%s: Will not low-pass filter signal\n", __func__);
        return 0;
    }
    tapsLen = winLen;
    // Make filter have an odd length for my sanity - phase shift is tapsLen/2
    if (winLen%2 == 0){tapsLen = tapsLen + 1;}
    if (tapsLen > xcfft->lxc)
    {
        fprintf(stderr, "%s: Error winLen=%d > lxc=%d\n",
                __func__, tapsLen, xcfft->lxc);
        return -1;
    }
    if (winLen < 1)
    {
        fprintf(stdout, "%s: Will not low-pass filter signal\n", __func__);
        return 0;
    }
    // Create the taps
    pTaps = ippsMalloc_32f(tapsLen);
    avg = (float) (1.0/(double) tapsLen); // Normalize FIR summation by 1/N
    ippsSet_32f(avg, pTaps, tapsLen);
    // Get the space
    status = ippsFIRSRGetSize(tapsLen, ipp32f, &specSize, &bufferSize);
    if (status != ippStsNoErr)
    {
        fprintf(stderr, "%s: Error getting state size\n", __func__);
        return -1;
    }
    pSpec = (IppsFIRSpec_32f *) ippsMalloc_8u(specSize);
    pBuf = ippsMalloc_8u(bufferSize); 
    status = ippsFIRSRInit_32f(pTaps, tapsLen, ippAlgAuto, pSpec);
    if (status != ippStsNoErr)
    {
        fprintf(stderr, "%s: Error initializing FIR filter\n", __func__);
        return -1;
    }
    // Set it on xcfft
    xcfft->fir.pSpec = (void *) pSpec;
    xcfft->fir.pBuf = (void *) pBuf;
    xcfft->fir.pTaps = pTaps;
    xcfft->fir.tapsLen = tapsLen;
    xcfft->fir.bufferSize = bufferSize;
    xcfft->fir.specSize = specSize;
    xcfft->fir.lfilter = true;
    return 0;
}
//============================================================================//
/*!
 * @brief Creates the FFTw forward and inverse transform plans.
 *
 * @param[in,out] xcfft  On input holds the number of signals to forward
 *                       and inverse transform, the precision, the length of the
 *                       cross-correlation signals, and the leading dimensions
 *                       of the time domain and frequency domain data. \n
 *                       On output holds the  FFTw plans.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_makeFFTWDescriptors(struct xcfft_struct *xcfft)
{
    fftw_plan  forwardPlan64; 
    fftw_plan  inversePlan64;
    fftwf_plan forwardPlan32;
    fftwf_plan inversePlan32;
    fftw_complex  *ini64, *outf64;
    fftwf_complex *ini32, *outf32;
    double *inf64, *outi64;
    float *inf32, *outi32;
    int nf[1], ni[1], idistf, idisti, howmanyf, howmanyi, odistf, odisti;
    const int *inembed = NULL;
    const int *onembed = NULL;
    const int rank = 1;    // Computing multiple 1D transforms
    const int istride = 1; // Distance between two elements in same input column
    const int ostride = 1; // Distance between two elements in same outpt column
    static bool linitFFTw = false;
    if (!linitFFTw)
    {
#ifdef _OPENMP
        int fftwSuccess;
        int nthreads = omp_get_num_threads();
        fftwSuccess = fftw_init_threads();
        if (fftwSuccess == 0)
        {
            fprintf(stderr, "%s: Error in fftw_init_threads\n", __func__);
        } 
        fftw_plan_with_nthreads(nthreads);
#endif
        linitFFTw = true;
    }
    nf[0]    = xcfft->lxc;        // 1D real transform length
    howmanyf = xcfft->nsignals;   // Number of transforms
    idistf = xcfft->dataOffset;   // Distance between start of k'th input
    odistf = xcfft->ftOffset;     // Distance between start of k'th output 
    ni[0]    = xcfft->lxc;        // Length of time-domain data to inverse
                                  // transform.  This will be the length of
                                  // the cross-correlation.
    howmanyi = xcfft->ntfSignals; // Number of inverse transforms
    idisti = xcfft->ftOffset;     // Distance between start of k'th input
    odisti = xcfft->dataOffset;   // Distance between start of k'th output 
    if (xcfft->precision == XCLOC_SINGLE_PRECISION)
    {
        inf32  = (float *) xcfft->x;
        outf32 = (fftwf_complex *) xcfft->fts;
        forwardPlan32 = fftwf_plan_many_dft_r2c(rank, nf, howmanyf,
                                                inf32, inembed,
                                                istride, idistf,
                                                outf32, onembed,
                                                ostride, odistf,
                                                FFTW_PATIENT);
        ini32  = (fftwf_complex *) xcfft->xcfts;
        outi32 = (float *) xcfft->y;
        inversePlan32 = fftwf_plan_many_dft_c2r(rank, ni, howmanyi,
                                                ini32, inembed,
                                                istride, idisti,
                                                outi32, onembed,
                                                ostride, odisti,
                                                FFTW_PATIENT);
        xcfft->forwardPlan = (void *) forwardPlan32;
        xcfft->inversePlan = (void *) inversePlan32;
    }
    else if (xcfft->precision == XCLOC_DOUBLE_PRECISION)
    {
        inf64  = (double *) xcfft->x;
        outf64 = (fftw_complex *) xcfft->fts;
        forwardPlan64 = fftw_plan_many_dft_r2c(rank, nf, howmanyf,
                                               inf64, inembed,
                                               istride, idistf,
                                               outf64, onembed,
                                               ostride, odistf,
                                               FFTW_PATIENT);
        ini64  = (fftw_complex *) xcfft->xcfts;
        outi64 = (double *) xcfft->y;
        inversePlan64 = fftw_plan_many_dft_c2r(rank, ni, howmanyi,
                                               ini64, inembed,
                                               istride, idisti,
                                               outi64, onembed,
                                               ostride, odisti,
                                               FFTW_PATIENT);
        xcfft->forwardPlan = (void *) forwardPlan64;
        xcfft->inversePlan = (void *) inversePlan64;
    }
    else
    {
        fprintf(stderr, "%s: Invalid precision\n", __func__);
        return -1;
    }
    return 0; 
}
//============================================================================//
/*!
 * @brief Sets the float data for the signal'th signal on the xcfft structure.
 *
 * @param[in] npts       Number of points in signal.
 * @param[in] signal     Signal index.  This is C indexed and in the range
 *                       [0,xcfft->nsignals-1].
 * @param[in] x          Signal to set on x.  This is an array of dimension
 *                       [npts].
 *
 * @param[in,out] xcfft  On input contains space to hold the signal. \n
 *                       On output the signal'th signal has been set on the
 *                       structure.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_setData32f(const int npts, const int signal,
                           const float *__restrict__ x,
                           struct xcfft_struct *xcfft)
{
    int indx;
    Ipp64f *pDst64;
    Ipp32f *pDst32;
    IppStatus status;
    indx = signal*xcfft->dataOffset;
    if (xcfft->precision == XCLOC_SINGLE_PRECISION)
    {
        pDst32 = (Ipp32f *) &xcfft->x[indx];
        status = ippsCopy_32f(x, pDst32, npts);
        if (status != ippStsNoErr)
        {
            fprintf(stderr, "%s: Error copying signal\n", __func__);
            return -1;
        }
        pDst32 = NULL;
    }
    else
    {
        pDst64 = (Ipp64f *) &xcfft->x[indx];
        status = ippsConvert_32f64f(x, pDst64, npts);
        if (status != ippStsNoErr)
        {
            fprintf(stderr, "%s: Error converting signal\n", __func__);
            return -1;
        }
        pDst64 = NULL;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets double data for the signal'th signal on the xcfft structure.
 *
 *
 * @param[in] npts       Number of points in signal.
 * @param[in] signal     Signal index.  This is C indexed and in the range
 *                       [0,xcfft->nsignals-1].
 * @param[in] x          Signal to set on x.  This is an array of dimension
 *                       [npts].
 *
 * @param[in,out] xcfft  On input contains space to hold the signal. \n
 *                       On output the signal'th signal has been set on the
 *                       structure.
 *
 * @result 0 indicates success. 
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_setData64f(const int npts, const int signal,
                           const double *__restrict__ x,
                           struct xcfft_struct *xcfft)
{
    int indx;
    Ipp64f *pDst64;
    Ipp32f *pDst32;
    IppStatus status;
    indx = signal*xcfft->dataOffset;
    if (xcfft->precision == XCLOC_SINGLE_PRECISION)
    {
        pDst32 = (Ipp32f *) &xcfft->x[indx];
        status = ippsConvert_64f32f(x, pDst32, npts);
        if (status != ippStsNoErr)
        {
            fprintf(stderr, "%s: Error converting signal\n", __func__);
            return -1;
        }
    }
    else
    {
        pDst64 = (Ipp64f *) &xcfft->x[indx];
        status = ippsCopy_64f(x, pDst64, npts);
        if (status != ippStsNoErr)
        {
            fprintf(stderr, "%s: Error copying signal\n", __func__);
            return -1;
        }
    }
    return 0;
}
//============================================================================//
/*! 
 * @brief Gets a pointer to the float data for the ixc'th cross-correlation
 *        on the xcfft structure.
 *
 * @param[in] lxc    Length of cross-correlation.
 * @param[in] ixc    Cross-correlation index.  This is in the range
 *                   [0,xcfft->ntfSignals-1].
 * @param[in] xcfft  Structure from which to get data pointer.
 * 
 * @param[out] ierr  0 indicates success.
 *
 * @result Pointer to cross-correlation.  
 * 
 * @author Ben Baker distributed under the MIT license.
 *
 */
const float *dales_xcfft_getXCDataPointer32f(const int lxc, const int ixc,
                                             const struct xcfft_struct xcfft,
                                             int *ierr)
{
    const float *xc = NULL;;
    int indx;
    *ierr = 0;
    if (lxc != xcfft.lxc || ixc < 0 || ixc >= xcfft.ntfSignals)
    {
        *ierr = 1;
        if (lxc != xcfft.lxc)
        {
            fprintf(stderr, "%s: lxc=%d != xcfft.lxc=%d\n",
                    __func__, lxc, xcfft.lxc);
        }
        if (ixc < 0 || ixc >= xcfft.ntfSignals)
        {
            fprintf(stderr, "%s: XC signal %d must be in bounds [%d,%d]\n",
                    __func__, ixc, 0, xcfft.ntfSignals-1);
        }
        return xc;
    }
    indx = ixc*xcfft.dataOffset;
    xc = (const float *) &xcfft.y[indx]; 
    return xc;
}
//============================================================================//
/*!
 * @brief Gets all the time-domain cross-correlations from the xcfft structure.
 *
 * @param[in] nxcs       Number of cross-correlations.  This must match
 *                       xcfft.ntfSignals.
 * @param[in] ldxc       Leading dimension of xcs.  This must be >= lxc.
 * @param[in] lxc        Length of cross-correlation.  This must equal
 *                       xcfft.lxc.
 * @param[in] precision  Precision of xcs.  This can be float or double. 
 * @param[in] xcfft      Structure with the time-domain cross-correlograms.
 *
 * @param[out] xcs       Time domain cross-correlograms.  This is an array
 *                       of dimension [ldxc x nxcs] with leading dimension
 *                       ldxc.  It's precision can be float or double and is
 *                       dictated by the precision.
 * 
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_getAllXCData(const int ntfSignals,
                             const int ldxc, const int lxc,
                             enum xclocPrecision_enum precision,
                             const struct xcfft_struct xcfft,
                             void *__restrict__ xcs)
{
    int ixc;
    size_t offset;
    if (ntfSignals == 0){return 0;} // Nothing to do
    if (ntfSignals != xcfft.ntfSignals || lxc != xcfft.lxc ||
        ldxc < lxc || xcs == NULL)
    {
        if (ntfSignals != xcfft.ntfSignals)
        {
            fprintf(stderr, "%s: ntfSignals=%d != xcfft.ntfSignals=%d\n",
                    __func__, ntfSignals, xcfft.ntfSignals);
        }
        if (lxc != xcfft.lxc)
        {
            fprintf(stderr, "%s: lxc=%d != xcfft.lxc=%d\n",
                    __func__, lxc, xcfft.lxc);
        }
        if (ldxc < lxc)
        {
            fprintf(stderr, "%s: ldxc=%d < lxc=%d\n", __func__, ldxc, lxc);
        }
        if (xcs == NULL){fprintf(stderr, "%s: xcs is NULL\n", __func__);}
        return -1;
    }
    if (precision != XCLOC_SINGLE_PRECISION &&
        precision != XCLOC_DOUBLE_PRECISION)
    {
        fprintf(stderr, "%s: Invalid precision\n", __func__);
        return -1;
    }
    for (ixc=0; ixc<xcfft.ntfSignals; ixc++)
    {
        if (precision == XCLOC_SINGLE_PRECISION)
        {
            offset = (size_t) (ixc*ldxc)*sizeof(float);
            xcloc_xcfft_getXCData(lxc, ixc, XCLOC_SINGLE_PRECISION,
                                  xcfft, &xcs[offset]);
        }
        else
        { 
            offset = (size_t) (ixc*ldxc)*sizeof(double);
            xcloc_xcfft_getXCData(lxc, ixc, XCLOC_DOUBLE_PRECISION,
                                  xcfft, &xcs[offset]);
        }
    }
    return 0;
} 
//============================================================================//
/*!
 * @brief Gets the cross-correlated data for the ixc'th cross-correlation from
 *        the data structure.
 *
 * @param[in] lxc        Length of the cross-correlation. 
 * @param[in] ixc        Cross-correlation index.  This is in the range of
 *                       [0,xcfft.ntfSignals-1].
 * @param[in] precision  Precision of data to get.
 * @param[in] xcfft      Structure containing cross-correlations.
 *
 * @param[out] xc        Array with the ixc'th cross-correlation.  This has 
 *                       dimension [lxc] and corresponds to a float or
 *                       double pointer depending on precision.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 * 
 */
int xcloc_xcfft_getXCData(const int lxc, const int ixc,
                          enum xclocPrecision_enum precision,
                          const struct xcfft_struct xcfft,
                          void *__restrict__ xc)
{
    Ipp64f *src64, *dst64;
    Ipp32f *src32, *dst32;
    int ierr, indx;
    ierr = 0;
    if (lxc != xcfft.lxc || ixc < 0 || ixc >= xcfft.ntfSignals)
    {
        if (lxc != xcfft.lxc)
        {
            fprintf(stderr, "%s: lxc=%d != xcfft.lxc=%d\n",
                    __func__, lxc, xcfft.lxc);
        }
        if (ixc < 0 || ixc >= xcfft.ntfSignals)
        {
            fprintf(stderr, "%s: XC signal %d must be in bounds [%d,%d]\n",
                    __func__, ixc, 0, xcfft.ntfSignals-1);
        }
        return -1;
    }
    if (xcfft.precision == XCLOC_SINGLE_PRECISION)
    {
        indx = ixc*xcfft.dataOffset;
        src32 = (Ipp32f *) &xcfft.y[indx];
        if (precision == XCLOC_SINGLE_PRECISION)
        {
            dst32 = (Ipp32f *) xc;
            ippsCopy_32f(src32, dst32, lxc);
        }
        else
        {
            dst64 = (Ipp64f *) xc;
            ippsConvert_32f64f(src32, dst64, lxc);
        }
    }
    else if (xcfft.precision == XCLOC_DOUBLE_PRECISION)
    {
        indx = ixc*xcfft.dataOffset;
        src64 = (Ipp64f *) &xcfft.y[indx];
        if (precision == XCLOC_SINGLE_PRECISION)
        {
            dst32 = (Ipp32f *) xc;
            ippsConvert_64f32f(src64, dst32, lxc);
        }
        else
        {
            dst64 = (Ipp64f *) xc;
            ippsCopy_64f(src64, xc, lxc);
        }
    }
    else
    {
        fprintf(stderr, "%s: Invalid precision\n", __func__);
        return -1;
    }
    return ierr;
}
//============================================================================//
const float *dales_xcfft_getDataPointer32f(const int npts, const int signal,
                                           const struct xcfft_struct xcfft,
                                           int *ierr)
{
    const float *x = NULL;
    int indx;
    *ierr = 0;
    if (npts != xcfft.npts || signal < 0 || signal >= xcfft.nsignals)
    {
        *ierr = 1;
        if (npts != xcfft.npts)
        {
            fprintf(stderr, "%s: npts=%d != xcfft.npts=%d\n",
                    __func__, npts, xcfft.npts);
        }
        if (signal < 0 || signal >= xcfft.nsignals)
        {
            fprintf(stderr, "%s: Error signal %d must be in bounds [%d,%d]\n",
                    __func__, signal, 0, xcfft.nsignals-1); 
        }
        return x;
    }
    indx = signal*xcfft.dataOffset;
    x = (const float *) &xcfft.x[indx];
    return x; 
}
//============================================================================//
/*!
 * @brief Sets the float data for the signal'th signal on the xcfft structure.
 *
 * @param[in] npts       Number of points in signal.
 * @param[in] signal     Signal index.  This is C indexed and in the range
 *                       [0,xcfft->nsignals-1].
 * @param[in] precision  Precision of the data t oset.
 * @param[in] x          Signal to set on x.  This is an array of dimension
 *                       [npts].  x's type is float or double depending on
 *                       the precision.
 *
 * @param[in,out] xcfft  On input contains space to hold the signal. \n
 *                       On output the signal'th signal has been set on the
 *                       structure.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */

int xcloc_xcfft_setData(const int npts, const int is,
                        const enum xclocPrecision_enum precision,
                        const void *__restrict__ x,
                        struct xcfft_struct *xcfft)
{
    Ipp64f *pSrc64, *pDst64;
    Ipp32f *pSrc32, *pDst32;
    int indx;
    indx = is*xcfft->dataOffset;
    if (precision == XCLOC_SINGLE_PRECISION)
    {
        pSrc32 = (Ipp32f *) x;
        if (xcfft->precision == XCLOC_SINGLE_PRECISION)
        {
            pDst32 = (Ipp32f *) &xcfft->x[indx]; 
            ippsCopy_32f(pSrc32, pDst32, npts);
        }
        else
        {
            pDst64 = (Ipp64f *) &xcfft->x[indx]; 
            ippsConvert_32f64f(pSrc32, pDst64, npts);
        }
    }
    else if (precision == XCLOC_DOUBLE_PRECISION)
    {
        pSrc64 = (Ipp64f *) x;
        if (xcfft->precision == XCLOC_SINGLE_PRECISION)
        {
            pDst32 = (Ipp32f *) &xcfft->x[indx]; 
            ippsConvert_64f32f(pSrc64, pDst32, npts);
        }
        else
        {
            pDst64 = (Ipp64f *) &xcfft->x[indx]; 
            ippsCopy_64f(pSrc64, pDst64, npts);
        }
    }
    else
    {
        fprintf(stderr, "%s: Invalid precision\n", __func__);
        return -1;
    }
    return 0;
}
//============================================================================//
/*! 
 * @brief Convenience function to set all signals on the data structure.
 *
 * @param[in] nsignals   Number of signals to set.  This should equal
 *                       xcfft->nsignals.
 * @param[in] lds        Leading dimension of x.  This should be >= npts.
 * @param[in] npts       Number of points in each signal.  This should
 *                       equal xcfft->npts.
 * @param[in] precision  Defines the precision of x.
 * @param[in] x          Input signals.  This is an array of dimension
 *                       [lds x nsignals] with leading dimension lds.
 *                       x's type is float or double as dictated by the 
 *                       precision variable.
 *
 * @param[in,out] xcfft  On input contains space to hold the input signal,
 *                       the number of signals, and the number of points in 
 *                       each signal. \n
 *                       On output contains the input signals.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_setAllData(const int nsignals, const int lds,
                           const int npts,
                           const enum xclocPrecision_enum precision,
                           const void *__restrict__ x,
                           struct xcfft_struct *xcfft)
{
    size_t offset;
    int ierr, is;
    if (nsignals != xcfft->nsignals || npts != xcfft->npts || lds < npts)
    {
        if (nsignals != xcfft->nsignals)
        {
            fprintf(stderr, "%s: Input %d signals but expecting %d signals\n",
                    __func__, nsignals, xcfft->nsignals);
        }
        if (npts != xcfft->npts)
        {
            fprintf(stderr, "%s: Input %d poitns but expecting %d points\n",
                    __func__, npts, xcfft->npts);
        }
        if (lds < npts)
        {
            fprintf(stderr, "%s: Error lds = %d < npts = %d\n",
                     __func__, lds, npts);
        }
    }
    if (precision != XCLOC_SINGLE_PRECISION &&
        precision != XCLOC_DOUBLE_PRECISION)
    {
        fprintf(stderr, "%s: Invalid precision: %d\n",
                __func__, (int) precision);
        return -1;
    }
    for (is=0; is<nsignals; is++)
    {
        if (precision == XCLOC_SINGLE_PRECISION)
        {
            offset = (size_t) (is*lds)*sizeof(float);
            ierr = xcloc_xcfft_setData(npts, is, XCLOC_SINGLE_PRECISION,
                                       &x[offset], xcfft);
        }
        else
        {
            offset = (size_t) (is*lds)*sizeof(double);
            ierr = xcloc_xcfft_setData(npts, is, XCLOC_DOUBLE_PRECISION,
                                       &x[offset], xcfft);
        }
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error setting data\n", __func__);
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Creates the cross-correlation table.  For this application the 
 *        auto-correlations are not required.
 *
 * @param[in,out] xcfft    On input contains the number of transformed signals
 *                         which is a triangular number. \n
 *                         On output contains the xcfft->ntfSignals transform
 *                         pairs.  The pairs describe an upper triangular 
 *                         matrix with nSignals rows and columns.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int dales_xcfft_createXCTable(struct xcfft_struct *xcfft)
{
    size_t nbytes;
    int ierr; //, i, indx, j;
    const bool lwantDiag = false;
    if (xcfft->nsignals < 2)
    {
        fprintf(stderr, "%s: Error at least 2 signals required\n", __func__);
        return -1;
    }
    if (xcfft->ntfSignals < 1)
    {
        fprintf(stderr, "%s: Error no transform signals\n", __func__);
        return -1;
    }
    if (xcfft->xcPairs != NULL){free(xcfft->xcPairs);}
    nbytes = (size_t) (2*xcfft->ntfSignals)*sizeof(int);
    xcfft->xcPairs = aligned_alloc(ALIGNMENT, nbytes);
    memset(xcfft->xcPairs, 0, nbytes);
    ierr = xcloc_xcfft_computeXCTable(lwantDiag, xcfft->nsignals,
                                      xcfft->ntfSignals, xcfft->xcPairs);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing XC pairs\n", __func__);
        return -1;
    }
    return 0; 
}
//============================================================================//
/*!
 * @brief Computes the cross-correlation pairs for the nsignals input signals.
 *
 * @param[in] lwantDiag   If true then the auto-correlations are desired
 *                        in the correlation pairs.  Otherwise, only the
 *                        super diagonal of the correlation matrix is computed.
 * @param[in] nsignals    Number of signals to cross-correlate.
 * @param[in] npairs      Number of transform pairs.  If lwantDiag is false
 *                        then this should be (nsignals*(nsignals-1))/2.  \n
 *                        Otherwise, this shoudl be (nsignals*(nsignals+1))/2.
 *
 * @param[out] xcPairs    Cross-correlation pairs.  This is an array of
 *                        dimension [2 x npairs] with leading dimension 2.
 *
 * @result 0 indicates success. 
 *
 * @author Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_computeXCTable(const bool lwantDiag,
                               const int nsignals, const int npairs, 
                               int *__restrict__ xcPairs)
{
    int i, j, indx, triNumber;
    // Compute the correlation pairs for the super-diagonal only
    if (!lwantDiag)
    {
        triNumber = (nsignals*(nsignals-1))/2;
        if (npairs != triNumber)
        {
            fprintf(stdout, "%s: Warning npairs=%d != triNumber=%d\n",
                    __func__, npairs, triNumber);
        }
        for (i=0; i<nsignals; i++) 
        {
            for (j=i+1; j<nsignals; j++) 
            {
                indx = nsignals*i - ((i+1)*(i+2))/2 + j; 
                xcPairs[2*indx]   = i;
                xcPairs[2*indx+1] = j; 
            }
        }    
        if (triNumber != indx + 1)
        {
            fprintf(stderr, "%s: Some pairs not computed\n", __func__);
            return -1;
        }
    }
    // Compute the correlation pairs for the upper triangle (includes auto-xc)
    else
    {
        triNumber = (nsignals*(nsignals+1))/2;
        if (npairs != triNumber)
        {
            fprintf(stdout, "%s: Warning npairs=%d != triNumber=%d\n",
                    __func__, npairs, triNumber);
        }
        for (i=0; i<nsignals; i++)
        {
            for (j=i; j<nsignals; j++)
            {
                indx = nsignals*i - (i*(i+1))/2 + j;
                xcPairs[2*indx]   = i;
                xcPairs[2*indx+1] = j;
            }
        }
        if (triNumber != indx + 1)
        {
            fprintf(stderr, "%s: Some pairs not computed\n", __func__);
            return -1;
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the phase correlation.
 *
 * @param[in,out] xcfft  On input contains the FFT descriptors and data. \n
 *                       On output contains the corresponding phase
 *                       correlations.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_computePhaseCorrelation(struct xcfft_struct *xcfft)
{
    int ierr;
    ierr = xcloc_xcfft_apply(true, xcfft);
    return ierr;
}
//============================================================================//
/*!
 * @brief Computes cross-correlation using the FFT.
 *
 * @param[in,out] xcfft  On input contains the FFT descriptors and data. \n
 *                       On output contains the corresponding cross
 *                       correlations.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker distributed under the MIT license.
 *
 */ 
int xcloc_xcfft_computeFFTCrossCorrelation(struct xcfft_struct *xcfft)
{
    int ierr;
    ierr = xcloc_xcfft_apply(false, xcfft);
    return ierr;
}
//============================================================================//
/*!
 * @brief Computes the Fourier transform of the input signals.
 *
 * @param[in,out] xcfft   On input contains the nsignals time domain data 
 *                        to transform. \n
 *                        On output contains the Fourier transforms of
 *                        the input data.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_forwardTransform(struct xcfft_struct *xcfft)
{
    fftw_plan  forwardPlan64 = 0;
    fftwf_plan forwardPlan32 = 0;
    if (xcfft->nsignals == 0){return 0;} // Nothing to do
    if (xcfft->precision == XCLOC_SINGLE_PRECISION)
    {
        forwardPlan32 = (fftwf_plan) xcfft->forwardPlan;
        fftwf_execute(forwardPlan32);
    }
    else
    {
        forwardPlan64 = (fftw_plan) xcfft->forwardPlan;
        fftw_execute(forwardPlan64);
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the inverse Fourier transform of the cross-correlation
 *        pairs then reorders the signals such that the time domain 
 *        signals go from -lxc/2 to +lxc/2 samples.
 *
 * @param[in,out] xcfft  On input contains the frequency domain (phase)
 *                       correlated signals. \n
 *                       On output contains the corresonding time domain
 *                       correlations. 
 *                         
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_inverseTransformXC(struct xcfft_struct *xcfft)
{
    Ipp64f *y64, *ywork64;
    Ipp32f *y32, *ywork32;
    fftw_plan  inversePlan64;
    fftwf_plan inversePlan32;
    double scal64;
    float scal32;
    int indx, is, jndx, ncopy1, ncopy2;
    if (xcfft->ntfSignals == 0){return 0;} // Nothing to do
    if (xcfft->precision == XCLOC_SINGLE_PRECISION)
    {
        inversePlan32 = (fftwf_plan) xcfft->inversePlan;
        fftwf_execute(inversePlan32);
        // Unshuffle signals and apply scaling
        #pragma omp parallel default(none) \
         shared(xcfft) \
         private(indx, is, jndx, ncopy1, ncopy2, scal32, y32, ywork32)
        {
        scal32 = (float) xcfft->lxc;
        ywork32 = ippsMalloc_32f(xcfft->lxc);
        y32 = (Ipp32f *) xcfft->y;
        #pragma omp for
        for (is=0; is<xcfft->ntfSignals; is++)
        {
            ncopy1 = xcfft->lxc/2;
            ncopy2 = xcfft->lxc/2 + 1; 
            indx = is*xcfft->dataOffset + ncopy1 + 1;
            jndx = is*xcfft->dataOffset;
            //ippsCopy_32f(&y32[jndx], ywork32, ncopy2);
            //ippsCopy_32f(&y32[indx], &y32[jndx], ncopy1);
            ippsDivC_32f(&y32[jndx], scal32, ywork32, ncopy2);
            ippsDivC_32f(&y32[indx], scal32, &y32[jndx], ncopy1);
            ippsCopy_32f(ywork32, &y32[indx-1], ncopy2);
        }
        ippsFree(ywork32);
        } // End parallel
    }
    else
    {
        inversePlan64 = (fftw_plan) xcfft->inversePlan;
        fftw_execute(inversePlan64);
        // Unshuffle signals and apply scaling
        #pragma omp parallel default(none) \
         shared(xcfft) \
         private(indx, is, jndx, ncopy1, ncopy2, scal64, y64, ywork64)
        {
        scal64 = (double) xcfft->lxc;
        ywork64 = ippsMalloc_64f(xcfft->lxc);
        y64 = (Ipp64f *) xcfft->y;
        #pragma omp for
        for (is=0; is<xcfft->ntfSignals; is++)
        {
            ncopy1 = xcfft->lxc/2;
            ncopy2 = xcfft->lxc/2 + 1;
            indx = is*xcfft->dataOffset + ncopy1 + 1;
            jndx = is*xcfft->dataOffset;
            //ippsCopy_64f(&y64[jndx], ywork64, ncopy2);
            //ippsCopy_64f(&y64[indx], &y64[jndx], ncopy1);
            ippsDivC_64f(&y64[jndx], scal64, ywork64, ncopy2);
            ippsDivC_64f(&y64[indx], scal64, &y64[jndx], ncopy1);
            ippsCopy_64f(ywork64, &y64[indx-1], ncopy2);
        }
        ippsFree(ywork64);
        } // End parallel
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Given the Fourier transformed data the corresponding frequency 
 *        domain cross-correlations.
 *
 * @param[in] lphaseOnly   If true then the phase-correlations are desired. \n
 *                         Otherwise, the cross-correlations will be computed.
 * @param[in] ntfSignals   Number of Fourier transformed input signals.
 * @param[in] ntfPts       Number of points in Fourier transforms.
 * @param[in] ftOffset     This is the leading dimension of ftsIn and xcftsOut.
 *                         This necessarily must be greater than or equal to 
 *                         ntfPts.
 * @param[in] precision    Precision indicating whter ftsIn and xcftsOut are
 *                         float complex or double complex.
 * @param[in] xcPairs      Maps from the is'th transform signal to the 
 *                         indices of the input signals constituting this 
 *                         correlation pair.  This is an array of dimension
 *                         [2 x ntfSignals] with leading dimension 2.
 * @param[in] ftsIn        Fourier transforms of input signals.  This is 
 *                         an array of dimension [ftOffset x nsignals] 
 *                         with leading dimension ftOffset.  This has type
 *                         float complex or double complex depending on
 *                         the precision.
 * @param[out] xcftsOut    Frequency domain representation of the correlations.
 *                         This is an array of dimension [ftOffset x ntfSignals]
 *                         with leading dimension [ftOffset].  This has type
 *                         float complex or double complex depending on 
 *                         the precision. 
 *
 * @result 0 indicates success. 
 *
 * @author Ben Baker distributed under the MIT license.
 * 
 */
int xcloc_xcfft_computeCorrelationsWithFFTData(
    const bool lphaseOnly,
    const int ntfSignals,
    const int ntfPts,
    const int ftOffset,
    const enum xclocPrecision_enum precision,
    const int *__restrict__ xcPairs,
    const void *__restrict__ ftsIn, //const float complex *__restrict__ ftsIn,
    const void *__restrict__ xcftsOut) //float complex *__restrict__ xcftsOut)
{
    //enum xclocPrecision_enum precision = XCLOC_SINGLE_PRECISION;
    MKL_Complex16 *fts64, *xcfts64, *xc64;
    MKL_Complex8 *fts32, *xcfts32, *xc32;
    double complex *xcorr64;
    float complex *xcorr32;
    Ipp64f *mag64;
    Ipp32f *mag32;
    MKL_INT n;
    int iw, i, is, j;
    const double tol64 = DBL_EPSILON*100.0;
    const float tol32 = FLT_EPSILON*10.0;
    // Correlate in the frequency domain
    if (precision == XCLOC_SINGLE_PRECISION)
    {
        #pragma omp parallel default(none) \
         shared(lphaseOnly, tol32, fts32, ftsIn, xcPairs, xcfts32, xcftsOut) \
         private(i, is, iw, j, n, mag32, xc32, xcorr32)
        {
        n = (MKL_INT) ntfPts;
        mag32 = (Ipp32f *) ippsMalloc_32f(ntfPts);
        fts32 = (MKL_Complex8 *) ftsIn;
        xcfts32 = (MKL_Complex8 *) xcftsOut;
        #pragma omp for
        for (is=0; is<ntfSignals; is++)
        {
            // Get the transform pairs
            i = xcPairs[2*is];
            j = xcPairs[2*is+1];
            // Point to workspace in output
            xc32 = (MKL_Complex8 *) &xcfts32[is*ftOffset];
            // Compute correlation: u1*conj(u2)
            vcMulByConj(n, &fts32[i*ftOffset], &fts32[j*ftOffset], xc32);
            // Normalize by absolute value of frequency
            if (lphaseOnly)
            {
                // Compute magnitude of correlation
                vcAbs(n, xc32, mag32);
                // Avoid division by zero
                #pragma omp simd aligned(mag32: ALIGNMENT)
                for (iw=0; iw<n; iw++)
                {
                    mag32[iw] = MAX(mag32[iw], tol32);
                }
                // Apply normalization to all frequencies
                xcorr32 = (float complex *) xc32;
                #pragma omp simd aligned(xcorr32, mag32: ALIGNMENT)
                for (iw=0; iw<n; iw++)
                {
                    xcorr32[iw] = xcorr32[iw]/mag32[iw];
                }
                xcorr32 = NULL;
            }
            // Dereference pointers
            xc32 = NULL;
        } // Loop on transform signals
        ippsFree(mag32);
        fts32 = NULL;
        xcfts32 = NULL;
        } // End parallel
    }
    else
    {
        #pragma omp parallel default(none) \
         shared(lphaseOnly, tol64, fts64, ftsIn, xcPairs, xcfts64, xcftsOut) \
         private(i, is, iw, j, n, mag64, xc64, xcorr64)
        {
        n = (MKL_INT) ntfPts;
        mag64 = (Ipp64f *) ippsMalloc_64f(ntfPts);
        fts64 = (MKL_Complex16 *) ftsIn;
        xcfts64 = (MKL_Complex16 *) xcftsOut;
        #pragma omp for
        for (is=0; is<ntfSignals; is++)
        {
            // Get the transform pairs
            i = xcPairs[2*is];
            j = xcPairs[2*is+1];
            // Point to workspace in output
            xc64 = (MKL_Complex16 *) &xcfts64[is*ftOffset];
            // Compute correlation: u1*conj(u2)
            vzMulByConj(n, &fts64[i*ftOffset], &fts64[j*ftOffset], xc64);
            // Normalize by absolute value of frequency
            if (lphaseOnly)
            {
                // Compute magnitude of correlation
                vzAbs(n, xc64, mag64);
                // Avoid division by zero
                #pragma omp simd aligned(mag64: ALIGNMENT)
                for (iw=0; iw<n; iw++)
                {
                    mag64[iw] = MAX(mag64[iw], tol64);
                }
                // Apply normalization to all frequencies
                xcorr64 = (double complex *) xc64;
                #pragma omp simd aligned(xcorr64, mag64: ALIGNMENT)
                for (iw=0; iw<n; iw++)
                {
                    xcorr64[iw] = xcorr64[iw]/mag64[iw];
                }
                xcorr64 = NULL;
            }
            // Dereference pointers
            xc64 = NULL;
        } // Loop on transform signals
        ippsFree(mag64);
        fts64 = NULL;
        xcfts64 = NULL;
        } // End parallel
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the cross-correlation or phase correlation of many signals.
 *        The algorithm is to: \n
 *         (1) Fourier transform all nsignal input signals. \n
 *         (2) Compute the (nsignal*(nsignal-1))/2 cross-correlations. \n
 *         (3) Optionally, apply the amplitude normalization to the
 *              cross-correlations. \n
 *         (4) Inverse transform the (nsignal*(nsignal-1))/2 cross-correlations.
 *         (5) Unshuffle the time axis in the cross-correlations so that the
 *             time axis monotonic increasing.
 *
 * @param[in] lphaseOnly  If true then then this performs the phase 
 *                        correlations. \n 
 *                        Otherwise, this computes cross-correlations.
 *
 * @param[in,out] xcfft   On input contains the workspace and input data. \n
 *                        On input contains the phase or cross-correlations.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright Ben Baker distributed under the the MIT license. 
 *
 */
int xcloc_xcfft_apply(const bool lphaseOnly, struct xcfft_struct *xcfft)
{
    int ierr;
    // Compute the Fourier transform
    ierr = xcloc_xcfft_forwardTransform(xcfft); 
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing forward transforms\n", __func__);
        return -1;
    }
    // Compute the correlations
    ierr = xcloc_xcfft_computeCorrelationsWithFFTData(lphaseOnly,
                                                      xcfft->ntfSignals,
                                                      xcfft->ntfPts,
                                                      xcfft->ftOffset,
                                                      xcfft->precision,
                                                      xcfft->xcPairs,
                                                      xcfft->fts,
                                                      xcfft->xcfts);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing correlations\n", __func__);
        return -1;
    }
    // Inverse transform the signals
    ierr = xcloc_xcfft_inverseTransformXC(xcfft);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing inverse transforms\n", __func__);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the sliding RMS averaging filter of the cross-correlations.
 *
 * @author Ben Baker
 *
 */
int dales_xcfft_computeXCRMSWindow(struct xcfft_struct *xcfft)
{
    IppsFIRSpec_32f *pSpec = NULL;
    Ipp8u *pBuf = NULL;
    float *xcSqr, *xcFiltered, *xwork, *y, *yfilt, pMaxAbs;
    int dataOffset, filterLen, indx, ixc, ntfSignals, winLen2;
    size_t nbytes;
    bool lfilter;
    MKL_INT lxc;
    winLen2 = 0;
    lxc = xcfft->lxc;
    // Extract the filter information
    lfilter = xcfft->fir.lfilter;
    if (lfilter)
    {
        pSpec = (IppsFIRSpec_32f *) xcfft->fir.pSpec; 
        pBuf = (Ipp8u *) xcfft->fir.pBuf;
        winLen2 = xcfft->fir.tapsLen/2;
    }
    y = xcfft->y;
    yfilt = xcfft->yfilt;
    ntfSignals = xcfft->ntfSignals;
    dataOffset = xcfft->dataOffset;
/*
    #pragma omp parallel default(none) \
     private(indx, ixc, filterLen, nbytes, pMaxAbs, xcFiltered, xcSqr, xwork) \
     shared(dataOffset, lfilter, lxc, ntfSignals, winLen2, y, yfilt) \
     firstprivate(pSpec, pBuf)
*/
    {
    filterLen = lxc + winLen2; // Length of filtered signal
    nbytes = (size_t) filterLen*sizeof(float); 
    xcSqr = (float *) aligned_alloc(ALIGNMENT, nbytes);
    xcFiltered = NULL;
    if (lfilter)
    {
        xcFiltered = (float *) aligned_alloc(ALIGNMENT, nbytes);
        xwork = (float *) aligned_alloc(ALIGNMENT, nbytes);
    }
    memset(xcSqr, 0, nbytes);  // Pre-pad signal 
    // Loop on the transform signals and compute RMS
    for (ixc=0; ixc<ntfSignals; ixc++)
    {
        // Index of cross-correlation pair 
        indx = ixc*dataOffset;
        // Compute RMS over sliding window
        if (lfilter)
        {
             // Deal with potential precision problem
             ippsMaxAbs_32f(&y[indx], lxc, &pMaxAbs);
             //ippsMean_32f(&y[indx], lxc, &pMaxAbs, ippAlgHintFast);
             if (pMaxAbs == 0.0){continue;} // Dead trace
             pMaxAbs = (float) (1.0/(double) pMaxAbs); // Precision issue
             if (pMaxAbs == 0.0){continue;} // Dead trace
             ippsMulC_32f(&y[indx], pMaxAbs, xwork, lxc);
             // Square all elements of cross-correlation
             vsSqr(lxc, xwork, xcSqr);
             // Low pass filter
             ippsFIRSR_32f(xcSqr, xcFiltered, filterLen, pSpec,
                           NULL, NULL, pBuf);
/*
for (int i=0; i<lxc; i++) 
{
if (xwork[i] > 0.0f){printf("%e %e\n", yfilt[i], sqrtf(xcFiltered[i]));}//xcSqr[i]);
}
getchar();
*/
             // Take sqrt understanding xcFiltered is delayed by winLen2 samples
             // Note, the importance of winLen2 here.  The trace has been
             // delayed grpdelay(b) = nb/2 samples where nb is odd.  
             // Consequently, the phase distortion is shifted by ignoring
             // the first winLen2 elements.
             vsSqrt(lxc, &xcFiltered[winLen2], xwork);
             // Undo the scaling
             pMaxAbs = (float) (1.0/(double) pMaxAbs); // Precision issue
             ippsMulC_32f(xwork, pMaxAbs, &yfilt[indx], lxc); 
        }
        // Window is length 1 so square then square root means absolute value
        else
        {
             vsAbs(lxc, &y[indx], &yfilt[indx]);
        }
    } // Loop on cross-correlation pairs
    if (lfilter)
    {
        free(xcFiltered);
        free(xwork);
    }
    free(xcSqr);
    } // End parallel
    return 0;
}
//============================================================================//
/*!
 * @brief Frees memory on the xcorr data structure and sets all variables
 *        to 0.
 *
 * @param[out] xcfft    On exit all memory has been freed and all variables
 *                      set to 0.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker distributed under the MIT license.
 *
 */
int xcloc_xcfft_finalize(struct xcfft_struct *xcfft)
{
    fftw_plan  plan64;
    fftwf_plan plan32;
    if (xcfft->x != NULL){free(xcfft->x);}
    if (xcfft->y != NULL){free(xcfft->y);}
    if (xcfft->yfilt != NULL){free(xcfft->yfilt);}
    if (xcfft->fts != NULL){free(xcfft->fts);}
    if (xcfft->xcfts != NULL){free(xcfft->xcfts);}
    if (xcfft->xcPairs != NULL){free(xcfft->xcPairs);}
    if (xcfft->precision == XCLOC_SINGLE_PRECISION)
    {
        if (xcfft->forwardPlan != NULL)
        {
            plan32 = (fftwf_plan) xcfft->forwardPlan; 
            fftwf_free(plan32);
        }
        if (xcfft->inversePlan != NULL)
        {
            plan32 = (fftwf_plan) xcfft->inversePlan;
            fftwf_free(plan32);
        }
        fftwf_cleanup(); // Pass through in MKL
    }
    else
    {
        if (xcfft->forwardPlan != NULL)
        {
            plan64 = (fftw_plan) xcfft->forwardPlan;
            fftw_free(plan64);
        }
        if (xcfft->inversePlan != NULL)
        {
            plan64 = (fftw_plan) xcfft->inversePlan;
            fftw_free(plan64);
        }
        fftw_cleanup(); // Pass through in MKL
    }
    memset(xcfft, 0, sizeof(struct xcfft_struct));
    return 0;
}
//============================================================================//
//                         Begin the static functions                         //
//============================================================================//
/*!
 * @brief Utility funciton to pad a 64 bit float (or 32 bit complex) array
 *        to the desired PADDING.  The resulting array length should be 
 *        n + computePadding64f(n).
 */
static int computePadding64f(const int n)
{
    size_t mod, pad;
    int ipad;
    // Set space and make G matrix
    pad = 0; 
    mod = ((size_t) n*sizeof(double))%PADDING;
    if (mod != 0)
    {
        pad = (PADDING - mod)/sizeof(double);
    }
    ipad = (int) pad;
    return ipad;
}
//============================================================================//
/*!
 * @brief Utility function to pad a 32 bit float array to the desired PADDING.
 *        The resulting array length should be n + computePadding32f(n).
 */
static int computePadding32f(const int n)
{
    size_t mod, pad;
    int ipad;
    pad = 0; 
    mod = ((size_t) n*sizeof(float))%PADDING;
    if (mod != 0)
    {    
        pad = (PADDING - mod)/sizeof(float);
    }
    ipad = (int) pad;
    // Intel lecture at 3:47 does this:
    // https://software.intel.com/en-us/videos/episode-59-elimination-of-false-cache-line-sharing
    // const int paddingBytes = 64;
    // const int paddingEements = paddingBytes/sizeof(float);
    // const int ipad = (paddingElements - n%paddingElements);
    // Then pad to n + ipad;
    return ipad;
}
