#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "xcloc_envelope.h"
#ifdef XCLOC_USE_MPI
#include <mpi.h>
#endif
#include <mkl.h>
#include <ipps.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


static void sinc(const int n,
                 const double *__restrict__ x,
                 double *__restrict__ sinc);

/*!
 * @brief Performs the FIR Hilbert transform design.  It is a complex valued
 *        Kaiser window filter design.  For more information see 12.4 of
 *        Oppenheim and Schafer - Discrete Time Signals Processing 3rd ed.
 *
 * @param[in] n          Number of samples in filter.
 * @param[in] precision  Controls whether the filter will be for single or
 *                       double precision signals.
 * @param[in] accuracy   Controls the accuracy when computing the magnitude
 *                       of the Hilbert transform.  Lower accuracy results
 *                       in better performance.
 *
 * @param[out] envelope  Envelope structure with the FIR Hilbert filter design
 *                       and work-space sizes.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_envelope_initialize(const int n,
                              const enum xclocPrecision_enum precision,
                              const enum xclocAccuracy_enum accuracy,
                              struct xclocEnvelope_struct *envelope)
{
    const Ipp64f beta = 8.0;
    const double fc = 1.0;
    double di, dn, fc2, gain, xfact;
    double *hfiltR, *hfiltI, *kaiser, *sinct, *t;
    Ipp64f *hfiltR64, *hfiltI64;
    Ipp32f *hfiltR32, *hfiltI32;
    int i;
    memset(envelope, 0, sizeof(struct xclocEnvelope_struct));
    if (n < 1)
    {
        fprintf(stderr, "%s: Error n must be positive\n", __func__);
        return -1;
    }
    if (precision != XCLOC_SINGLE_PRECISION && 
        precision != XCLOC_DOUBLE_PRECISION)
    {
        fprintf(stderr, "%s: Invalid precision\n", __func__);
        return -1;
    }
    if (accuracy != XCLOC_HIGH_ACCURACY &&
        accuracy != XCLOC_MEDIUM_ACCURACY &&
        accuracy != XCLOC_EXTENDED_ACCURACY)
    {
        fprintf(stderr, "%s: Invalid accuracy\n", __func__);
        return -1;
    }
    // Create the kaiser window
    t      = (double *) ippsMalloc_64f(n);
    hfiltR = (double *) ippsMalloc_64f(n);
    hfiltI = (double *) ippsMalloc_64f(n);
    sinct  = (double *) ippsMalloc_64f(n);
    kaiser = (double *) ippsMalloc_64f(n);
    ippsSet_64f(1.0, kaiser, n);
    if (n > 1)
    {
        xfact = 2.0/(double) (n - 1)*beta; 
        ippsWinKaiser_64f_I(kaiser, n, xfact);
    }
    // Create the ideal filter 
    fc2 = 0.5*fc; 
    dn = (double) n;
    for (i=0; i<n; i++)
    {
        di = (double) i;
        t[i] = fc2*(0.5*(1.0 - dn) + di);
    }
    sinc(n, t, sinct);
    // Create the complex valued FIR coefficients with 12.66 of O & S:
    // hfilt = sinc(t)*exp(i*pi*t) is what matlab uses.  however
    // openheimer has sin*sinc in 12.67.  i'll stick with matlab for
    // consistency.
    gain = 0.0;
    for (i=0; i<n; i++)
    {
        //printf("%e\n", kaiser[i]);
        hfiltR[i] = (kaiser[i]*sinct[i])*cos(M_PI*t[i]);
        hfiltI[i] = (kaiser[i]*sinct[i])*sin(M_PI*t[i]);
        gain = gain + hfiltR[i];
    }
    ippsDivC_64f_I(gain, hfiltR, n);
    ippsDivC_64f_I(gain, hfiltI, n);
    // Fix the Type IV case
    if (n%2 == 1)
    {
        for (i=0; i<n; i++)
        {
            hfiltR[i] = 0.0;
            if (i == n/2){hfiltR[i] = 1.0;}
            if (i%2 == 1){hfiltI[i] = 0.0;}
        }
    }
    envelope->tapsLen = n;
    if (envelope->tapsLen%2 == 1){envelope->ltype4 = true;}
    // Copy filter coefficients
    if (precision == XCLOC_SINGLE_PRECISION)
    {
        hfiltR32 = ippsMalloc_32f(n);
        hfiltI32 = ippsMalloc_32f(n);
        ippsConvert_64f32f(hfiltR, hfiltR32, n);
        ippsConvert_64f32f(hfiltI, hfiltI32, n);
        envelope->hfiltR = (void *) hfiltR32;
        envelope->hfiltI = (void *) hfiltI32;
        ippsFIRSRGetSize(envelope->tapsLen, ipp32f,
                         &envelope->specSize, &envelope->bufferSize);
    }
    else
    {
        hfiltR64 = ippsMalloc_64f(n);
        hfiltI64 = ippsMalloc_64f(n);
        ippsCopy_64f(hfiltR, hfiltR64, n);
        ippsCopy_64f(hfiltI, hfiltI64, n);
        envelope->hfiltR = (void *) hfiltR64;
        envelope->hfiltI = (void *) hfiltI64;
        ippsFIRSRGetSize(envelope->tapsLen, ipp32f,
                         &envelope->specSize, &envelope->bufferSize);
    }
    envelope->linit = true;
/*
for (i=0; i<n; i++)
{
printf("%e %e\n", hfiltR[i], hfiltI[i]);
}
*/ 
    // Release memory and copy to handles
    ippsFree(t);
    ippsFree(kaiser);
    ippsFree(sinct);
    ippsFree(hfiltR);
    ippsFree(hfiltI);
    return 0;
}
//============================================================================//
/*!
 * @brief Finalizes the FIR filter structure.
 *
 * @param[out] envelope   On exit all memory has been released and all 
 *                        variables set to 0.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_envelope_finalize(struct xclocEnvelope_struct *envelope)
{
    if (!envelope->linit){return 0;}
    if (envelope->hfiltR != NULL){ippsFree(envelope->hfiltR);}
    if (envelope->hfiltI != NULL){ippsFree(envelope->hfiltI);}
    memset(envelope, 0, sizeof(struct xclocEnvelope_struct));
    return 0;
}
#ifdef XCLOC_USE_MPI
int xcloc_envelope_applyMPI(const MPI_Comm comm, const int root,
                            const int nsignalsIn, const int ldsIn,
                            const int nptsIn,
                            const MPI_Datatype dataType,
                            struct xclocEnvelope_struct *envelope,
                            const void *__restrict__ x,
                            void *__restrict__ xfilt)
{
    void *xwork1, *xwork2;
    double dPart;
    int *displs, *sendCounts, *sendPtr, ip, ierr, lds, myid, nprocs,
        npts, nsignals, nsignalsLoc, recvCount;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myid);
    if (nprocs == 1)
    {
        if (dataType == MPI_DOUBLE)
        {
            ierr = xcloc_envelope_apply(nsignalsIn, ldsIn, nptsIn,
                                        XCLOC_DOUBLE_PRECISION, envelope,
                                        x, xfilt);
        }
        else
        {
            ierr = xcloc_envelope_apply(nsignalsIn, ldsIn, nptsIn,
                                        XCLOC_SINGLE_PRECISION, envelope,
                                        x, xfilt);
        }
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error applying envelope\n", __func__);
            ierr = 1;
        }
        return ierr;
    }
    if (myid == root)
    {
        nsignals = nsignalsIn;
        lds = ldsIn;
        npts = nptsIn;
    } 
    MPI_Bcast(&nsignals, 1, MPI_INT, root, comm);
    MPI_Bcast(&lds,      1, MPI_INT, root, comm);
    MPI_Bcast(&npts,     1, MPI_INT, root, comm);
    // Split the grid
    dPart = (double) nsignals/(double) nprocs;
    displs     = (int *) calloc((size_t) nprocs, sizeof(int));
    sendCounts = (int *) calloc((size_t) nprocs, sizeof(int));
    sendPtr    = (int *) calloc((size_t) nprocs + 1, sizeof(int));
    for (ip=0; ip<nprocs; ip++)
    {
        sendPtr[ip] = (int) round((double) ip*dPart);
        if (ip == 0){sendPtr[ip] = 0;}
        sendPtr[ip] = MIN(sendPtr[ip], nsignals);
    }
    sendPtr[nprocs] = nsignals;
    nsignalsLoc = sendPtr[myid+1] - sendPtr[myid]; 
    recvCount = nsignalsLoc*lds;
    for (ip=0; ip<nprocs; ip++)
    {
        sendCounts[ip] = (sendPtr[ip+1] - sendPtr[ip])*lds;
        displs[ip] = sendPtr[ip]*lds; 
    }
    if (dataType == MPI_DOUBLE)
    {
        xwork1 = aligned_alloc(XCLOC_MEM_ALIGNMENT,
                               (size_t) (nsignalsLoc*lds)*sizeof(double));
        xwork2 = aligned_alloc(XCLOC_MEM_ALIGNMENT,
                               (size_t) (nsignalsLoc*lds)*sizeof(double));
        MPI_Scatterv(x, sendCounts, displs,
                     MPI_DOUBLE, xwork1, recvCount,
                     MPI_DOUBLE, root, comm);
        ierr = xcloc_envelope_apply(nsignalsLoc,
                                    lds, npts,
                                    XCLOC_DOUBLE_PRECISION,
                                    envelope,
                                    xwork1, xwork2);
        MPI_Allgatherv(xwork2, recvCount, MPI_DOUBLE,
                       xfilt, sendCounts, displs, MPI_DOUBLE, comm);
    }
    else
    {
        xwork1 = aligned_alloc(XCLOC_MEM_ALIGNMENT,
                               (size_t) (nsignalsLoc*lds)*sizeof(float));
        xwork2 = aligned_alloc(XCLOC_MEM_ALIGNMENT,
                               (size_t) (nsignalsLoc*lds)*sizeof(float));
        MPI_Scatterv(x, sendCounts, displs,
                     MPI_FLOAT, xwork1, recvCount,
                     MPI_FLOAT, root, comm);
        ierr = xcloc_envelope_apply(nsignalsLoc,
                                    lds, npts,
                                    XCLOC_SINGLE_PRECISION,
                                    envelope,
                                    xwork1, xwork2);
        MPI_Allgatherv(xwork2, recvCount, MPI_FLOAT,
                       xfilt, sendCounts, displs, MPI_FLOAT, comm);
    }
    free(xwork1);
    free(xwork2);
    return 0;
}
#endif
//============================================================================//
/*!
 * @brief Computes the envelope of a group of signals.
 *
 * @param[in] nsignals       Number of signals.
 * @param[in] lds            Leading dimension of x and xfilt.
 * @param[in] npts           Number of points in signals.
 * @param[in] precision      Precision of the input signals.
 * @param[in,out] envelope   On input contains the initialized envelope
 *                           structure. \n
 *                           On exit some workspace variables may be updated.
 * @param[in] x              Contains the signals of which to compute the
 *                           envelope.  This is an array of dimension 
 *                           [lds x nsignals] whose precision, float or
 *                           double, is dictated by the precision variable.
 * @param[out] xfilt         The envelope of the input signals.  This is an
 *                           array of dimension [lds x nsignals] whose,
 *                           single or double, is dicated by the precision
 *                           variable. 
 *                           
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_envelope_apply(const int nsignals,
                         const int lds, const int npts, 
                         const enum xclocPrecision_enum precision,
                         struct xclocEnvelope_struct *envelope,
                         const void *__restrict__ x,
                         void *__restrict__ xfilt)
{
    IppsFIRSpec_64f *pSpec64R, *pSpec64I;
    IppsFIRSpec_32f *pSpec32R, *pSpec32I;
    Ipp8u *pBuf;
    Ipp64f *hfiltR64, *hfiltI64, *ptr64R, *ptr64I, *xmean64,
           *xwork64R, *xwork64I, pMean64;
    Ipp32f *hfiltR32, *hfiltI32, *ptr32R, *ptr32I, *xmean32,
           *xwork32R, *xwork32I, pMean32;
    int bufferSize, filterLen, is, specSize, tapsLen, winLen2;
    bool ltype4;
    size_t indx;
    enum xclocAccuracy_enum accuracy;
    if (!envelope->linit)
    {
        fprintf(stderr, "%s: envelope not initialized\n", __func__);
        return -1;
    }
    if (precision != envelope->precision)
    {
        fprintf(stderr, "%s: Mixed precision not yet handled\n", __func__);
        return -1;
    }
    ltype4 = envelope->ltype4;
    winLen2 = envelope->tapsLen/2;
    filterLen = npts + winLen2;
    tapsLen = envelope->tapsLen;
    accuracy = envelope->accuracy;
    specSize = envelope->specSize;
    bufferSize = envelope->bufferSize;
    if (envelope->precision == XCLOC_SINGLE_PRECISION)
    {
        hfiltR32 = (Ipp32f *) envelope->hfiltR;
        hfiltI32 = (Ipp32f *) envelope->hfiltI;
        #pragma omp parallel default(none) \
         shared(accuracy, bufferSize, hfiltR32, hfiltI32, filterLen, ltype4, \
                npts, specSize, tapsLen, x, xfilt, winLen2) \
         private(indx, is, pBuf, pMean32, pSpec32R, pSpec32I, ptr32R, ptr32I, \
                 xmean32, xwork32R, xwork32I)
        {
        pSpec32R = (IppsFIRSpec_32f *) ippsMalloc_8u(specSize);
        pSpec32I = (IppsFIRSpec_32f *) ippsMalloc_8u(specSize);
        pBuf = ippsMalloc_8u(bufferSize);
        ippsFIRSRInit_32f(hfiltR32, tapsLen, ippAlgAuto, pSpec32R);
        ippsFIRSRInit_32f(hfiltI32, tapsLen, ippAlgAuto, pSpec32I);
        xwork32I = ippsMalloc_32f(filterLen);
        xwork32R = ippsMalloc_32f(filterLen);
        xmean32  = ippsMalloc_32f(filterLen);
        ippsZero_32f(xwork32R, filterLen);
        ippsZero_32f(xwork32I, filterLen);
        ippsZero_32f(xmean32,  filterLen);
        // Loop on the transforms signals and compute envelope
        #pragma omp for
        for (is=0; is<nsignals; is++)
        {
            // Demean the signal
            indx = (size_t) (is*lds)*sizeof(float);
            ippsMean_32f(&x[indx], npts, &pMean32, ippAlgHintFast);
            ippsSubC_32f(&x[indx], pMean32, xmean32, npts);
            if (ltype4)
            {
                // The real part of the FIR filter is just a delay of 
                // winLen/2 samples  
                //ippsCopy_32f(xmean32, &xwork32R[winLen2], npts);
                // TODO - This is actually a sparse filter
                ippsFIRSR_32f(xmean32, xwork32I, filterLen, pSpec32I,
                              NULL, NULL, pBuf);
                // Handle the phase delay
                ptr32R = (Ipp32f *) &xmean32[0];
                ptr32I = (Ipp32f *) &xwork32I[winLen2];
            }
            else
            {
                // Apply the complex filter
                ippsFIRSR_32f(xmean32, xwork32R, filterLen, pSpec32R,
                              NULL, NULL, pBuf);
                ippsFIRSR_32f(xmean32, xwork32I, filterLen, pSpec32I,
                              NULL, NULL, pBuf);
                ptr32R = (Ipp32f *) &xwork32R[winLen2];
                ptr32I = (Ipp32f *) &xwork32I[winLen2];
            }
            // Compute the absolute value of the Hilbert transform.  Note,
            // winLen2 this removes the phase delay.
            if (accuracy == XCLOC_MEDIUM_ACCURACY)
            {
                vmsHypot(npts, ptr32R, ptr32I, xmean32, VML_LA);
            }
            else if (accuracy == XCLOC_HIGH_ACCURACY)
            {
                vmsHypot(npts, ptr32R, ptr32I, xmean32, VML_HA);
            }
            else
            {
                vmsHypot(npts, ptr32R, ptr32I, xmean32, VML_EP);
            }
            // Reincorporate the mean into the signal
            ippsAddC_32f(xmean32, pMean32, &xfilt[indx], npts);
        } // Loop on signals
        ippsFree(pSpec32R);
        ippsFree(pSpec32I);
        ippsFree(pBuf);
        ippsFree(xwork32I);
        ippsFree(xwork32R);
        ippsFree(xmean32);
        } // End parallel
    }
    else
    {
        hfiltR64 = (Ipp64f *) envelope->hfiltR;
        hfiltI64 = (Ipp64f *) envelope->hfiltI;
        #pragma omp parallel default(none) \
         shared(accuracy, bufferSize, hfiltR64, hfiltI64, filterLen, ltype4, \
                npts, specSize, tapsLen, x, xfilt, winLen2) \
         private(indx, is, pBuf, pMean64, pSpec64R, pSpec64I, ptr64R, ptr64I, \
                 xmean64, xwork64R, xwork64I)
        {
        pSpec64R = (IppsFIRSpec_64f *) ippsMalloc_8u(specSize);
        pSpec64I = (IppsFIRSpec_64f *) ippsMalloc_8u(specSize);
        pBuf = ippsMalloc_8u(bufferSize);
        ippsFIRSRInit_64f(hfiltR64, tapsLen, ippAlgAuto, pSpec64R);
        ippsFIRSRInit_64f(hfiltI64, tapsLen, ippAlgAuto, pSpec64I);
        xwork64I = ippsMalloc_64f(filterLen);
        xwork64R = ippsMalloc_64f(filterLen);
        xmean64  = ippsMalloc_64f(filterLen);
        ippsZero_64f(xwork64R, filterLen);
        ippsZero_64f(xwork64I, filterLen);
        ippsZero_64f(xmean64,  filterLen);
        // Loop on the transforms signals and compute envelope
        #pragma omp for
        for (is=0; is<nsignals; is++)
        {
            // Demean the signal
            indx = (size_t) (is*lds)*sizeof(float);
            ippsMean_64f(&x[indx], npts, &pMean64);
            ippsSubC_64f(&x[indx], pMean64, xmean64, npts);
            if (ltype4)
            {
                // The real part of the FIR filter is just a delay of 
                // winLen/2 samples  
                // TODO - This is actually a sparse filter
                ippsFIRSR_64f(xmean64, xwork64I, filterLen, pSpec64I,
                              NULL, NULL, pBuf);
                // Handle the phase delay
                ptr64R = (Ipp64f *) &xmean64[0];
                ptr64I = (Ipp64f *) &xwork64I[winLen2];
            }
            else
            {
                // Apply the complex filter
                ippsFIRSR_64f(xmean64, xwork64R, filterLen, pSpec64R,
                              NULL, NULL, pBuf);
                ippsFIRSR_64f(xmean64, xwork64I, filterLen, pSpec64I,
                              NULL, NULL, pBuf);
                ptr64R = (Ipp64f *) &xwork64R[winLen2];
                ptr64I = (Ipp64f *) &xwork64I[winLen2];
            }
            // Compute the absolute value of the Hilbert transform.  Note,
            // winLen2 this removes the phase delay.
            if (accuracy == XCLOC_MEDIUM_ACCURACY)
            {
                vmdHypot(npts, ptr64R, ptr64I, xmean64, VML_LA);
            }
            else if (accuracy == XCLOC_HIGH_ACCURACY)
            {
                vmdHypot(npts, ptr64R, ptr64I, xmean64, VML_HA);
            }
            else
            {
                vmdHypot(npts, ptr64R, ptr64I, xmean64, VML_EP);
            }
            // Reincorporate the mean into the signal
            ippsAddC_64f(xmean64, pMean64, &xfilt[indx], npts);
        } // Loop on signals
        ippsFree(pSpec64R);
        ippsFree(pSpec64I);
        ippsFree(pBuf);
        ippsFree(xwork64I);
        ippsFree(xwork64R);
        ippsFree(xmean64);
        } // End parallel
    }
    return 0;
}

static void sinc(const int n,
                 const double *__restrict__ x,
                 double *__restrict__ sinc)
{
    double pix;
    int i;
    #pragma omp simd
    for (i=0; i<n; i++)
    {   
        pix = M_PI*x[i];
        sinc[i] = 1.0; // sinc will evaluate to 1 if x is 0
        if (x[i] != 0.0){sinc[i] = sin(pix)/pix;}
    }
    return;
}
