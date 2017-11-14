#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "xcloc_envelope.h"
#include <ipps.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


static void sinc(const int n,
                 const double *__restrict__ x,
                 double *__restrict__ sinc);

int xcloc_envelope_initialize(const int n,
                              const enum xclocPrecision_enum precision,
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
    // Fix the Type III case
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

int xcloc_envelope_finalize(struct xclocEnvelope_struct *envelope)
{
    memset(envelope, 0, sizeof(struct xclocEnvelope_struct));
    return 0;
}
int xcloc_envelope_apply(const int nsignals,
                          const int lds, const int npts, 
                          const enum xclocPrecision_enum precision,
                          struct xclocEnvelope_struct envelope,
                          const void *__restrict__ x,
                          void *__restrict__ xfilt)
{
    IppsFIRSpec_32f *pSpec32R, *pSpec32I;
    Ipp8u *pBuf;
    Ipp32f *xmean32, *xwork32R, *xwork32I, pMean32;
    int filterLen, is, winLen2;
    size_t indx;
    if (!envelope.linit)
    {
        fprintf(stderr, "%s: envelope not initialized\n", __func__);
        return -1;
    }
    if (precision != envelope.precision)
    {
        fprintf(stderr, "%s: Mixed precision not yet handled\n", __func__);
        return -1;
    }
    winLen2 = envelope.tapsLen/2;
    filterLen = npts + winLen2;
    if (envelope.precision == XCLOC_SINGLE_PRECISION)
    {
        pSpec32R = (IppsFIRSpec_32f *) ippsMalloc_8u(envelope.specSize);
        pSpec32I = (IppsFIRSpec_32f *) ippsMalloc_8u(envelope.specSize);
        pBuf = ippsMalloc_8u(envelope.bufferSize);
        ippsFIRSRInit_32f(envelope.hfiltR, envelope.tapsLen,
                          ippAlgAuto, pSpec32R);
        ippsFIRSRInit_32f(envelope.hfiltI, envelope.tapsLen,
                          ippAlgAuto, pSpec32I);
        xwork32I = ippsMalloc_32f(filterLen);
        xwork32R = ippsMalloc_32f(filterLen);
        xmean32  = ippsMalloc_32f(filterLen);
        ippsZero_32f(xmean32, filterLen);
        // Loop on the transforms signals and compute envelope
        for (is=0; is<nsignals; is++)
        {
            // Demean the signal
            ippsMean_32f(&x[indx], npts, &pMean32, ippAlgHintFast);
            ippsSubC_32f(&x[indx], pMean32, xmean32, npts);
            indx = (size_t) (is*lds)*sizeof(float);
            ippsFIRSR_32f(xmean32, xwork32R, filterLen, pSpec32R,
                          NULL, NULL, pBuf);
            ippsFIRSR_32f(xmean32, xwork32I, filterLen, pSpec32I,
                          NULL, NULL, pBuf); 
            // Compute the absolute value of the Hilbert transform.  Note,
            // winLen2 this removes the phase delay.
            ippsMagnitude_32f(&xwork32R[winLen2], &xwork32I[winLen2],
                              xmean32, npts);
            // Re-instate the mean into the signal
            ippsAddC_32f(xmean32, pMean32, &xfilt[indx], npts); 
        }
        ippsFree(pSpec32R);
        ippsFree(pSpec32I);
        ippsFree(pBuf);
        ippsFree(xwork32I);
        ippsFree(xwork32R);
        ippsFree(xmean32);
    }
    else
    {

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
