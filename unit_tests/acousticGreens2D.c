#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include <float.h>
#include "xcloc_config.h"
#include "acoustic_greens.h"
#include "iscl/fft/fft.h"
#ifdef XCLOC_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#ifndef DCMPLX
#define DCMPLX(r,i) ((double) (r) + (double) (i)*I)
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*!
 * @brief Computes a ricker wavelet.
 *
 * @param[in] npts      Number of points in signal. 
 * @param[in] dt        Sampling period (seconds).
 * @param[in] peakFreq  Peak frequency of Ricker wavelet (Hz).
 * @param[in] lnorm     If true then normalize the Ricker wavelet
 *                      by the square root of the energy in the signal.
 * @param[in] lshift    If true then shift the Ricker wavelet to the
 *                      start of the trace.
 * @param[out] ricker   Ricker wavelet.  This is an array of dimension
 *                      [npts].
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int dales_unitTests_computeRickerWavelet(const int npts,
                                         const double dt, 
                                         const double peakFreq,
                                         const bool lnorm,
                                         const bool lshift,
                                         double *__restrict__ ricker)
{
    double *work, area, pi2f2, pi2f2t2, t, t2, xmax;
    int i, ncopy, npts2;
    const double tol = 1.0 - 0.995; //1.e-2; //FLT_EPSILON*1000.0; //1.e-6;
    xmax = 0.0;
    pi2f2 = pow(M_PI*peakFreq, 2); 
    npts2 = npts/2;
    for (i=0; i<npts; i++)
    {   
        t = (double) (-npts2 + i)*dt;
        t2 = t*t;
        pi2f2t2 = pi2f2*t2;
        ricker[i] = (1.0 - 2.0*pi2f2t2)*exp(-pi2f2t2);
        //printf("%e %e\n", t, ricker[i]);
        xmax = fmax(fabs(ricker[i]), xmax);
    }   
    // Shift the wavelet to the start of the trace - easier to work in
    // pct of 1 b/c the wavelet has max value of 1
    if (lshift)
    {   
        work = (double *) calloc((size_t) npts, sizeof(double));
        for (i=0; i<npts; i++)
        {
            if (fabs(ricker[i]) > tol*xmax)
            {
                ncopy = npts - i;
                cblas_dcopy(ncopy, &ricker[i], 1, work, 1); 
                break;
            }
        }
        memset(ricker, 0, (size_t) npts*sizeof(double));
        cblas_dcopy(npts, work, 1, ricker, 1); 
        free(work);
    }
    // Normalize the area energy in the signal 
    if (lnorm)
    {   
        area = cblas_dnrm2(npts, ricker, 1); 
        cblas_dscal(npts, 1.0/area, ricker, 1); 
    }
    return EXIT_SUCCESS;
}
/*!
 * @brief Computes the time domain Green's functions in a half-space where
 *        the source is assumed to extend infinitely in the downward 
 *        direction.  This implements the time domain version of Eqn 11 
 *        of Fichtner et al., 2016.
 *
 * @param[in] nrec     Number of receivers to at which to compute Greens
 *                     functions.
 * @param[in] vel      Velocity of medium (m/s).
 * @param[in] rho      Density in (kg/m**3)
 * @param[in] Q        Quality factor.
 * @param[in] xr       Position (x, y, z) of receiver (meters).
 * @param[in] xs       Position (x, y, z) of source (meters).
 * @param[in] stf      If NULL then a Dirac delta function will be assumed. \n
 *                     Otherwise, this is an array of dimension [npts] 
 *                     that represents the time domain source time function
 *                     to be convolved with the Greens function.
 * @param[out] G       Time domain Greens function.  This is an array of 
 *                     dimension [npts x nrec] with leading dimension nrec.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int dales_unitTests_acousticGreensLineSource(
    const int nrec,
    const double vel,
    const double rho,
    const double Q,
    const int npts,
    const double dt,
    const double *__restrict__ xs,
    const double *__restrict__ xr,
    const double *__restrict__ stf,
    double *__restrict__ G)
{
    double complex *Gw = NULL;
    double complex *stfF = NULL;
    double *omega = NULL;
    int i, ierr, irec, nomega;
    const bool lverb = false; 
    // Set space
    nomega = npts/2 + 1; // Number of fft samples
    omega = (double *) calloc((size_t) nomega, sizeof(double));
    Gw = (double complex *) calloc((size_t) nomega, sizeof(double complex));
    stfF  =(double complex *) calloc((size_t) nomega, sizeof(double complex));
    // Handle the source time function
    if (stf == NULL)
    {
        for (i=0; i<nomega; i++)
        {
            stfF[i] = DCMPLX(1.0, 0.0);
        }
    }
    else
    {
        ierr = fft_rfft64f_work(npts, stf, npts, nomega, stfF);
        if (ierr != EXIT_SUCCESS)
        {
            fprintf(stderr, "%s: Failed to fourier transform stf\n", __func__);
            return EXIT_FAILURE;
        }
    }
    // Fix the stf at the zero and Nyquist frequencies
    stfF[0] = DCMPLX(0.0, 0.0);
    stfF[nomega-1] = DCMPLX(0.0, 0.0);
    // Compute FFT frequencies at which to evaluate Greens function
    ierr = fft_rfftfreqs64f_work(npts, dt, nomega, omega);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "%s: Error computing rfft freqs\n", __func__);
        return EXIT_FAILURE;
    }
    cblas_dscal(nomega, 2.0*M_PI, omega, 1); // Convert to angular frequency
    // Evaluate Greens functions for the different receivers 
    for (irec=0; irec<nrec; irec++)
    {
        ierr = dales_unitTests_acousticGreensLineSourceFD(lverb, vel, rho, Q,
                                                          xs,
                                                          &xr[3*irec],
                                                          nomega, omega,
                                                          stfF, Gw);
        if (ierr != EXIT_SUCCESS)
        {
            fprintf(stderr, "%s: Error computing Greens function", __func__);
            return EXIT_FAILURE;
        }
        // Inverse transform Gw
        ierr = fft_irfft64z_work(nomega, Gw, npts, &G[irec*npts]);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error inverse transforming Gw\n", __func__);
            return EXIT_FAILURE;
        }
    } // Loop on receivers
    // Release memory
    free(Gw);
    free(omega);
    return EXIT_SUCCESS;
}
//============================================================================//
/*!
 * @brief Computes the Green's functions in a half-space where the source
 *        is assumed to extend infinitely in the downward z direction.
 *        This implements Equation 11 of Fichtner et al., 2016. 
 *
 * @param[in] lverb    If true then should the zero-frequency be encountered
 *                     then a warning message will be issued that the 
 *                     function is setting G to 0.
 * @param[in] vel      Velocity of medium (m/s).
 * @param[in] rho      Density in (kg/m**3)
 * @param[in] Q        Quality factor.
 * @param[in] xr       Position (x, y, z) of receiver (meters).
 * @param[in] xs       Position (x, y, z) of source (meters).
 * @param[in] nomega   Number of frequencies 
 * @param[in] omega    Angular frequencies (rad/s) at which to evaluate Green's
 *                     functions.  This is an array of dimension [nomega].
 * @param[in] stf      If not NULL then this is the frequency domain 
 *                     repesentation of the source time function to convolve
 *                     into the Green's functions.  If used then it must be
 *                     an array of dimension [nomega].
 *
 * @param[out] G       Green's functions evaluated at each frequency.  This
 *                     is an array of dimension [nomega].
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int dales_unitTests_acousticGreensLineSourceFD(
    const bool lverb,
    const double vel,
    const double rho,
    const double Q,
    const double xs[3],
    const double xr[3],
    const int nomega,
    const double *__restrict__ omega,
    const double complex *__restrict__ stf,
    double complex *__restrict__ G)
{
    double complex carg, coeff;
    double omegar, r, xscal;
    int i;
    const double complex expipi4 = DCMPLX(sqrt(2.0)/2.0, sqrt(2.0)/2.0);
    // Compute the distance - the distance in z is excluded b/c the line
    // source extends infinitely in z-
    r = pow(xs[0] - xr[0], 2) + pow(xs[1] - xr[1], 2);
    r = sqrt(r);
    coeff = DCMPLX(0.0, -1.0/(4.0*rho*vel*vel))*expipi4;
    for (i=0; i<nomega; i++)
    {
        xscal = sqrt(2.0*vel/(M_PI*omega[i]*r));
        if (omega[i] == 0.0)
        {
            if (lverb)
            {
                fprintf(stderr, "%s: Warning - setting 0 frequency to 0\n",
                       __func__);
            }
            xscal = 0.0;
        }
        omegar = omega[i]*r;
        carg  = DCMPLX(-omegar/(2.0*vel*Q), -omegar/vel);
        G[i] = coeff*(xscal*cexp(carg));
    }
    // Convolve source time function
    if (stf != NULL)
    {
        for (i=0; i<nomega; i++){G[i] = stf[i]*G[i];}
    }
    return EXIT_SUCCESS;
}
