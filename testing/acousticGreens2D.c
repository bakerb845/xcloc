#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include <float.h>
#include "xcloc_config.h"
#include "acousticGreens2D.h"
#include "iscl_hack.h" //fft/fft.h"
#include <ipps.h>

#ifndef DCMPLX
#define DCMPLX(r,i) ((double) (r) + (double) (i)*I)
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*!
 * @brief Driver routine for computing Greens functions that have been
 *        convolved with a Ricker wavelet.  If multiple sources are
 *        present then they will be stacked into the signal.
 * @param[in] nsrc      Number of sources.
 * @param[in] nrec      Number of receivers.
 * @param[in] nptsSig   Number of points in signal.
 * @param[in] fcent     Center frequency (Hz) of Ricker wavelet.
 * @param[in] dt        Sampling period.
 * @param[in] lnorm     If true then normalize the Ricker wavelet by the
 *                      square root of the energy in the signal.
 * @param[in] lshift    If true then shift the Ricker wavelet to the start
 *                      of the trace.
 * @param[in] vel       Model velocity (m/s).
 * @param[in] rho       Model mensity (kg/m**3).
 * @param[in] Q         Quality factor of model.
 * @param[in] srcScale  Magnitude of each source.  This is an array of
 *                      dimension [nsrc].
 * @param[in] xs        Source locations (m).  This is a [nsrc x 3] matrix
 *                      stroed in row major format.  Each row is a (x,y,z)
 *                      triplet.
 * @param[in] xr        Receiver locations (m).  This is a [nrec x 3] matrix
 *                      stored in row major format.  Each row is a (x,y,z)
 *                      triplet.
 * @param[out] obsOut   Synthetic seismograms.  This is a [nrec x nptsSig]
 *                      matrix stored in row major order.
 * @result 0 indicates success.
 */
int acousticGreens2D_computeGreensFunctions(
    const int nsrc, const int nrec, const int nptsSig,
    const double fcent, const double dt,
    const bool lnorm, const bool lshift,
    const double vel, const double rho,
    const double Q,
    const double srcScale[],
    const double xs[],
    const double xr[],
    double **obsOut)
{
    double *obs, *obsTemp, *stf;
    int ierr, isrc, k;
    // Compute the Green's functions using a Ricker wavelet
    fprintf(stdout, "%s: Computing Ricker wavelet...\n", __func__);
    stf = (double *) calloc((size_t) nptsSig, sizeof(double));
    ierr = acousticGreens2D_computeRickerWavelet(nptsSig, dt, fcent,
                                                 lnorm, lshift, stf);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to compute ricker wavelet\n", __func__);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "%s: Computing synthetics...\n", __func__);
    obsTemp = (double *) calloc((size_t) (nptsSig*nrec), sizeof(double));
    obs  = (double *) calloc((size_t) (nptsSig*nrec), sizeof(double));
    for (isrc=0; isrc<nsrc; isrc++)
    {
        ierr = acousticGreens2D_computeGreensLineSource(nrec, vel, rho, Q,
                                                        nptsSig, dt,
                                                        &xs[3*isrc], xr,
                                                        stf, obsTemp);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing line source greens fns\n",
                    __func__);
            return -1;
        }
        for (k=0; k<nptsSig*nrec; k++)
        {
            obs[k] = obs[k] + srcScale[isrc]*obsTemp[k];
        }
        //ippsAddProductC_64f(obsTemp, srcScale[isrc], obs, nptsSig*nrec);
    }
    free(stf);
    free(obsTemp);
    *obsOut = obs;
    return EXIT_SUCCESS;
}
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
int acousticGreens2D_computeRickerWavelet(const int npts,
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
                ippsCopy_64f(&ricker[i], work, ncopy);
                break;
            }
        }
        memset(ricker, 0, (size_t) npts*sizeof(double));
        ippsCopy_64f(work, ricker, npts);
        free(work);
    }
    // Normalize the area energy in the signal 
    if (lnorm)
    {
        //area = cblas_dnrm2(npts, ricker, 1);
        //cblas_dscal(npts, 1.0/area, ricker, 1);
        ippsNorm_L2_64f(ricker, npts, &area);
        ippsDivC_64f_I(area, ricker, npts);
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
int acousticGreens2D_computeGreensLineSource(
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
    //const bool lverb = false; 
    // Set space
    nomega = npts/2 + 1; // Number of fft samples
    omega = (double *) calloc((size_t) (nomega+8), sizeof(double));
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
        ierr = fft_rfft64f_work_hack(npts, stf, npts, nomega, stfF);
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
    ierr = fft_rfftfreqs64f_work_hack(npts, dt, nomega, omega);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "%s: Error computing rfft freqs\n", __func__);
        return EXIT_FAILURE;
    }
    //cblas_dscal(nomega, 2.0*M_PI, omega, 1); // Convert to angular frequency
    const Ipp64f twopi = 2.0*M_PI;
    ippsMulC_64f_I(twopi, omega, nomega); // Convert to angular frequency
    // Evaluate Greens functions for the different receivers 
    for (irec=0; irec<nrec; irec++)
    {
        ierr = acousticGreens2D_computeLineSourceFD(//lverb,
                                                    vel, rho, Q,
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
        ierr = fft_irfft64z_work_hack(nomega, Gw, npts, &G[irec*npts]);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error inverse transforming Gw\n", __func__);
            return EXIT_FAILURE;
        }
    } // Loop on receivers
    // Release memory
    free(stfF);
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
int acousticGreens2D_computeLineSourceFD(
    //const bool lverb,
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
    #pragma omp simd
    for (i=0; i<nomega; i++)
    {
        xscal = 0.0;
        if (fabs(omega[i]) > 0.0){xscal = sqrt(2.0*vel/(M_PI*omega[i]*r));}
/*
        if (omega[i] == 0.0)
        {
            if (lverb)
            {
                fprintf(stderr, "%s: Warning - setting 0 frequency to 0\n",
                       __func__);
            }
            xscal = 0.0;
        }
*/
        omegar = omega[i]*r;
        carg  = DCMPLX(-omegar/(2.0*vel*Q), -omegar/vel);
        G[i] = coeff*(xscal*cexp(carg));
    }
    // Convolve source time function
    if (stf != NULL)
    {
        const Ipp64fc* pSrc = (const Ipp64fc *) stf;
        Ipp64fc *pSrcDst = (Ipp64fc *) G;
        ippsMul_64fc_I(pSrc, pSrcDst, nomega);
        //for (i=0; i<nomega; i++){G[i] = stf[i]*G[i];}
    }
    return EXIT_SUCCESS;
}
/*!
 * @brief Randomly places receivers in a 3D model for use by the 
 *        acoustic2D Green's functions calculators.
 * @param[in] nrec  Number of receivers.
 * @param[in] x0    Origin of model in x (m).
 * @param[in] x1    End of model in x (m).
 * @param[in] y0    Origin of model in y (m).
 * @param[in] y1    End of model in y (m).
 * @param[in] z0    Origin of model in z (m).
 * @param[in] z1    End of model in z (m).
 * @param[out] xr   The randomly placed receiver positions.  This is 
 *                  a matrix of dimension [nrec x 3] in row major format
 *                  where each row is an (x,y,z) triplet.
 * @result 0 indicates success.
 */
int acousticGreens2D_computeRandomReceiverLocations(
    const int nrec,
    const double x0, const double y0, const double z0,
    const double x1, const double y1, const double z1,
    double *xr)
{
    double dx, dy, dz;
    int irec;
    dx = 0.0;
    dy = 0.0;
    dz = 0.0;
    if (fabs(x1 - x0) > 1.e-8){dx = fabs(x1 - x0);}
    if (fabs(y1 - y0) > 1.e-8){dy = fabs(y1 - y0);}
    if (fabs(z1 - z0) > 1.e-8){dz = fabs(z1 - z0);}
    for (irec=0; irec<nrec; irec++)
    {
        xr[3*irec+0] = x0;
        xr[3*irec+1] = y0;
        xr[3*irec+2] = z0;
        xr[3*irec+0] = ((double) rand()/RAND_MAX)*dx;
        xr[3*irec+1] = ((double) rand()/RAND_MAX)*dy;
        xr[3*irec+2] = ((double) rand()/RAND_MAX)*dz;
    }
    return 0;
}

/*!
 * @brief Computes the reciprocity travel time field in homogeneous medium.
 * @param[in] nx       Number of x grid points.
 * @param[in] ny       Number of y grid points.
 * @param[in] nz       Number of z grid points.
 * @param[in] vel      Constant velocity of model (m/s).
 * @param[in] x0       x origin (m).
 * @param[in] y0       y origin (m).
 * @param[in] z0       z origin (m).
 * @param[in] dx       x grid spacing (m).
 * @param[in] dy       y grid spacing (m).
 * @param[in] dz       z grid spacing (m).
 * @param[in] xr       Receiver position in x (m).
 * @param[in] yr       Receiver position in y (m).
 * @param[in] zr       Receiver position in z (m).
 * @param[out] ttable  Travel times (s) from the receiver at (xr,yr,zr) to
 *                     all poitns in the medium.  This is an array of
 *                     dimension [nz x ny x nx] with leading dimension nx.
 * @result 0 indicates success.
 */
int acousticGreens2D_computeTravelTimeTable(
    const int nx, const int ny, const int nz,
    const double vel,
    const double x0, const double y0, const double z0,
    const double dx, const double dy, const double dz,
    double xr, double yr, double zr,
    double ttable[])
{
    double diff_x, diff_y, diff_z, dist2, slow, x, y, z;
    int igrd, ix, iy, iz;
    slow = 1.0/vel;
    for (iz=0; iz<nz; iz++)
    {
        for (iy=0; iy<ny; iy++)
        {
            #pragma omp simd
            for (ix=0; ix<nx; ix++)
            {
                x = x0 + (double) ix*dx;
                y = y0 + (double) iy*dy;
                z = z0 + (double) iz*dz;
                diff_x = xr - x;
                diff_y = yr - y;
                diff_z = zr - z;
                dist2 = diff_x*diff_x + diff_y*diff_y + diff_z*diff_z;
                igrd = iz*nx*ny + iy*nx + ix;
                ttable[igrd] = sqrt(dist2)*slow;
            }
        }
    }
    return 0;
}
