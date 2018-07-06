#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "xcloc_finter.h"
#include "xcloc_enum.h"
#include "test_suite.h"
#include "iscl_hack.h"
#ifdef XCLOC_PROFILE
#include "advisor-annotate.h"
#endif


void getReferenceSoln(int *npts, double *xc, double *phaseXC);
int test_fdxc(void);
int test_serial_fdxc_hardwired(void);
int test_serial_fdxc_random(const int precision);
#ifndef CHKERR
#define CHKERR(ierr, msg) \
{ \
   if (ierr != EXIT_SUCCESS) \
   { \
       fprintf(stderr, "ERROR calling %s: %s line %d\n", msg, __func__, __LINE__); \
       return EXIT_FAILURE; \
   } \
};
#endif

#define NPTS 1201   /*!< 0.2 s window w/ rate of 6000 Hz */
#define NSIGNALS 45 /*!< make my computer sweat but don't tank the debugger */

/*!
 * @brief Serial tests for the frequency domain cross correlation.
 */
int test_serial_fdxc(void)
{
    int ierr;
    fprintf(stdout, "%s: Performing hardwired tests...\n", __func__);
    ierr = test_serial_fdxc_hardwired();
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed hardwired tests", __func__);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "%s: Performing float random tests...\n", __func__);
    ierr = test_serial_fdxc_random((int) XCLOC_SINGLE_PRECISION);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed float random tests", __func__);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "%s: Performing double random tests...\n", __func__);
    ierr = test_serial_fdxc_random((int) XCLOC_DOUBLE_PRECISION);
    if (ierr != 0)
    {   
        fprintf(stderr, "%s: Failed float random tests", __func__);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int test_serial_fdxc_hardwired(void)
{
    const int npts = 6;
    const int nptsPad = 10;
    const int nsignals = 5;
    const int verbose = 0;
    const int precision = (int) XCLOC_SINGLE_PRECISION;
    const int xcPairRef1[20] = {1, 2, 1, 3, 1, 4, 1, 5,
                                2, 3, 2, 4, 2, 5,
                                3, 4, 3, 5,
                                4, 5};
    const int xcPairRef2[30] = {1, 1, 1, 2, 1, 3, 1, 4, 1, 5,
                                2, 2, 2, 3, 2, 4, 2, 5,
                                3, 3, 3, 4, 3, 5,
                                4, 4, 4, 5,
                                5, 5};
    float *xcs, res;
    // Get reference solutions
    double xcRef[190], phaseRef[190]; 
    int npref;
    getReferenceSoln(&npref, xcRef, phaseRef);
    // Create some bogus time series each with signal length npts=6 and
    // padded to nptsPad=10.
    const float xall[50] = {1, 2, 3, 4, 5, -1,  0, 0,  0,  0,  // signal 1
                            3, 2, 1, 3, 2,  1,  0, 0,  0,  0,  // signal 2
                           -1, 2,-1,-2, 1, -1,  0, 0,  0,  0,  // signal 3
                            4,-2, 3,-1, 5, -1,  0, 0,  0,  0,  // signal 4
                            0,-3,-4, 5,-2,  3,  0, 0,  0,  0}; // signal 5
    int *xcPairs, i, ierr, ixc, nwork, nxcs;
    bool ldoAutoCorrs;
    int accuracy = XCLOC_HIGH_ACCURACY;
    xcPairs = NULL;
    // Compute a cross correlation table with autocorrelations
    ldoAutoCorrs = true;
    nwork =-1; // Workspace query
    xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                      XCLOC_FORTRAN_NUMBERING, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable workspace query");
    nwork = 2*nxcs;
    xcPairs = calloc((size_t) nwork, sizeof(int));
    xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                      XCLOC_FORTRAN_NUMBERING, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable");
    if (nxcs != 15)
    {
        fprintf(stderr, "%s: all nxcs is wrong\n", __func__);
        return EXIT_FAILURE;
    }
    for (i=0; i<2*nxcs; i++)
    {
        if (xcPairRef2[i] != xcPairs[i])
        {
            fprintf(stderr, "%s: error computing pairs diag\n", __func__);
            return EXIT_FAILURE;
        }
    }
    free(xcPairs); xcPairs = NULL;
    // Compute the cross-correlation table
    nwork =-1; // Workspace query
    ldoAutoCorrs = false;
    xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                      XCLOC_FORTRAN_NUMBERING, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable workspace query");
    nwork = 2*nxcs;
    xcPairs = calloc((size_t) nwork, sizeof(int));
    xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                      XCLOC_FORTRAN_NUMBERING, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable");
    if (nxcs != 10)
    {
        fprintf(stderr, "%s: nxcs is wrong\n", __func__);
        return EXIT_FAILURE;
    }
    for (i=0; i<2*nxcs; i++)
    {
        if (xcPairRef1[i] != xcPairs[i])
        {
            fprintf(stderr, "%s: error computing pairs no diag\n", __func__);
            return EXIT_FAILURE;
        }
    }


    // Initialize
    xcloc_fdxc_initialize(npts, nptsPad, 
                          nxcs, xcPairs,
                          verbose, precision, accuracy, &ierr); 
    CHKERR(ierr, "initialize");
    // Compute a cross correlation table with autocorrelations
/*
    nwork =-1; // Workspace query
    xcloc_fdxc_computeDefaultXCTableF(true, nwork, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable workspace query");
    nwork = 2*nxcs;
    xcPairs = calloc((size_t) nwork, sizeof(int));
    ldoAutoCorrs = true;
    xcloc_fdxc_computeDefaultXCTableF(ldoAutoCorrs, nwork, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable");
    if (nxcs != 15)
    {
        fprintf(stderr, "%s: all nxcs is wrong\n", __func__); 
        return EXIT_FAILURE; 
    }   
    for (i=0; i<2*nxcs; i++)
    {
        if (xcPairRef2[i] != xcPairs[i])
        {
            fprintf(stderr, "%s: error computing pairs diag\n", __func__);
            return EXIT_FAILURE;
        }   
    }
    free(xcPairs); xcPairs = NULL;
    // Compute the cross-correlation table
    nwork =-1; // Workspace query
    ldoAutoCorrs = false;
    xcloc_fdxc_computeDefaultXCTableF(ldoAutoCorrs, nwork, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable workspace query");
    nwork = 2*nxcs;
    xcPairs = calloc((size_t) nwork, sizeof(int));
    xcloc_fdxc_computeDefaultXCTableF(false, nwork, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable");
    if (nxcs != 10)
    {
        fprintf(stderr, "%s: nxcs is wrong\n", __func__); 
        return EXIT_FAILURE; 
    }
    for (i=0; i<2*nxcs; i++)
    {
        if (xcPairRef1[i] != xcPairs[i])
        {
            fprintf(stderr, "%s: error computing pairs no diag\n", __func__);
            return EXIT_FAILURE;
        }
    }
    // Set the table
    xcloc_fdxc_setXCTableF(nxcs, xcPairs, &ierr);
    CHKERR(ierr, "xcloc_fdxc_setXCTableF");
*/
    free(xcPairs); xcPairs = NULL;
    // Set the signals
    int ldx = nptsPad;
    xcloc_fdxc_setSignals32f(ldx, npts, nsignals, xall, &ierr);
    CHKERR(ierr, "xcloc_fdxc_setSignals32f");
    // Compute the phase correlograms
    xcloc_fdxc_computePhaseCorrelograms(&ierr);
    CHKERR(ierr, "xcloc_fdxc_computePhaseCorrelograms");
    // Get the result
    int nptsInXC = 2*nptsPad - 1;
    xcs = calloc((size_t) (nptsInXC*nxcs), sizeof(float));
    xcloc_fdxc_getCorrelograms32f(nptsInXC, nxcs, xcs, &ierr);
    CHKERR(ierr, "xcloc_fdxc_getCorrelograms32f");
    for (ixc=0; ixc<nxcs; ixc++)
    {
        for (i=0; i<nptsInXC; i++) //i<nptsInXC*nxcs; i++)
        {
            res = fabsf(xcs[nptsInXC*ixc+i] - (float) phaseRef[nptsInXC*ixc+i]);
            if (res > 90.0*FLT_EPSILON)
            {
                fprintf(stderr, "Failed in xc phase test: %d %d %f %e %e\n",
                        ixc, i, xcs[nptsInXC*ixc+i], phaseRef[nptsInXC*ixc+i], res);
                return EXIT_FAILURE;
            }
        }
    }
    free(xcs);
    // Repeat for cross-correlograms
    xcloc_fdxc_computeCrossCorrelograms(&ierr);
    xcs = calloc((size_t) (nptsInXC*nxcs), sizeof(float));
    xcloc_fdxc_getCorrelograms32f(nptsInXC, nxcs, xcs, &ierr);
    CHKERR(ierr, "xcloc_fdxc_getCorrelograms32f");
    for (ixc=0; ixc<nxcs; ixc++)
    {   
        for (i=0; i<nptsInXC; i++) //i<nptsInXC*nxcs; i++)
        {
            res = fabsf(xcs[nptsInXC*ixc+i] - (float) xcRef[nptsInXC*ixc+i]);
            if (res > 90.0*FLT_EPSILON)
            {
                fprintf(stderr, "Failed in xc test: %d %d %f %e %e\n",
                        ixc, i, xcs[nptsInXC*ixc+i], xcRef[nptsInXC*ixc+i], res);
                return EXIT_FAILURE;
            }
        }
    }
    free(xcs);
    // Finalize
    xcloc_fdxc_finalize(); 
    return EXIT_SUCCESS;
}
//============================================================================//
int test_serial_fdxc_random(const int precision)
{
    int npts = NPTS; // 0.2 s window w/ rate of 6000 Hz
    int nptsPad = npts;
    int lxc = 2*npts - 1;
    int nsignals = NSIGNALS; // make it work but don't brick my computer in debug mode
    const float *xc;
    float *phaseXcsAll, *xcsAll;
    double diffTime, diffTimeRef, res;
    double *xrand;
    int i, ierr, indx, ixc, nxcs;
    clock_t start;
    int accuracy = XCLOC_HIGH_ACCURACY;
    const int verbose = 0;
    //const int precision = (int) XCLOC_DOUBLE_PRECISION; //(int) XCLOC_SINGLE_PRECISION; //0;
    // Compute the cross-correlation table
    int nwork =-1; // Workspace query
    bool ldoAutoCorrs = false;
    int *xcPairs = NULL;
    xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                      XCLOC_FORTRAN_NUMBERING, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable workspace query");
    nwork = 2*nxcs;
    xcPairs = calloc((size_t) nwork, sizeof(int));
    xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                      XCLOC_FORTRAN_NUMBERING, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable");
    // Compute some random signals
    fprintf(stdout, "%s: Generating random numbers...\n", __func__);
    xrand = xcfft_createRandomSignals(NULL, nsignals, npts, &ierr);
    // Initialize
    xcloc_fdxc_initialize(npts, nptsPad,
                          nxcs, xcPairs,
                          verbose, precision, accuracy, &ierr); 
    CHKERR(ierr, "initialize");
/*
    // Compute the cross-correlation table
    int nwork =-1; // Workspace query
    bool ldoAutoCorrs = false;
    int *xcPairs = NULL;
    xcloc_fdxc_computeDefaultXCTableF(ldoAutoCorrs, nwork, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable workspace query");
    nwork = 2*nxcs;
    xcPairs = calloc((size_t) nwork, sizeof(int));
    xcloc_fdxc_computeDefaultXCTableF(false, nwork, &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable");
    // Set the table
    xcloc_fdxc_setXCTableF(nxcs, xcPairs, &ierr);
    CHKERR(ierr, "xcloc_fdxc_setXCTableF");
*/
    free(xcPairs); xcPairs = NULL;
    // Set the signals
    start = clock();
#ifdef XCLOC_PROFILE
//ANNOTATE_SITE_BEGIN("fdxc: random serial");
#endif
    int ldx = nptsPad;
    xcloc_fdxc_setSignals64f(ldx, npts, nsignals, xrand, &ierr);
    CHKERR(ierr, "xcloc_fdxc_setSignals32f");
    // Compute the phase correlograms
    xcloc_fdxc_computePhaseCorrelograms(&ierr);
    CHKERR(ierr, "xcloc_fdxc_computePhaseCorrelograms");
    // Get the result
    int nptsInXC = 2*nptsPad - 1;
    phaseXcsAll = calloc((size_t) (nptsInXC*nxcs), sizeof(float));
    xcloc_fdxc_getCorrelograms32f(nptsInXC, nxcs, phaseXcsAll, &ierr);
    // Compute the cross-correlograms
    xcsAll = calloc((size_t) (nptsInXC*nxcs), sizeof(float));
    // Compute the cross correlograms
    xcloc_fdxc_computeCrossCorrelograms(&ierr);
    CHKERR(ierr, "xcloc_fdxc_computeCrossCorrelograms");
    xcloc_fdxc_getCorrelograms32f(nptsInXC, nxcs, xcsAll, &ierr); 
#ifdef XCLOC_PROFILE
//ANNOTATE_SITE_END("fdxc: random serial");
#endif
    CHKERR(ierr, "xcloc_fdxc_getCorrelograms32f");
    diffTime = (double) (clock() - start)/CLOCKS_PER_SEC;
    // Do it the straight forward wa 
    start = clock();
    double *phaseXcs = (double *) calloc((size_t) (nxcs*lxc)+32, sizeof(double));
    double *xcs = (double *) calloc((size_t) (nxcs*lxc)+32, sizeof(double));
    bool ldoPhase = true;
    ierr = xcfft_computeXCsWithISCL(ldoPhase, nsignals, nxcs,
                                    npts, lxc, xrand, phaseXcs, &diffTimeRef);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "%s: Failed to compute reference signals\n", __func__);
        return EXIT_FAILURE;
    }
    ldoPhase = false;
    ierr = xcfft_computeXCsWithISCL(ldoPhase, nsignals, nxcs,
                                    npts, lxc, xrand, xcs, &diffTimeRef);
    if (ierr != EXIT_SUCCESS)
    {   
        fprintf(stderr, "%s: Failed to compute reference signals\n", __func__);
        return EXIT_FAILURE;
    }
    diffTimeRef = (double) (clock() - start)/CLOCKS_PER_SEC;
    // check phase xcs
    double resMax = 0;
    for (ixc=0; ixc<nxcs; ixc++)
    {
        xc = &phaseXcsAll[ixc*nptsInXC];
        for (i=0; i<lxc; i++)
        {   
            indx = ixc*lxc + i; 
            res = fabs((double) xc[i] - phaseXcs[indx]);
            if (res > 10.0*FLT_EPSILON)
            {   
                fprintf(stderr, "Failed on phase test: %d %d %f %e %e\n",
                        ixc, i, xc[i], phaseXcs[indx], res);
                return EXIT_FAILURE;
            }
            resMax = fmax(resMax, res);
        }
    }
    fprintf(stdout, "%s: phaseXC L1 error=%e\n", __func__, resMax);
    // check xc's
    resMax = 0;
    for (ixc=0; ixc<nxcs; ixc++)
    {
        xc = &xcsAll[ixc*nptsInXC];
        for (i=0; i<lxc; i++)
        {
            indx = ixc*lxc + i;  
            res = fabs((double) xc[i] - xcs[indx]);
            if (res > 100.0*FLT_EPSILON)
            {
                fprintf(stderr, "Failed on xcs test: %d %d %f %e %e\n",
                        ixc, i, xc[i], xcs[indx], res);
                return EXIT_FAILURE;
            }
            resMax = fmax(resMax, res);
        }
    }   
    fprintf(stdout, "%s: XC L1 error=%e\n", __func__, resMax);
    fprintf(stdout, "%s: Processing time %e (s)\n", __func__, diffTime);
    fprintf(stdout, "%s: Reference processing time %e (s)\n", __func__, diffTimeRef);
    fprintf(stdout, "%s: Performance improvement %lf pct\n",
            __func__, (diffTimeRef - diffTime)/diffTimeRef*100.0);
    // Finalize
    xcloc_fdxc_finalize();
    free(xrand);
    free(phaseXcsAll);
    free(xcsAll);
    free(xcs);
    free(phaseXcs);
    return EXIT_SUCCESS;
}
//============================================================================//
/*!
 * @brief Creates nsignals random signals each of npts each.
 *
 * @param[in] seed      Optional random seed.  If NULL then SEED will be used.
 * @param[in] nsignals  Number of signals.
 * @param[in] npts      Number of points in each signal.
 *
 * @result An array of random signals whose dimension is [npts x nsignals]
 *         with leading dimension npts.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
double *xcfft_createRandomSignals(int *seed, const int nsignals, const int npts,
                                  int *ierr)
{
    double *xrand = NULL;
    const int pad = 16; // Pad a little bit so cblas doesn't complain
    int i, is; 
    *ierr = 0;
    if (nsignals < 1 || npts < 1)
    {
        fprintf(stderr, "%s: Invalid number of signals or points\n", __func__);
        *ierr = 1;
        return xrand;
    }
    if (seed != NULL){srand(*seed);}
    // Set space
    xrand = (double *) calloc((size_t) (nsignals*npts+pad), sizeof(double));
    // Compute random numbers from [-1,1]
    for (is=0; is<nsignals; is++)
    {
        // Generate a random signal 
        for (i=0; i<npts; i++)
        {
            xrand[is*npts+i] = ((double) rand()/RAND_MAX - 0.5)*2.0;
        }
    }
    return xrand;
}
//============================================================================//
/*!
 * @brief Computes the phase correlations slowly but correctly.
 *
 * @param[in] ldoPhase     If true then compute phase correlograms.
 * @param[in] nsignals     Number of signals to transform.
 * @param[in] ntfSignals   Number of cross-correlations.
 * @param[in] npts         Number of points in signals.
 * @param[in] lxc          Length of the cross-correlations.
 * @param[in] x            Signals to transform.  This is an array of dimension
 *                         [npts x nsignals] with leading dimension npts.
 *
 * @param[out] xcs         The cross-correlations.  This is an array of
 *                         dimension [lxc x ntfSignals] with leading dimension
 *                         lxc.
 * @param[out] diffTime    Time to compute the reference solution (seconds).
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcfft_computeXCsWithISCL(const bool ldoPhase,
                             const int nsignals, const int ntfSignals,
                             const int npts, const int lxc,
                             const double x[],
                             double xcs[],
                             double *diffTime)
{
    double complex *tforms, *xcfd, xnum;
    double *work, xden;
    int i, ierr, indx, is, j, jndx, k, kndx, ntf;
    clock_t start;
    *diffTime = 0.0;
    ierr = EXIT_SUCCESS;
    // Fourier transform the signals
    ntf = lxc/2 + 1;
    tforms = (double complex *)
             calloc((size_t) (ntf*nsignals)+16, sizeof(double complex));
    start = clock();
    for (is=0; is<nsignals; is++)
    {
        indx = is*ntf;
        jndx = is*npts;
        fft_rfft64f_work_hack(npts, &x[jndx], lxc, ntf, &tforms[indx]);
    }
    // Now let's do it the dumb way
    xcfd = (double complex *) calloc((size_t) ntf+32, sizeof(double complex));
    work = (double *) calloc((size_t) lxc+32, sizeof(double));
    for (i=0; i<nsignals; i++)
    {
        for (j=i+1; j<nsignals; j++)
        {
            kndx = nsignals*i - ((i+1)*(i+2))/2 + j;
            if (kndx < 0 || kndx >= ntfSignals)
            {
                fprintf(stderr, "%s: kndx=%d out of bounds [0,%d]\n", __func__,
                        kndx, ntfSignals-1);
                return EXIT_FAILURE;
            }
            if (ldoPhase)
            {
                for (k=0; k<ntf; k++)
                {
                    xnum = tforms[i*ntf+k]*conj(tforms[j*ntf+k]);
                    xden = fmax(DBL_EPSILON, cabs(xnum));
                    xcfd[k] = xnum/xden;
                }
            }
            else
            {
                for (k=0; k<ntf; k++)
                {
                    xcfd[k] = tforms[i*ntf+k]*conj(tforms[j*ntf+k]);
                }
            }
            // Inverse fourier transform
            fft_irfft64z_work_hack(ntf, xcfd, lxc, work);
            // Need to shuffle the transform
            fft_fftshift64f_work_hack(lxc, work, &xcs[kndx*lxc]);
        }
    }
    free(tforms);
    free(xcfd);
    free(work);
    *diffTime = (double) (clock() - start)/CLOCKS_PER_SEC;
    return ierr;
}

void getReferenceSoln(int *npts, double *xc, double *phaseXC)
{
    *npts = 190;
    const double xcRef[190]={0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         1.000000000000e+00, 4.000000000000e+00, 1.000000000000e+01, 
         1.700000000000e+01, 2.600000000000e+01, 3.100000000000e+01, 
         2.900000000000e+01, 1.900000000000e+01, 2.100000000000e+01, 
         1.300000000000e+01, -3.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, -1.000000000000e+00, -1.000000000000e+00, 
         -3.000000000000e+00, -6.000000000000e+00, -7.000000000000e+00, 
         -2.000000000000e+00, -1.100000000000e+01, 2.000000000000e+00, 
         7.000000000000e+00, -7.000000000000e+00, 1.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, -1.000000000000e+00, 
         3.000000000000e+00, 6.000000000000e+00, 1.200000000000e+01, 
         1.600000000000e+01, 3.100000000000e+01, 4.000000000000e+00, 
         2.000000000000e+01, 3.000000000000e+00, 2.200000000000e+01, 
         -4.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         3.000000000000e+00, 4.000000000000e+00, 1.000000000000e+01, 
         1.200000000000e+01, 1.100000000000e+01, -1.100000000000e+01, 
         2.000000000000e+00, -3.700000000000e+01, -1.100000000000e+01, 
         3.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, -3.000000000000e+00, 1.000000000000e+00, 
         -5.000000000000e+00, -9.000000000000e+00, 3.000000000000e+00, 
         -5.000000000000e+00, -6.000000000000e+00, 1.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, -1.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, -3.000000000000e+00, 
         1.300000000000e+01, 6.000000000000e+00, 9.000000000000e+00, 
         1.200000000000e+01, 1.700000000000e+01, 1.800000000000e+01, 
         3.000000000000e+00, 1.100000000000e+01, 6.000000000000e+00, 
         4.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         9.000000000000e+00, 0.000000000000e+00, 1.400000000000e+01, 
         5.000000000000e+00, -1.200000000000e+01, 4.000000000000e+00, 
         -7.000000000000e+00, -1.200000000000e+01, -1.000000000000e+01, 
         -3.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 1.000000000000e+00, -7.000000000000e+00, 
         1.200000000000e+01, -8.000000000000e+00, -2.000000000000e+00, 
         -3.000000000000e+00, -2.000000000000e+00, 4.000000000000e+00, 
         -1.300000000000e+01, 6.000000000000e+00, -4.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, -3.000000000000e+00, 
         8.000000000000e+00, -1.200000000000e+01, 1.000000000000e+01, 
         -3.000000000000e+00, -1.700000000000e+01, 1.800000000000e+01, 
         -3.000000000000e+00, 1.000000000000e+00, 3.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         1.200000000000e+01, -1.400000000000e+01, 3.300000000000e+01, 
         -3.500000000000e+01, 2.800000000000e+01, -2.400000000000e+01, 
         2.200000000000e+01, -2.200000000000e+01, -1.100000000000e+01, 
         3.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00};
    const double phaseRef[190]={-6.151232118130e-02, 
         -1.271463664953e-01, 9.234959206218e-02, 5.188365449900e-02, 
         1.007946926863e-01, -1.384503278758e-01, -4.599154463037e-02, 
         -6.692471528557e-02, 3.462207254361e-01, 3.382015905608e-01, 
         4.211022201744e-01, -1.708913738597e-01, 3.141261296725e-01, 
         2.171428871850e-01, -5.416230630190e-01, 1.130543416856e-01, 
         1.694096565823e-02, 1.927866849194e-01, -5.206377219263e-02, 
         -6.110640595977e-02, 6.078960614197e-02, -3.636048250895e-02, 
         7.143334282381e-02, -5.332972772089e-02, 5.446901187884e-02, 
         -1.684042674602e-01, -1.622948985139e-01, -4.044897573956e-01, 
         -1.501133580202e-01, -5.443235724940e-01, 2.201560630917e-01, 
         4.897520977562e-01, -3.524816664157e-01, 9.411259099467e-03, 
         -9.483152694585e-02, 1.368126230662e-01, -4.441472720788e-02, 
         2.932638678475e-02, 1.513611694217e-01, -8.707414296398e-03, 
         4.476859850668e-02, 6.066889533386e-02, -2.312901472205e-01, 
         -1.121518234196e-01, 8.754046594817e-02, 1.274588810001e-02, 
         3.712140822885e-01, 6.228319826878e-01, -1.835418800844e-01, 
         1.252372535053e-01, -5.147736475058e-02, 4.936927552075e-01, 
         -1.558544249262e-01, -1.796303761686e-01, 7.367991703163e-02, 
         -1.181569610040e-01, -2.930616160822e-03, -4.657576245393e-02, 
         -3.346773480180e-02, -9.391716602589e-03, -3.150363464305e-03, 
         1.719709124728e-02, 6.837526324578e-03, 9.243102368635e-02, 
         1.206705026536e-01, 1.504667227417e-01, -2.162812312978e-01, 
         1.558790515985e-01, -8.723456120237e-01, -2.205387761067e-01, 
         1.788609656633e-01, -3.777589919814e-02, -1.540324956250e-01, 
         -1.117870940386e-01, -4.039521432802e-03, -1.295667687013e-02, 
         -4.610713963271e-02, -1.039259432700e-01, -1.496421396623e-02, 
         1.105277447623e-01, 8.872632026263e-02, 2.485513591691e-02, 
         -2.575439670438e-01, -7.296426217411e-01, 2.644206684547e-01, 
         -4.018669769521e-01, -1.927627669846e-01, 6.646648979167e-02, 
         2.133413804166e-01, 3.405787599792e-02, -1.344335879828e-01, 
         -1.150400574661e-01, 1.913617169847e-02, 1.077268828733e-01, 
         6.702860486500e-02, 1.475951151632e-01, 6.455309958769e-02, 
         -4.838415458234e-03, -1.600762141834e-01, -4.171380751621e-01, 
         3.854786165132e-01, 2.620456517960e-02, 1.320522810091e-01, 
         2.039503717482e-01, 3.606465838675e-01, 3.711038052805e-01, 
         -3.282224918508e-01, 2.356607492812e-01, 1.242011226350e-01, 
         3.584132309845e-02, -6.641256589359e-02, -8.087082308615e-02, 
         -2.460889460540e-01, 2.163598983249e-01, 1.322476437923e-01, 
         -1.039455309052e-02, -1.003739578050e-03, -2.234628828427e-01, 
         1.492312591243e-01, -1.002433712327e-01, 4.766514668362e-01, 
         1.419307847166e-01, -4.781172639093e-01, 1.588244163980e-01, 
         -2.570642993161e-01, -3.540899709393e-01, -3.792448224629e-01, 
         6.728977844981e-02, 3.736448003084e-02, 1.053783452716e-02, 
         -2.125248286460e-01, -8.238288729399e-02, -7.554904456361e-02, 
         1.597250253161e-01, -1.332154261830e-01, -1.194161018407e-01, 
         -2.116534279576e-02, 1.314730421445e-02, -2.572205078330e-02, 
         5.814538103175e-01, -2.357972765031e-01, -3.709761209623e-01, 
         -3.022458526333e-01, -9.454845429065e-02, 1.743675308363e-02, 
         -5.013379702591e-01, 8.370718296356e-02, -1.379571026681e-01, 
         -4.836907080211e-02, 9.174166169559e-02, 1.222881665516e-01, 
         -7.874913442103e-02, -1.309352236855e-02, 9.409001912605e-02, 
         -5.381068165009e-02, 7.185474923593e-02, -9.770780415090e-02, 
         1.677586032335e-01, -1.729897898263e-01, 1.732511518470e-01, 
         -2.963523833114e-01, -5.123603647543e-01, 5.618361459272e-01, 
         2.734104780777e-01, 2.607572880995e-01, 9.586954870677e-02, 
         2.115251540230e-01, 4.106353961454e-02, 1.036793646667e-01, 
         -9.210729235147e-03, 1.004292327388e-01, 2.220107895525e-03, 
         -1.109707703659e-01, -6.928752962202e-02, 1.892797336907e-01, 
         5.044686352879e-02, -2.024350617909e-02, 3.270074815130e-01, 
         -3.664379123018e-01, 1.114898153236e-01, -9.370640422693e-02, 
         1.633438905570e-01, -4.985933611715e-01, -5.989413859685e-01, 
         3.595926973357e-02, 3.197731547029e-02, -5.634579947733e-02, 
         1.188118241113e-01, -5.935997412746e-02, -1.566496583831e-01}; 
    memcpy(xc, xcRef, 190*sizeof(double));
    memcpy(phaseXC, phaseRef, 190*sizeof(double));
    return;
}
