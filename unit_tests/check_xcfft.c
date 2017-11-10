#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "xcloc_config.h"
#ifdef XCLOC_USE_MPI
#include "xcloc_xcfftMPI.h"
#endif
#include "xcloc_xcfft.h"
#ifdef XCLOC_USE_INTEL
#include <mkl_cblas.h>
#include <fftw/fftw3.h>
#else
#include <cblas.h>
#include <fftw3.h>
#endif
#include "iscl/fft/fft.h"
#include "iscl/iscl/iscl.h"

#define SEED 40935  // random seed - but predictable results
#define NPTS 1201   // 0.2 s window w/ rate of 6000 Hz
#define NSIGNALS 45 // make my computer sweat but don't tank the debugger 

int xcfft_hardwiredSerialTest(void);
int xcfft_randomSerialTest(void);
int xcfft_randomSerialTestMPI(MPI_Comm comm);
int xcfft_computeXCsWithISCL(const int nsignals, const int ntfSignals,
                             const int npts, const int lxc,
                             const double *__restrict__ x,
                             double *__restrict__ xcs,
                             double *diffTime);
double *xcfft_createRandomSignals(int *seed, const int nsignals, const int npts,
                                  int *ierr);


int main(int argc, char *argv[])
{
    int provided, ierr, ierrAll, myid, nprocs;
    const int master = 0;
    // Start MPI
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    iscl_init();
    // Have the master run the serial test
    if (myid == master)
    {
        ierr = xcfft_hardwiredSerialTest();
        if (ierr != EXIT_SUCCESS)
        {
            fprintf(stderr, "%s: Error in hardwiredSerialTest\n\n", __func__);
        }
        else
        {
            fprintf(stdout, "%s: Passed hardwired serial test\n\n", __func__);
        }
    }
    MPI_Bcast(&ierr, 1, MPI_INT, master, MPI_COMM_WORLD);
    if (ierr != EXIT_SUCCESS){goto QUIT;}
    // Have the master run the more impressive serial; this will create 
    // a reference solution for the MPI test 
    if (myid == master)
    {
        ierr = xcfft_randomSerialTest();
        if (ierr == EXIT_SUCCESS)
        {
            fprintf(stdout, "%s: Passed random serial test\n\n", __func__);
        }
        else
        {
            fprintf(stdout, "%s: Failed random serial test\n\n", __func__);
        }
    }
    MPI_Bcast(&ierr, 1, MPI_INT, master, MPI_COMM_WORLD);
    if (ierr != EXIT_SUCCESS){goto QUIT;}
    // Redo the random test but in parallel
    ierr = xcfft_randomSerialTestMPI(MPI_COMM_WORLD);
    ierr = abs(ierr);
    if (ierr != 0)
    {
        fprintf(stdout, "%s: Error on process: %d\n", __func__, myid);
        MPI_Abort(MPI_COMM_WORLD, 20);
    }
    //MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    //if (ierrAll != 0){goto QUIT;}

    // All done
    if (myid == master){fprintf(stdout, "\n%s: Passed XC tests!\n", __func__);}
QUIT:;
    iscl_finalize();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
//============================================================================//
int xcfft_randomSerialTestMPI(MPI_Comm comm)
{
    struct xcfftMPI_struct xcfftMPI;
    double *xcsRef = NULL;
    double *xrand = NULL;
    float *xcs = NULL;
    double diffTime, diffTimeRef, res, t0;
    const int master = 0;
    int i, ierr, ixc, lxc, ldxc, myid, nprocs, npts,
        nptsPad, nsignals, ntfSignals;
    MPI_Barrier(comm);
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &nprocs);
    // Always initialize the structure to nothing
    lxc = 0;
    memset(&xcfftMPI, 0, sizeof(struct xcfftMPI_struct));
    // Have master set some variables
    nsignals = 0;
    npts = 0;
    if (myid == master)
    {
        npts = NPTS;
        nptsPad = npts;
        lxc = 2*npts - 1;
        nsignals = NSIGNALS; 
        fprintf(stdout, "%s: Initializing...\n", __func__); 
    }
    // Initialize
    ierr = xcloc_xcfftMPI_initialize(npts, nptsPad, nsignals,
                                     comm, master, &xcfftMPI);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error initializing structure on rank %d\n",
                __func__, myid);
        return EXIT_FAILURE;
    }
    if (myid == master)
    {
        fprintf(stdout, "%s: Generating random numbers...\n" ,__func__);
        xrand = xcfft_createRandomSignals(NULL, nsignals, npts, &ierr);
    }
    MPI_Barrier(comm);
    // Scatter the data
    if (myid == master){fprintf(stdout, "%s: Scattering data...\n", __func__);}
    t0 = MPI_Wtime();
    ierr = xcloc_xcfftMPI_scatterData(master, nsignals, npts, npts,
                                      MPI_DOUBLE, xrand, &xcfftMPI);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error scattering data on rank %d\n",
                __func__, myid);
        return EXIT_FAILURE;
    }
    // Transform it
    if (myid == master){fprintf(stdout, "%s: Transforming data...\n", __func__);}
    ierr = xcloc_xcfftMPI_computePhaseCorrelation(&xcfftMPI);
    diffTime = MPI_Wtime() - t0;
    if (myid == master)
    {
        fprintf(stdout, "%s: Processing time %e (s)\n", __func__, diffTime);
    }
    // Gather the answer
    ldxc = xcfftMPI.lxc;
    lxc = xcfftMPI.lxc;
    ntfSignals = xcfftMPI.ntfSignals;
    if (myid == master)
    {
        xcs = (float *) calloc((size_t) (ldxc*ntfSignals), sizeof(float));
    }
    if (myid == master){fprintf(stdout, "%s: Gathering data...\n", __func__);}
    ierr = xcloc_xcfftMPI_gatherXCs(master, ntfSignals,
                                    ldxc, lxc,
                                    MPI_FLOAT, xcfftMPI, xcs);
    // Check it
    if (myid == master)
    {
/*
float *xcs2 = xcfftMPI.xcInv.y;
float *data = xcfftMPI.sigFwd.x;
printf("%e %e\n", data[2*], xrand[2*npts]);
*/
        xcsRef = (double *)
                 calloc((size_t) (xcfftMPI.ntfSignals*lxc), sizeof(double));
        fprintf(stdout, "%s: Computing reference solution...\n", __func__);
        ierr = xcfft_computeXCsWithISCL(nsignals, xcfftMPI.ntfSignals,
                                        npts, lxc, xrand, xcsRef, &diffTimeRef);
        if (ierr != EXIT_SUCCESS)
        {
            fprintf(stderr, "%s: Error computing ref XCs\n", __func__); 
            goto BCAST_ERROR;
        }
        for (ixc=0; ixc<xcfftMPI.ntfSignals; ixc++)
        {
            for (i=0; i<lxc; i++)
            {
                res = fabs((double) xcs[ldxc*ixc+i] - xcsRef[ixc*lxc+i]);
//fprintf(stdout, "%d %e %e %e\n", i, xcs[ldxc*ixc+i], xcs[ixc*lxc+i], res);
                if (res > 10.0*FLT_EPSILON)
                {
                    fprintf(stderr, "Failed on phase test: %d %d %f %e %e\n",
                            ixc, i, xcs[ldxc*ixc+i], xcsRef[ixc*lxc+i], res);
                    ierr = EXIT_FAILURE;
                    goto BCAST_ERROR;
                }
            }
        }
        free(xcsRef);
        fprintf(stdout, "%s: Reference time %e (s)\n", __func__, diffTimeRef);
        fprintf(stdout, "%s: Improvement: %4.2f pct\n", __func__,
                (diffTimeRef - diffTime)/diffTimeRef*100.0);
        fprintf(stdout, "\n");
    }
    // Finalize
BCAST_ERROR:;
    MPI_Bcast(&ierr, 1, MPI_INT, master, xcfftMPI.comm);
    if (myid == master){fprintf(stdout, "%s: Finalizing...\n", __func__);}
    ierr = xcloc_xcfftMPI_finalize(&xcfftMPI);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error finalizing on rank %d\n", __func__, myid);
        return EXIT_FAILURE;
    }
    if (xrand != NULL){free(xrand);}
    if (xcs != NULL){free(xcs);}
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
    // Set the seed - can be consistent for different tests or not
    srand(SEED);
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
 * @brief This is a more substantial test.  This performs the XC's on a 
 *        `random' sequence.  This does the computation carefully and
 *        quickly. 
 *
 * @result EXIT_SUCCESS indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcfft_randomSerialTest(void)
{
    struct xcfft_struct xcfft;
    int npts = NPTS; // 0.2 s window w/ rate of 6000 Hz
    int nptsPad = npts;
    int lxc = 2*npts - 1;
    int nsignals = NSIGNALS; // make it work but don't brick my computer in debug mode
    double *xcs;
    const float *xc;
    float *xcsAll;
    double diffTime, diffTimeRef, res;
    double *xrand;
    int i, ierr, indx, ixc;
    clock_t start;
    // Compute some random signals
    fprintf(stdout, "%s: Generating random numbers...\n", __func__);
    xrand = xcfft_createRandomSignals(NULL, nsignals, npts, &ierr);
    // Okay - let's run my stuff 
    memset(&xcfft, 0, sizeof(struct xcfft_struct));
    fprintf(stdout, "%s: Initializing...\n", __func__);
    ierr = xcloc_xcfft_initialize(npts, nptsPad, nsignals, &xcfft);
    if (ierr != 0)
    {   
        fprintf(stderr, "%s: Error initializing xcfft structure\n", __func__);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "%s: Setting data...\n", __func__);
    start = clock();
    ierr = xcloc_xcfft_setAllData(nsignals, npts, npts,
                                  XCLOC_DOUBLE_PRECISION, xrand,
                                  &xcfft);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error setting data\n", __func__);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "%s: Computing phase-correlations...\n", __func__);
    ierr = dales_xcfft_computePhaseCorrelation(&xcfft);
    if (ierr != 0)
    {   
        fprintf(stderr, "%s: Error computing phase correlations\n", __func__);
        return EXIT_FAILURE;
    }
    diffTime = (double) (clock() - start)/CLOCKS_PER_SEC;
    fprintf(stdout, "%s: Processing time %e (s)\n", __func__, diffTime);
    // Do it the straight forward wa 
    xcs = (double *) calloc((size_t) (xcfft.ntfSignals*lxc)+32, sizeof(double));
    ierr = xcfft_computeXCsWithISCL(nsignals, xcfft.ntfSignals,
                                    npts, lxc, xrand, xcs, &diffTimeRef);
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "%s: Failed to compute reference signals\n", __func__);
        return EXIT_FAILURE;
    } 
    // Now compare
//FILE *f = fopen("debug.txt", "w");
    xcsAll = (float *) calloc((size_t) (lxc*xcfft.ntfSignals), sizeof(float));
    xcloc_xcfft_getAllXCData(xcfft.ntfSignals, 
                                      lxc, lxc, XCLOC_SINGLE_PRECISION,
                                      xcfft, xcsAll); 
    for (ixc=0; ixc<xcfft.ntfSignals; ixc++)
    {
        xc = dales_xcfft_getXCDataPointer32f(lxc, ixc, xcfft, &ierr);
        for (i=0; i<lxc; i++)
        {
            indx = ixc*lxc + i;
            res = fabs((double) xc[i] - xcs[indx]);
//fprintf(f, "%d %e %e %e\n", i, xc[i], xcs[ixc*lxc+i], res);
            if (res > 10.0*FLT_EPSILON)
            {
                fprintf(stderr, "Failed on phase test: %d %d %f %e %e\n",
                        ixc, i, xc[i], xcs[indx], res);
                return EXIT_FAILURE;
            }
            res = fabsf(xc[i] - xcsAll[indx]);
            if (res > FLT_EPSILON)
            {
                fprintf(stderr, "Failed on getAllXC: %d %d %f %e %e\n",
                        ixc, i, xc[i], xcs[indx], res);
                return EXIT_FAILURE;
            }
        }
        xc = NULL;
    }
//fclose(f);
    fprintf(stdout, "%s: Reference time %e (s)\n", __func__, diffTimeRef);
    fprintf(stdout, "%s: Improvement: %4.2f pct\n", __func__,
            (diffTimeRef - diffTime)/diffTimeRef*100.0);
    fprintf(stdout, "\n");
    free(xcs);
    free(xcsAll);
    // Free it up 
    xcloc_xcfft_finalize(&xcfft);
    free(xrand);
    return EXIT_SUCCESS;
}
//============================================================================//
/*!
 * @brief Computes the phase correlations slowly but correctly.
 *
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
int xcfft_computeXCsWithISCL(const int nsignals, const int ntfSignals,
                             const int npts, const int lxc,
                             const double *__restrict__ x,
                             double *__restrict__ xcs,
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
        fft_rfft64f_work(npts, &x[jndx], lxc, ntf, &tforms[indx]); 
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
            for (k=0; k<ntf; k++)
            {
                xnum = tforms[i*ntf+k]*conj(tforms[j*ntf+k]);
                xden = fmax(DBL_EPSILON, cabs(xnum));
                xcfd[k] = xnum/xden;
            }
            // Inverse fourier transform
            fft_irfft64z_work(ntf, xcfd, lxc, work);
            // Need to shuffle the transform
            fft_fftshift64f_work(lxc, work, &xcs[kndx*lxc]);
        }
    }   
    free(tforms);
    free(xcfd);
    free(work);
    *diffTime = (double) (clock() - start)/CLOCKS_PER_SEC;
    return ierr;
}
//============================================================================//
/*!
 * @brief This is the most basic serial test with hard-wired values computed
 *        by Python.  If you fail this test then you should stop immediately.
 * @result EXIT_SUCCESS indicates a success.
 *
 */
int xcfft_hardwiredSerialTest(void)
{
    int npts = 6;
    int nptsPad = 10; 
    int nsignals = 5;
    double res;
    struct xcfft_struct xcfft;
    int ierr, ixc;
    const int lxc=19;
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
    // Initialize the cross-correlation structure
    memset(&xcfft, 0, sizeof(struct xcfft_struct));
    printf("Initializing...\n");
    ierr = xcloc_xcfft_initialize(npts, nptsPad, nsignals, &xcfft);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error initializing xcfft structure\n", __func__);
        return EXIT_FAILURE;
    }   
    // Create some bogus time series
    const float a[10] = { 1, 2, 3, 4, 5, -1,  0, 0,  0,  0}; 
    const float b[10] = { 3, 2, 1, 3, 2,  1,  0, 0,  0,  0}; 
    const float c[10] = {-1, 2,-1,-2, 1, -1,  0, 0,  0,  0}; 
    const float d[10] = { 4,-2, 3,-1, 5, -1,  0, 0,  0,  0}; 
    const float e[10] = { 0,-3,-4, 5,-2,  3,  0, 0,  0,  0}; 
    const float *xc;
    float *xall = (float *) calloc((size_t) (npts*nsignals)+16, sizeof(float));
    int i; 
    for (i=0; i<nsignals; i++)
    {   
        if (i == 0){cblas_scopy(npts, a, 1, &xall[npts*i], 1);}
        if (i == 1){cblas_scopy(npts, b, 1, &xall[npts*i], 1);}
        if (i == 2){cblas_scopy(npts, c, 1, &xall[npts*i], 1);}
        if (i == 3){cblas_scopy(npts, d, 1, &xall[npts*i], 1);}
        if (i == 4){cblas_scopy(npts, e, 1, &xall[npts*i], 1);}
    }
    // Set it on the data
    printf("%s: Setting data...\n", __func__);
    ierr = xcloc_xcfft_setAllData(nsignals, npts, npts,
                                  XCLOC_SINGLE_PRECISION, xall,
                                  &xcfft);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error setting signal\n", __func__);
        return EXIT_FAILURE;
    }
    // Compute the cross-correlations
    fprintf(stdout, "%s: Computing cross-correlations...\n", __func__);
    ierr = dales_xcfft_computeFFTCrossCorrelation(&xcfft);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing cross-correlations\n", __func__);
        return EXIT_FAILURE;
    }
    // Check it
    for (ixc=0; ixc<xcfft.ntfSignals; ixc++)
    {
        xc = dales_xcfft_getXCDataPointer32f(lxc, ixc, xcfft, &ierr);
        for (i=0; i<lxc; i++)
        {
            res = fabs((double) xc[i] - xcRef[ixc*lxc+i]);
            if (res > 90.0*FLT_EPSILON)
            {
                fprintf(stderr, "Failed on xctest: %d %d %f %e %e\n",
                        ixc, i, xc[i], xcRef[ixc*lxc+i], res);
                return EXIT_FAILURE;
            }
        } 
        xc = NULL;
    } 
    // Compute the phase-correlations
    fprintf(stdout, "%s: Computing phase-correlations...\n", __func__);
    ierr = dales_xcfft_computePhaseCorrelation(&xcfft);
    if (ierr != 0)
    {   
        fprintf(stderr, "%s: Error computing phase correlations\n", __func__);
        return EXIT_FAILURE;
    }
    // Check it
    for (ixc=0; ixc<xcfft.ntfSignals; ixc++)
    {   
        xc = dales_xcfft_getXCDataPointer32f(lxc, ixc, xcfft, &ierr);
        for (i=0; i<lxc; i++)
        {
            res = fabs((double) xc[i] - phaseRef[ixc*lxc+i]);
            if (res > 10.0*FLT_EPSILON)
            {
                fprintf(stderr, "Failed on phase test: %d %d %f %e %e\n",
                        ixc, i, xc[i], phaseRef[ixc*lxc+i], res);
                return EXIT_FAILURE;
            }
        }
        xc = NULL;
    }
    // Free the memory
    xcloc_xcfft_finalize(&xcfft);
    free(xall);
    return EXIT_SUCCESS;
}


