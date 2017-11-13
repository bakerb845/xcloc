#define _BSD_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "xcloc_xcfft.h"
#include "xcloc_rmsFilter.h"
#include "xcloc_xcfftMPI.h"
#include "xcloc_migrate.h"
#include "xcloc_hdf5.h"
#include "acoustic_greens.h"
#include <ipps.h>
#include "sacio.h" 
#include "iscl/signal/window.h"
#include "iscl/array/array.h"
#include "iscl/fft/fft.h"
#include "iscl/time/time.h"
#include "iscl/os/os.h"
#include <mkl_vsl.h>
#include <mkl_cblas.h>
#ifdef XCLOC_PROFILE
#include "advisor-annotate.h"
#endif

#define SEED 518815
#define DTOUT 0.001
#define TIME0 1509409617
#define SET_DIRECTORY(outdir, dirRoot) \
    memset(dirRoot, 0, PATH_MAX*sizeof(char)); \
    if (outdir == NULL) \
    {   \
        strcpy(dirRoot, "./\0"); \
    }   \
    else \
    {   \
        if (strlen(outdir) == 0) \
        { \
            strcpy(dirRoot, "./\0"); \
        } \
        else \
        { \
            strcpy(dirRoot, outdir); \
            if (dirRoot[strlen(dirRoot)] != '/') \
            { \
                strcat(dirRoot, "/\0"); \
            } \
        } \
    }  \
    if (!os_path_isdir(dirRoot)) \
    {   \
        os_makedirs(dirRoot); \
    }

#define LNORM false     // Don't normalize ricker wavelet (max will be 1)
#define LSHIFT = true;     // Make wavelet start at time 0 

int create2DReceiversAndTravelTimes(const int seed, const int nrec,
                                    const double x0, const double x1,
                                    const double y0, const double y1,
                                    const double z0, const double vel,
                                    const double xs[3],
                                    double *__restrict__ xr,
                                    double *__restrict__ ttimes);
int compute2DGreensFunctions(const int nrec, const int nptsSig,
                             const double fcent, const double dt, 
                             const bool lnorm, const bool lshift,
                             const double vel, const double rho,
                             const double Q, const double snr,
                             const double pct,
                             const double xs[6], 
                             const double *__restrict__ xr, 
                             double **obsOut,
                             double **obs2,
                             double **obsNoisyOut);
int compute2DTravelTimes(const int nx, const int ny, const int nz, 
                         const double x0, const double y0, const double z0, 
                         const double dx, const double dy, const double dz, 
                         const double vel, const double xs[3],
                         float *__restrict__ ttimes);
int writeSynthetics(const char *outdir,
                    const int nrec, const int npts,
                    const double dt, 
                    const double xs[3],
                    const double *__restrict__ xr, 
                    const double *__restrict__ ttimes,
                    const double *obs);
int writeCrossCorrelations(const bool lfiltered, const char *outdir,
                           struct xcfft_struct xcfft);
int addGaussianNoise(const int nrec, const int npts, const double snr,
                     const double *__restrict__ obs,
                     double *__restrict__ obsNoisy);

int main(int argc, char *argv[])
{
    struct xclocHDF5Grid_struct h5io;
    struct xcfft_struct xcfft;
    struct xcfftMPI_struct xcfftMPI;
    struct xcfftRMSFilter_struct rms;
    struct migrate_struct migrate;
    const double snr = -14.0; //10.0; //-14.0 works
    const double pct = 8.0;  // 5 pct taper
    int chunkSize = 2048;
    int nrec = 30; // TODO make 100 or 200
    double dt = 1.0/6000.0; // Sampling rate is 6000 Hz
    double twin = 1.0;      // Window is 0.2 seconds
    double fcent = 800.0;  // Dominant resolution is vmin/fcent ~ 5.0m (for plotting)
    bool lnorm = false;     // Don't normalize ricker wavelet (max will be 1)
    bool lshift = true;     // Make wavelet start at time 0 
    int winLen = (int) ceil((1.0/fcent)/dt); // make the window about 3 cycles 
    if (winLen%2 == 0){winLen = winLen + 1;} // Make it odd
    //int npts = (int) (round(twin/dt)) + 1; // Number of points in time series
    double x0 = 0.0;    // Model origin in x is 0 km 
    double x1 = 1000.0; // Model extent in x is 1 km
    double y0 = 0.0;    // Model origin in y is 0 km
    double y1 = 1000.0; // Model extent in y is 1 km
    double z0 = 0.0;    // Problem is 2d - make z equal to 0.0
    double z1 = 0.0;    // Problem is 2d - make z equal to 0.0
    double vel = 3100.0; // constant velocity (m/s)
    double rho = 2700.0; // constant density (kg/m**3)
    double Q = 1.e5;     // effectively remove damping
    double tmodel = 0.5; // max modeling time is the traveltime from the
                         // furthest point in the medium to the reciever
                         // plus some
    int nptsSig = (int) (round(tmodel/dt)) + 1;
    int nx = 512;
    int ny = 512;
    int nz = 1;
    int ngrd = nx*ny*nz;
    double dx = (x1 - x0)/(double) (nx - 1); // should be less than 2m
    double dy = (y1 - y0)/(double) (ny - 1); // should be less than 2m
    double dz = 0.0;
    // Set the receiver location to the center of the model
    double xs[6] = {x0 + (double) (nx/2)*dx,
                    y0 + (double) (ny/2)*dy,
                    z0,// Source has to be at zero for the 2d example
                    x0 + (double) (5*nx/8)*dx,
                    y0 + (double) (3*nx/8)*dy,
                    z0};
    // Randomly create the receiver locations and theoretical picks
    int ierr, irec, myid, nprocs, provided;
    double *obsPtr, *obs, *obs2, *obsNoisy, *xr, *ttimes;
    float *tTable;
    const int master = 0;
    // Initialize MPI 
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 
    // Initialize some things
    memset(&xcfft,   0, sizeof(struct xcfft_struct));
    memset(&migrate, 0, sizeof(struct migrate_struct));
    memset(&h5io,    0, sizeof(struct xclocHDF5Grid_struct));
    // Have the master compute the problem set-up 
    obs = NULL;
    obs2 = NULL;
    obsNoisy = NULL;
    xr = NULL;
    xr = (double *) calloc((size_t) (3*nrec), sizeof(double));
    ttimes = (double *) calloc((size_t) nrec, sizeof(double));
    if (myid == master)
    {
        ierr = create2DReceiversAndTravelTimes(SEED, nrec,
                                               x0, x1, y0,  y1, z0, vel,
                                               xs, xr, ttimes);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error creating recv locs and ttimes\n",
                    __func__);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        fprintf(stdout, "%s: Setting cross-correlations...\n", __func__);
        // Compute the Green's functions using a Ricker wavelet
        fprintf(stdout, "%s: Computing Green's functions...\n", __func__);
        ierr = compute2DGreensFunctions(nrec, nptsSig, 
                                        fcent, dt,
                                        lnorm, lshift,
                                        vel, rho,
                                        Q, snr,
                                        pct,
                                        xs, xr, &obs, &obs2, &obsNoisy);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to compute synthetics\n", __func__);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        // Write the Green's functions
        fprintf(stdout, "%s: Writing Green's functions...\n", __func__);
        ierr = writeSynthetics("./greens", nrec, nptsSig, dt, xs, xr,
                               ttimes, obs);
        ierr = writeSynthetics("./greensNoisy", nrec, nptsSig, dt, xs, xr,
                               ttimes, obsNoisy);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Cross-correlate the signals
    fprintf(stdout, "%s: Computing cross-correlations...\n", __func__);
    // Initialize
    int nsignals = nrec;
    int npts = nptsSig;
    int nptsPad = nptsSig;
    ierr = xcloc_xcfft_initialize(npts, nptsPad, nsignals, &xcfft);
    ierr = xcloc_xcfftMPI_initialize(npts, nptsPad, nsignals, MPI_COMM_WORLD,
                                     master, NULL, &xcfftMPI);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error initializing xcfft structure\n", __func__);
        return EXIT_FAILURE;
    }
    ierr = xcloc_rmsFilter_initialize(winLen, xcfftMPI.xcInv.lxc,
                                      xcfftMPI.xcInv.precision,
                                      &rms); 
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error creating RMS filter\n", __func__);
        return EXIT_FAILURE;
    }
    ierr = dales_xcfft_createRMSFilter(winLen, &xcfft);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error creating RMS filter\n", __func__);
        return EXIT_FAILURE;
    }
    // Loop on perfect and noisy signals
    double tall = 0.0; 
    int i, ncopy;
    ncopy = xcfft.ntfSignals*xcfft.dataOffset;
    float *yf1, *yf2;
    yf1 = (float *) calloc((size_t) ncopy, sizeof(float));
    yf2 = (float *) calloc((size_t) ncopy, sizeof(float));
    for (i=0; i<2; i++)
    {
        time_tic();
        // Get a pointer to the signal
        obsPtr = NULL;
        if (i == 0)
        {
            obsPtr = obs;
        }
        else if (i == 1)
        {
            obsPtr = obsNoisy;
        }
        else
        {
            obsPtr = obs2;
        } 
        if (myid == master)
        {
            ierr = xcloc_xcfft_setAllData(nsignals, nptsSig, nptsSig,
                                          XCLOC_DOUBLE_PRECISION, obsPtr,
                                          &xcfft);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error setting data\n", __func__);
                MPI_Abort(MPI_COMM_WORLD, 30);
            }
        }
        ierr = xcloc_xcfftMPI_scatterData(master, nsignals, nptsSig, nptsSig,
                                          MPI_DOUBLE, obsPtr, &xcfftMPI);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error scattering data\n", __func__);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        ierr = xcloc_xcfft_computePhaseCorrelation(&xcfft);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing phase correlations\n",
                    __func__);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        ierr = xcloc_xcfftMPI_computePhaseCorrelation(&xcfftMPI);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing phase correlations\n",
                    __func__);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        ierr = xcloc_rmsFilter_apply(xcfft.ntfSignals,
                                     xcfft.dataOffset,
                                     xcfft.lxc,
                                     xcfft.precision,
                                     &rms,
                                     xcfft.y, xcfft.yfilt);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing RMS window\n", __func__);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        tall = tall + time_toc();
        fprintf(stdout, "%s: Writing correlations...\n", __func__);
        if (i == 0)
        {
            ippsCopy_32f(xcfft.yfilt, yf1, ncopy);
            writeCrossCorrelations(false, "./xcorr", xcfft);
            writeCrossCorrelations(true,  "./xcorrProc", xcfft);
        }
        else
        {
            ippsCopy_32f(xcfft.yfilt, yf2, ncopy);
            writeCrossCorrelations(false, "./xcorrNoise", xcfft);
            writeCrossCorrelations(true,  "./xcorrNoiseProc", xcfft);
        }
    } // Loop on cross-correlograms
    if (myid == master)
    {
        fprintf(stdout, "%s: Average cross-correlation time: %e\n",
                __func__, tall*0.5);
        fprintf(stdout, "%s: Initializing H5 archive...\n", __func__);
        ierr = xcloc_h5ioGrid_open("./migrate.h5", "./migrate.xdmf",
                                   "/migrate",
                                   nx, ny, nz,
                                   dx, dy, dz,
                                   x0, y0, z0,
                                   &h5io);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error initializing H5 file\n", __func__);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
    }
    // Migrate the cross-correlations
    int ldxc = xcfft.dataOffset;
    int lxc = xcfft.lxc;
    int nxc = xcfft.ntfSignals;
    ierr = xcloc_migrate_initialize(nrec, ngrd, nxc, lxc, chunkSize,
                                    dt, xcfft.xcPairs, &migrate); 
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error initializing migration struct\n", __func__);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    // Set the travel time tables
    if (myid == master)
    {
        fprintf(stdout, "%s: Computing travel time tables...\n", __func__);
    }
    tTable = (float *) calloc((size_t) ngrd, sizeof(float));
    for (irec=0; irec<nrec; irec++)
    {
        if (myid == master)
        {
            compute2DTravelTimes(nx, ny, nz, 
                                 x0, y0, z0, 
                                 dx, dy, dz, 
                                 vel, &xr[3*irec], tTable);
        }
        ierr = xcloc_migrate_setTable32f(irec, ngrd, tTable, &migrate);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to set it'th=%d table\n",
                    __func__, irec);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
    }
    fprintf(stdout, "%s: Setting cross-correlations...\n", __func__);
    for (i=0; i<2; i++)
    {
        if (i == 0)
        {
            ierr = xcloc_migrate_setCrossCorrelations(ldxc, lxc, nxc, yf1,
                                                      XCLOC_SINGLE_PRECISION,
                                                      &migrate);
        }
        else
        {
            ierr = xcloc_migrate_setCrossCorrelations(ldxc, lxc, nxc, yf2,
                                                      XCLOC_SINGLE_PRECISION, 
                                                      &migrate);
        }
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to set xc %d\n", __func__, i+1);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        fprintf(stdout, "%s: Computing image %d...\n", __func__, i+1);
        time_tic();
        ANNOTATE_SITE_BEGIN("blocking");
        xcloc_migrate_computeMigrationImage(&migrate);
        ANNOTATE_SITE_END();
        fprintf(stdout, "Migration time %e (s)\n", time_toc());
        const float *image = xcloc_migrate_getImagePointer(ngrd, migrate, &ierr);
        fprintf(stdout, "writing h5 file...\n");
        if (i == 0)
        {
           ierr = xcloc_h5ioGrid_writeDataSet32f("image", nx, ny, nz,
                                                 image, &h5io);
        }
        else
        {
           ierr = xcloc_h5ioGrid_writeDataSet32f("imageNoise", nx, ny, nz,
                                                 image, &h5io);
        }
        image = NULL;
    }
    ierr = xcloc_h5ioGrid_close(&h5io);
    ierr = xcloc_xcfft_finalize(&xcfft);
    ierr = xcloc_xcfftMPI_finalize(&xcfftMPI);
    ierr = xcloc_migrate_finalize(&migrate);
    ierr = xcloc_rmsFilter_finalize(&rms);
    free(yf1);
    free(yf2);
    free(xr);
    free(obs);
    free(obs2);
    free(obsNoisy);
    free(tTable);
    free(ttimes);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
//============================================================================//
int create2DReceiversAndTravelTimes(const int seed, const int nrec,
                                    const double x0, const double x1,
                                    const double y0, const double y1,
                                    const double z0, const double vel,
                                    const double xs[3],
                                    double *__restrict__ xr,
                                    double *__restrict__ ttimes)
 
{
    double dist;
    int irec; 
    // Randomly create the receiver locations
    srand(seed);
    for (irec=0; irec<nrec; irec++)
    {
        xr[3*irec+0] = ((double) rand()/RAND_MAX)*(x1 - x0);
        xr[3*irec+1] = ((double) rand()/RAND_MAX)*(y1 - y0);
        xr[3*irec+2] = z0; // Receiver has to be zero for the 2d example
        if (xr[3*irec+0] < x0 || xr[3*irec+0] > x1)
        {   
            fprintf(stderr, "%s: Invalid x location\n", __func__);
            return EXIT_FAILURE;
        }
        if (xr[3*irec+1] < y0 || xr[3*irec+1] > y1)
        {   
            fprintf(stderr, "%s: Invalid y location\n", __func__);
            return EXIT_FAILURE;
        }
        dist = sqrt( pow(xs[0] - xr[3*irec+0], 2)
                   + pow(xs[1] - xr[3*irec+1], 2)
                   + pow(xs[2] - xr[3*irec+2], 2));
        ttimes[irec] = dist/vel;
        //printf("%d %f %f %f\n", irec, xr[3*irec+0], xr[3*irec+1], xr[3*irec+2]);
    }
    return 0;
}
//============================================================================//
int compute2DGreensFunctions(const int nrec, const int nptsSig,
                             const double fcent, const double dt,
                             const bool lnorm, const bool lshift,
                             const double vel, const double rho,
                             const double Q, const double snr,
                             const double pct,
                             const double xs[6],
                             const double *__restrict__ xr,
                             double **obsOut,
                             double **obsOut2,
                             double **obsNoisyOut)
{
    double *obs, *obs2, *obsNoisy, *stf, *window, *taper;
    int i, ierr, irec;
    const bool lsym = true; // Create symmetric windows
    const int m = (int) MIN(nptsSig - 2, MAX(0, round(pct/100.0*nptsSig)));
    const int mp12 = (m+1)/2;
    // Compute the Green's functions using a Ricker wavelet
    fprintf(stdout, "%s: Computing synthetics...\n", __func__);
    stf = (double *) calloc((size_t) nptsSig, sizeof(double));
    ierr = dales_unitTests_computeRickerWavelet(nptsSig, dt, fcent,
                                                lnorm, lshift, stf);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to compute wavelet\n", __func__);
        return EXIT_FAILURE;
    }
    obs  = (double *) calloc((size_t) (nptsSig*nrec), sizeof(double));
    obs2 = (double *) calloc((size_t) (nptsSig*nrec), sizeof(double));
    // Compute the Green's functions
    ierr = dales_unitTests_acousticGreensLineSource(nrec, vel, rho, Q,
                                                    nptsSig, dt, xs, xr,
                                                    stf, obs);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to compute Greens functions", __func__);
        return EXIT_FAILURE;
    }
    // Compute the Green's functions
    ierr = dales_unitTests_acousticGreensLineSource(nrec, vel, rho, Q,
                                                    nptsSig, dt, &xs[3], xr, 
                                                    stf, obs2);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to compute Greens functions", __func__);
        return EXIT_FAILURE;
    }
    cblas_daxpy(nptsSig*nrec, 1.1, obs, 1, obs2, 1);
    // Make the noisy signals
    obsNoisy = (double *) calloc((size_t) (nptsSig*nrec), sizeof(double));
    ierr = addGaussianNoise(nrec, nptsSig, snr, obs, obsNoisy);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to make noisy signals\n", __func__);
        return EXIT_FAILURE;
    }
    // Taper ends of noisy signals
    window = (double *) calloc((size_t) (m+2), sizeof(double));
    signal_hammingWindow64f_work(m, lsym, &window[1]);
    taper = (double *) calloc((size_t) nptsSig, sizeof(double));
    for (irec=0; irec<nrec; irec++)
    {
        for (i=0; i<mp12; i++)
        {
            obsNoisy[irec*nptsSig+i]  = obsNoisy[irec*nptsSig+i]*window[i];
        }
        for (i=0; i<mp12; i++)
        {
            obsNoisy[(irec+1)*nptsSig-1-i]
               = obsNoisy[(irec+1)*nptsSig-1-i]*window[m+1-i];
        }
    }
    *obsOut = obs;
    *obsOut2 = obs2;
    *obsNoisyOut = obsNoisy;
    free(taper);
    free(window);
    free(stf);
    return 0;
}
//============================================================================//
int writeCrossCorrelations(const bool lfiltered, const char *outdir,
                           struct xcfft_struct xcfft)
{
    char dirRoot[PATH_MAX], filename[PATH_MAX];
    const char *khole = "00\0";
    const char *knetwk = "GRNS\0";
    const char *kcmpnm = "HHZ\0";
    char kstnm[8];
    const float *xc;
    int ierr, ixc, npts, st1, st2;
    struct sacData_struct sac;
    SET_DIRECTORY(outdir, dirRoot);
    memset(&sac, 0, sizeof(struct sacData_struct));
    sacio_setDefaultHeader(&sac.header);
    npts = xcfft.lxc;
    sacio_setCharacterHeader(SAC_CHAR_KNETWK, knetwk, &sac.header);
    sacio_setCharacterHeader(SAC_CHAR_KCMPNM, kcmpnm, &sac.header);
    sacio_setCharacterHeader(SAC_CHAR_KHOLE,  khole,  &sac.header);
    sacio_setIntegerHeader(SAC_INT_NPTS, npts, &sac.header);
    sacio_setFloatHeader(SAC_FLOAT_DELTA, DTOUT, &sac.header);
    sacio_setEpochalStartTime(TIME0-(double) npts*DTOUT/2, &sac.header);
    sac.data = sacio_malloc64f(npts);
    sac.npts = npts;
    for (ixc=0; ixc<xcfft.ntfSignals; ixc++)
    {
        st1 = xcfft.xcPairs[2*ixc+0];
        st2 = xcfft.xcPairs[2*ixc+1];
        sprintf(kstnm, "S%03d%03d", st1, st2);
        sacio_setCharacterHeader(SAC_CHAR_KSTNM, kstnm, &sac.header);
        if (lfiltered)
        {
            xc = (float *) &xcfft.yfilt[ixc*xcfft.dataOffset];
        }
        else
        {
            xc = dales_xcfft_getXCDataPointer32f(xcfft.lxc, ixc, xcfft, &ierr);
        }
        ippsConvert_32f64f(xc, sac.data, npts);
        xc = NULL;
        memset(filename, 0, PATH_MAX*sizeof(char));
        sprintf(filename, "%s/%s.%s.%s.%s.SAC", dirRoot,
                sac.header.knetwk, sac.header.kstnm, sac.header.kcmpnm,
                sac.header.khole);
        sacio_writeTimeSeriesFile(filename, sac);
    }
    sacio_free(&sac);
    return 0;
}
//============================================================================//
int writeSynthetics(const char *outdir,
                    const int nrec, const int npts,
                    const double dt,
                    const double xs[3],
                    const double *__restrict__ xr,
                    const double *__restrict__ ttimes,
                    const double *obs)
{
    char dirRoot[PATH_MAX], filename[PATH_MAX];
    const char *khole = "00\0";
    const char *knetwk = "GRNS\0"; 
    const char *kcmpnm = "HHZ\0";
    char kstnm[8];
    double dist, pick;
    int i;
    struct sacData_struct sac;
    SET_DIRECTORY(outdir, dirRoot);
    memset(&sac, 0, sizeof(struct sacData_struct));
    sacio_setDefaultHeader(&sac.header);
    sacio_setCharacterHeader(SAC_CHAR_KNETWK, knetwk, &sac.header);
    sacio_setCharacterHeader(SAC_CHAR_KCMPNM, kcmpnm, &sac.header);
    sacio_setCharacterHeader(SAC_CHAR_KHOLE,  khole,  &sac.header);
    sacio_setCharacterHeader(SAC_CHAR_KT0,    "P\0",  &sac.header);
    sacio_setIntegerHeader(SAC_INT_NPTS, npts, &sac.header);
    sacio_setFloatHeader(SAC_FLOAT_DELTA, DTOUT, &sac.header);
    sacio_setEpochalStartTime(TIME0, &sac.header);
    sac.data = sacio_malloc64f(npts); 
    sac.npts = npts;
    for (i=0; i<nrec; i++)
    {
        memset(kstnm, 0, 8*sizeof(char));
        sprintf(kstnm, "S%03d", i);
        dist = sqrt( pow(xs[0] - xr[3*i+0], 2)
                   + pow(xs[1] - xr[3*i+1], 2) 
                   + pow(xs[2] - xr[3*i+2], 2) );
        sacio_setCharacterHeader(SAC_CHAR_KSTNM, kstnm, &sac.header);
        pick = ttimes[i]/dt*DTOUT;
        sacio_setFloatHeader(SAC_FLOAT_T0,    pick, &sac.header);
        sacio_setFloatHeader(SAC_FLOAT_DIST,  dist, &sac.header);
        sacio_setFloatHeader(SAC_FLOAT_GCARC, dist, &sac.header);
        ippsCopy_64f(&obs[i*npts], sac.data, sac.npts); 
        memset(filename, 0, PATH_MAX*sizeof(char));
        sprintf(filename, "%s/%s.%s.%s.%s.SAC", dirRoot,
                sac.header.knetwk, sac.header.kstnm, sac.header.kcmpnm,
                sac.header.khole);
        sacio_writeTimeSeriesFile(filename, sac);
    }
    sacio_free(&sac);
    return EXIT_SUCCESS;
}  
//============================================================================//
int compute2DTravelTimes(const int nx, const int ny, const int nz,
                         const double x0, const double y0, const double z0,
                         const double dx, const double dy, const double dz,
                         const double vel, const double xs[3],
                         float *__restrict__ ttimes)
{
    double dist, dist2, slow, x, y, z;
    int indx, ix, iy, iz;
    slow = 1.0/vel;
    for (iz=0; iz<nz; iz++)
    {
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                x = (double) ix*dx + x0;
                y = (double) iy*dy + y0;
                z = (double) iz*dz + z0;
                dist2 = pow(x - xs[0], 2)
                      + pow(y - xs[1], 2)
                      + pow(z - xs[2], 2);
                dist = sqrt(dist2);
                indx = iz*nx*ny + nx*iy + ix;
                ttimes[indx] = (float) (dist*slow);
            }
        }
    }
    return EXIT_SUCCESS;
}
//============================================================================//
int addGaussianNoise(const int nrec, const int npts, const double snr,
                     const double *__restrict__ obs,
                     double *__restrict__ obsNoisy)
{
    VSLStreamStatePtr stream;
    double en, es, xfact;
    double *noise;
    const double xmean = 0.0;
    const double xstd = 1.0;
    int irec;
    noise = (double *) calloc((size_t) npts, sizeof(double));
    cblas_dcopy(nrec*npts, obs, 1, obsNoisy, 1); // Copy obs -> obsNoisy
    vslNewStream(&stream, VSL_BRNG_MT19937, 777);
    // Compute the max signal energy
    es = 0.0;
    for (irec=0; irec<nrec; irec++)
    {
        es = fmax(es, cblas_dnrm2(npts, &obs[irec*npts], 1));
    }
    // Add noise
    for (irec=0; irec<nrec; irec++)
    {
        // Compute the energy in the signal
        vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, npts, noise,
                      xmean, xstd);
        en = cblas_dnrm2(npts, noise, 1);
        xfact = (1.0/en)*(es/pow(10.0, 0.05*snr));
        // obsNoisy = obsNoisy + xfact*noise
        cblas_daxpy(npts, xfact, noise, 1, &obsNoisy[irec*npts], 1);  
    }
    vslDeleteStream(&stream);
    return 0;
}
