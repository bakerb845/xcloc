#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ipps.h>
#include "xcloc_finter.h"
#include "acousticGreens2D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static
int computeTravelTimeTable(const int nx, const int ny, const int nz, 
                           const double vel,
                           const double x0, const double y0, const double z0, 
                           const double dx, const double dy, const double dz, 
                           double xr, double yr, double zr, 
                           double ttable[]);
static
int computeRandomReceiverLocations(const int nrec,
                                   const double x0, const double y0, const double z0, 
                                   const double x1, const double y1, const double z1, 
                                   double *xr);

static
int compute2DGreensFunctions(const int nsrc, const int nrec, const int nptsSig,
                             const double fcent, const double dt, 
                             const bool lnorm, const bool lshift,
                             const double vel, const double rho,
                             const double Q,
                             const double pct,
                             const double srcScale[],
                             const double xs[],
                             const double xr[],
                             double **obsOut);
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

int test_serial_dsmLocation(void)
{
    const int nrec = 20;
    double dt = 1.0/6000.0; // Sampling rate is 6000 Hz
    double fcent = 400.0;  // Dominant resolution is vmin/fcent ~ 5.0m (for plotting)
    bool lnorm = false;     // Don't normalize ricker wavelet (max will be 1)
    bool lshift = true;     // Make wavelet start at time 0 
    const double pct = 8.0;  // 8 pct taper
    double x0 = 0.0;    // Model origin in x is 0 km 
    double x1 = 1000.0; // Model extent in x is 1 km
    double y0 = 0.0;    // Model origin in y is 0 km
    double y1 = 1000.0; // Model extent in y is 1 km
    double z0 = 0.0;    // Problem is 2d - make z equal to 0.0
    double z1 = 0.0;    // Problem is 2d - make z equal to 0.0
    double vel = 3100.0; // constant velocity (m/s)
    double rho = 2700.0; // constant density (kg/m**3)
    double Q = 9.e2;     // add some damping
    double tmodel = 0.5; // max modeling time is the traveltime from the
                         // furthest point in the medium to the reciever
                         // plus some
    int nptsSig = (int) (round(tmodel/dt)) + 1;
    // rest of geometry
    int nx = 512;
    int ny = 512;
    int nz = 1;
    int ngrd = nx*ny*nz;
    double dx = (x1 - x0)/(double) (nx - 1); // should be less than 2m
    double dy = (y1 - y0)/(double) (ny - 1); // should be less than 2m
    double dz = 0.0;
    int i, ierr, it;
    // Set the receiver location to the center of the model
    const int nsrc = 2; 
    const double srcScale[2] = {1, 1.1};
    double xs[6] = {x0 + (double) (2*nx/7)*dx,
                    y0 + (double) (6*ny/7)*dy,
                    z0,// Source has to be at zero for the 2d example
                    x0 + (double) (5*nx/8)*dx,
                    y0 + (double) (3*ny/8)*dy,
                    z0};
    double *xr = NULL;
    double *obs = NULL;
    xr = (double *) calloc((size_t) (3*nrec), sizeof(double));
    ierr = computeRandomReceiverLocations(nrec,
                                          x0, y0, z0,
                                          x1, y1, z1,
                                          xr);
    CHKERR(ierr, "failed making receiver locations");
    // Scatter the receivers
    ierr = compute2DGreensFunctions(nsrc, nrec, nptsSig,
                                    fcent, dt, 
                                    lnorm, lshift,
                                    vel, rho,
                                    Q, pct,
                                    srcScale, xs, xr,
                                    &obs);
    CHKERR(ierr, "failed computing acoustic greens fns");
    //------------------------Compute the correlograms -----------------------//
    fprintf(stdout, "%s: Computing correlograms...\n", __func__); 
    int nsignals = nrec;
    int verbose = 0;
    int prec = XCLOC_SINGLE_PRECISION; //0;
    int accuracy = XCLOC_HIGH_ACCURACY; //0;
    bool ldoAutoCorrs = false;
    int nwork =-1; 
    int nxcs;
    int *xcPairs = NULL;
    xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                      XCLOC_FORTRAN_NUMBERING,
                                      &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable workspace query");
    nwork = 2*nxcs;
    xcPairs = calloc((size_t) nwork, sizeof(int));
    xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                      XCLOC_FORTRAN_NUMBERING,
                                      &nxcs, xcPairs, &ierr);
    CHKERR(ierr, "computeDefaultXCTable");
    // Initialize the cross-correlation table
    xcloc_fdxc_initialize(nptsSig, nptsSig, 
                          nxcs, xcPairs,
                          verbose, prec, accuracy, &ierr);
    CHKERR(ierr, "failed initializing fdxc");
    xcloc_fdxc_setSignals64f(nptsSig, nptsSig, nsignals, obs, &ierr);
    CHKERR(ierr, "failed to set signals");
    xcloc_fdxc_computePhaseCorrelograms(&ierr);
    CHKERR(ierr, "failed to compute phase correlograms");
    int nptsInXCs;
    xcloc_fdxc_getCorrelogramLength(&nptsInXCs, &ierr);
    CHKERR(ierr, "failed to get number of points in xcs");
    double *xcs = (double *) calloc((size_t) (nxcs*nptsInXCs), sizeof(double));
float *xcs32 = (float *) calloc((size_t) (nxcs*nptsInXCs), sizeof(float));
    xcloc_fdxc_getCorrelograms64f(nptsInXCs, nxcs, xcs, &ierr);
xcloc_fdxc_getCorrelograms32f(nptsInXCs, nxcs, xcs32, &ierr);
    //------------------------Filter the Correlograms-------------------------//
    int nTaps = 301;
    int ftype = XCLOC_SPXC_ENVELOPE_FILTER;
    xcloc_spxc_initialize(nTaps, ftype, &ierr);
double *xcsFilt = (double *) calloc((size_t) (nxcs*nptsInXCs), sizeof(double));
    xcloc_spxc_filterXCsOutOfPlace64f(nptsInXCs, nptsInXCs, nxcs,
                                      xcs, xcsFilt, &ierr);
float *xcsFilt32 = (float *) calloc((size_t) (nxcs*nptsInXCs), sizeof(float));
    xcloc_spxc_filterXCsOutOfPlace32f(nptsInXCs, nptsInXCs, nxcs,
                                      xcs32, xcsFilt32, &ierr);
FILE *fl = fopen("envelope.txt", "w");
for (int i=0; i<nptsInXCs; i++)
{
 int j = i;
 fprintf(fl, "%f %e %e %e\n", (i-nptsInXCs/2)*dt, xcs32[j], xcsFilt32[j], xcsFilt[j]);
}
fclose(fl);
    xcloc_spxc_finalize();
//  return 0;
    //---------------------------Compute the DSM------------------------------//
    fprintf(stdout, "%s: Initializing DSM...\n", __func__);
    int nxcPairs = nxcs;
    xcloc_dsmxc_initialize(ngrd, nxcPairs, nptsInXCs,
                           dt, xcPairs, verbose, &ierr);
    CHKERR(ierr, "failed to initialize dsm");
    fprintf(stdout, "%s: Setting travel time tables...\n", __func__);
    double *ttable = (double *) calloc((size_t) ngrd, sizeof(double));
    for (i=0; i<nrec; i++)
    {
        computeTravelTimeTable(nx, ny, nz, vel, x0, y0, z0, dx, dy, dz,
                               xr[3*i], xr[3*i+1], xr[3*i+2], ttable);
        xcloc_dsmxc_signalToTableIndex(i+1, &it, &ierr);
        xcloc_dsmxc_setTable64f(it, ngrd, ttable, &ierr);
    } 
    free(ttable);
    fprintf(stdout, "%s: Setting observations...\n", __func__);
    xcloc_dsmxc_setCorrelograms64f(nptsInXCs, nptsInXCs, nxcs, xcsFilt, &ierr);
    fprintf(stdout, "%s: Computing dsm...\n", __func__);
    xcloc_dsmxc_compute(&ierr);
    CHKERR(ierr, "failed to compute dsm");
    float *image = (float *) calloc((size_t) ngrd, sizeof(float));
    xcloc_dsmxc_getImage32f(ngrd, image, &ierr);
    CHKERR(ierr, "failed to get image");
printf("src1: %f %f\n", xs[0], xs[1]);
printf("src2: %f %f\n", xs[3], xs[4]);
FILE *ftemp = fopen("dsm2d.txt", "w");
for (int iy=0; iy<ny; iy++)
{
 for (int ix=0; ix<nx; ix++)
 {
  fprintf(ftemp, "%e %e %e\n", x0+ix*dx, y0+iy*dy, image[iy*nx+ix]);
 }
 fprintf(ftemp, "\n");
}
fclose(ftemp);
    free(image);
    // Free memory
    xcloc_fdxc_finalize();
    xcloc_dsmxc_finalize();
    free(xcPairs);
    free(xcs);
    if (xcs32 != NULL){free(xcs32);}
    free(obs);
    free(xr);
    free(xcsFilt);
    free(xcsFilt32);
    return EXIT_SUCCESS;
}

int computeTravelTimeTable(const int nx, const int ny, const int nz,
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

double computeTravelTimesToReceivers(const int nrec,
                                     const double vel,
                                     const double xsrc[],
                                     const double xrec[],
                                     double ttimes[])
{
    double dx, dy, dz, slow;
    int i;
    if (vel <= 0.0)
    {
        fprintf(stderr, "%s: velocity must be positive\n", __func__);
        return -1;
    }
    slow = 1.0/vel;
    #pragma omp simd
    for (i=0; i<nrec; i++)
    {
        dx = xrec[3*i+0] - xsrc[0];
        dy = xrec[3*i+1] - xsrc[1];
        dz = xrec[3*i+2] - xsrc[2];
        ttimes[i] = sqrt(dx*dx + dy*dy + dz*dz)*slow;
    }
    return 0; 
}

double gaussianCorrelogram(const int nptsInXC,
                           const double samplingPeriod,
                           const double deltaT,
                           const double t0,
                           const double sigma1,
                           const double sigma2,
                           double xc[])
{
    int i;
    double arg, res, res2, t, two_sigma12_p_sigma22, xden, xfact;
    two_sigma12_p_sigma22 = 2.0*(sigma1*sigma1 + sigma2*sigma2);
    xfact = 1.0/sqrt(M_PI*two_sigma12_p_sigma22);
    xden = 1.0/two_sigma12_p_sigma22;
    #pragma omp simd
    for (i=0; i<nptsInXC; i++)
    {
        t = t0 + (double) i*samplingPeriod;
        res = t - deltaT;
        res2 = res*res;
        arg =-res2*xden;
        xc[i] = xfact*exp(arg);
    }
    return 0;
}

int computeRandomReceiverLocations(const int nrec,
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

int compute2DGreensFunctions(const int nsrc, const int nrec, const int nptsSig,
                             const double fcent, const double dt,
                             const bool lnorm, const bool lshift,
                             const double vel, const double rho,
                             const double Q,
                             const double pct,
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
