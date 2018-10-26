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

int test_serial_xcloc(void)
{
    const int nrec = 20; 
    int ntaps = 301; // Number of FIR filter taps
    double dt = 1.0/6000.0; // Sampling rate is 6000 Hz
    double fcent = 400.0;  // Dominant resolution is vmin/fcent ~ 5.0m (for plotting)
    bool lnorm = false;     // Don't normalize ricker wavelet (max will be 1)
    bool lshift = true;     // Make wavelet start at time 0 
    //const double pct = 8.0;  // 8 pct taper
    double x0 = 0.0;    // Model origin in x is 0 km 
    double x1 = 2600.0; // Model extent in x is 1 km
    double y0 = 0.0;    // Model origin in y is 0 km
    double y1 = 2800.0; // Model extent in y is 1 km
    double z0 = 0.0;    // Problem is 2d - make z equal to 0.0
    double z1 = 0.0;    // Problem is 2d - make z equal to 0.0
    double vel = 3100.0; // constant velocity (m/s)
    double rho = 2700.0; // constant density (kg/m**3)
    double Q = 9.e2;     // add some damping
    double tmodel = 1.5; // max modeling time is the traveltime from the
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
    int i, ierr, is, it;
    // Set the receiver location to the center of the model
    const int nsrc = 2;  
    const double srcScale[2] = {1, 1.1};
    double xs[6] = {x0 + (double) (1*nx/4)*dx,
                    y0 + (double) (5*ny/6)*dy,
                    z0,// Source has to be at zero for the 2d example
                    x0 + (double) (7*nx/9)*dx,
                    y0 + (double) (2*ny/9)*dy,
                    z0};
    double *xr = NULL;
    double *obs = NULL;
    float *image = NULL;
    bool lfound;
    // Scatter the receivers
    xr = (double *) calloc((size_t) (3*nrec), sizeof(double));
    ierr = acousticGreens2D_computeRandomReceiverLocations(nrec,
                                                           x0, y0, z0, 
                                                           x1, y1, z1, 
                                                           xr);
    CHKERR(ierr, "failed making receiver locations");
    // Compute the greens functions
    fprintf(stdout, "%s: Computing Green's functions...\n", __func__);
    ierr = acousticGreens2D_computeGreensFunctions(nsrc, nrec, nptsSig,
                                                   fcent, dt,
                                                   lnorm, lshift,
                                                   vel, rho,
                                                   Q,
                                                   srcScale, xs, xr,
                                                   &obs);
    CHKERR(ierr, "failed computing acoustic greens fns");
    //-------------------------- Initialize xcloc ----------------------------//
    fprintf(stdout, "%s: Creating correlation table...\n", __func__); 
    int nsignals = nrec;
    int verbose = XCLOC_PRINT_ERRORS;
    int prec = XCLOC_SINGLE_PRECISION; //0;
    int accuracy = XCLOC_HIGH_ACCURACY; //0;
    bool ldoAutoCorrs = false;
    int nwork =-1; 
    int maxIndex, nxcs;
    int *xcPairs = NULL;
    float maxValue;
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
    xcloc_initialize(nptsSig, nptsSig, nxcs,
                     XCLOC_MIGRATE_PHASE_XCS, dt, ngrd,
                     ntaps, XCLOC_SPXC_ENVELOPE_FILTER,
                     xcPairs,
                     verbose, prec, accuracy, &ierr);
    CHKERR(ierr, "failed to initialize xcloc");
    // Set the travel time tables
    fprintf(stdout, "%s: Computing and setting tables...\n", __func__);
    double *ttable = (double *) calloc((size_t) ngrd, sizeof(double));
    for (i=0; i<nrec; i++)
    {
        acousticGreens2D_computeTravelTimeTable(nx, ny, nz,
                                                vel, x0, y0, z0, dx, dy, dz,
                                                xr[3*i], xr[3*i+1], xr[3*i+2],
                                                ttable);
        xcloc_signalToTableIndex(i+1, &it, &ierr);
        CHKERR(ierr, "error getting table index");
        xcloc_setTable64f(it, ngrd, ttable, &ierr);
        CHKERR(ierr, "error settint table");
    }
    free(ttable);
    // Set the signals
    fprintf(stdout, "%s: Setting signals...\n", __func__);
    xcloc_setSignals64f(nptsSig, nptsSig, nsignals, obs, &ierr);
    CHKERR(ierr, "failed to set signals");
    // Do the heavy lifting
    fprintf(stdout, "%s: Computing...\n", __func__);
    xcloc_compute(&ierr);
    // Get the max of image and th image itself
    xcloc_getImageMax(&maxIndex, &maxValue, &ierr);
    CHKERR(ierr, "failed to get image max");
    maxIndex = maxIndex - 1; // Fortran to C
    image = (float *) calloc((size_t) ngrd, sizeof(float));
    xcloc_getImage32f(ngrd, image, &ierr);
    CHKERR(ierr, "failed to get image");
    // Search
    lfound = false;
    for (is=0; is<nsrc; is++)
    {
        int ixs = (int) ((xs[3*is+0] - x0)/dx + 0.5);
        int iys = (int) ((xs[3*is+1] - y0)/dy + 0.5);
        if (iys*nx + ixs == maxIndex){lfound = true;}
        fprintf(stdout, "%s: (maxIndex,trueIndex)=(%d,%d) has value %f\n",
                __func__, maxIndex, iys*nx + ixs, image[iys*nx+ixs]);
    }
    if (!lfound)
    {
        fprintf(stderr, "%s: Failed to find an event\n", __func__);
        return EXIT_FAILURE;
    }
    // Free xcloc
    xcloc_finalize();
    if (xr != NULL){free(xr);}
    if (obs != NULL){free(obs);}
    if (xcPairs != NULL){free(xcPairs);}
    if (image != NULL){free(image);}
    return EXIT_SUCCESS;
}
