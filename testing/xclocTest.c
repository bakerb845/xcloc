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
    // Scatter the receivers
    xr = (double *) calloc((size_t) (3*nrec), sizeof(double));
    ierr = acousticGreens2D_computeRandomReceiverLocations(nrec,
                                                           x0, y0, z0, 
                                                           x1, y1, z1, 
                                                           xr);
    CHKERR(ierr, "failed making receiver locations");
    // Compute the greens functions
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
printf("hey\n");

    if (xr != NULL){free(xr);}
    if (obs != NULL){free(obs);}
    if (xcPairs != NULL){free(xcPairs);}
    return EXIT_SUCCESS;
}
