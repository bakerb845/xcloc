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

#define NPTS 1201   /*!< 0.2 s window w/ rate of 6000 Hz */
#define NSIGNALS 45 /*!< make my computer sweat but don't tank the debugger */

int test_parallel_xcloc(const MPI_Comm comm, const int root)
{
    const int nrec = 20; //5; //20;
    int nsignals = nrec;
    int nfcoeffs = 301; // Number of FIR filter taps
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
    int npts = (int) (round(tmodel/dt)) + 1;
    // rest of geometry
    int nx = 512;
    int ny = 512;
    int nz = 1;
    int ngrd = nx*ny*nz;
    double dx = (x1 - x0)/(double) (nx - 1); // should be less than 2m
    double dy = (y1 - y0)/(double) (ny - 1); // should be less than 2m
    double dz = 0.0;
    int *xcPairs = NULL;
    int ierr, ierrLoc, myid, nprocs, nxcs, nwork;
    int accuracy  = (int) XCLOC_HIGH_ACCURACY;
    int verbose   = (int) XCLOC_PRINT_INFO;
    int precision = (int) XCLOC_SINGLE_PRECISION;
    int ftype     = (int) XCLOC_SPXC_ENVELOPE_FILTER;
    int s2m       = (int) XCLOC_MIGRATE_PHASE_XCS;
    // Set the receiver location to the center of the model
    const int nsrc = 1; // TODO change here
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
    float *imageRef = NULL;
    bool lfound;
    int nxcdsmGroups = 1;
    int i, it;
    ierr = EXIT_SUCCESS;
    MPI_Fint fcomm = MPI_Comm_c2f(comm); // why is this necessary openMPI?
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &nprocs);
    srand(4096);
    if (myid == root)
    {
        // Define the number of xc/dsm groups
        nxcdsmGroups = 1;
        if (nprocs == 4)
        {
            nxcdsmGroups = 2;
        }
        else if (nprocs == 2)
        {
            nxcdsmGroups = 1;
        }
        else if (nprocs == 6)
        {
            nxcdsmGroups = 3;
        }
        else if (nprocs == 8)
        {
            nxcdsmGroups = 4;
        }
        // Scatter the receivers
        xr = (double *) calloc((size_t) (3*nrec), sizeof(double));
        ierr = acousticGreens2D_computeRandomReceiverLocations(nrec,
                                                               x0, y0, z0,
                                                               x1, y1, z1,
                                                               xr);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to compute recv locs", __func__);
            goto ERROR1;
        }
        // Compute the greens functions
        fprintf(stdout, "%s: Computing Green's functions...\n", __func__);
        ierr = acousticGreens2D_computeGreensFunctions(nsrc, nrec, npts,
                                                       fcent, dt,
                                                       lnorm, lshift,
                                                       vel, rho,
                                                       Q,
                                                       srcScale, xs, xr,
                                                       &obs);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to compute Green's funcctions\n",
                    __func__);
            goto ERROR1;
        }
        // Define the cross-correlation pairs
        bool ldoAutoCorrs = false;
        nwork =-1; // space inquiry
        xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork, 
                                           XCLOC_FORTRAN_NUMBERING, &nxcs,
                                           xcPairs, &ierr);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to compute XC table workspace",
                    __func__);
            goto ERROR1;
        }
        nwork = 2*nxcs;
        xcPairs = calloc((size_t) nwork, sizeof(int));
        xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                          XCLOC_FORTRAN_NUMBERING, &nxcs,
                                          xcPairs, &ierr);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to compute XC table", __func__);
            goto ERROR1;
        }
    }
ERROR1:;
    MPI_Bcast(&ierr, 1, MPI_INTEGER, root, comm);
    if (ierr != EXIT_SUCCESS){return EXIT_FAILURE;}
    if (myid == root)
    {
        fprintf(stdout, "%s: Initializing xclocMPI...\n", __func__);
    }
    // Initialize
    int nptsPad = npts;
    xclocMPI_initialize((int64_t) fcomm, root, 
                        nxcdsmGroups,
                        npts, nptsPad, nxcs,
                        s2m, dt, ngrd,
                        nfcoeffs, ftype,
                        xcPairs,
                        verbose, precision, accuracy,
                        &ierrLoc);
    if (ierrLoc != 0){ierrLoc = 1;}
    MPI_Allreduce(&ierrLoc, &ierr, 1, MPI_INTEGER, MPI_MAX, comm);
    if (ierr != 0)
    {
        if (myid == root)
        {
            fprintf(stderr, "%s: Failed to initialize xclocMPI\n", __func__);
        }
        return EXIT_FAILURE;
    }
    // Set the tables
    double *ttable = NULL;
    if (myid == root)
    {
        fprintf(stdout, "%s: Making and setting tables...\n", __func__);
        ttable = (double *) calloc((size_t) ngrd, sizeof(double));
    }
    for (i=0; i<nrec; i++)
    {
        if (myid == root)
        {
            acousticGreens2D_computeTravelTimeTable(nx, ny, nz,
                                                    vel, x0, y0, z0, dx, dy, dz,
                                                    xr[3*i],xr[3*i+1],xr[3*i+2],
                                                    ttable);
        }
        xclocMPI_signalToTableIndex(i+1, root, &it, &ierr);
        if (ierr != 0 || i + 1 != it)
        {
            fprintf(stderr, "%s: Error getting table index\n", __func__);
            return EXIT_FAILURE;
        }
        xclocMPI_setTable64f(it, ngrd, root, ttable, &ierr);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error setting table %d on process %d\n",
                    __func__, i+1, myid);
            return EXIT_FAILURE;
        }
    }
    MPI_Barrier(comm);
    // Set the signals
    if (myid == root){fprintf(stdout, "%s: Setting signals..\n", __func__);}
    xclocMPI_setSignals64f(npts, npts, nsignals, root, obs, &ierr);
    if (ierr != 0){return EXIT_FAILURE;}
    // Do the heavy lifting
    if (myid == root){fprintf(stdout, "%s: Computing...\n", __func__);}
    xclocMPI_compute(&ierr);
    if (ierr != 0){return EXIT_FAILURE;}
    // Get the results - verify correlogram lengths
    int nptsInXCs, nxcsTotal;
    xclocMPI_getCorrelogramLength(&nptsInXCs, &ierr);
    if (ierr != 0){return EXIT_FAILURE;}
    xclocMPI_getNumberOfCorrelograms(&nxcsTotal, &ierr);
    if (ierr != 0){return EXIT_FAILURE;}
    if (nptsInXCs != 2*npts - 1)
    {
        fprintf(stderr, "%s: Incorrect number of pts in xcs\n", __func__);
        return EXIT_FAILURE;
    }
    if (nxcsTotal != (nrec*(nrec-1))/2)
    {
        fprintf(stderr, "%s: Incorrect number of xcs\n", __func__);
        return EXIT_FAILURE;
    }
    // Get the results - the xcs
    double *xcs = NULL;
    if (myid == root)
    {
        xcs = (double *) calloc((size_t) (nptsInXCs*nxcsTotal), sizeof(double));
        image = (float *) calloc((size_t) ngrd, sizeof(float));
        imageRef = (float *) calloc((size_t) ngrd, sizeof(float));
    }
    xclocMPI_getCorrelograms64f(nptsInXCs, nxcsTotal, root, xcs, &ierr);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to get correlograms\n", __func__);
        return EXIT_FAILURE;
    }
    int maxIndex;
    float maxValue;
    xclocMPI_getImage32f(ngrd, root, image, &ierr); 
    xclocMPI_getImageMax(&maxIndex, &maxValue, &ierr);
    if (ierr != 0){return EXIT_FAILURE;}
    if (myid == root)
    {
        fprintf(stdout, "Est max: %d %.8f\n", maxIndex, maxValue);
    }
    // Free space
    xclocMPI_finalize();
    // Do it again
    if (myid == root)
    {
        fprintf(stdout, "%s: Generating serial reference xcloc...\n", __func__);
        double *xcsRef = (double *)
                         calloc((size_t) (nptsInXCs*nxcsTotal), sizeof(double));
        // Create a reference
        xcloc_initialize(npts, npts, nxcs,
                         s2m, dt, ngrd,
                         nfcoeffs, ftype,
                         xcPairs,
                         verbose, precision, accuracy, &ierr);
        for (i=0; i<nrec; i++)
        {
            acousticGreens2D_computeTravelTimeTable(nx, ny, nz,
                                                    vel, x0, y0, z0, dx, dy, dz,
                                                    xr[3*i],xr[3*i+1],xr[3*i+2],
                                                    ttable);
            xcloc_signalToTableIndex(i+1, &it, &ierr);
            xcloc_setTable64f(it, ngrd, ttable, &ierr);
        }
        xcloc_setSignals64f(npts, npts, nsignals, obs, &ierr);
        xcloc_compute(&ierr); if (ierr != 0){goto EXIT;}
        xcloc_getCorrelograms64f(nptsInXCs, nxcsTotal, xcsRef, &ierr);
        xcloc_getImageMax(&maxIndex, &maxValue, &ierr);
        xcloc_getImage32f(ngrd, imageRef, &ierr);
        printf("True max %d %.8f\n", maxIndex, maxValue);
        double resMax = 0;
        ierr = 0;
        ippsNormDiff_L1_64f(xcs, xcsRef, nptsInXCs*nxcsTotal, &resMax);
        if (resMax > 1.e-10)
        {
            fprintf(stderr, "%s: Correlograms inconsistent with error %lf\n",
                    __func__, resMax);
            ierr = 1;
        }
        xcloc_finalize();
        free(xcsRef);
        FILE *fwork = fopen("image_mpi.txt", "w");
        for (int iy=0; iy<ny; iy++)
        {
            for (int ix=0; ix<nx; ix++)
            {
                fprintf(fwork, "%lf %lf %f %f %f\n",
                        x0 + (double) ix*dx,
                        y0 + (double) iy*dy,
                        image[iy*nx+ix], imageRef[iy*nx+ix],
                        fabsf(image[iy*nx+ix] - imageRef[iy*nx+ix]));
            }
            fprintf(fwork, "\n");
        }
        fclose(fwork);
    }
EXIT:;
    MPI_Bcast(&ierr, 1, MPI_INTEGER, root, comm);
    if (xcPairs != NULL){free(xcPairs);}
    if (xr != NULL){free(xr);}
    if (image != NULL){free(image);}
    if (imageRef != NULL){free(imageRef);}
    if (obs != NULL){free(obs);}
    if (xcs != NULL){free(xcs);}
    if (ttable != NULL){free(ttable);}
    if (ierr != 0){return EXIT_FAILURE;}
    return EXIT_SUCCESS;
}
