#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "xcloc.h"
#include "iscl/iscl/iscl.h"
#include "acoustic_greens.h"
#ifdef XCLOC_PROFILE
#include "advisor-annotate.h"
#endif
#include <ipps.h>

#define SEED 518815
#define N_P_SIGNALS 15
#define N_S_SIGNALS 15
#define NSIGNALS 30

int compute3DTravelTimes(const int nx, const int ny, const int nz,
                         const double x0, const double y0, const double z0,
                         const double dx, const double dy, const double dz,
                         const double vel, const double xs[3],
                         float *__restrict__ ttimes);
int create2DReceivers(const int seed, const int nrec,
                      const double x0, const double x1,
                      const double y0, const double y1,
                      const double z0,
                      double *__restrict__ xr);
int computeTravelTimes(const int nrec,
                       const double vel,
                       const double xs[3],
                       const double *__restrict__ xr,
                       double *__restrict__ ttimes);

int main(int argc, char *argv[])
{
    struct xclocParms_struct xclocParms;
    struct xcloc_struct xcloc;
    int i, ierr, ierrAll, it, myid, nprocs, provided;
    int envFIRLen = 251;  // Number of points in Hilbert FIR filter
    double fcent = 800.0; // Dominant frequency (Hz) in Ricker wavelet
    bool lnorm = false;   // Don't normalize Ricker wavelet
    bool lshift = true;   // Do shift Ricker wavelte
    double x0 = 0.0;    // Model origin in x is 0 km 
    double x1 = 1000.0; // Model extent in x is 1 km
    double y0 = 0.0;    // Model origin in y is 0 km
    double y1 = 1000.0; // Model extent in y is 1 km
    double z0 = 0.0;    // Problem is 2d - make z equal to 0.0
    double z1 = 0.0;    // Problem is 2d - make z equal to 0.0
    double dt = 1.0/6000.0;       // Sampling rate is 6000 Hz
    double vp = 3100.0;           // P velocity
    double vs = 3100.0/sqrt(3.0); // S velocity
    double rho = 2700.0;          // A nice crustal density
    double Qp = 600.0;            // P quality factor
    double Qs = 400.0;            // Q quality factor
    double tmodel = 0.5;          // Max modeling is the traveltime from the
                                  // furthest point in the medium to the
                                  // receiver plus some.
    int nx = 512;
    int ny = 512;
    int nz = 1;
    int ngrd = nx*ny*nz;
    double dx = (x1 - x0)/(double) (nx - 1); // should be less than 2m
    double dy = (y1 - y0)/(double) (ny - 1); // should be less than 2m
    double dz = 0.0;
    // Set the receiver location to the center of the model
    double xs[3] = {x0 + (double) (nx/2)*dx,
                    y0 + (double) (ny/2)*dy,
                    z0};
    int nptsSig = (int) (round(tmodel/dt)) + 1;
    double *obs = NULL;
    double *stf = NULL;
    double *ttimes = NULL;
    double *xr = NULL;
    float *tt = NULL;
    const int master = 0;
    double t0, tAll;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    iscl_init();
    // Initialize xcloc
    memset(&xclocParms, 0, sizeof(struct xclocParms_struct));
    memset(&xcloc, 0, sizeof(struct xcloc_struct));
    if (myid == master)
    {
        fprintf(stdout, "%s: Generating test...\n", __func__);
        xclocParms.nfftProcs = 1;
        xclocParms.ngridProcs = nprocs;
        if (nprocs == 4)
        {
            xclocParms.nfftProcs = 2;
            xclocParms.ngridProcs = 2;
        }
        if (nprocs == 6)
        {
            xclocParms.nfftProcs = 3;
            xclocParms.ngridProcs = 2;
        }
        if (nprocs == 8)
        {
            xclocParms.nfftProcs = 2;
            xclocParms.ngridProcs = 4;
        }
        xclocParms.dt = dt;
        xclocParms.lphaseXCs = true;
        xclocParms.chunkSize = 2048;
        xclocParms.envFIRLen = envFIRLen;
        xclocParms.npts    = nptsSig;
        xclocParms.nptsPad = nptsSig;
        xclocParms.ngrd    = ngrd;
        xclocParms.nsignals = NSIGNALS;
        xclocParms.signalGroup
            = (int *) calloc((size_t) xclocParms.nsignals, sizeof(int));
        for (i=0; i<xclocParms.nsignals; i++)
        {
            xclocParms.signalGroup[i] = 0;
            if (i + 1 > N_P_SIGNALS){xclocParms.signalGroup[i] = 1;}
        }
        // Compute the receiver locations 
        xr = (double *) calloc((size_t) (3*N_P_SIGNALS), sizeof(double));
        ttimes = (double *) calloc((size_t) NSIGNALS, sizeof(double));
        fprintf(stdout, "%s: Computing synthetics...\n", __func__);
        create2DReceivers(SEED, N_P_SIGNALS, x0, x1, y0, y1, z0, xr);
        // Compute the theoretical P travel-times
        computeTravelTimes(N_P_SIGNALS, vp, xs, xr, ttimes);
        // Compute the theoretical S travel-times
        computeTravelTimes(N_S_SIGNALS, vs, xs, xr,
                           &ttimes[N_P_SIGNALS]);
        stf = (double *) calloc((size_t) nptsSig, sizeof(double));
        obs  = (double *) calloc((size_t) (nptsSig*(NSIGNALS)), sizeof(double));
        ierr = dales_unitTests_computeRickerWavelet(nptsSig, dt, fcent,
                                                    lnorm, lshift, stf);
        // Compute theoretical P seismograms
        ierr = dales_unitTests_acousticGreensLineSource(N_P_SIGNALS, vp,
                                                        rho, Qp, nptsSig, dt,
                                                        xs, xr, stf,
                                                        obs);
        // Compute theoretical S seismograms
        ierr = dales_unitTests_acousticGreensLineSource(N_S_SIGNALS, vs,
                                                     rho, Qs, nptsSig, dt,
                                                     xs, xr, stf,
                                                     &obs[N_P_SIGNALS*nptsSig]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == master)
    {
        fprintf(stdout, "%s: Initializing xcloc...\n", __func__);
    }
    t0 = MPI_Wtime();
    ierr = xcloc_initialize(MPI_COMM_WORLD, xclocParms, &xcloc);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error initializing\n", __func__);
        ierr = 1;
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (ierrAll != 0){goto END;}
    // Set the travel time tables
    if (myid == master)
    {
        fprintf(stdout, "%s: Setting travel time tables...\n", __func__);
        tt = (float *) calloc((size_t) ngrd, sizeof(float));
    }
    for (it=0; it<xcloc.nTotalSignals; it++)
    {
        if (myid == master)
        {
            if (it < N_P_SIGNALS) 
            {
                 compute3DTravelTimes(nx, ny, nz,
                                      x0, y0, z0,
                                      dx, dy, dz,
                                      vp, xs,
                                      tt);
            }
            else
            {
                 compute3DTravelTimes(nx, ny, nz,
                                      x0, y0, z0,
                                      dx, dy, dz,
                                      vs, xs,
                                      tt);
            } 
        }
        ierr = xcloc_setTableFromRoot(it, ngrd, 
                                      XCLOC_SINGLE_PRECISION,
                                      tt, &xcloc);
        MPI_Barrier(MPI_COMM_WORLD);
    } 
    if (tt != NULL){free(tt); tt = NULL;}
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == master)
    {
        fprintf(stdout, "%s: Initializing time %5.2f (s)\n",
                __func__, MPI_Wtime() - t0);
    }
    // Loop so that my scaling test computes an average
    tAll = 0.0;
    for (it=0; it<1; it++)
    {
        // Scatter the data
        t0 = MPI_Wtime();
        ierr = xcloc_scatterDataFromRoot(NSIGNALS, nptsSig, nptsSig,
                                         MPI_DOUBLE, obs, &xcloc);
        if (ierr != 0)
        {
            fprintf(stdout, "%s: Error scattering data\n", __func__);
            ierr = 1;
        }
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (ierrAll != 0){goto END;}
        // Apply
        ierr = xcloc_apply(&xcloc);
        if (ierr != 0)
        {
            fprintf(stdout, "%s: Error applying migration\n", __func__);
            ierr = 1;
        }
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (ierrAll != 0){goto END;} 
        // Update the timing
        tAll = tAll + (MPI_Wtime() - t0);
    }
    MPI_Barrier(MPI_COMM_WORLD); 
    if (myid == master)
    {
        fprintf(stdout, "%s: Average executing time: %6.3f (s)\n",
                __func__, tAll/1.0);
    }
END:;
    // Finalize xcloc
    xcloc_finalize(&xcloc);
    if (xr != NULL){free(xr);}
    if (tt != NULL){free(tt);}
    if (ttimes != NULL){free(ttimes);}
    if (stf != NULL){free(stf);}
    if (obs != NULL){free(obs);}
    iscl_finalize();
    MPI_Finalize();
    return EXIT_SUCCESS;
}
//============================================================================//
int create2DReceivers(const int seed, const int nrec,
                      const double x0, const double x1,
                      const double y0, const double y1,
                      const double z0,
                      double *__restrict__ xr)
{
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
    }
    return 0;
}

int computeTravelTimes(const int nrec,
                       const double vel,
                       const double xs[3],
                       const double *__restrict__ xr,
                       double *__restrict__ ttimes)
{
    double dist;
    int irec;
    for (irec=0; irec<nrec; irec++)
    {
        dist = sqrt( pow(xs[0] - xr[3*irec+0], 2)
                   + pow(xs[1] - xr[3*irec+1], 2)
                   + pow(xs[2] - xr[3*irec+2], 2));
        ttimes[irec] = dist/vel;
    }
    return 0;
} 

int compute3DTravelTimes(const int nx, const int ny, const int nz,
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
