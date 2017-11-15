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
    int i, ierr, myid, nprocs, provided;
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
    const int master = 0;
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
        xclocParms.npts    = nptsSig;
        xclocParms.nptsPad = nptsSig;
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
printf("%d %d\n", xclocParms.nfftProcs, xclocParms.ngridProcs);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ierr = xcloc_initialize(MPI_COMM_WORLD, xclocParms, &xcloc);
    // Finalize xcloc
    xcloc_finalize(&xcloc);
    if (xr != NULL){free(xr);}
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

