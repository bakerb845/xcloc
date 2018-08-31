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

#define NPTS 1201   /*!< 0.2 s window w/ rate of 6000 Hz */
#define NSIGNALS 45 /*!< make my computer sweat but don't tank the debugger */

int test_parallel_fdxc(const MPI_Comm comm, const int root)
{
    int npts = NPTS;
    int nptsPad = npts;
    int nsignals = NSIGNALS;
    int nptsInXC;
    int *xcPairs = NULL;
    int i, ierr, indx, ixc, myid, nxcs, nwork;
    int accuracy = (int) XCLOC_HIGH_ACCURACY;
    int verbose = (int) XCLOC_PRINT_INFO;
    int precision = (int) XCLOC_SINGLE_PRECISION;
    double diffTimeRef1, diffTimeRef2, res;
    double *phaseXCs = NULL;
    double *phaseXCsAll = NULL;
    double *xcsAll = NULL;
    double *xrand = NULL;
    ierr = EXIT_SUCCESS;
    MPI_Fint fcomm = MPI_Comm_c2f(comm); // why is this necessary openMPI?
    MPI_Comm_rank(comm, &myid);
    if (myid == root)
    {
        bool ldoAutoCorrs = false;
        nwork =-1; // space inquiry
        xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork, 
                                           XCLOC_FORTRAN_NUMBERING, &nxcs,
                                           xcPairs, &ierr);
        nwork = 2*nxcs;
        xcPairs = calloc((size_t) nwork, sizeof(int));
        xcloc_utils_computeDefaultXCTable(ldoAutoCorrs, nsignals, nwork,
                                          XCLOC_FORTRAN_NUMBERING, &nxcs,
                                          xcPairs, &ierr);
        // Compute some random signals
        fprintf(stdout, "%s: Generating random numbers...\n", __func__);
        xrand = xcfft_createRandomSignals(NULL, nsignals, npts, &ierr);
        fprintf(stdout, "%s: Initializing...\n", __func__);
    }
    else
    {
        xcPairs = (int *) calloc(1, sizeof(int));
    }
    xcloc_fdxcMPI_initialize(//comm, root,
                             fcomm, root,
                             npts, nptsPad,
                             nxcs, xcPairs,
                             verbose, precision, accuracy, &ierr);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Initialization failed\n", __func__);
        return EXIT_FAILURE;
    }
    nptsInXC = 2*nptsPad - 1;
    if (myid == root)
    {   
        phaseXCsAll = calloc((size_t) (nptsInXC*nxcs), sizeof(double));
        xcsAll = calloc((size_t) (nptsInXC*nxcs), sizeof(double));
    } 
#ifdef XCLOC_PROFILE
//ANNOTATE_SITE_BEGIN("fdxcMPI: random serial");
#endif
    // Set the signals
    double tbeg = MPI_Wtime();
    xcloc_fdxcMPI_setSignals64f(npts, npts, nsignals, root, xrand, &ierr);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error setting signals\n", __func__);
        return -1;
    }
    xcloc_fdxcMPI_computePhaseCorrelograms(&ierr);
    xcloc_fdxcMPI_gatherCorrelograms64f(nptsInXC, nxcs, root,
                                        phaseXCsAll, &ierr);
    double tend1 = MPI_Wtime() - tbeg;
    tbeg = MPI_Wtime();
    xcloc_fdxcMPI_setSignals64f(npts, npts, nsignals, root, xrand, &ierr);
    xcloc_fdxcMPI_computeCrossCorrelograms(&ierr);
    xcloc_fdxcMPI_gatherCorrelograms64f(nptsInXC, nxcs, root,
                                        xcsAll, &ierr);
    double tend2 = MPI_Wtime() - tbeg;
#ifdef XCLOC_PROFILE
//ANNOTATE_SITE_END("fdxcMPI: random serial");
#endif
    // Compare
    if (myid == root)
    {
        ierr = EXIT_SUCCESS;
        double resMax;
        int lxc = nptsInXC; 
        // Compute a reference answer 
        bool ldoPhase = true;
        double *xcref = (double *) calloc((size_t) (nxcs*lxc)+32, sizeof(double));
        tbeg = MPI_Wtime();
        ierr = xcfft_computeXCsWithISCL(ldoPhase, nsignals, nxcs,
                                        npts, lxc, xrand,
                                        xcref, &diffTimeRef1);
        diffTimeRef1 = MPI_Wtime() - tbeg;
        // Check phase xcs
        resMax = 0;
        for (ixc=0; ixc<nxcs; ixc++)
        {
            double *xc = &phaseXCsAll[ixc*nptsInXC];
            for (i=0; i<lxc; i++)
            {
                indx = ixc*lxc + i;
                res = fabs(xc[i] - xcref[indx]);
                if (res > 10.0*FLT_EPSILON)
                {
                    fprintf(stderr, "Failed on phase test: %d %d %f %e %e\n",
                            ixc, i, xc[i], xcref[indx], res);
                    //return EXIT_FAILURE;
                }
                resMax = fmax(resMax, res);
            }
        }
        if (resMax > 200.*FLT_EPSILON)
        {
            fprintf(stderr, "%s: Failed to compute phase XCs; %lf\n",
                      __func__, resMax);
            ierr = EXIT_FAILURE;
        }
        // Repeat for the regular correlograms
        ldoPhase = false;
        tbeg = MPI_Wtime();
        ierr = xcfft_computeXCsWithISCL(ldoPhase, nsignals, nxcs,
                                        npts, lxc, xrand,
                                        xcref, &diffTimeRef2);
        diffTimeRef2 = MPI_Wtime() - tbeg;
        // Check phase xcs
        resMax = 0;
        for (ixc=0; ixc<nxcs; ixc++)
        {
            double *xc = &xcsAll[ixc*nptsInXC];
            for (i=0; i<lxc; i++)
            {
                indx = ixc*lxc + i;
                res = fabs(xc[i] - xcref[indx]);
                if (res > 2.e-5) //100.0*FLT_EPSILON)
                {
                    fprintf(stderr, "Failed on phase test: %d %d %f %e %e\n",
                            ixc, i, xc[i], xcref[indx], res);
                    //return EXIT_FAILURE;
                }
                resMax = fmax(resMax, res);
            }
        }
        if (resMax > 200.*FLT_EPSILON)
        {
            fprintf(stderr, "%s: Failed to compute XCs\n; %lf", __func__, resMax);
            ierr = EXIT_FAILURE;
        }
        free(xcref);
        fprintf(stdout, "%s: Phase XC Reference time %lf; MPI time %lf\n",
                 __func__, diffTimeRef1, tend1);
        fprintf(stdout, "%s: XC Reference time %lf; MPI time %lf\n",
                 __func__, diffTimeRef2, tend2);
printf("%lf\n", resMax);
    }
    MPI_Bcast(&ierr, 1, MPI_INTEGER, root, comm); 
    // Clean up
    if (xcPairs != NULL){free(xcPairs);}
    xcloc_fdxcMPI_finalize();
    if (xrand != NULL){free(xrand);}
    if (phaseXCsAll != NULL){free(phaseXCsAll);}
    if (xcsAll != NULL){free(xcsAll);}
    if (phaseXCs != NULL){free(phaseXCs);}
    return ierr; //EXIT_SUCCESS;
}
