#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "xcloc_finter.h"
#include "test_suite.h"
#include "mkl.h"
//#include "iscl/iscl/iscl.h"

#define SEED 40935  /*!< random seed - but predictable results. */

int main(int argc, char *argv[])
{
    int ierr, myid, provided;
    ierr = EXIT_SUCCESS;
    const int root = 0;
#ifdef XCLOC_USE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    //MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &myid);
#else
    myid = root;
#endif 
    ierr = EXIT_SUCCESS;
    // Initialize ISCL to make reference solution
    srand(SEED);
    if (myid == root)
    {
        fprintf(stdout, "%s: Performing serial tests...\n", __func__);
        ierr = test_serial_fdxc();
        if (ierr != EXIT_SUCCESS)
        {
            fprintf(stderr, "%s: Failed serial tests!\n", __func__);
            goto BCAST_ERROR_1;
        }
        fprintf(stdout, "%s: Performing serial DSM tests...\n", __func__);
        ierr = test_serial_dsmLocation();
        if (ierr != EXIT_SUCCESS)
        {
            fprintf(stderr, "%s: Failed serial dsm test!\n", __func__);
            goto BCAST_ERROR_1;
        }
    }
BCAST_ERROR_1:;
#ifdef XCLOC_USE_MPI
    MPI_Bcast(&ierr, 1, MPI_INTEGER, root, MPI_COMM_WORLD);
#endif
    if (ierr !=0 ){goto END;}
#ifdef XCLOC_USE_MPI
    // Start the parallel tests
    ierr = test_parallel_fdxc(comm, root);
#endif
END:;
    mkl_finalize();
#ifdef XCLOC_USE_MPI
    MPI_Finalize();
#endif
    return ierr;
}
