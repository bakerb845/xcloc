#ifndef TEST_SUITE
#define TEST_SUITE 1
#include "xcloc_finter.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif


int test_serial_fdxc(void);
int test_serial_dsmLocation(void);
int xcfft_computeXCsWithISCL(const bool ldoPhase,
                             const int nsignals, const int ntfSignals,
                             const int npts, const int lxc,
                             const double x[],
                             double xcs[],
                             double *diffTime);
double *xcfft_createRandomSignals(int *seed, const int nsignals, const int npts,
                                  int *ierr);

#ifdef XCLOC_USE_MPI
int test_parallel_fdxc(const MPI_Comm comm, const int root);
#endif

#ifdef __cplusplus
}
#endif
#endif
