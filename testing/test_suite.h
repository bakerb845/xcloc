#ifndef TEST_SUITE
#define TEST_SUITE 1

#ifdef __cplusplus
extern "C"
{
#endif


int test_serial_fdxc(void);
int xcfft_computeXCsWithISCL(const int nsignals, const int ntfSignals,
                             const int npts, const int lxc,
                             const double x[],
                             double xcs[],
                             double *diffTime);
double *xcfft_createRandomSignals(int *seed, const int nsignals, const int npts,
                                  int *ierr);

#ifdef __cplusplus
}
#endif
#endif
