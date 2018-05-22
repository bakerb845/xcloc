#ifndef XCLOC_FINTER_H__
#define XCLOC_FINTER_H__
#include "xcloc_config.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

/* Initialize the frequency domain cross-correlation. */
void xcloc_fdxc_initialize(const int npts,
                           const int nsignals,
                           const int nptsPad,
                           const int verbose,
                           const int precision,
                           int *ierr);
/* Finalize the cross-correlation. */
void xcloc_fdxc_finalize(void);
/* Set the Fortran indexed XC table. */
void xcloc_fdxc_setXCTableF(const int nxcs, const int xcPairs[], int *ierr);
/* Computes a default cross-correlation XC table. */
void xcloc_fdxc_computeDefaultXCTableF(const bool ldoAutoCorrs,
                                       const int nwork,
                                       int *nxcs, int xcPairs[], int *ierr);
/* Sets a float signal. */ 
void xcloc_fdxc_setSignal32fF(const int signalNumber,
                              const int npts,
                              const float x[],
                              int *ierr);


#ifdef __cplusplus
}
#endif
#endif
