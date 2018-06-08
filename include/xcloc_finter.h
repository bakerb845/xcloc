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
                           const int accuracy,
                           int *ierr);
/* Finalize the cross-correlation. */
void xcloc_fdxc_finalize(void);
/* Set the Fortran indexed XC table. */
void xcloc_fdxc_setXCTableF(const int nxcs, const int xcPairs[], int *ierr);
/* Computes a default cross-correlation XC table. */
void xcloc_fdxc_computeDefaultXCTableF(const bool ldoAutoCorrs,
                                       const int nwork,
                                       int *nxcs, int xcPairs[], int *ierr);
/* Sets many signals. */
void xcloc_fdxc_setSignals64f(const int ldx, const int npts, const int nsignals,
                              const double x[], int *ierr);
void xcloc_fdxc_setSignals32f(const int ldx, const int npts, const int nsignals,
                              const float x[], int *ierr);
/* Sets a float signal. */ 
void xcloc_fdxc_setSignal64fF(const int signalNumber,
                              const int npts,
                              const double x[],
                              int *ierr);
void xcloc_fdxc_setSignal32fF(const int signalNumber,
                              const int npts,
                              const float x[],
                              int *ierr);
/* Gets all the cross-correlograms. */
void xcloc_fdxc_getCorrelograms32f(const int ldxc, const int nxcs,
                                   float xcs[], int *ierr); 
/* Gets a cross-correlogram. */
void xcloc_fdxc_getCorrelogram32fF(const int corrNumber, const int lwork,
                                   float xc[], int *ierr);
/* Computes the phase correlograms. */
void xcloc_fdxc_computePhaseCorrelograms(int *ierr);
/* Computes the cross-correlograms. */
void xcloc_fdxc_computeCrossCorrelograms(int *ierr);


#ifdef __cplusplus
}
#endif
#endif
