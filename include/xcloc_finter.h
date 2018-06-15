#ifndef XCLOC_FINTER_H__
#define XCLOC_FINTER_H__
#include "xcloc_config.h"
#include "xcloc_enum.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

/*----------------------------------------------------------------------------*/
/*                            Convenience Utilities                           */
/*----------------------------------------------------------------------------*/
/* Make a full table of correlation pairs. */
void xcloc_utils_computeDefaultXCTable(const bool ldoAutoCorrs,
                                       const int nsignals,
                                       const int nwork,
                                       const int numbering,
                                       int *nxcs,
                                       int xcPairs[],
                                       int *ierr); 
/*----------------------------------------------------------------------------*/
/*                Diffraction Stack Migration of Correlograms                 */
/*----------------------------------------------------------------------------*/
/* Initialize the DSM. */
void xcloc_dsmxc_initialize(const int ntables, const int ngrd,
                            const int nxcPairs, const int nptsInXCs,
                            const double dt, const int xcPairs[],
                            int *ierr);
/* Computes the DSM image. */
void xcloc_dsmxc_compute(int *ierr);
/* Sets correlograms. */
void xcloc_dsmxc_setCorrelograms64f(const int ldxc,
                                    const int nptsInXCs,
                                    const int nxcPairs,
                                    const double xcs[], int *ierr);
void xcloc_dsmxc_setCorrelograms32f(const int ldxc,
                                    const int nptsInXCs,
                                    const int nxcPairs,
                                    const float xcs[], int *ierr);
/* Gets the migration image. */
void xcloc_dsmxc_getImage64f(const int nwork,
                             double image[],
                             int *ierr);
void xcloc_dsmxc_getImage32f(const int nwork, 
                             float image[],
                             int *ierr);
/* Sets a travel time table. */
void xcloc_dsmxc_setTable64fF(const int tableNumber,
                              const int ngrd,
                              const double table[],
                              int *ierr);
void xcloc_dsmxc_setTable32fF(const int tableNumber,
                              const int ngrd,
                              const float table[],
                              int *ierr);
/* Release the memory on the DSM structure. */
void xcloc_dsmxc_finalize(void);
/*----------------------------------------------------------------------------*/
/*                 Cross-Correlation in the Frequency Domain                  */
/*----------------------------------------------------------------------------*/
/* Initialize the frequency domain cross-correlation. */
void xcloc_fdxc_initialize(const int npts,
                           const int nsignals,
                           const int nptsPad,
                           const int nxcs,
                           const int xcPairs[],
                           const int verbose,
                           const int precision,
                           const int accuracy,
                           int *ierr);
/* Finalize the cross-correlation. */
void xcloc_fdxc_finalize(void);
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
/* Gets the number of points in the correlograms. */
void xcloc_fdxc_getCorrelogramLength(int *nptsInXCs, int *ierr);
/* Gets the number of time domain signals. */
void xcloc_fdxc_getNumberOfSignals(int *nsignals, int *ierr);
/* Gets all the cross-correlograms. */
void xcloc_fdxc_getCorrelograms64f(const int ldxc, const int nxcs,
                                   double xcs[], int *ierr);
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
