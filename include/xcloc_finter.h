#ifndef XCLOC_FINTER_H__
#define XCLOC_FINTER_H__
#include <stdbool.h>
#include "xcloc_config.h"
#include "xcloc_enum.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*----------------------------------------------------------------------------*/
/*                                    xcloc                                   */
/*----------------------------------------------------------------------------*/
/* Initializes xcloc */
void xcloc_initialize(const int npts, const int nptsPad, const int nxcs,
                      const enum xclocMigrateXC_enum s2m,
                      const double dt, const int ngrd,
                      const int nfcoeffs,
                      const enum xclocXCFilterType_enum ftype,
                      const int xcPairs[],
                      const enum xclocXCVerbose_enum verbose,
                      const enum xclocPrecision_enum precision,
                      const enum xclocAccuracy_enum accuracy,
                      int *ierr);
/* Finalizes the module. */
void xcloc_finalize(void);
/* Gets the travel time table corresponding to the signal index. */
void xcloc_signalToTableIndex(const int is, int *it, int *ierr);
/* Sets the travel time table. */
void xcloc_setTable64f(const int tableNumber, const int ngrd,
                       const double table[], int *ierr);
void xcloc_setTable32f(const int tableNumber, const int ngrd,
                       const float table[], int *ierr);
/* Sets the signals to correlate. */
void xcloc_setSignals64f(const int ldx, const int npts, const int nsignals,
                         const double x[], int *ierr);
void xcloc_setSignals32f(const int ldx, const int npts, const int nsignals,
                         const float x[], int *ierr);
/* Gets the migrated image. */
void xcloc_getImage64f(const int nwork, const double image[], int *ierr);
void xcloc_getImage32f(const int nwork, const float image[], int *ierr);
/* Convenience utility to get the max of the image. */
void xcloc_getImageMax(int *imageMax, float *maxValue, int *ierr);
/* Sets the filter to apply to the correlograms. */
void xcloc_setXCFilter(const int nfcoeffs,
                       const enum xclocXCFilterType_enum ftype,
                       int *ierr);
/* Sets the type of correlogram to migrate. */
void xcloc_setXCTypeToMigrate(const  enum xclocMigrateXC_enum,
                              int *ierr);
/* Performs cross-correlations and computes the xcloc image. */
void xcloc_compute(int *ierr);
/* Gets the number of grid points in the image. */
void xcloc_getNumberOfGridPointsInImage(const int ngrd, int *ierr);


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
void xcloc_dsmxc_initialize(const int ngrd,
                            const int nxcPairs, const int nptsInXCs,
                            const double dt, const int xcPairs[],
                            const int verbose, int *ierr);
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
/* Get maximum of image. */
void xcloc_dsmxc_getImageMax(int *maxIndex, float *maxValue, int *ierr);
/* Maps from the signal in the xcPairs table to the */
void xcloc_dsmxc_signalToTableIndex(const int is, int *it, int *ierr);
/* Sets a travel time table. */
void xcloc_dsmxc_setTable64f(const int tableNumber, /* Fortran indexed. */
                             const int ngrd,
                             const double table[],
                             int *ierr);
void xcloc_dsmxc_setTable32f(const int tableNumber, /* Fortran indexed. */
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

/*----------------------------------------------------------------------------*/
/*                     Signals Processing of Correlograms                     */
/*----------------------------------------------------------------------------*/
/* Initializes the correlogram filtering. */
void xcloc_spxc_initialize(const int n, const int ftype, int *ierr);
/* Filters signals out of place. */
void xcloc_spxc_filterXCsOutOfPlace64f(const int ldxc,
                                       const int nptsInXCs,
                                       const int nxcs,
                                       const double xcs[], double xcsFilt[],
                                       int *ierr);
void xcloc_spxc_filterXCsOutOfPlace32f(const int ldxc,
                                       const int nptsInXCs,
                                       const int nxcs,
                                       const float xcs[], float xcsFilt[],
                                       int *ierr);
/* Filters signals in place. */
void xcloc_spxc_filterXCsInPlace64f(const int ldxc,
                                    const int nptsInXCs,
                                    const int nxcs,
                                    double xcs[],
                                    int *ierr);
void xcloc_spxc_filterXCsInPlace32f(const int ldxc,
                                    const int nptsInXCs,
                                    const int nxcs,
                                    float xcs[],
                                    int *ierr);
/* Finalizes the correlogram filtering. */
void xcloc_spxc_finalize(void);

#ifdef XCLOC_USE_MPI
#include <mpi.h>
/*----------------------------------------------------------------------------*/
/*                              Parallel xcloc                                */
/*----------------------------------------------------------------------------*/
/* Initializes the parallel xcloc */
void xclocMPI_initialize(const MPI_Fint comm, //const MPI_Comm comm,
                         const int root,
                         const int dsmGroupSize,
                         const int npts,
                         const int nptsPad,
                         const int nxcs,
                         const enum xclocMigrateXC_enum s2m,
                         const double dt,
                         const int ngrd,
                         const int nfCoeffs,
                         const enum xclocXCFilterType_enum ftype,
                         const int xcPairs[],
                         const enum xclocXCVerbose_enum verbose,
                         const enum xclocPrecision_enum precision,
                         const enum xclocAccuracy_enum accuracy,
                         int *ierr);
/* Releases memory on the parallel xcloc module */
void xclocMPI_finalize(void);
                        
/*----------------------------------------------------------------------------*/
/*                 Parallel Frequency Domain Cross-Correlation                */
/*----------------------------------------------------------------------------*/
/* Initialize parallel frequency domain cross correlation. */
void xcloc_fdxcMPI_initialize(const MPI_Fint comm, //const MPI_Comm comm,
                              const int root, 
                              const int npts,
                              const int nptsPad,
                              const int nxcs,
                              const int xcPairs[],
                              const int verbose,
                              const int precision,
                              const int accuracy,
                              int *ierr);
/* Set many signals. */
void xcloc_fdxcMPI_setSignals64f(const int ldx,
                                 const int npts,
                                 const int nsignals,
                                 const int root,
                                 const double x[], int *ierr);
void xcloc_fdxcMPI_setSignals32f(const int ldx,
                                 const int npts,
                                 const int nsignals,
                                 const int root,
                                 const float x[], int *ierr);
/* Gathers the correlogram onto the root process. */
void xcloc_fdxcMPI_gatherCorrelograms64f(const int ldxc, const int nxcs,
                                         const int root,
                                         double xcs[], int *ierr);
void xcloc_fdxcMPI_gatherCorrelograms32f(const int ldxc, const int nxcs,
                                         const int root,
                                         float xcs[], int *ierr);
/* Compute the cross-correlograms. */
void xcloc_fdxcMPI_computeCrossCorrelograms(int *ierr);
/* Compute the phase-correlograms. */
void xcloc_fdxcMPI_computePhaseCorrelograms(int *ierr);
/* Finalize the module. */
void xcloc_fdxcMPI_finalize(void);
/*----------------------------------------------------------------------------*/
/*                      Parallel DSM of Correlograms                          */
/*----------------------------------------------------------------------------*/
/* Initialize parallel DSM. */
void xcloc_dsmxcMPI_initialize(const MPI_Fint comm, //const MPI_Comm comm,
                               const int root,
                               const int ngrd, const int nxcPairs,
                               const int nptsInXCs, 
                               const double dt,
                               const int xcPairs[],
                               const int verbose,
                               int *ierr);
/* Set a table. */
//void xcloc_dsmxcMPI_setTable64f(const int tableNumber, 
#endif /* MPI */

#ifdef __cplusplus
}
#endif
#endif
