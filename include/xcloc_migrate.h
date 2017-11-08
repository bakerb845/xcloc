#ifndef XCLOC_MIGRATE__
#define XCLOC_MIGRATE__ 1
#include "xcloc_config.h"
#include "xcloc_enum.h"

struct migrate_struct
{
    int *ttimesInt; /*!< Travel-times table. This is computed by dividing 
                         the travel-times by the sampling period.  This
                         is an array of dimension [mgrd x ntables] where
                         mgrd is the leading dimension. */
    int *xcPairs;   /*!< Maps from the ixc'th cross-correlation to the 
                         travel-time indices.  This is an array of dimension
                         [2 x nxc] with leading dimension 2. */
    float *migrate; /*!< The migrated image.  This is an array of dimension
                         [ngrd]. */
    float *xcs;     /*!< These are the cross-correlograms.  This is an 
                         array of dimension [ldxc x nxc] with leading
                         dimension ldxc. */
    double dt;      /*!< Sampling period (seconds) of signals to migrate. */
    int ntables;    /*!< Number of travel-time tables. */
    int nxc;        /*!< Number of cross-correlation pairs to migrate. */
    int lxc;        /*!< Length of the cross-correlations. */
    int ldxc;       /*!< Leading dimension of the cross-correlations. */
    int ngrd;       /*!< Number of grid points. */
    int mgrd;       /*!< Leading dimension of ttimes. */
    int chunkSize;  /*!< Chunksize in migration loop.  This is a tuning
                         parameter. */
};

#ifdef __cplusplus
extern "C"
{
#endif

/* Set the chunkSize */
int xcloc_migrate_setChunkSize(const int chunkSize,
                               struct migrate_struct *migrate);
/* Initializes the migration structure */
int xcloc_migrate_initialize(const int ntables, const int ngrd,
                             const int nxc, const int lxc,
                             const int chunkSize,
                             const double dt,
                             const int *__restrict__ xcPairs,
                             struct migrate_struct *migrate);
/* Deallocates/frees the migration structure */
int xcloc_migrate_finalize(struct migrate_struct *migrate);
/* Sets the it'th travel time table */
int xcloc_migrate_setTable64f(const int it, const int ngrd,
                              const double *__restrict__ ttimes, 
                              struct migrate_struct *migrate);
int xcloc_migrate_setTable32f(const int it, const int ngrd,
                              const float *__restrict__ ttimes,
                              struct migrate_struct *migrate);
/* Sets the ixc'th cross-correlation. */
int xcloc_migrate_setCrossCorrelations(const int ldxc, const int lxc,
                                       const int nxc, 
                                       const void *__restrict__ xcs,
                                       const enum xclocPrecision_enum precision,
                                       struct migrate_struct *migrate);
/* Sets the cross-correlations. */
int xcloc_migrate_setTable64f(const int it, const int ngrd,
                              const double *__restrict__ ttimes,
                              struct migrate_struct *migrate);
/* Sets the diffraction stack migration image to 0 */
int xcloc_migrate_setImageToZero(struct migrate_struct *migrate);
/* Convenience function to compute the stack of migration images */
int xcloc_migrate_computeMigrationImage(struct migrate_struct *migrate);
int xcloc_migrate_computeXCDSMImage(struct migrate_struct *migrate);
/* Stacks the cross-correlation signal into the migration image */
int xcloc_migrate_updateXCDSMImage(const int it1, const int it2,
                                   const int ixc,
                                   struct migrate_struct *migrate);
/* Gets the image */
int xcloc_migrate_getImage64f(const int ngrd,
                              const struct migrate_struct migrate,
                              double *__restrict__ image);
int xcloc_migrate_getImage32f(const int ngrd,
                              const struct migrate_struct migrate,
                              float *__restrict__ image);
const float *xcloc_migrate_getImagePointer(const int ngrd,
                                           const struct migrate_struct migrate,
                                           int *ierr);


int xcloc_migrate_migrateDifferentialTimesOnGrid(
    const int ld1, const int ld2,
    const int n1, const int n2, const int n3,
    const int lxc, const float dt,
    const float *__restrict__ xc,
    const float *__restrict__ ttimes1,
    const float *__restrict__ ttimes2,
    float *__restrict__ migrate);
int xcloc_migrate_migrateSignalsOnGrid(
    const int ld1, const int ld2,
    const int n1, const int n2, const int n3,
    const int npts, const float dt,
    const float *__restrict__ signal,
    const float *__restrict__ ttimes,
    float *__restrict__ migrate);

#ifdef __cplusplus
}
#endif
#endif
