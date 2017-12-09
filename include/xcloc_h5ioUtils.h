#ifndef XCLOC_H5IOUTILS_H__
#define XCLOC_H5IOUTILS_H__ 1
#include "xcloc_config.h"
#include <hdf5.h>
#include "sacio.h"

enum h5Model_enum
{
    XCLOC_MODEL_CARTESIAN_GRIDDED = 0,
    XCLOC_MODEL_SPHERICAL_GRIDDED = 1,
    XCLOC_MODEL_CARTESIAN_UNSTRUCTURED = 2,
    XCLOC_MODEL_SPHERICAL_UNSTRUTCURED = 3
};

struct h5TravelTimeModel_struct
{
    float *ttimes;     /*!< Travel-times (seconds) from receiver to all
                            points in model.  It is an array of dimension
                            [ngrd]. */
    float *locs;       /*!< For unstructured models these are the node
                            locations in (x, y, z) (meters) or 
                            (longitude (deg), latitude (deg), radius (m)).
                            This is an array of dimension [3 x ngrd] with
                            leading dimension 3. */
    int *ieng;         /*!< For unstructured models this the node
                            connectivity map.  This is an array of dimension
                            [iengPtr[nelem]]. */
    int *iengPtr;      /*!< Maps from elem'th element ot start index of
                            iengv.  This is an array of dimension 
                            [nelem+1]. */
    float grid0[3];    /*!< Model origin {z0 (m), y0 (m), x0 (m)} or 
                            {r0 (m),  theta0 (deg), phi0} of
                            gridded model. */
    float dGrid[3];    /*!< Grid spacing {dz (m), dy (m), dx (m)} or
                            {dr (m), dtheta (deg), dphi (deg)} of
                             gridded model. */
    int nGrid[3];      /*!< {nz, ny, nx} or {nr, ntheta, nphi} of
                            gridded model. */
    int ngrd;          /*!< Number of grid points in model. */
    int nelem;         /*!< For unstructured models this is the number
                            of elements in the model. */
    int signalNumber;  /*!< This corresponds to the signal number. */
    enum h5Model_enum 
           modelType;  /*!< Defines the model type. */
};

#ifdef __cplusplus
extern "C"
{
#endif

int xcloc_h5ioUtils_data_initialize(const char *archiveName,
                                    const int nwin,
                                    const int nsgroups,
                                    hid_t *h5fl);
int xcloc_h5ioUtils_data_writeWindow(const int iwin, const int ig, 
                                     const int *compress,
                                     const size_t *traceChunkSz,
                                     const size_t *dataChunkSz,
                                     const hid_t fileID,
                                     const int ntrace,
                                     const struct sacData_struct *sac);
int xcloc_h5ioUtils_data_finalize(const hid_t fileID);


int xcloc_h5ioUtils_data_readWindowAndGroup(
    const hid_t fileID,
    const int iwin, const int ig,
    int *npts, double *dt, int *ntrace,
    double **data);
int xcloc_h5ioUtils_data_getNumberOfWindowsAndGroups(
    const hid_t fileID,
    int *nwin, int *nsgroups);

/*----------------------------------------------------------------------------*/
int xcloc_h5ioUtils_initTravelTimeArchive(
    const char *archiveName,
    const int nmodels,
    const struct h5TravelTimeModel_struct ttimes,
    hid_t *h5fl);
int xcloc_h5ioUtils_closeTravelTimeArchive(const hid_t h5fl);
int xcloc_h5ioUtils_writeTravelTimeModel(
    const hid_t fileID,
    const struct h5TravelTimeModel_struct ttimes);

int xcloc_h5ioUtils_openTravelTimeModelForReading(
    const char *archiveName,
    struct h5TravelTimeModel_struct *ttimes,
    hid_t *fileID);
int xcloc_h5ioUtils_readTravelTimeModel(
    const hid_t fileID,
    struct h5TravelTimeModel_struct *ttimes);


#ifdef __cplusplus
}
#endif
#endif
