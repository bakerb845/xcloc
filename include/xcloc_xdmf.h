#ifndef INTERLOC_XDMF_H__
#define INTERLOC_XDMF_H__ 1
#include <stdbool.h>
#include "xcloc_config.h"
#include "xcloc_enum.h"

struct xclocXDMFGrid_struct
{
    enum xclocPrecision_enum *precision; /*!< Precision of datsets. */
    char **h5flNames; /*!< Number of HDF5 file names. */
    char **dataSets;  /*!< Names of datasets. */
    double x0;        /*!< x origin (meters). */
    double y0;        /*!< y origin (meters). */
    double z0;        /*!< z origin (meters). */
    double dx;        /*!< Grid spacing in x (meters). */
    double dy;        /*!< Grid spacing in y (meters). */
    double dz;        /*!< Grid spacing in z (meters). */
    int nx;           /*!< Number of x grid points. */
    int ny;           /*!< Number of y grid points. */
    int nz;           /*!< Number of z grid points. */
    size_t nDataSets; /*!< Number of models in collection. */
    size_t mDataSets; /*!< Max number of datasets. */
    bool lis2d;       /*!< Flag indicating the model is 2D. */
    bool linit;       /*!< Flag indicating whether or not the structure
                           is initialized. */
};

#ifdef __cplusplus
extern "C"
{
#endif

int xcloc_xdmfGrid_free(struct xclocXDMFGrid_struct *xdmf);
int xcloc_xdmfGrid_initialize(
    const int nx, const int ny, const int nz,
    const double dx, const double dy, const double dz,
    const double x0, const double y0, const double z0,
    struct xclocXDMFGrid_struct *xdmf);
int xcloc_xdmfGrid_add(const char *h5flName, const char *dataSet,
                       const enum xclocPrecision_enum precision,
                       struct xclocXDMFGrid_struct *xdmf);
int xcloc_xdmf_writeGrid(const char *xdmfFile,
                         const struct xclocXDMFGrid_struct xdmf);

#ifdef __cplusplus
}
#endif
#endif
