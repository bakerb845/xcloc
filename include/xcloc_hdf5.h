#ifndef XCLOC_HDF5_H__
#define XCLOC_HDF5_H__ 1
#include <limits.h>
#include <hdf5.h>
#include "xcloc_xdmf.h"

struct xclocHDF5Grid_struct
{
    struct xclocXDMFGrid_struct xdmf; /*!< XDMF description of HDF5 file. */
    char outDir[PATH_MAX];    /*!< Output directory. */
    char groupName[PATH_MAX]; /*!< Name of group where datasets are added. */
    char h5flName[PATH_MAX];  /*!< Name of HDF5 file. */
    char xdmfFile[PATH_MAX];  /*!< Name of HDF5 file. */
    hid_t dataSpaceID;/*!< Data space descriptor for HDF5. */
    hid_t h5fl;       /*!< HDF5 file handle. */
    hid_t groupID;    /*!< Group where datasets will be written. */
    double x0;        /*!< x origin (meters). */
    double y0;        /*!< y origin (meters). */
    double z0;        /*!< z origin (meters). */
    double dx;        /*!< Grid spacing in x (meters). */
    double dy;        /*!< Grid spacing in y (meters). */
    double dz;        /*!< Grid spacing in z (meters). */
    int nx;           /*!< Number of x grid points. */
    int ny;           /*!< Number of y grid points. */
    int nz;           /*!< Number of z grid points. */
    bool linit;       /*!< Flag indicating the structure was initialized. */
};

#ifdef __cplusplus
extern "C"
{
#endif

/* opens an h5 grid archive */
int xcloc_h5ioGrid_open(
    const char *h5flName,
    const char *xdmfFile,
    const char *groupName,
    const int nx, const int ny, const int nz,
    const double dx, const double dy, const double dz,
    const double x0, const double y0, const double z0,
    struct xclocHDF5Grid_struct *h5);

/* writes a dataset to the archive */
int xcloc_h5ioGrid_writeDataSet32f(
    const char *dataSet,
    const int nx, const int ny, const int nz,
    const float *__restrict__ grid,
    struct xclocHDF5Grid_struct *h5);

/* open the hdf5 grid archive */
int xcloc_h5ioGrid_close(struct xclocHDF5Grid_struct *h5);

/* private function - don't call this */
int xcloc_h5ioGrid_free(struct xclocHDF5Grid_struct *h5);

#ifdef __cplusplus
}
#endif

#endif
