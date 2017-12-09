#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <libgen.h>
#include "xcloc_hdf5.h"
#include "iscl/os/os.h"

static char *strdup(const char *string)
{
    char *copy;
    size_t lenos;
    lenos = strlen(string);
    copy = (char *) calloc(lenos+1, sizeof(char));
    strcpy(copy, string);
    return copy;
}

/*!
 * @brief Deallocates memory on the HDF5 structure.
 *
 * @param[in,out] h5    On input contains an intialized h5io structure. \n
 *                      On exit all memory has been released and variables
 *                      set to 0.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_h5ioGrid_free(struct xclocHDF5Grid_struct *h5)
{
    int ierr;
    if (!h5->linit){return 0;}
    H5Gclose(h5->groupID);
    H5Sclose(h5->dataSpaceID);
    H5Fclose(h5->h5fl);
    ierr = xcloc_xdmfGrid_free(&h5->xdmf);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error freeing xdmfGrid\n", __func__);
    }
    memset(h5, 0, sizeof(struct xclocHDF5Grid_struct));
    return ierr;
}
//============================================================================//
/*!
 * @brief Closes the HDF5 file and writes the corresponding XDMF file.
 *        Because the H5 file has been closed this will also finalize the 
 *        data structures.
 *
 * @param[in,out] h5    On input contains the HDF5 structure. \n
 *                      On exit the corresponding XDMF file has been written,
 *                      the HDF5 file has been closed, all memory has been
 *                      released, and all variables set to 0.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int xcloc_h5ioGrid_close(struct xclocHDF5Grid_struct *h5)
{
    int ierr;
    if (!h5->linit){return 0;}
    ierr = xcloc_xdmf_writeGrid(h5->xdmfFile, h5->xdmf);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error writing XDMF grid\n", __func__);
    }
    ierr = xcloc_h5ioGrid_free(h5);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error freeing H5\n", __func__);
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Writes a floating dataset to the HDF5 file.
 *
 * @param[in] dataSet   Name of dataset.  This should not include a path.
 * @param[in] nx        Number of grid points in x.  This is the fastest
 *                      changing dimension in the model.
 * @param[in] ny        Number of grid points in y.
 * @param[in] nz        Number of grid points in z.  This is the slowest
 *                      changing dimension in the model.
 * @param[in] grid      Values to write.  This is an array of dimension
 *                      [nx x ny x nz] stored in column major format.
 *
 * @param[in,out] h5    On input this contains the initialized H5 structure. \n
 *                      On exit the dataSet has been added to the XDMF struct.
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distributed under MIT.
 *
 */
int xcloc_h5ioGrid_writeDataSet32f(
    const char *dataSet,
    const int nx, const int ny, const int nz,
    const float *__restrict__ grid,
    struct xclocHDF5Grid_struct *h5)
{
    char path[PATH_MAX];
    hid_t dataSetID;
    int ierr;
    herr_t status;
    if (!h5->linit)
    {
        fprintf(stderr, "%s: h5 structure not initialized\n", __func__);
        return -1;
    }
    if (h5->nx != nx || h5->ny != ny || h5->nz != nz)
    {
        fprintf(stderr, "%s: Inconsistent sizes\n", __func__);
        return -1;
    }
    if (dataSet == NULL)
    {
        fprintf(stderr, "%s: dataSet is NULL\n", __func__);
        return -1;
    } 
    // Create the dataset
    dataSetID = H5Dcreate2(h5->groupID, dataSet, H5T_NATIVE_FLOAT,
                           h5->dataSpaceID,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Write it
    status = H5Dwrite(dataSetID, H5T_NATIVE_FLOAT, H5S_ALL,
                      H5S_ALL, H5P_DEFAULT, grid);
    if (status < 0)
    {
        fprintf(stderr, "%s: Error writing %s\n", __func__, dataSet);
        return -1;
    }
    H5Dclose(dataSetID);
    // Update the XDMF
    memset(path, 0, PATH_MAX*sizeof(char)); 
    strcpy(path, h5->groupName);
    strcat(path, dataSet);
    ierr = xcloc_xdmfGrid_add(h5->h5flName, path, //dataSet,
                                 XCLOC_SINGLE_PRECISION, &h5->xdmf);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error appending %s to xdmf\n",
                __func__, dataSet);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Creates an H5 archive and associated with XDMF file.
 *
 * @param[in] h5flName    Name of HDF5 file archive to open.
 * @param[in] xdmfFile    Name of the corresponding XDMF file.  If NULL then
 *                        this will be named h5flName.xdmf.
 * @param[in] groupName   Name of the group in H5 file where datasets will
 *                        be written.  If NULL then this will be "/Images". 
 * @param[in] nx          Number of x grid points in model.  This will be
 *                        the fastest changing dimension in the model in 3D.
 * @param[in] ny          Number of y grid points in model.
 * @param[in] nz          Number of z grid points in model.  This will be 
 *                        the slowest changing dimension in the model in 3D.
 * @param[in] dx          Grid spacing in x (meters).  If nx == 1 then this
 *                        can be 0.
 * @param[in] dy          Grid spacing in y (meters).  If ny == 1 then this
 *                        can be 0.
 * @param[in] dz          Grid spacing in z (meters).  If nz == 1 then this
 *                        can be 0.
 * @param[in] x0          x model origin (meters).
 * @param[in] y0          y model origin (meters).
 * @param[in] z0          z model origin (meters).
 *
 * @result 0 indicates success.
 *
 * @copyright Ben Baker distirbuted under the MIT license.
 *
 */ 
int xcloc_h5ioGrid_open(
    const char *h5flName,
    const char *xdmfFile,
    const char *groupName,
    const int nx, const int ny, const int nz, 
    const double dx, const double dy, const double dz, 
    const double x0, const double y0, const double z0,
    struct xclocHDF5Grid_struct *h5)
{
    char *bname, *dname, *temp;
    int ierr;
    size_t lenos;
    hsize_t dims[1];
    memset(h5, 0, sizeof(struct xclocHDF5Grid_struct));
    if (nx < 1 || ny < 1 || nz < 1 || 
        (nx > 1 && dx == 0.0) || (ny > 1 && dy == 0.0) || (nz > 1 && dz == 0.0))
    {
        if (nx < 1){fprintf(stderr, "%s: nx must be positive\n", __func__);}
        if (ny < 1){fprintf(stderr, "%s: ny must be positive\n", __func__);}
        if (nz < 1){fprintf(stderr, "%s: nz must be positive\n", __func__);}
        if (nx > 1 && dx == 0.0)
        {
            fprintf(stderr, "%s: dx cannot be 0\n", __func__);
        }
        if (ny > 1 && dy == 0.0)
        { 
            fprintf(stderr, "%s: dy cannot be 0\n", __func__);
        }
        if (nz > 1 && dz == 0.0)
        {
            fprintf(stderr, "%s: dz cannot be 0\n", __func__);
        }
        return -1; 
    }
    if (h5flName == NULL)
    {
        fprintf(stderr, "%s: h5flName is NULL\n", __func__);
        return -1;
    }
    if (strlen(h5flName) < 1)
    {
        fprintf(stderr, "%s: h5flName is emtpy\n", __func__);
        return -1;
    }
    // Figure out the output directory
    temp = strdup(h5flName);
    dname = dirname(temp);
    if (dname != NULL){strcpy(h5->outDir, dname);}
    free(temp); 
    if (!os_path_isdir(h5->outDir))
    {
        if (os_makedirs(h5->outDir) != ISCL_SUCCESS)
        {
            fprintf(stderr, "%s: Failed to make directory %s\n",
                     __func__, h5->outDir);
            return -1;
        }
    }
    // Get the XDMF file name (without directory) b/c it will exist next to the
    // HDF5 archive
    temp = strdup(h5flName);
    bname = basename(temp);
    strcpy(h5->h5flName, bname);
    free(temp);
    // Set the xdmf file name
    if (xdmfFile == NULL)
    {
        strcpy(h5->xdmfFile, h5flName);
        strcat(h5->xdmfFile, ".xdmf");
    }
    else
    {
        if (strlen(xdmfFile) == 0)
        {
            strcpy(h5->xdmfFile, h5flName);
            strcat(h5->xdmfFile, ".xdmf");
        }
        else
        {
            strcpy(h5->xdmfFile, xdmfFile);
        }
    }
    // Set the mesh constants
    h5->nx = nx;
    h5->ny = ny; 
    h5->nz = nz; 
    h5->x0 = x0; 
    h5->y0 = y0; 
    h5->z0 = z0; 
    h5->dx = dx; 
    h5->dy = dy; 
    h5->dz = dz;
    // Initialize the XDMF grid
    ierr = xcloc_xdmfGrid_initialize(nx, ny, nz,
                                     dx, dy, dz,
                                     x0, y0, z0, &h5->xdmf);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to initialize XDMF grid\n", __func__);
        return -1;
    }
    // Open the HDF5 file
    h5->h5fl = H5Fcreate(h5flName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (groupName == NULL)
    {
        strcpy(h5->groupName, "/Images/");
    }
    else
    {
        lenos = strlen(groupName);
        if (lenos == 0)
        {
            strcpy(h5->groupName, "/Images/");
        }
        else
        {
            strcpy(h5->groupName, groupName);
            if (h5->groupName[lenos-1] != '/'){h5->groupName[lenos] = '/';}
        }
    }
    h5->groupID = H5Gcreate2(h5->h5fl, h5->groupName, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = (hsize_t) (nx*ny*nz);
    h5->dataSpaceID = H5Screate_simple(1, dims, NULL);
    h5->linit = true;
    return 0;
}

