#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "xcloc_h5ioUtils.h"
#include "sacioh5.h"
#include "iscl/os/os.h"

#define DATA_WINDOWS "/DataWindows"
#define WINDOW_GROUP "Window"
#define DATA_NAME "Data"
#define NUMBER_OF_WINDOWS "NumberOfWindows"
#define NUMBER_OF_SIGNAL_GROUPS "NumberOfSignalGroups"
#define MODEL_GROUP "/TravelTimeModels"
#define NUMBER_OF_MODELS "NumberOfModels"
#define MODEL_TYPE "ModelType"
#define NUMBER_OF_GRIDPOINTS "NumberOfGridPoints"
#define GEOMETRY_GROUP "/ModelGeometry"

int xcloc_h5ioUtils_data_initialize(const char *archiveName,
                                    const int nwin,
                                    const int nsgroups,
                                    hid_t *h5fl)
{
    char dirname[PATH_MAX], window[128];
    hid_t attrID, dataSpace, groupID, windowGroupID;
    hsize_t dims[1] = {1};
    enum isclError_enum isclError;
    int i;
    if (nwin < 1)
    {   
        fprintf(stderr, "%s: No data windows\n", __func__);
        return -1;
    }
    if (nsgroups < 1)
    {   
        fprintf(stderr, "%s: No signal groups\n", __func__);
        return -1;
    }
    os_dirname_work(archiveName, dirname); 
    if (!os_path_isdir(dirname))
    {
        isclError = os_makedirs(dirname);
        if (isclError != ISCL_SUCCESS)
        {
            fprintf(stderr, "%s: Failed to make directory %s\n",
                    __func__, dirname);
            return -1;
        }
    }
    if (os_path_isfile(archiveName))
    {
        fprintf(stdout, "%s: Clobbering %s\n", __func__, archiveName);
    }
    // Create the archive
    *h5fl = H5Fcreate(archiveName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    groupID = H5Gcreate(*h5fl, DATA_WINDOWS, //"/DataWindows\0", 
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (i=0; i<nwin; i++)
    {
        memset(window, 0, 128*sizeof(char));
        sprintf(window, "%s_%d", WINDOW_GROUP, i+1); 
        windowGroupID = H5Gcreate(groupID, window, 
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
        H5Gclose(windowGroupID);
    }
    // Aadd some attributes
    dataSpace = H5Screate_simple(1, dims, NULL);
    // Number of windows
    attrID = H5Acreate2(groupID, NUMBER_OF_WINDOWS,
                        H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT); 
    H5Awrite(attrID, H5T_NATIVE_INT, &nwin); 
    H5Aclose(attrID);
    attrID = H5Acreate2(groupID, NUMBER_OF_SIGNAL_GROUPS,
                        H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT); 
    H5Awrite(attrID, H5T_NATIVE_INT, &nsgroups); 
    H5Aclose(attrID);
    H5Sclose(dataSpace);
    H5Gclose(groupID);
    return 0;
}
//============================================================================//
/*!
 */
int xcloc_h5ioUtils_data_writeWindow(const int iwin, const int ig,
                                     const int *compress,
                                     const size_t *traceChunkSz,
                                     const size_t *dataChunkSz,
                                     const hid_t fileID,
                                     const int ntrace,
                                     const struct sacData_struct *sac)
{
    char groupName[128], dataName[128];
    hid_t groupID;
    int ierr;
    int c = 0;
    if (compress != NULL){c = MAX(0, MIN(9, *compress));}
    memset(groupName, 0, 128*sizeof(char));
    sprintf(groupName, "%s/%s_%d", DATA_WINDOWS, WINDOW_GROUP, iwin);
    if (H5Lexists(fileID, groupName, H5P_DEFAULT) < 0)
    {
        fprintf(stderr, "%s: Group %s doesn't exist\n", __func__, groupName);
        return -1;
    }
    groupID = H5Gopen2(fileID, groupName, H5P_DEFAULT);
    memset(dataName, 0, 128*sizeof(char));
    sprintf(dataName, "%s_%d", DATA_NAME, ig);
    ierr = sacioh5_writeTracesToChunk(dataName, ntrace, c,
                                      groupID, traceChunkSz, dataChunkSz, sac);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error writing data\n", __func__);
    }
    H5Gclose(groupID); 
    return ierr;
}
//============================================================================//
int xcloc_h5ioUtils_data_finalize(const hid_t fileID)
{
    H5Fclose(fileID);
    return 0;
}
//============================================================================//
int xcloc_h5ioUtils_data_getNumberOfWindowsAndGroups(
    const hid_t fileID,
    int *nwin, int *nsgroups)
{
    hid_t attrID, groupID;
    herr_t status;
    *nwin = 0;
    *nsgroups = 0;
    if (H5Lexists(fileID, DATA_WINDOWS, H5P_DEFAULT) < 0)
    {
        fprintf(stderr, "%s: %s doesn't exist\n", __func__, DATA_WINDOWS);
        return -1;
    }
    groupID = H5Gopen2(fileID, DATA_WINDOWS, H5P_DEFAULT);

    attrID = H5Aopen_name(groupID, NUMBER_OF_WINDOWS);
    status = H5Aread(attrID, H5T_NATIVE_INT, nwin);
    if (status < 0)
    {
        fprintf(stderr, "%s: Error reading %s\n", __func__, NUMBER_OF_WINDOWS);
    }
    H5Aclose(attrID);

    attrID = H5Aopen_name(groupID, NUMBER_OF_SIGNAL_GROUPS);
    status = H5Aread(attrID, H5T_NATIVE_INT, nsgroups);
    if (status < 0)
    {   
        fprintf(stderr, "%s: Error reading %s\n",
                __func__, NUMBER_OF_SIGNAL_GROUPS);
    }   
    H5Aclose(attrID);

    H5Gclose(groupID);
    return 0;
}
//============================================================================//
int xcloc_h5ioUtils_data_readWindowAndGroup(
    const hid_t fileID,
    const int iwin, const int ig,
    int *npts, double *dt, int *ntrace,
    double **data)
{
    char groupName[128], dataName[128];
    struct sacData_struct *sac = NULL;
    int i, ierr;
    size_t nbytes;
    double *ts = NULL;
    hid_t groupID;
    *dt = 0.0;
    *npts = 0;
    memset(groupName, 0, 128*sizeof(char));
    sprintf(groupName, "%s/%s_%d", DATA_WINDOWS, WINDOW_GROUP, iwin);
    groupID = H5Gopen2(fileID, groupName, H5P_DEFAULT);
    memset(dataName, 0, 128*sizeof(char));
    sprintf(dataName, "%s_%d", DATA_NAME, ig);
    ierr = sacioh5_readChunkedTraces(dataName,
                                     groupID,
                                     ntrace, &sac);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error reading %s\n", __func__, dataName);
        return -1;
    }
    *npts = sac[0].header.npts;
    *dt = sac[0].header.delta;
    nbytes = (size_t) (*ntrace * *npts)*sizeof(double);
    ts = aligned_alloc(XCLOC_MEM_ALIGNMENT, nbytes);
    for (i=0; i<*ntrace; i++)
    {
        memcpy(&ts[*npts*i], sac[i].data, (size_t) *npts*sizeof(double));
    }
    H5Gclose(groupID);
    if (ierr != 0){return -1;}
    if (sac != NULL)
    {
        for (i=0; i<*ntrace; i++)
        {
            sacio_free(&sac[i]);
        }
        free(sac);
    }
    *data = ts;
    return 0;
}
//============================================================================//
//                           Write/Read Model                                 //
//============================================================================//
int xcloc_h5ioUtils_initTravelTimeArchive(
    const char *archiveName,
    const int nmodels,
    const struct h5TravelTimeModel_struct ttimes,
    hid_t *h5fl)
{
    char dirname[PATH_MAX]; //, groupName[128];
    hid_t attrID, dataSet, dataSpace, groupID; //, modelID;
    hsize_t dims[1] = {1};
    int type;
    enum isclError_enum isclError;
    os_dirname_work(archiveName, dirname); 
    if (!os_path_isdir(dirname))
    {
        isclError = os_makedirs(dirname);
        if (isclError != ISCL_SUCCESS)
        {
            fprintf(stderr, "%s: Failed to make directory %s\n",
                    __func__, dirname);
            return -1; 
        }
    }
    if (os_path_isfile(archiveName))
    {
        fprintf(stdout, "%s: Clobbering %s\n", __func__, archiveName);
    }
    if (ttimes.modelType != XCLOC_MODEL_CARTESIAN_GRIDDED)
    {
        fprintf(stderr, "%s: Only gridded cartesian models are done\n",
                __func__);
        return -1;
    }
    // Create the archive
    *h5fl = H5Fcreate(archiveName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    groupID = H5Gcreate(*h5fl, MODEL_GROUP, //"/TravelTimeModels\0",
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataSpace = H5Screate_simple(1, dims, NULL);
    attrID = H5Acreate2(groupID, NUMBER_OF_MODELS,
                        H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrID, H5T_NATIVE_INT, &nmodels);
    H5Aclose(attrID);

    type = (int) ttimes.modelType;
    attrID = H5Acreate2(groupID, MODEL_TYPE,
                        H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrID, H5T_NATIVE_INT, &type);
    H5Aclose(attrID);

    attrID = H5Acreate2(groupID, NUMBER_OF_GRIDPOINTS,
                        H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrID, H5T_NATIVE_INT, &ttimes.ngrd);
    H5Aclose(attrID);

    H5Sclose(dataSpace);
/*
    for (i=0; i<nmodels; i++)
    {
        memset(groupName, 0, 128*sizeof(char));
        sprintf(groupName, "Model_%d", i+1);
        modelID = H5Gcreate(groupID, groupName,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Gclose(modelID);
    }
*/
    H5Gclose(groupID);

    groupID = H5Gcreate(*h5fl, GEOMETRY_GROUP, //"/ModelGeometry"
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ttimes.modelType == XCLOC_MODEL_CARTESIAN_GRIDDED)
    {
        dims[0] = 3;
        dataSpace = H5Screate_simple(1, dims, NULL);

        dataSet = H5Dcreate2(groupID, "GridSpacing\0", H5T_NATIVE_FLOAT,
                             dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 ttimes.dGrid);
        H5Dclose(dataSet);

        dataSet = H5Dcreate2(groupID, "GridOrigin\0", H5T_NATIVE_FLOAT,
                             dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 ttimes.grid0);
        H5Dclose(dataSet);

        dataSet = H5Dcreate2(groupID, "GridSize\0", H5T_NATIVE_INT,
                             dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 ttimes.nGrid);
        H5Dclose(dataSet);

        H5Sclose(dataSpace);
    }
    H5Gclose(groupID);
    return 0;
}

int xcloc_h5ioUtils_closeTravelTimeArchive(const hid_t h5fl)
{
    H5Fclose(h5fl);
    return 0;
}
//============================================================================//
int xcloc_h5ioUtils_openTravelTimeModelForReading(
    const char *archiveName,
    struct h5TravelTimeModel_struct *ttimes,
    hid_t *fileID)
{
    int rank, type;
    hid_t attrID, dataSpace, dataSet, groupID, modelGroup, memSpace;
    hsize_t dims[1];
    herr_t status;
    memset(ttimes, 0, sizeof(struct h5TravelTimeModel_struct));
    if (!os_path_isfile(archiveName))
    {
        fprintf(stderr, "%s: Error travel time archive %s doesn't exist\n",
                __func__, archiveName);
        return -1;
    }
    *fileID = H5Fopen(archiveName, H5F_ACC_RDONLY, H5P_DEFAULT);
    // Lift some information
    groupID = H5Gopen2(*fileID, MODEL_GROUP, H5P_DEFAULT); 

    attrID = H5Aopen_name(groupID, MODEL_TYPE);
    status = H5Aread(attrID, H5T_NATIVE_INT, &type);
    if (status < 0)
    {   
        fprintf(stderr, "%s: Error reading %s\n", __func__, MODEL_TYPE);
    }   
    H5Aclose(attrID);
    ttimes->modelType = type;

    attrID = H5Aopen_name(groupID, NUMBER_OF_GRIDPOINTS);
    status = H5Aread(attrID, H5T_NATIVE_INT, &ttimes->ngrd);
    if (status < 0)
    {   
        fprintf(stderr, "%s: Error reading %s\n",
                __func__, NUMBER_OF_GRIDPOINTS);
    }   
    H5Aclose(attrID);

    // Model information
    if (ttimes->modelType == XCLOC_MODEL_CARTESIAN_GRIDDED)
    {
        modelGroup = H5Gopen2(*fileID, GEOMETRY_GROUP, H5P_DEFAULT); 
        // Grid origin
        dataSet = H5Dopen(modelGroup, "GridOrigin", H5P_DEFAULT);
        dataSpace = H5Dget_space(dataSet);
        rank = H5Sget_simple_extent_ndims(dataSpace);
        H5Sget_simple_extent_dims(dataSpace, dims, NULL);
        memSpace = H5Screate_simple(rank, dims, NULL);
        ttimes->ttimes = (float *) calloc((size_t) dims[0], sizeof(float));
        status = H5Dread(dataSet, H5T_NATIVE_FLOAT, memSpace, dataSpace,
                         H5P_DEFAULT, ttimes->grid0);
        H5Sclose(memSpace);
        H5Sclose(dataSpace);
        H5Dclose(dataSet);

        // Grid spacing
        dataSet = H5Dopen(modelGroup, "GridSpacing", H5P_DEFAULT);
        dataSpace = H5Dget_space(dataSet);
        rank = H5Sget_simple_extent_ndims(dataSpace);
        H5Sget_simple_extent_dims(dataSpace, dims, NULL);
        memSpace = H5Screate_simple(rank, dims, NULL);
        ttimes->ttimes = (float *) calloc((size_t) dims[0], sizeof(float));
        status = H5Dread(dataSet, H5T_NATIVE_FLOAT, memSpace, dataSpace,
                         H5P_DEFAULT, ttimes->dGrid);
        H5Sclose(memSpace);
        H5Sclose(dataSpace);
        H5Dclose(dataSet);

        // Grid size
        dataSet = H5Dopen(modelGroup, "GridSize", H5P_DEFAULT);
        dataSpace = H5Dget_space(dataSet);
        rank = H5Sget_simple_extent_ndims(dataSpace);
        H5Sget_simple_extent_dims(dataSpace, dims, NULL);
        memSpace = H5Screate_simple(rank, dims, NULL);
        ttimes->ttimes = (float *) calloc((size_t) dims[0], sizeof(float));
        status = H5Dread(dataSet, H5T_NATIVE_INT, memSpace, dataSpace,
                         H5P_DEFAULT, ttimes->nGrid);
        H5Sclose(memSpace);
        H5Sclose(dataSpace);
        H5Dclose(dataSet);
        // Close model group
        H5Gclose(modelGroup);
    } 
    else
    {
        fprintf(stderr, "%s: Non cartesian gridded model not yet done\n",
                __func__);
        return -1;
    }

    H5Gclose(groupID);
    return 0; 
}
//============================================================================//
/*
int xcloc_h5ioUtils_closeTravelTimeModelForReading(
    const hid_t fileID,
    struct h5TravelTimeModel_struct *ttimes)
{
    if (ttimes->ttimes != NULL){free(ttimes->ttimes);}
}
*/
//============================================================================//
int xcloc_h5ioUtils_readTravelTimeModel(
    const hid_t fileID,
    struct h5TravelTimeModel_struct *ttimes)
{
    char dataName[128];
    hid_t dataSet, dataSpace, memSpace;
    hsize_t dims[1];
    herr_t status;
    int rank;
    if (ttimes->signalNumber < 1)
    {
        fprintf(stderr, "%s: Error signalNumber must be positive\n",
                __func__);
        return -1;
    }
    memset(dataName, 0, 128*sizeof(char));
    sprintf(dataName, "%s/Model_%d", MODEL_GROUP, ttimes->signalNumber); 
    if (H5Lexists(fileID, dataName, H5P_DEFAULT) < 0)
    {
        fprintf(stderr, "%s: Error dataSet %s doesn't exist\n", __func__,
                dataName);
    }
    dataSet = H5Dopen(fileID, dataName, H5P_DEFAULT);
    dataSpace = H5Dget_space(dataSet);
    rank = H5Sget_simple_extent_ndims(dataSpace);
    if (rank != 1)
    {
        fprintf(stderr, "%s: Only 1 dimensional arrays are done\n", __func__);
        return -1;
    }
    H5Sget_simple_extent_dims(dataSpace, dims, NULL);
    memSpace = H5Screate_simple(rank, dims, NULL);
    ttimes->ttimes = (float *) calloc((size_t) dims[0], sizeof(float));
    status = H5Dread(dataSet, H5T_NATIVE_FLOAT, memSpace, dataSpace,
                     H5P_DEFAULT, ttimes->ttimes);
    if (status < 0)
    {
        fprintf(stderr, "%s: Error reading data\n", __func__);
        return -1;
    }
    H5Sclose(memSpace);
    H5Sclose(dataSpace);
    H5Dclose(dataSet);

    return 0; 
}

//============================================================================//
int xcloc_h5ioUtils_writeTravelTimeModel(
    const hid_t fileID,
    const struct h5TravelTimeModel_struct ttimes)
{
    char groupName[128], dataName[128];
    hsize_t dims[1];
    hid_t dataSet, dataSpace, groupID;
    if (ttimes.modelType != XCLOC_MODEL_CARTESIAN_GRIDDED)
    {
        fprintf(stderr, "%s: Only gridded cartesian models are done\n",
                __func__);
        return -1;
    }
    sprintf(groupName, "%s", MODEL_GROUP);
    if (H5Lexists(fileID, groupName, H5P_DEFAULT) < 0)
    {
        fprintf(stderr, "%s: Group %s doesn't exist\n", __func__, groupName);
        return -1;
    }
    if (ttimes.signalNumber < 1)
    {
        fprintf(stderr, "%s: Error signalNumber must be positive\n",
                __func__);
        return -1;
    }
    memset(dataName, 0, 128*sizeof(char));
    sprintf(dataName, "%s_%d", "Model", ttimes.signalNumber);
    groupID = H5Gopen2(fileID, groupName, H5P_DEFAULT);
    if (ttimes.modelType == XCLOC_MODEL_CARTESIAN_GRIDDED)
    {
        dims[0] = ttimes.ngrd;
        dataSpace = H5Screate_simple(1, dims, NULL);
        dataSet = H5Dcreate2(groupID, dataName, H5T_NATIVE_FLOAT,
                             dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 ttimes.ttimes);
        H5Dclose(dataSet);
        H5Sclose(dataSpace);
    }
    H5Gclose(groupID);
    return 0;
}

