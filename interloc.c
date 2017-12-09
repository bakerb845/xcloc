#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <getopt.h>
#include <limits.h>
#include "xcloc.h"
#include "xcloc_h5ioUtils.h"
#include <mpi.h>
#include "iscl/os/os.h"
#include <iniparser.h>
#include <ipps.h>
#include "iscl/time/time.h"


int interloc_readIni(const char *iniFile,
                     int *chunkSize,
                     int *envFIRlen,
                     bool *lphaseXCs,
                     char ttFile[PATH_MAX],
                     char obsFile[PATH_MAX]);
static void printUsage(void);
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX]);

int main(int argc, char *argv[])
{
    char iniFile[PATH_MAX], ttFile[PATH_MAX], obsFile[PATH_MAX];
    struct xclocParms_struct xclocParms;
    struct xcloc_struct xcloc;
    struct h5TravelTimeModel_struct ttimes;
    struct xclocHDF5Grid_struct h5io;
    double *data, dt, t0;
    float *image;
    int *tpg, chunkSize, envFIRlen,
        i, ierr, ierrAll, iwin, j, k, myid, npts, nprocs, nsgroups, ntrace, 
        ntraceAll, nwin, provided;
    bool lphaseXCs;
    hid_t dataFID, ttFID;
    const int master = 0;
    // Initialize MPI 
    nwin = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    // Null out structures
    image = NULL;
    memset(&xclocParms, 0, sizeof(struct xclocParms_struct));
    memset(&xcloc,      0, sizeof(struct xcloc_struct));
    memset(&h5io,       0, sizeof(struct xclocHDF5Grid_struct));
    memset(&ttimes,     0, sizeof(struct h5TravelTimeModel_struct));
    // Get the command line parameters
    ierr = 0;
    if (myid == master)
    {
        ierr = parseArguments(argc, argv, iniFile);
        if (ierr != 0)
        {
            if (ierr !=-2)
            {
                fprintf(stderr, "%s: Invalid input arguments\n", __func__);
            }
        }
    }
    MPI_Bcast(&ierr, 1, MPI_INT, master, MPI_COMM_WORLD); 
    if (ierr != 0)
    {
        if (ierr ==-2){ierr = 0;}
        goto EXIT_MPI;
    }
    // Have the master read the ini file
    if (myid == master)
    {
        // Do something until user gives me a better idea
        xclocParms.nfftProcs = 1;
        xclocParms.ngridProcs = nprocs;
        // Read parameters from ini file
        fprintf(stdout, "%s: Reading ini file...\n", __func__);
        ierr = interloc_readIni(iniFile,
                                &chunkSize,
                                &envFIRlen,
                                &lphaseXCs,
                                ttFile, obsFile);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error reading ini file\n", __func__);
            goto BCAST_ERROR;
        }
        fprintf(stdout, "%s: Getting information from observation file...\n",
                __func__);
        dataFID = H5Fopen(obsFile, H5F_ACC_RDONLY, H5P_DEFAULT);
        ierr = xcloc_h5ioUtils_data_getNumberOfWindowsAndGroups(
                        dataFID, &nwin, &nsgroups);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error getting number of windows/groups\n",
                    __func__);
            goto BCAST_ERROR;
        } 
        fprintf(stdout, "%s: Number of signal groups: %d\n",
                __func__, nsgroups);
        fprintf(stdout, "%s: Number of windows: %d\n", __func__, nwin);
        ntraceAll = 0;
        iwin = 1;
        tpg = (int *) calloc((size_t) nsgroups, sizeof(int));
        for (i=1; i<=nsgroups; i++)
        {
            ierr = xcloc_h5ioUtils_data_readWindowAndGroup(dataFID,
                                                           iwin, i,
                                                           &npts, &dt, &ntrace,
                                                           &data);
            tpg[i-1] = ntrace;
            ntraceAll = ntraceAll + ntrace;
            free(data);
        }
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error getting dt/npts/ntrace\n", __func__);
            goto BCAST_ERROR;
        }
        fprintf(stdout, "%s: Total number of traces: %d\n",
                __func__, ntraceAll);
        fprintf(stdout, "%s: Number of points in traces: %d\n", 
                __func__, npts);
        fprintf(stdout, "%s: Sampling period %e (s)\n", __func__, dt);
        H5Fclose(dataFID);
        // Read the travel time info 
        ierr = xcloc_h5ioUtils_openTravelTimeModelForReading(ttFile,
                                                             &ttimes,
                                                             &ttFID);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error lifting travel time model info\n",
                    __func__);
            goto BCAST_ERROR;
        }
        H5Fclose(ttFID);
        fprintf(stdout, "%s: Number of grid points: %d\n",
                __func__, ttimes.ngrd);
        // Initialize the rest of the parameters 

        xclocParms.dt = dt; 
        xclocParms.lphaseXCs = lphaseXCs;
        xclocParms.chunkSize = chunkSize;
        xclocParms.envFIRLen = envFIRlen;
        xclocParms.npts    = npts;
        xclocParms.nptsPad = npts;
        xclocParms.ngrd    = ttimes.ngrd;
        xclocParms.nsignals = ntraceAll;
        xclocParms.signalGroup 
            = (int *) calloc((size_t) xclocParms.nsignals, sizeof(int));
        k = 0;
        for (i=0; i<nsgroups; i++)
        {
            for (j=0; j<tpg[i]; j++)
            {
                xclocParms.signalGroup[k] = i;
                k = k + 1;
            }
        }
        free(tpg);
    }
BCAST_ERROR:;
    MPI_Bcast(&ierr, 1, MPI_INT, master, MPI_COMM_WORLD); 
    if (ierr != 0){goto EXIT_MPI;}
    // Initialize
    if (myid == master)
    {
        fprintf(stdout, "%s: Initializing xcloc...\n", __func__);
    }
    t0 = MPI_Wtime();
    ierr = xcloc_initialize(MPI_COMM_WORLD, xclocParms, &xcloc);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error initializing on rank %d\n", __func__, myid);
        ierr = 1;
    }
    MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (ierr != 0){goto EXIT_MPI;}
    if (myid == master)
    {
        fprintf(stdout, "%s: Setting travel time tables....\n", __func__);
        ttFID = H5Fopen(ttFile, H5F_ACC_RDONLY, H5P_DEFAULT);
    } 
    for (i=0; i<xcloc.nTotalSignals; i++)
    {
        if (myid == master)
        {
            ttimes.signalNumber = i + 1;
            xcloc_h5ioUtils_readTravelTimeModel(ttFID, &ttimes);
        }
        ierr = xcloc_setTableFromRoot(i, ttimes.ngrd,
                                      XCLOC_SINGLE_PRECISION,
                                      ttimes.ttimes, &xcloc);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error setting ttimes on rank %d\n",
                    __func__, myid); 
            ierr = 1;
        }
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (ierr != 0){goto EXIT_MPI;}
        if (myid == master && ttimes.ttimes != NULL)
        {
            free(ttimes.ttimes);
        }
    }
    if (myid == master){H5Fclose(ttFID);}
    // Open the data window
    MPI_Bcast(&nwin, 1, MPI_INT, master, MPI_COMM_WORLD);
    if (myid == master)
    {
        dataFID = H5Fopen(obsFile, H5F_ACC_RDONLY, H5P_DEFAULT);
        fprintf(stdout, "%s: Beginning the processing...\n", __func__);
    }
    // Loop on the windows in the observation file
    for (iwin=0; iwin<nwin; iwin++)
    {
        if (myid == master)
        {
            fprintf(stdout, "%s: Loading observation window: %d...\n",
                    __func__, iwin + 1);
            printf("fix here\n");
            double *data1, *data2;
            time_tic();
            ierr = xcloc_h5ioUtils_data_readWindowAndGroup(dataFID,
                                                           iwin+1, 1,
                                                           &npts, &dt,
                                                           &ntrace,
                                                           &data1);
            ierr = xcloc_h5ioUtils_data_readWindowAndGroup(dataFID,
                                                           iwin+1, 1,
                                                           &npts, &dt,
                                                           &ntrace,
                                                           &data2); 
            data = (double *) calloc(2*ntrace*npts, sizeof(double));
            ippsCopy_64f(data1,  data,              ntrace*npts); 
            ippsCopy_64f(data2, &data[ntrace*npts], ntrace*npts);
            free(data1);
            free(data2);
            fprintf(stdout, "%s: File load time: %f\n", __func__, time_toc());
        } 
        // Scatter the observations
        ierr = xcloc_scatterDataFromRoot(2*ntrace, npts, npts,
                                         MPI_DOUBLE, data, &xcloc);
        if (myid == master){free(data);}
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error scattering data on rank %d\n",
                    __func__, myid+1);
            ierr = 1;
        }
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (ierr != 0){goto EXIT_MPI;}
        // Apply
        t0 = MPI_Wtime();
        ierr = xcloc_apply(&xcloc); 
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Application error on rank %d\n",
                    __func__, myid+1);
            ierr = 1;
        }
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (ierr != 0){goto EXIT_MPI;}
        if (myid == master)
        {
            fprintf(stdout, "%s: Application time: %f\n",
                    __func__, MPI_Wtime() - t0);
        } 
    }
    if (myid == master){H5Fclose(dataFID);}

    // Get the image
    if (myid == master)
    {
        fprintf(stdout, "%s: Initializing H5 archive\n", __func__);
        ierr = xcloc_h5ioGrid_open("./puget_image.h5", "./puget_image.xdmf",
                           "/migrate",
                           ttimes.nGrid[2], ttimes.nGrid[1], ttimes.nGrid[0], // nx, ny, nz, 
                           ttimes.dGrid[2], ttimes.dGrid[1], ttimes.dGrid[0], // dx, dy, dz, 
                           ttimes.grid0[2], ttimes.grid0[1], ttimes.grid0[0], //x0, y0, z0, 
                           &h5io);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error initializing H5 file\n", __func__);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        fprintf(stdout, "%s: Gathering image...\n", __func__);
        image = (float *) calloc((size_t) ttimes.ngrd, sizeof(float));
    }   
    ierr = xcloc_gatherMigrationImage(xcloc.ngrd, xcloc, image);
    if (myid == master)
    {
        if (ierr == 0)
        {
            int optIndx, iz, ix, iy;
            float xopt, yopt, zopt, xmax;
            int nz = ttimes.nGrid[0];
            int ny = ttimes.nGrid[1];
            int nx = ttimes.nGrid[2];
            ippsMaxIndx_32f(image, xcloc.ngrd, &xmax, &optIndx); 
            for (iz=0; iz<nz; iz++)
            {
                for (iy=0; iy<ny; iy++)
                {
                    for (ix=0; ix<nx; ix++)
                    {
                        if (iz*nx*ny + iy*nx + ix == optIndx)
                        {
                            zopt = ttimes.grid0[0] + iz*ttimes.dGrid[0];
                            yopt = ttimes.grid0[1] - iy*ttimes.dGrid[1];
                            xopt = ttimes.grid0[2] + ix*ttimes.dGrid[2];
                            goto EXIT_ME;
                        }
                    }
                }
            }
            EXIT_ME:;
            fprintf(stdout, "%s: Optimium location: %f %f %f\n",
                    __func__, zopt, yopt, xopt);
            fprintf(stdout, "%s: Writing image...\n", __func__);
            ierr = xcloc_h5ioGrid_writeDataSet32f("vpvs", nx, ny, nz,
                                                  image, &h5io);
            ierr = xcloc_h5ioGrid_close(&h5io);
        }
        else
        {
            fprintf(stderr, "%s: Error gathering migration image\n", __func__);
        }
    }
    if (image != NULL){free(image);}

EXIT_MPI:;
    MPI_Finalize();
    if (ierr != 0){return EXIT_FAILURE;}
    return EXIT_SUCCESS;
}
//============================================================================//
int interloc_readIni(const char *iniFile,
                     int *chunkSize,
                     int *envFIRlen,
                     bool *lphaseXCs,
                     char ttFile[PATH_MAX],
                     char obsFile[PATH_MAX])
{
    dictionary *ini;
    const char *s;
    int ierr;
    // Load the ini file
    ierr = 1;
    *chunkSize = 2048;
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: Error ini file %s doesn't exist\n",
                 __func__, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    // Get the chunk size
    *chunkSize = iniparser_getint(ini, "interloc:chunkSize\0", 2048);
    if (*chunkSize%XCLOC_MEM_ALIGNMENT != 0)
    {
        fprintf(stderr, "%s: chunkSize = %d not divisible by %d\n",
                __func__, *chunkSize, XCLOC_MEM_ALIGNMENT);
        goto ERROR;
    }
    // Get the FIR filter length for the envelope computation
    *envFIRlen = iniparser_getint(ini, "interloc:envelopeFIRLength\0", 251);
    if (*envFIRlen < 1)
    {
        fprintf(stdout, "%s: Setting envFIRlen to 251\n", __func__);
        *envFIRlen = 251;
    }
    if (*envFIRlen%2 != 1){*envFIRlen = *envFIRlen + 1;}
    // Using phase correlations or cross-correlations?
    *lphaseXCs = iniparser_getboolean(ini, "interloc:usePhaseXC\0", true);
    // Get the travel time table archive
    s = iniparser_getstring(ini, "interloc:travelTimeFile\0", NULL);
    if (s != NULL)
    {
        strcpy(ttFile, s);
    }
    else
    {
        fprintf(stderr, "%s: Error travelTimeFile not defined\n", __func__);
        goto ERROR;
    }
    // Get the observations archive
    s = iniparser_getstring(ini, "interloc:observationFile\0", NULL);
    if (s != NULL)
    {
        strcpy(obsFile, s);
    }
    else
    {
        fprintf(stderr, "%s: Error observationFile not defined\n", __func__);
        goto ERROR;
    }
    ierr = 0;
ERROR:;
    iniparser_freedict(ini);
    return ierr;
}
//============================================================================//
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX])
{
    bool liniFile;
    memset(iniFile, 0, PATH_MAX*sizeof(char));
    liniFile = false;
    while (true)
    {
        static struct option longOptions[] =
        {
            {"help", no_argument, 0, '?'}, 
            {"help", no_argument, 0, 'h'}, 
            {"iniFile", required_argument, 0, 'i'},
            {0, 0, 0, 0}
        };
        int c, optionIndex;
        c = getopt_long(argc, argv, "?hi:",
                        longOptions, &optionIndex);
        if (c ==-1){break;}
        if (c == 'i')
        {
            strcpy(iniFile, (const char *) optarg);
            liniFile = true;
        }
        else if (c == 'h' || c == '?')
        {
            printUsage();
            return -2;
        }
        else
        {
            fprintf(stderr, "%s: Unknown options: %s\n",
                    __func__, argv[optionIndex]);
        }
    }
    if (liniFile)
    {
        if (!os_path_isfile(iniFile))
        {
            fprintf(stderr, "%s: ini file %s does not exist\n",
                    __func__, iniFile);
            return -1;
        }
    }
    else
    {
        fprintf(stderr, "%s: Invalid arguments\n", __func__);
        printUsage();
        return -1;
    }
    return 0;
}
//============================================================================//
static void printUsage(void)
{
    fprintf(stdout, "Usage:\n    interloc -i iniFile\n\n");
    fprintf(stdout, "Require arguments:\n");
    fprintf(stdout, "   -i iniFile is the initialization file\n\n");
    fprintf(stdout, "Optional arguments:\n");
    fprintf(stdout, "   -h displays this message\n\n");
    return; 
}
