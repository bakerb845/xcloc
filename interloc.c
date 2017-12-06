#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <getopt.h>
#include <limits.h>
#include "xcloc.h"
#include <mpi.h>
#include "iscl/os/os.h"
#include <iniparser.h>


static void printUsage(void);
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX]);

int main(int argc, char *argv[])
{
    char iniFile[PATH_MAX], ttFile[PATH_MAX], obsFile[PATH_MAX];
    struct xclocParms_struct xclocParms;
    struct xcloc_struct xcloc;
    double t0;
    int ierr, ierrAll, iwin, myid, nprocs, nwin, provided;
    const int master = 0;
    // Initialize MPI 
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    // Null out structures
    memset(&xclocParms, 0, sizeof(struct xclocParms_struct));
    memset(&xcloc,      0, sizeof(struct xcloc_struct));
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
nwin = 1;
        fprintf(stdout, "%s: Reading ini file...\n", __func__);
        ierr = interloc_readIni(iniFile, &xclocParms.chunkSize, ttFile, obsFile);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error reading ini file\n", __func__);
            goto BCAST_ERROR;
        }
        fprintf(stdout, "%s: Getting information from observation file...\n",
                __func__);
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
    //ierr = xcloc_initialize(MPI_COMM_WORLD, xclocParms, &xcloc);
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
    } 
    MPI_Bcast(&nwin, 1, MPI_INT, master, MPI_COMM_WORLD);
    // Loop on the windows in the observation file
    for (iwin=0; iwin<nwin; iwin++)
    {
        if (myid == master)
        {
            fprintf(stdout, "%s: Loading observation window: %d...\n",
                    __func__, iwin + 1);
        } 
        // Scatter the observations
        //ierr = xcloc_scatterDataFromRoot(nsignals, nptsSig, nptsSig,
        //                                 MPI_DOUBLE, obs, &xcloc);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error scattering data on rank %d\n",
                    __func__, myid+1);
            ierr = 1;
        }
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (ierr != 0){goto EXIT_MPI;}
        // Apply
        //ierr = xcloc_apply(&xcloc); 
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Application error on rank %d\n",
                    __func__, myid+1);
            ierr = 1;
        }
        MPI_Allreduce(&ierr, &ierrAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (ierr != 0){goto EXIT_MPI;}
    }
EXIT_MPI:;
    MPI_Finalize();
    if (ierr != 0){return EXIT_FAILURE;}
    return EXIT_SUCCESS;
}
//============================================================================//
int interloc_readIni(const char *iniFile,
                     int *chunkSize,
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
