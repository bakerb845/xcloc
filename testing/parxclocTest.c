#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ipps.h>
#include "xcloc_finter.h"
#include "acousticGreens2D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef CHKERR
#define CHKERR(ierr, msg) \
{ \
   if (ierr != EXIT_SUCCESS) \
   { \
       fprintf(stderr, "ERROR calling %s: %s line %d\n", msg, __func__, __LINE__); \
       return EXIT_FAILURE; \
   } \
};
#endif

int test_parallel_xcloc(const MPI_Comm comm, const int root)
{

    return EXIT_SUCCESS;
}
