#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "test_suite.h"
#include "iscl/iscl/iscl.h"

#define SEED 40935  /*!< random seed - but predictable results. */

int main()
{
    int ierr;
    // Initialize ISCL to make reference solution
    iscl_init();
    srand(SEED);

    fprintf(stdout, "%s: Performing serial tests...\n", __func__);
    ierr = test_serial_fdxc();
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "%s: Failed serial tests!\n", __func__);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "%s: Performing serial DSM tests...\n", __func__);
    ierr = test_serial_dsmLocation();
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "%s: Failed serial dsm test!\n", __func__);
        return EXIT_FAILURE;
    }

    iscl_finalize();
    return EXIT_SUCCESS;
}
