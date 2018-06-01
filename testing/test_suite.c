#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "test_suite.h"

int main()
{
    int ierr;
    fprintf(stdout, "%s: Performing serial tests...\n", __func__);
    ierr = test_serial_fdxc();
    if (ierr != EXIT_SUCCESS)
    {
        fprintf(stderr, "%s: Failed serial tests!\n", __func__);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
