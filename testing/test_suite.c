#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "test_suite.h"


int main()
{
    int ierr;
    printf("%s: Performing serial tests...\n", __func__);
    test_serial_fdxc();
    return EXIT_SUCCESS;
}
