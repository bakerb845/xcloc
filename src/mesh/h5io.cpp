#include <cstdio>
#include <cstring>
#include <H5Cpp.h>


int test(const std::string &fileName)
{
    H5::H5File *file = new H5::H5File(fileName, H5F_ACC_TRUNC);
}
