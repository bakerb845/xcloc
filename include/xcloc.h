#ifndef XCLOC_H__
#define XCLOC_H__ 1
#include "xcloc_config.h"
#include "xcloc_enum.h"
#include "xcloc_xcfft.h"
#include "xcloc_migrate.h"
#include "xcloc_hdf5.h"
#include "xcloc_xdmf.h"
#include "xcloc_rmsFilter.h"
#ifdef XCLOC_USE_MPI
#include "xcloc_xcfftMPI.h"
//#include "xcloc_migrateMPI.h"
#endif


#ifdef __cplusplus
extern "C"
{
#endif

/*
int dales_compute_travelTimesInHomogenousMedia(
    const int nx, const int ny, const int nz, 
    const double dx, const double dy, const double dz, 
    const double x0, const double y0, const double z0, 
    const double xs, const double ys, const double zs, 
    const double vel, double *__restrict__ ttimes);
*/

#ifdef __cplusplus
}
#endif
#endif
