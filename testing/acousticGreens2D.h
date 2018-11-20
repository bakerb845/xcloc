#ifndef XCLOC_ACOUSTIC_GREENS_H__
#define XCLOC_ACOUSTIC_GREENS_H__
#include <stdbool.h>
#include <complex.h>

#ifdef __cplusplus
extern "C"
{
#endif

int acousticGreens2D_computeGreensFunctions(
    const int nsrc, const int nrec, const int nptsSig,
    const double fcent, const double dt, 
    const bool lnorm, const bool lshift,
    const double vel, const double rho,
    const double Q,
    const double srcScale[],
    const double xs[],
    const double xr[],
    double **obsOut);
int acousticGreens2D_computeRandomReceiverLocations(
    const int nrec,
    const double x0, const double y0, const double z0, 
    const double x1, const double y1, const double z1, 
    double *xr);
int acousticGreens2D_computeRickerWavelet(const int npts,
                                          const double dt,
                                          const double peakFreq,
                                          const bool lnorm,
                                          const bool lshift,
                                          double ricker[]);
int acousticGreens2D_computeGreensLineSource(
    const int nrec,
    const double vel,
    const double rho,
    const double Q,
    const int npts,
    const double dt, 
    const double xs[],
    const double xr[],
    const double stf[],
    double G[]);
int acousticGreens2D_computeLineSourceFD(
    //const bool lverb,
    const double vel,
    const double rho,
    const double Q,
    const double xs[3],
    const double xr[3],
    const int nomega,
    const double omega[],
    const double complex stf[],
    double complex G[]);
int acousticGreens2D_computeTravelTimeTable(
    const int nx, const int ny, const int nz, 
    const double vel,
    const double x0, const double y0, const double z0, 
    const double dx, const double dy, const double dz, 
    double xr, double yr, double zr, 
    double ttable[]);

#ifdef __cplusplus
}
#endif
#endif
