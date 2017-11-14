#ifndef XCLOC_ACOUSTIC_GREENS_H__
#define XCLOC_ACOUSTIC_GREENS_H__

#ifdef __cplusplus
extern "C"
{
#endif

// Greens functions stuff
int dales_unitTests_computeRickerWavelet(const int npts,
                                         const double dt,
                                         const double peakFreq,
                                         const bool lnorm,
                                         const bool lshift,
                                         double *__restrict__ ricker);
int dales_unitTests_acousticGreensLineSource(
    const int nrec,
    const double vel,
    const double rho,
    const double Q,
    const int npts,
    const double dt, 
    const double *__restrict__ xs, 
    const double *__restrict__ xr, 
    const double *__restrict__ stf,
    double *__restrict__ G);
int dales_unitTests_acousticGreensLineSourceFD(
    const bool lverb,
    const double vel,
    const double rho,
    const double Q,
    const double xs[3],
    const double xr[3],
    const int nomega,
    const double *__restrict__ omega,
    const double complex *__restrict__ stf,
    double complex *__restrict__ G);

#ifdef __cplusplus
}
#endif
#endif