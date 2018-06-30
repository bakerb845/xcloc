#ifndef ISCL_HACK_H__
#define ISCL_HACK_H__

#ifdef __cplusplus
extern "C"
{
#endif
#include <complex.h>

int fft_rfftfreqs64f_work_hack(const int n, const double dt,
                          const int lenf, double *__restrict__ freqs);
int fft_irfft64z_work_hack(const int nx, const double complex *__restrict__ x,
                      const int n, double *__restrict__ y);
int fft_fftshift64f_work_hack(const int n, const double *__restrict__ x,
                         double *__restrict__ xshift);
int fft_rfft64f_work_hack(const int nx, const double *__restrict__ x, const int n,
                     const int ny, double complex *__restrict__ y);

#ifdef __cplusplus
}
#endif

#endif
