#ifndef XCLOC_XCFFT_H__
#define XCLOC_XCFFT_H__ 1
#include "xcloc_config.h"
#include "xcloc_enum.h"
#include <stdbool.h>
#include <complex.h>
/*
#include <fftw/fftw3.h>
#ifdef XLOC_USE_MPI
#include <mpi.h>
#endif
*/

struct xcfftFIR_struct
{
    void *pSpec;  /*!< IppsFIRSpec_32f FIR filter structure.  This has
                       dimension [pSpecSize]. */
    void *pBuf;   /*!< Ipp8u workspace buffer. This has dimension
                       [bufferSize]. */
    float *pTaps; /*!< Averaging coefficients in FIR filter. */
    int tapsLen;  /*!< Number of taps. */
    int bufferSize; /*!< Size of buffer. */
    int specSize;  /*!< Size of pSpec. */
    bool lfilter; /*!< If true then the filter will be applied. */
    char pad[3];
};

struct xcfft_struct
{
    struct xcfftFIR_struct fir; /*!< FIR filter for RMS window averaging. */
    void *forwardPlan;      /*!< Pointer to FFTW forward transform. */
    void *inversePlan;      /*!< Pointer to FFTW inverse transform. */
    //void *forwardHandlePtr; /*! Handle describing forward transform. */
    //void *inverseHandlePtr; /*! Handle describing inverse transform. */
    float complex *fts;   /*!< Fourier transforms of input signals.  This is an
                               array of dimension [ntfPts x nsignals] with
                               leading dimension ntfPts.  */
    float complex *xcfts; /*!< Fourier transforms of the cross-correlations. */
    float *x;       /*!< Signals to transform.  This is an array of dimension
                         [dataOffset x nsignals] with leading dimension
                         dataOffset. */
    float *y;       /*!< Holds the correlated data.  This is an array of 
                         dimension [dataOffset x ntfSignals] with leading 
                         dimension dataOffset. */
    float *yfilt;   /*!< Holds the filtered correlated data.  This is an array
                         of dimension [dataOffset x ntfSignals] with leading
                         dimension dataOffset. */ 
    int *xcPairs;   /*!< This is the table of unique cross-correlation pairs; e.g.,
                         for two stations this would like: {0, 1, 0, 2, 1, 2}
                         where station 0 is correleated with 1 (0, 1), and 
                         station 2 (0, 2), and station 1 is correlated with
                         station 2 (1, 2).  This is an array of dimension 
                         [2 x ntfSignals] with leading dimension. */
    int ntfSignals; /*!< Number of transformed signals.  This is a trianglular
                         number equal to (nsignals*(nsignals+1))/2 which holds
                         all of the correlations. */
    int nsignals;   /*!< Number of signals to transform. */
    int dataOffset; /*!< Data offset of input time domain and output time
                         domain cross-correlation signals. */
    int ftOffset;   /*!< Fourier transform offset. */
    int ntfPts;     /*!< Number of points in transform. */
    int npts;       /*!< Number of points in input signals. */
    int lxc;        /*!< Length of the cross-correlation.  This is an array
                         of dimension [2*npts-1]. */
    enum xclocPrecision_enum
         precision; /*!< Precision - FLOAT or DOUBLE. */
    enum xclocAccuracy_enum
         accuracy;  /*!< Controls the accuracy in the MKL computations. */
};

#ifdef __cplusplus
extern "C"
{
#endif

/* Applies the cross-correlation. */
int xcloc_xcfft_apply(const bool lphaseOnly, struct xcfft_struct *xcfft);
/* Compute the phase only cross-correlation. */
int dales_xcfft_computePhaseCorrelation(struct xcfft_struct *xcfft);
/* Compute the cross-correlation with the FFT. */
int dales_xcfft_computeFFTCrossCorrelation(struct xcfft_struct *xcfft);
/* Compute the transform pairs */
int dales_xcfft_computeXCTable(const bool lwantDiag,
                               const int nsignals, const int npairs,
                               int *__restrict__ xcPairs);
/* Filters the signal */
int dales_xcfft_createRMSFilter(const int winLen, struct xcfft_struct *xcfft);
 int dales_xcfft_computeXCRMSWindow(struct xcfft_struct *xcfft);
/* Frees memory on the cross-correlation structure. */
int xcloc_xcfft_finalize(struct xcfft_struct *xcfft);
/* Gets the signal'th signal from xcfft. */
int dales_xcfft_getData32f(const int npts, const int signal,
                           const struct xcfft_struct xcfft,
                           float *__restrict__ x);
/* Fourier transforms the input signals. */
int dales_xcfft_forwardTransform(struct xcfft_struct *xcfft);
int xcloc_xcfft_forwardTransform(struct xcfft_struct *xcfft);
/* Inverse Fourier transforms the cross-correlations. */
int dales_xcfft_inverseTransformXC(struct xcfft_struct *xcfft);
int xcloc_xcfft_inverseTransformXC(struct xcfft_struct *xcfft);
/* Compute frequency domain correlations. */
int xcloc_xcfft_computeCorrelationsWithFFTData(
    const bool lphaseOnly,
    const int ntfSignals,
    const int ntfPts,
    const int ftOffset,
    const enum xclocPrecision_enum precision,
    const int *__restrict__ xcPairs,
    const void *__restrict__ ftsIn, //const float complex *__restrict__ ftsIn,
    const void *__restrict__ xcftsOut); //float complex *__restrict__ xcftsOut);
/* Gets the cross-correlation. */
const float *dales_xcfft_getXCDataPointer32f(const int lxc, const int ixc,
                                             const struct xcfft_struct xcfft,
                                             int *ierr);
int xcloc_xcfft_getAllXCData(const int ntfSignals,
                             const int ldxc, const int lxc,
                             enum xclocPrecision_enum precision,
                             const struct xcfft_struct xcfft,
                             void *__restrict__ xcs);
int xcloc_xcfft_getXCData(const int lxc, const int ixc,
                          enum xclocPrecision_enum precision,
                          const struct xcfft_struct xcfft,
                          void *__restrict__ xcs);
int xcloc_xcfft_getXCData64f(const int lxc, const int ixc,
                             const struct xcfft_struct xcfft,
                             double *__restrict__ xc);
/* Initializes the cross-correlation. */
int xcloc_xcfft_initialize(const int npts, const int nptsPad,
                           const int nsignals, struct xcfft_struct *xcfft);
/* Sets all signals on xcfft. */
int xcloc_xcfft_setAllData(const int nsignals, const int lds,
                           const int npts,
                           const enum xclocPrecision_enum precision,
                           const void *__restrict__ x,
                           struct xcfft_struct *xcfft);
int dales_xcfft_setAllData32f(const int nsignals, const int lds,
                              const int npts, const float *__restrict__ x,
                              struct xcfft_struct *xcfft);
int dales_xcfft_setAllData64f(const int nsignals, const int lds,
                              const int npts, const double *__restrict__ x,
                              struct xcfft_struct *xcfft);
/* Sets the signal'th signal on xcfft. */
int xcloc_xcfft_setData(const int npts, const int is,
                        const enum xclocPrecision_enum precision,
                        const void *__restrict__ x,
                        struct xcfft_struct *xcfft);
int xcloc_xcfft_setData32f(const int npts, const int signal,
                              const float *__restrict__ x,
                              struct xcfft_struct *xcfft);
int xclo_xcfft_setData64f(const int npts, const int signal,
                          const double *__restrict__ x,
                              struct xcfft_struct *xcfft);


/* Internal function for creating an XC table. */
int dales_xcfft_createXCTable(struct xcfft_struct *xcfft);
/* Internal function to create Dfti FFT descriptors. */
int xcloc_xcfft_makeFFTWDescriptors(struct xcfft_struct *xcfft);
/* Internal function to allocate space. */
int xcloc_xcfft_setSpace(struct xcfft_struct *xcfft);

#ifdef __cplusplus
}
#endif
#endif /* XCLOC_XCFFT_H__ */
