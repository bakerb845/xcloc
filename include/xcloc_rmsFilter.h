#ifndef XCLOC_RMS_FILTER_H__
#define XCLOC_RMS_FILTER_H__
#include "xcloc_config.h"
#include "xcloc_enum.h"

struct xcfftRMSFilter_struct
{
    void *pSpec;    /*!< IppsFIRSpec_64f or IppsFIRSpec_32f FIR filter
                         structure.  This has dimension [pSpecSize]. */
    void *pBuf;     /*!< Ipp8u workspace buffer. This has dimension
                        [bufferSize]. */
    void *pTaps;    /*!< Ipp64f or Ipp32f Averaging coefficients in FIR
                         filter. */
    int tapsLen;    /*!< Number of taps. */
    int bufferSize; /*!< Size of buffer. */
    int specSize;   /*!< Size of pSpec. */
    int npts;       /*!< Number of points to filter. */
    enum xclocPrecision_enum precision; /*!< Precision of filter. */
    bool lfilter;   /*!< If true then the filter will be applied. */
    char pad[3];
};

#ifdef __cplusplus
extern "C"
{
#endif
/* Initialize the RMS filter. */
int xcloc_rmsFilter_initialize(const int winLen, const int npts,
                               const enum xclocPrecision_enum precision,
                               struct xcfftRMSFilter_struct *rms);
/* Apply the RMS filter. */
int xcloc_rmsFilter_apply(const int nsignals,
                          const int lds, const int npts, 
                          const enum xclocPrecision_enum precision,
                          struct xcfftRMSFilter_struct *rms,
                          const void *__restrict__ x,
                          void *__restrict__ xfilt);
/* Finalize the RMS filter. */
int xcloc_rmsFilter_finalize(struct xcfftRMSFilter_struct *rms);

#ifdef __cplusplus
}
#endif
#endif
