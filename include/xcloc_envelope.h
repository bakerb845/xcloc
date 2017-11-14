#ifndef XCLOC_ENVELOPE_H__
#define XCLOC_ENVELOPE_H__ 1
#include <stdbool.h>
#include "xcloc_config.h"
#include "xcloc_enum.h"

struct xclocEnvelope_struct
{
    void *hfiltR;
    void *hfiltI;
    int tapsLen;
    int specSize;
    int bufferSize;
    enum xclocPrecision_enum precision;
    bool linit;
};

#ifdef __cplusplus
extern "C"
{
#endif

/*! Initialize the envelope structure. */
int xcloc_envelope_initialize(const int n,
                              const enum xclocPrecision_enum precision,
                              struct xclocEnvelope_struct *envelope);
/*! Apply the envelope. */
int xcloc_envelope_apply(const int nsignals,
                          const int lds, const int npts, 
                          const enum xclocPrecision_enum precision,
                          struct xclocEnvelope_struct envelope,
                          const void *__restrict__ x,
                          void *__restrict__ xfilt);

#ifdef __cplusplus
}
#endif 
#endif
