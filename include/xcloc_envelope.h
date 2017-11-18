#ifndef XCLOC_ENVELOPE_H__
#define XCLOC_ENVELOPE_H__ 1
#include <stdbool.h>
#include "xcloc_config.h"
#include "xcloc_enum.h"
#ifdef XCLOC_USE_MPI
#include <mpi.h>
#endif

struct xclocEnvelope_struct
{
    void *hfiltR;             /*!< Real coefficients of Hilbert FIR 
                                   transform.  This is an Ipps32f or
                                   Ipps64f array of dimension [tapsLne]. */
    void *hfiltI;             /*!< Imaginary coefficients of Hilbert FIR
                                   transform.  This is an Ipps32f or 
                                   Ipps64f array of dimension [tapsLen]. */
    int tapsLen;              /*!< Number of taps in FIR filters. */
    int specSize;             /*!< Size of workspace for FIR state. */
    int bufferSize;           /*!< Size of workspace for filter. */
    enum xclocAccuracy_enum
         accuracy;            /*!< Controls numerical precision in envelope
                                   computation. */
    enum xclocPrecision_enum
         precision;           /*!< Precision of filter.  This can be
                                    single or double precision. */
          
    bool ltype4;              /*!< If true then the FIR transform is of
                                   Type 4 (tapsLen is odd).  This results in
                                   a lot of sparseness. */
    bool linit;               /*!< If true then the filter was initialized. */
};

#ifdef __cplusplus
extern "C"
{
#endif

/*! Initialize the envelope structure. */
int xcloc_envelope_initialize(const int n,
                              const enum xclocPrecision_enum precision,
                              const enum xclocAccuracy_enum accuracy,
                              struct xclocEnvelope_struct *envelope);
/*! Apply the envelope. */
int xcloc_envelope_apply(const int nsignals,
                          const int lds, const int npts, 
                          const enum xclocPrecision_enum precision,
                          struct xclocEnvelope_struct *envelope,
                          const void *__restrict__ x,
                          void *__restrict__ xfilt);
#ifdef XCLOC_USE_MPI
int xcloc_envelope_applyMPI(const MPI_Comm comm, const int root,
                            const int nsignalsIn, const int ldsIn,
                            const int nptsIn,
                            const MPI_Datatype dataType,
                            struct xclocEnvelope_struct *envelope,
                            const void *__restrict__ x,
                            void *__restrict__ xfilt);
#endif
/*! Finalizes the envelope structure. */
int xcloc_envelope_finalize(struct xclocEnvelope_struct *envelope);
#ifdef __cplusplus
}
#endif 
#endif
