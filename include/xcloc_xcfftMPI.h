#ifndef XCLOC_XCFFTMPI_H__
#define XCLOC_XCFFTMPI_H__
#include "xcloc_config.h"
#include "xcloc_enum.h"
#include "xcloc_xcfft.h"
#ifdef XCLOC_USE_MPI
#include <mpi.h>

struct xcfftMPI_struct
{
    struct xcfft_struct sigFwd; /*!< Forward transforms of input signals. */
    struct xcfft_struct xcInv;  /*!< Fourier domain and inverse transforms
                                     of cross-correlated signals. */
    float complex *fts; /*!< All Fourier transforms of signals.  This is
                             an array of dimension [localSizeOfFts] and
                             whose memory has been allocated by MPI. */
    MPI_Comm comm;  /*!< MPI communicator. */
    MPI_Win tfWindow; /*!< RMA window defining the signal transforms
                           one-sided communication. */
//    int *tfRankID;
    /*!< Maps from the transform of the signal'th
                           signal to the target rank on the communicator.
                           This is an array of dimension
                           [tfRankIDPtr[nsignalsLoc+1]. */
//    int *tfDisp;
      /*!< Maps from the transform of the signal'th
                           signal to the target rank displacement.
                           This is an array of dimension
                           [tfRankIDPtr[nsignalsLoc+1]. */
//    int *tfRankIDPtr;
 /*!< Maps from the signal'th signal to the start index
                           of tfRankID.  This is an array of dimension
                           [nsignalsLoc+1]. */
    int *proc2ftPtr;    /*!< The ip'th process owns global signal Fourier
                             transform proc2tfPtr[ip]:proc2tfPtr[ip+1]-1.
                             This is an array of dimension [nprocs+1]. */
    int *proc2xcPtr;    /*!< The ip'th process owns global cross-correlations
                             proc2xcPtr[ip]:proc2xcptr[ip+1]-1.  This is
                             an array of dimension [nprocs+1]. */
    int *getFwdOffset;  /*!< When getting a signal from getFwdTarget this
                             corresponds to the local forward transform
                             number on getFwdTarget.  This is an array of
                             dimension [ngetFwdSignals]. */ 
    int *getFwdTarget;  /*!< Process containing the forward transform
                             signal to get.  This is an array of dimension
                             [ngetFwdSignals]. */
    int *xcrPairs;      /*!< Maps from the ixc'th transform pair to the
                             i'th and j'th global waveform that constitutes
                             this cross-correlation and the rank that has
                             compute this transform pair.  This is an array of
                             dimension [3 x ntfSignals] with leading dimension
                             3. */
    int ngetFwdSignals; /*!< Number of forward signals to get. */
    int nsignalsLoc;/*!< Number of local signals. */
    int nsignals;   /*!< Total number of signals. */
    int ntfSignalsLoc;
    int ntfSignals;
    int lxc;        /*!< Length of the cross-correlations. */
    int rank;       /*!< Rank of process. */
    int root;       /*!< Rank (ID) of root process.  This is likely 0. */
    int nprocs;     /*!< Number of processes on communicator. */
    //size_t localSizeOfFTs; /*!< Size of fts. */ 
    enum xclocPrecision_enum
         precision; /*!< Precision - FLOAT or DOUBLE. */
};

#ifdef __cplusplus
extern "C"
{
#endif

/* Initializes the memory on the xcfftMPI structure */
int xcloc_xcfftMPI_initialize(const int nptsIn, const int nptsPadIn,
                              const int nsignalsIn,
                              const int ntfSignalsIn,
                              const MPI_Comm comm, const int master,
                              const int xcPairsIn[],
                              struct xcfftMPI_struct *xcfftMPI);
/* Deallocates memory on the xcfftMPI structure */
int xcloc_xcfftMPI_finalize(struct xcfftMPI_struct *xcfftMPI);
/* Computes the phase correlations. */
int xcloc_xcfftMPI_computePhaseCorrelation(struct xcfftMPI_struct *xcfftMPI);
/* Gathers the data. */
int xcloc_xcfftMPI_gatherXCs(const int root, const int ntfSignals,
                             const int ldxc, const int lxc, 
                             const MPI_Datatype recvType,
                             const struct xcfftMPI_struct xcfftMPI,
                             void *__restrict__ x);
/* Scatters the data. */
int xcloc_xcfftMPI_scatterData(
    const int root, const int nsignals,
    const int lds, const int npts, 
    const MPI_Datatype sendType,
    const void *__restrict__ x,
    struct xcfftMPI_struct *xcfftMPI);

/*----------------------------------------------------------------------------*/
/*                          Internal Functions                                */
/*----------------------------------------------------------------------------*/
/* Gets the forward transforms onto the local cross-correlatoion FT
   structure. */
int xcloc_xcfftMPI_getForwardTransforms(struct xcfftMPI_struct *xcfftMPI);

#ifdef __cplusplus
}
#endif
#endif /* XCLOC_USE_MPI */
#endif /* XCLOC_XCFFTMPI_H__ */
