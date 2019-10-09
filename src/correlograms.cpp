#include <cstdio>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#else
namespace
{
int omp_get_thread_num(){return 0;}
int omp_get_num_threads(){return 1;}
}
#endif
#include <vector>
#include <algorithm>
#include <cstring>
#include <cfloat>
#include <cassert>
#include <iterator>
#include <complex> // Put this before fftw
#include <mkl.h>
#include <fftw/fftw3.h>
#include <ipps.h>
#include "xcloc/correlograms.hpp"
#include "xcloc/correlogramParameters.hpp"
#include "xcloc/waveformIdentifier.hpp"
#include "rtseis/utilities/transforms/firEnvelope.hpp"

using namespace XCLoc;

namespace
{

MKL_INT64 convertAccuracyToMKL(
    const MKLFloatingPointAccuracy accuracy)
{
    if (accuracy == MKLFloatingPointAccuracy::HIGH_ACCURACY)
    {
        return VML_HA;
    }
    else if (accuracy == MKLFloatingPointAccuracy::LOW_ACCURACY) 
    {
        return VML_LA;
    }
    else if (accuracy == MKLFloatingPointAccuracy::EP_ACCURACY)
    {
        return VML_EP;
    }
    return VML_HA;
}

inline int padLength(const int n, 
                     const size_t precisionSize = sizeof(double),
                     const int alignment=64)
{
    auto size = static_cast<int> (precisionSize); //sizeof(double));
    int padLength = 0;
    auto xmod = (n*size)%alignment;
    if (xmod != 0){padLength = (alignment - xmod)/size;}
    auto nptsPadded = n + padLength;
    return nptsPadded;
}

inline int getInputWaveformID(const int64_t waveID,
                              const std::vector<int64_t> &signalMap,
                              const bool luseFastSignalMap)
{
    if (luseFastSignalMap)
    {
        if (waveID < 0 || waveID >= static_cast<int> (signalMap.size()))
        {
            throw std::invalid_argument("waveID = " + std::to_string(waveID)
                                      + " must be in range [0," 
                                      + std::to_string(signalMap.size()-1)
                                      + "]\n"); 
        }
        return waveID;
    }
    // Hunt for it
    auto low = std::lower_bound(signalMap.begin(), signalMap.end(), waveID);
    if (*low == waveID)
    {
        return std::distance(signalMap.begin(), low);
    }
    else
    {
        fprintf(stderr, "using linear search\n");
        int index =-1;
        for (int i=0; i<static_cast<int> (signalMap.size()); ++i)
        {
            if (signalMap[i] == waveID){return index;}
        }
        throw std::invalid_argument("WaveID = " + std::to_string(waveID)
                                  + " not in signal map\n");
    }
}

bool useFastSignalMap(const std::vector<int64_t> &signalMap)
{
    int nsignals = static_cast<int> (signalMap.size());
    // If the unique elements are 0, 1, 2, ..., nsignals - 1 then I don't need
    // to use a bisection search, i.e., I can use a 'fastSignalMap'
    bool lUseFastSignalMap = true;
    if (nsignals > 0)
    {
        // Need to start at 0
        if (signalMap[0] != 0)
        {
            lUseFastSignalMap = false;
        }
        else
        {
            // Need to assert this counts 0, 1, 2, ..., nsignals - 1
            for (int i=0; i<nsignals-1; ++i)
            {
                if (signalMap[i+1] != signalMap[i] + 1)
                {
                    lUseFastSignalMap = false;
                    break;
                }
            }
        }
    }
    return lUseFastSignalMap;
}

std::vector<int64_t>
createSignalMap(
    const std::vector<std::pair<WaveformIdentifier, WaveformIdentifier>> &xcPairs)
{
    auto npairs = static_cast<int> (xcPairs.size());
    // Extract the indices
    std::vector<std::pair<int64_t, int64_t>> xcPairIDs(npairs);
    for (int i=0; i<npairs; ++i)
    {
        auto id1 = xcPairs[i].first.getWaveformIdentifier();
        auto id2 = xcPairs[i].second.getWaveformIdentifier();
        xcPairIDs[i] = std::make_pair(id1, id2);
    }
    // Now sort
    std::vector<int64_t> work(2*npairs+1);
    for (int i=0; i<npairs; ++i)
    {
        work[2*i]   = xcPairIDs[i].first;
        work[2*i+1] = xcPairIDs[i].second;
    } 
    std::sort(work.begin(), work.begin() + 2*npairs);
    work[2*npairs] = work[0] - 1; // Trick to handle edge
    // Now create the signal map 
    std::vector<int64_t> signalMapWork;
    signalMapWork.reserve(2*npairs); 
    int nsignals = 0;
    for (int i=0; i<2*npairs; ++i)
    {
        int j;
        for (j=i+1; j<2*npairs+1; ++j)
        {
            // Mismatch is another signal
            if (work[i] != work[j])
            {
                signalMapWork.push_back(work[i]);
                nsignals = nsignals + 1;
                break;
            }
        }
        i = j; // Start next iteration at mismatch
    }
    // Copy the elements I need from the signal map
    std::vector<int64_t> signalMap(signalMapWork.begin(),
                                   signalMapWork.begin() + nsignals);
#ifdef DEBUG
    assert(std::is_sorted(signalMap.begin(), signalMap.end()));
#endif
    return signalMap;
}

std::vector<std::pair<int, int>>
createCorrelationIndexPairs(
    const std::vector<std::pair<WaveformIdentifier, WaveformIdentifier>> &xcPairs,
    const std::vector<int64_t> &signalMap,
    const bool luseFastSignalMap)
{
    int nxcs = static_cast<int> (xcPairs.size()); 
    std::vector<std::pair<int, int>> correlationIndexPairs(nxcs);
    for (int ixc=0; ixc<nxcs; ixc++)
    {
        auto i1 = xcPairs[ixc].first.getWaveformIdentifier(); 
        auto i2 = xcPairs[ixc].second.getWaveformIdentifier();
        // From the XC pairs get the input signal indices
        int inputSignal1 = getInputWaveformID(i1, //xcPairs[ixc].first,
                                              signalMap,
                                              luseFastSignalMap);
        int inputSignal2 = getInputWaveformID(i2, //xcPairs[ixc].second,
                                              signalMap,
                                              luseFastSignalMap);
        // Make the input signal indices a pair for internal usage
        correlationIndexPairs[ixc] = std::make_pair(inputSignal1, inputSignal2);
    }
    return correlationIndexPairs;
}

template<typename T>
inline void checkInputSignalAndLength(const int nSamples,
                                      const int nSamplesRef,
                                      const T x[])
{
    if (nSamples != nSamplesRef)
    {
        throw std::invalid_argument("nSamples = " + std::to_string(nSamples)
                                  + " must equal + "
                                  + std::to_string(nSamplesRef) + "\n");
    }
    if (x == nullptr){throw std::invalid_argument("x is NULL\n");} 
}

}

//============================================================================//
//                              End Private Functions                         //
//============================================================================//

template<>
class Correlograms<double>::CorrelogramsImpl
{
public:
    /// Destructor
    ~CorrelogramsImpl()
    {
        clear();
    }
    /// Releases memory
    void clear() noexcept
    {
        mEnvelopes.clear();
        mParms.clear();
        mSignalMap.clear();
        mCorrelationIndexPairs.clear();
        if (mHavePlan){fftw_destroy_plan(mForwardPlan);}
        if (mHavePlan){fftw_destroy_plan(mInversePlan);}
        if (mInputData){MKL_free(mInputData);}
        if (mWork){MKL_free(mWork);}
        if (mProcessedCorrelograms && mFilterCorrelograms)
        {
            MKL_free(mProcessedCorrelograms);
        }
        if (mCorrelogramSpectra){MKL_free(mCorrelogramSpectra);}
        mInputSpectra = nullptr;
        mCorrelogramSpectra = nullptr;
        mProcessedCorrelograms = nullptr;
        mRawOutputCorrelograms = nullptr;
        mWork = nullptr;
        mAccuracyMode = VML_HA;
        mWorkSizeInBytes = 0;
        mNumberOfSignals = 0;
        mSamples = 0;
        mPaddedInputLength = 0;
        mNumberOfFrequencies = 0;
        mNumberOfCorrelograms = 0;
        mSamplesInCorrelogram = 0;
        mInputDataLeadingDimension = 0;
        mSpectraLeadingDimension = 0;
        mCorrelogramLeadingDimension = 0;
        mUseFastSignalMap = false;
        mHavePlan = false;
        mHaveCorrelograms = false;
        mFilterCorrelograms = false;
        mInitialized = false;
    }
    /// Initializes the class
    void initialize()
    {
        // Figure out the desired accuracy
        mAccuracyMode 
            = convertAccuracyToMKL(mParms.getMKLFloatingPointAccuracy());
        auto xcPairs = mParms.getCorrelationPairs();
        // Create the signal map
        mSignalMap = createSignalMap(xcPairs);
        mNumberOfSignals = static_cast<int> (mSignalMap.size());
        mUseFastSignalMap = useFastSignalMap(mSignalMap);
        mCorrelationIndexPairs = createCorrelationIndexPairs(xcPairs,
                                                             mSignalMap,
                                                             mUseFastSignalMap);
        // Get the number of samples and compute the correlogram length
        mSamples = mParms.getNumberOfSamples();
        mPaddedInputLength = mParms.getNumberOfPaddedSamples();
        mSamplesInCorrelogram = mParms.getCorrelogramLength(); //2*mPaddedInputLength - 1;
        mNumberOfFrequencies = mSamplesInCorrelogram/2 + 1;
        mNumberOfCorrelograms = mParms.getNumberOfCorrelationPairs();
        // Filter correlograms?
        mFilterCorrelograms = false;
        if (mParms.getFilteringType() != CorrelogramFilteringType::NO_FILTERING)
        {
            mFilterCorrelograms = true;
        }
        if (mParms.getFilteringType() ==
           CorrelogramFilteringType::FIR_ENVELOPE_FILTERING)
        {
            mEnvelopes.resize(omp_get_num_threads());
            for (int i=0; i<static_cast<int> (mEnvelopes.size()); ++i)
            {
                mEnvelopes[i].initialize(mParms.getFIREnvelopeFilterLength(),
                                         RTSeis::POST_PROCESSING);
            }
        }
        // Figure out the padding sizes.  Note, to perform the correlation we
        // need the input signals evaluated at the transform frequencies of
        // the correlograms.  Hence, the is significant padding.
        mInputDataLeadingDimension
            = padLength(mSamplesInCorrelogram, sizeof(double), 64); //PaddedInputLength, sizeof(double), 64);
        mCorrelogramLeadingDimension
            = padLength(mSamplesInCorrelogram, sizeof(double), 64);
        mSpectraLeadingDimension
            = padLength(mNumberOfFrequencies, sizeof(fftw_complex), 64);
        // Allocate space for the input signals
        size_t nbytes = static_cast<size_t> (mInputDataLeadingDimension)
                       *static_cast<size_t> (mNumberOfSignals)
                       *sizeof(double);
        // Apple sits on 1 trillion in market capital but cant implement C11.
/*
#ifdef __APPLE__
        void *mInputDataPtr = mInputData;
        posix_memalign(&mInputDataPtr, 64, nbytes);
#else
        mInputData = aligned_alloc(64, nbytes);
#endif
*/
        mInputData = static_cast<double *> (MKL_calloc(nbytes, 1, 64));
        // Allocate space for the workspace
        size_t work1 = static_cast<size_t> (mCorrelogramLeadingDimension)
                      *static_cast<size_t> (mNumberOfCorrelograms)
                      *sizeof(double);
        size_t work2 = static_cast<size_t> (mSpectraLeadingDimension)
                      *static_cast<size_t> (mNumberOfSignals)
                      *sizeof(std::complex<double>);
        mWorkSizeInBytes = std::max(work1, work2);
        mWork = MKL_calloc(mWorkSizeInBytes, 1, 64);
        // Set input spectra and output correlograms to point at workspace
        mInputSpectra = reinterpret_cast<fftw_complex *> (mWork);
        mRawOutputCorrelograms = reinterpret_cast<double *> (mWork);
        // Set the processed output correlograms
        if (mFilterCorrelograms)
        {
            size_t nbytes = work1;
            mProcessedCorrelograms
               = static_cast<double *> (MKL_calloc(nbytes, 1, 64)); 
        }
        else
        {
            mProcessedCorrelograms = mRawOutputCorrelograms;
        }
        // Allocate space for frequency domain correlograms
        nbytes = static_cast<size_t> (mNumberOfCorrelograms)
                *static_cast<size_t> (mSpectraLeadingDimension)
                *sizeof(fftw_complex);
        mCorrelogramSpectra
            = reinterpret_cast<fftw_complex *> (MKL_calloc(nbytes, 1, 64));
        // Initialize the FFTw
        constexpr int rank = 1;
        constexpr int istride = 1;
        constexpr int ostride = 1;
        constexpr int inembed[1] = {0};
        constexpr int onembed[1] = {0};
        // Forward transform input signals to correlogram frequencies
        int nf[1] = {mSamplesInCorrelogram};
        // Inverse transform correlogram frequencies to full correlograms
        int ni[1] = {mSamplesInCorrelogram};
        // Create the plans
        mForwardPlan = fftw_plan_many_dft_r2c(rank, nf, mNumberOfSignals,
                                              mInputData, inembed,
                                              istride, mInputDataLeadingDimension,
                                              mInputSpectra, onembed,
                                              ostride, mSpectraLeadingDimension,
                                              FFTW_PATIENT);
        mInversePlan = fftw_plan_many_dft_c2r(rank, ni, mNumberOfCorrelograms,
                                              mCorrelogramSpectra, inembed,
                                              istride, mSpectraLeadingDimension,
                                              mRawOutputCorrelograms, onembed,
                                              ostride, mCorrelogramLeadingDimension,
                                              FFTW_PATIENT);
        mHavePlan = true;
        mHaveCorrelograms = false;
        mInitialized = true;
    }
    // Computes the cross-correlgorams 
    void computeCorrelograms(const bool mCrossCorrelate)
    {
        // First Fourier transform the input signals
        fftw_execute_dft_r2c(mForwardPlan, mInputData, mInputSpectra);
        #pragma omp parallel default(none)
        {
        // Set workspace
        auto amplitudeSpectra
             = static_cast<double *> (ippsMalloc_64f(mNumberOfFrequencies));
        // Now perform the correlation 
        #pragma omp for
        for (int ixc=0; ixc<mNumberOfCorrelograms; ++ixc)
        {
            int signal1 = mCorrelationIndexPairs[ixc].first;
            int signal2 = mCorrelationIndexPairs[ixc].second;
            auto i1 = signal1*mSpectraLeadingDimension;
            auto i2 = signal2*mSpectraLeadingDimension;
            auto j1 = ixc*mSpectraLeadingDimension;
            auto ft1 = reinterpret_cast<MKL_Complex16 *> (mInputSpectra[i1]);
            auto ft2 = reinterpret_cast<MKL_Complex16 *> (mInputSpectra[i2]);
            auto __attribute__((aligned(64))) xcPtr
                = reinterpret_cast<MKL_Complex16 *> (mCorrelogramSpectra[j1]);
            // Correlate ft1*conj(ft2)
            vmzMulByConj(mNumberOfFrequencies, ft1, ft2, xcPtr, mAccuracyMode);
            // Do phase correlograms
            if (!mCrossCorrelate)
            {
                // Compute absolute value and avoid division by zero
                vmzAbs(mNumberOfFrequencies, xcPtr, amplitudeSpectra,
                       mAccuracyMode);
                // The numerator is close to zero so this will be 
                // substantially less than 1.
                constexpr double small = DBL_EPSILON*1.e2;
                ippsThreshold_LT_64f_I(amplitudeSpectra, mNumberOfFrequencies,
                                       small);
                // Now normalize
                auto __attribute__((aligned(64))) xcPtr2
                    = reinterpret_cast<std::complex<double> *> (xcPtr);
                #pragma omp simd
                for (int i=0; i<mNumberOfFrequencies; ++i)
                {
                    xcPtr2[i] = xcPtr2[i]/amplitudeSpectra[i];
                }
            }
        }
        ippsFree(amplitudeSpectra);
        } // End parallel
        // Inverse transform
        fftw_execute_dft_c2r(mInversePlan, mCorrelogramSpectra,
                             mRawOutputCorrelograms);
        #pragma omp parallel default(none)
        {
        auto xcTemp
            = static_cast<double *> (ippsMalloc_64f(mSamplesInCorrelogram));
        // Shuffle the correlograms to obtain causal and acausal part and 
        // normalize
        auto __attribute__((aligned(64))) xnorm
             = 1./static_cast<double> (mSamplesInCorrelogram);
        int ncopy1 = mSamplesInCorrelogram/2;
        int ncopy2 = mSamplesInCorrelogram - ncopy1;
        #pragma omp for
        for (int ixc=0; ixc<mNumberOfCorrelograms; ++ixc)
        {
            int j1 = ixc*mCorrelogramLeadingDimension;
            // Copy and scale
            auto __attribute__((aligned(64)))
                xcPtr = &mRawOutputCorrelograms[j1];
            ippsMulC_64f(xcPtr, xnorm, xcTemp, mSamplesInCorrelogram);
            // Shuffle
            //ippsCopy_64f(xcTemp, &xcPtr[ncopy2], ncopy1);
            //ippsCopy_64f(&xcTEmp[ncopy1], xcPtr, ncopy2);
            std::copy(xcTemp, xcTemp+ncopy2, &xcPtr[ncopy1]);
            std::copy(&xcTemp[ncopy2], &xcTemp[ncopy2]+ncopy1, xcPtr);
        }
        ippsFree(xcTemp);
        } // End parallel
        // Compute envelopes?
        if (mParms.getFilteringType() ==
            CorrelogramFilteringType::FIR_ENVELOPE_FILTERING)
        {
            #pragma omp parallel for
            for (int ixc=0; ixc<mNumberOfCorrelograms; ++ixc)
            {
                int threadID = omp_get_thread_num();
                int j1 = ixc*mCorrelogramLeadingDimension;
                auto xcPtr = &mRawOutputCorrelograms[j1];
                auto ycPtr = &mProcessedCorrelograms[j1];
                mEnvelopes[threadID].transform(mSamplesInCorrelogram,
                                               xcPtr, &ycPtr);
            }
        }
        mHaveCorrelograms = true;
    }
//private:
    /// FIR envelope filter
    std::vector<RTSeis::Utilities::Transforms::FIREnvelope<double>> mEnvelopes;
    /// Correlogram parameters
    CorrelogramParameters mParms;
    /// A pointer to the input data that was transformed to the frequency 
    /// domain This is a row major matrix whose dimensions are
    /// [mNumberOfSignals x mSpectraLeadingDimension].  Note, it shares the same
    /// space as mRawOutputCorrelgorams.
    fftw_complex *mInputSpectra = nullptr;
    /// Holds the spectra of the correlograms.  This is a row major matrix
    /// whose dimensions are [mNumberOfCorrelograms x mSpectraLeadingDimension]
    fftw_complex *mCorrelogramSpectra = nullptr;
    /// A pointer to the output correlograms.  This is a row major matrix
    /// whose dimensions are
    /// [mNumberOfCorrelograms x mCorrelogramLeadingDimension].
    /// Note, it shares the same space as mInputSpectra.
    double *mRawOutputCorrelograms = nullptr;    
    /// Holds the processed correlograms.  If no processing is required then
    /// this points back to the raw correlograms otherwise it is alllocated.
    /// Regardless, it has dimension 
    /// [mNumberOfCorrelograms x mCorrelogramLeadingDimension]. 
    double *mProcessedCorrelograms = nullptr;
    /// Holds the correlograms in the in the frequency domain.  These will
    /// be transformed back to the time domain.
    /// This performs double duty and is a chunk of memory that holds the
    /// time domain correlograms as well as the intermediate Fourier transforms
    /// of the input data.  
    void *mWork = nullptr;
    /// Holds the input data.  This is a row major matrix whose dimensions are 
    /// [mNumberOfSignals x mInputDataLeadingDimension].
    double *mInputData = nullptr;
    /// FFTw plan to go from input signals to frequency domain
    fftw_plan mForwardPlan;
    /// FFTw plan to go from frequency domain correlograms to time domain
    fftw_plan mInversePlan;
    /// Used in bisection search to map from a signal index to the
    /// local signal storage.
    std::vector<int64_t> mSignalMap;
    /// This is used internally to extract the (i,j)'th input signal 
    /// indices comprising a correlation pair.  This has dimension
    /// [mNumberOfCorrelograms].
    std::vector<std::pair<int, int>> mCorrelationIndexPairs;
    /// The desired accuracy for MKL
    MKL_INT16 mAccuracyMode = VML_HA;
    /// The size of workspace in bytes.
    size_t mWorkSizeInBytes = 0;
    /// Total number of input signals
    int mNumberOfSignals = 0;
    /// Number of samples in input signals
    int mSamples = 0;
    /// Number of 0 padded samples in input signals.
    int mPaddedInputLength = 0;
    /// The number of frequencies in the correlation.
    int mNumberOfFrequencies = 0;
    /// Total number of correlograms.
    int mNumberOfCorrelograms = 0;
    /// Number of samples in correlograms
    int mSamplesInCorrelogram = 0;
    /// Leading dimension of the input signals
    int mInputDataLeadingDimension = 0;
    /// Leading dimension of the signals in frequency domain signals.
    int mSpectraLeadingDimension = 0;
    /// Leading dimension of output time-domain correlograms
    int mCorrelogramLeadingDimension = 0;
    /// If true then I can use a fast signal map for mapping signal ID to 
    /// the local unique signal ID.
    bool mUseFastSignalMap = false;
    /// Flag indicating the plane was created
    bool mHavePlan = false;
    /// Flag indicating that the correlograms have been computed
    bool mHaveCorrelograms = false;
    /// Flag indicating that we are to filter the correlograms
    bool mFilterCorrelograms = false;
    /// Flag indicating that the class is initialized
    bool mInitialized = false;
};

template<>
class Correlograms<float>::CorrelogramsImpl
{
public:

//private:
};

//============================================================================//
//                              End the implementations                       //
//============================================================================//

/// Constructors
template<class T>
Correlograms<T>::Correlograms() :
    pImpl(std::make_unique<CorrelogramsImpl> ()) 
{
}

template<class T>
Correlograms<T>::Correlograms(
    const CorrelogramParameters &parameters) :
    pImpl(std::make_unique<CorrelogramsImpl> ())
{
    initialize(parameters);
}

/*
template<class T>
Correlograms<T>::Correlograms(
    const Correlograms &engine)
{
    *this = engine;
}
*/

/// Move constructor
template<class T>
Correlograms<T>::Correlograms(Correlograms &&engine) noexcept
{
    *this = std::move(engine);
}

/// Operators
template<class T> Correlograms<T>&
Correlograms<T>::operator=(Correlograms &&engine) noexcept
{
    if (&engine == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(engine.pImpl);
    return *this;
}

/// Destructor
template<class T>
Correlograms<T>::~Correlograms() = default;

template<class T>
void Correlograms<T>::clear() noexcept
{
    pImpl->clear();
}

/// Initializer
template<class T>
void Correlograms<T>::initialize(const CorrelogramParameters &parameters)
{
    if (!parameters.isValid())
    {
        throw std::invalid_argument("Parameters are invalid\n");
    }
    // Copy the parameters
    pImpl->mParms = parameters;
    // And initialize
    pImpl->initialize();
}

/// Set the input signal - double to double
template<>
void Correlograms<double>::setInputSignal(
    const int64_t waveID, const int nSamples, const double *__restrict__ x)
{
    pImpl->mHaveCorrelograms = false;
    // Verify class is initialized
    if (!isInitialized())
    {
        throw std::runtime_error("Class not yet initialized\n");
    }
    // Check nSamples and x are okay - this will throw
    checkInputSignalAndLength(nSamples, pImpl->mSamples, x);
    // Get the signal map - this will throw on waveID
    auto signalIndex = getInputWaveformID(waveID,
                                          pImpl->mSignalMap,
                                          pImpl->mUseFastSignalMap);
    // Copy
    auto index = static_cast<size_t> (signalIndex)
                *static_cast<size_t> (pImpl->mInputDataLeadingDimension);
    std::copy(x, x + nSamples, &pImpl->mInputData[index]);
}

/// Set the input signal - float to double
template<>
void Correlograms<double>::setInputSignal(
    const int64_t waveID, const int nSamples, const float x[])
{
    pImpl->mHaveCorrelograms = false;
    // Verify class is initialized
    if (!isInitialized())
    {
        throw std::runtime_error("Class not yet initialized\n");
    }
    // Check nSamples and x are okay - this will throw
    checkInputSignalAndLength(nSamples, pImpl->mSamples, x);
    // Get the signal map - this will throw on waveID
    auto signalIndex = getInputWaveformID(waveID,
                                          pImpl->mSignalMap,
                                          pImpl->mUseFastSignalMap);
    // Convert and copy
    auto index = static_cast<size_t> (signalIndex)
                *static_cast<size_t> (pImpl->mInputDataLeadingDimension);
    ippsConvert_32f64f(x, &pImpl->mInputData[index], pImpl->mSamples);
}

/// Set an input signal to zero
template<class T>
void Correlograms<T>::zeroInputSignal(const int waveID)
{
    pImpl->mHaveCorrelograms = false;
    // Verify class is initialized
    if (!isInitialized())
    {
        throw std::runtime_error("Class not yet initialized\n");
    }
    // Get the signal map - this will throw on waveID
    auto signalIndex = getInputWaveformID(waveID,
                                          pImpl->mSignalMap,
                                          pImpl->mUseFastSignalMap);
    // Copy
    auto index = static_cast<size_t> (signalIndex)
                *static_cast<size_t> (pImpl->mInputDataLeadingDimension);
    auto nzero = static_cast<size_t> (pImpl->mSamples)*sizeof(T);
    std::memset(&pImpl->mInputData[index], 0, nzero);
}

/// Compute cross correlograms
template<class T>
void Correlograms<T>::computeCrossCorrelograms()
{
    constexpr bool mCrossCorrelate = true;
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    pImpl->computeCorrelograms(mCrossCorrelate);
}

/// Compute phase correlograms
template<class T>
void Correlograms<T>::computePhaseCorrelograms()
{
    constexpr bool mCrossCorrelate = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    pImpl->computeCorrelograms(mCrossCorrelate);
}

/// Gets the signal pairs comprising the correlogram 
template<class T>
std::pair<WaveformIdentifier, WaveformIdentifier>
Correlograms<T>::getCorrelationPair(const int ixc) const
{
    auto nxc = getNumberOfCorrelograms(); // Throws on initialization
    // Check that ixc is in range
    if (ixc < 0 || ixc >= nxc)
    {
        throw std::invalid_argument("ixc = " + std::to_string(ixc)
                                  + " must be in range [0,"
                                  + std::to_string(nxc-1) + "]\n");
    }
    return pImpl->mParms.getCorrelationPair(ixc);
}

/// Get a correlogram
template<class T>
const T* Correlograms<T>::getRawCorrelogramPointer(const int ixc) const
{
    auto nxc = getNumberOfCorrelograms(); // Throws on initialization
    if (!haveCorrelograms())
    {
        throw std::runtime_error("Correlograms not yet computed\n");
    }
    // Check that ixc is in range
    if (ixc < 0 || ixc >= nxc)
    {
        throw std::invalid_argument("ixc = " + std::to_string(ixc)
                                  + " must be in range [0,"
                                  + std::to_string(nxc-1) + "]\n");
    }
    int index = static_cast<size_t> (ixc)
               *static_cast<size_t> (pImpl->mCorrelogramLeadingDimension);
    const T *xcPtr = &pImpl->mRawOutputCorrelograms[index];
    return xcPtr;
}

template<class T>
const T* Correlograms<T>::getProcessedCorrelogramPointer(
    const int ixc) const
{
    auto nxc = getNumberOfCorrelograms(); // Throws on initialization
    if (!haveCorrelograms())
    {
        throw std::runtime_error("Correlograms not yet computed\n");
    }
    // Check that ixc is in range
    if (ixc < 0 || ixc >= nxc)
    {
        throw std::invalid_argument("ixc = " + std::to_string(ixc)
                                  + " must be in range [0,"
                                  + std::to_string(nxc-1) + "]\n");
    }
    int index = static_cast<size_t> (ixc)
               *static_cast<size_t> (pImpl->mCorrelogramLeadingDimension);
    const T *xcPtr = &pImpl->mProcessedCorrelograms[index];
    return xcPtr;
}

template<class T>
void Correlograms<T>::getRawCorrelogram(
    const int ixc, const int nwork, T *xcIn[]) const
{
    // Get pointer to correlogram.  This throws an initialization error
    // and checks if ixc is valid
    auto xcPtr __attribute__((aligned(64))) = getRawCorrelogramPointer(ixc);
    // Get the correlogram length and check nwork and xcIn
    int lxc = getCorrelogramLength();
    T *xc = *xcIn;
    if (nwork < lxc || xc == nullptr)
    {
        if (nwork < lxc)
        {
            throw std::invalid_argument("nwork = " + std::to_string(nwork)
                                      + " must be at least = " 
                                      + std::to_string(lxc));
        }
        throw std::invalid_argument("xcIn is NULL\n"); 
    }
    // Finally perform the copy
    std::copy(xcPtr, xcPtr+nwork, xc);    
}

template<class T>
void Correlograms<T>::getProcessedCorrelogram(
    const int ixc, const int nwork, T *xcIn[]) const
{
    // Get pointer to correlogram.  This throws an initialization error
    // and checks if ixc is valid
    auto xcPtr __attribute__((aligned(64)))
        = getProcessedCorrelogramPointer(ixc);
    // Get the correlogram length and check nwork and xcIn
    int lxc = getCorrelogramLength();
    T *xc = *xcIn;
    if (nwork < lxc || xc == nullptr)
    {
        if (nwork < lxc)
        {
            throw std::invalid_argument("nwork = " + std::to_string(nwork)
                                      + " must be at least = "
                                      + std::to_string(lxc));
        }
        throw std::invalid_argument("xcIn is NULL\n");
    }
    // Finally perform the copy
    std::copy(xcPtr, xcPtr+nwork, xc);
}

/// Initialized?
template<class T>
bool Correlograms<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Checks if the correlograms have been computed 
template<class T>
bool Correlograms<T>::haveCorrelograms() const noexcept
{
    return pImpl->mHaveCorrelograms;
}

/// Number of samples in input signals
template<class T>
int Correlograms<T>::getInputSignalLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    return pImpl->mSamples;
}

/// Number of correlograms
template<class T>
int Correlograms<T>::getNumberOfCorrelograms() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    return pImpl->mNumberOfCorrelograms;
}

/// Number of samples in correlogram
template<class T>
int Correlograms<T>::getCorrelogramLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    return pImpl->mSamplesInCorrelogram;
}

/// Template instantiation
template class XCLoc::Correlograms<double>;
//template class XCLoc::Correlograms<float>;

