#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include <functional>
#include <algorithm>
#include <random>
#include <ipps.h>
#include <rtseis/utilities/transforms/dftRealToComplex.hpp>
#include <rtseis/utilities/transforms/utilities.hpp>
#include "xcloc/correlograms.hpp"
#include "xcloc/correlogramParameters.hpp"
#include <gtest/gtest.h>

namespace
{
std::vector<double>
computeReferenceCorrelogram(const int npts,
                            const double signal1[],
                            const double signal2[],
                            const bool doPhase = false);

using namespace XCLoc;

TEST(testCorrelograms, parameters)
{
    CorrelogramParameters parameters;

    EXPECT_NO_THROW(parameters.setNumberOfSamples(100));
    EXPECT_EQ(parameters.getNumberOfSamples(), 100);
    EXPECT_EQ(parameters.getNumberOfPaddedSamples(), 100);

    EXPECT_NO_THROW(parameters.setNumberOfPaddedSamples(120));
    EXPECT_EQ(parameters.getNumberOfPaddedSamples(), 120);
    EXPECT_EQ(parameters.getCorrelogramLength(), 2*120 - 1);

    parameters.setMKLFloatingPointAccuracy(MKLFloatingPointAccuracy::LOW_ACCURACY);
    EXPECT_EQ(parameters.getMKLFloatingPointAccuracy(),
              MKLFloatingPointAccuracy::LOW_ACCURACY);

    // Check default is set
    EXPECT_EQ(parameters.getFilteringType(),
              CorrelogramFilteringType::NO_FILTERING);
    // fir envelope
    EXPECT_NO_THROW(parameters.setFIREnvelopeFiltering(199));
    EXPECT_EQ(parameters.getFilteringType(),
              CorrelogramFilteringType::FIR_ENVELOPE_FILTERING);
    EXPECT_EQ(parameters.getFIREnvelopeFilterLength(), 199);
    // deal with odd number
    EXPECT_NO_THROW(parameters.setFIREnvelopeFiltering(200));
    EXPECT_EQ(parameters.getFIREnvelopeFilterLength(), 201);

    // still not initialized yet
    EXPECT_FALSE(parameters.isValid());

    std::vector<std::pair<int, int>> xcPairs;
    for (int job=0; job<2; ++job)
    {
        bool ldoAutoCorr = false;
        if (job == 1){ldoAutoCorr = true;}
        int iStart = 2;
        if (job == 1){iStart = 2;}
        for (int is=iStart; is<10; ++is)
        {
            EXPECT_NO_THROW(parameters.setCorrelationPairs(is,
                                                           ldoAutoCorr));
            EXPECT_TRUE(parameters.isValid());
            EXPECT_NO_THROW(xcPairs = parameters.getCorrelationPairs());
            int ixc = 0;
            for (int i=0; i<is; ++i) 
            {
                for (int j=i; j<is; ++j)
                {
                    if (!ldoAutoCorr && j == i){continue;}
                    if (!ldoAutoCorr && i == is - 1 && j == is - 1){continue;}
                    std::pair<int, int> xcPair;
                    std::pair<int, int> refPair(i, j);
                    EXPECT_NO_THROW(xcPair = parameters.getCorrelationPair(ixc));
                    EXPECT_EQ(xcPairs[ixc].first,  refPair.first);
                    EXPECT_EQ(xcPairs[ixc].second, refPair.second);
                    EXPECT_EQ(xcPair.first,  refPair.first);
                    EXPECT_EQ(xcPair.second, refPair.second);
                    ixc = ixc + 1;
                } // Loop on j
            } // Loop on i
            EXPECT_EQ(parameters.getNumberOfCorrelationPairs(), ixc);
        } // Loop on iStart
    }
    // Test the copy assignment and do it all over again
    CorrelogramParameters copyParms(parameters);
    EXPECT_EQ(copyParms.getNumberOfSamples(), 100);
    EXPECT_EQ(copyParms.getNumberOfPaddedSamples(), 120);
    EXPECT_EQ(copyParms.getMKLFloatingPointAccuracy(),
              MKLFloatingPointAccuracy::LOW_ACCURACY);
    EXPECT_EQ(copyParms.getFilteringType(),
              CorrelogramFilteringType::FIR_ENVELOPE_FILTERING);
    EXPECT_EQ(copyParms.getFIREnvelopeFilterLength(), 201);
    EXPECT_TRUE(copyParms.isValid());
    EXPECT_NO_THROW(xcPairs = parameters.getCorrelationPairs()); 
    EXPECT_EQ(static_cast<int> (xcPairs.size()),
              copyParms.getNumberOfCorrelationPairs()); 
    EXPECT_TRUE(xcPairs.size() > 0);
    for (int i=0; i<static_cast<int> (xcPairs.size()); ++i)
    {
        auto xcPair = copyParms.getCorrelationPair(i);
        EXPECT_EQ(xcPair.first,  xcPairs[i].first);
        EXPECT_EQ(xcPair.second, xcPairs[i].second);
    }

}

TEST(testCorrelograms, correlograms)
{
    // Set the parameters
    int nSignals = 5;
    int nSamples = 6;
    int nPaddedSamples = 10;
    bool ldoAutoCorrelograms = true;
    std::array<double, 50> signals({
        1, 2, 3, 4, 5, -1,  0, 0,  0,  0,  // signal 1
        3, 2, 1, 3, 2,  1,  0, 0,  0,  0,  // signal 2
       -1, 2,-1,-2, 1, -1,  0, 0,  0,  0,  // signal 3
        4,-2, 3,-1, 5, -1,  0, 0,  0,  0,  // signal 4
        0,-3,-4, 5,-2,  3,  0, 0,  0,  0   // signal 5
                                    });
    std::array<double, 19*15> xcsRef({
    0,  0,  0,  0, -1,  3, 11, 22, 35, 56, 35, 22, 11,  3, -1,  0,  0, 0,  0, // correlate(x1, x1)
    0,  0,  0,  0,  1,  4, 10, 17, 26, 31, 29, 19, 21, 13, -3,  0,  0, 0,  0, // correlate(x1, x2)
    0,  0,  0,  0, -1, -1, -3, -6, -7, -2, -11, 2,  7, -7,  1,  0,  0, 0,  0, // correlate(x1, x3)
    0,  0,  0,  0, -1,  3,  6, 12, 16, 31,  4, 20,  3, 22, -4,  0,  0, 0,  0, // correlate(x1, x4)
    0,  0,  0,  0,  3,  4, 10, 12, 11,-11,  2,-37,-11,  3,  0,  0,  0, 0,  0, // correlate(x1, x5)
    0,  0,  0,  0,  3,  8, 14, 14, 19, 28, 19, 14, 14,  8,  3,  0,  0, 0,  0, // correlate(x2, x2)
    0,  0,  0,  0, -3,  1, -5, -9,  3, -5, -6,  1,  0,  0, -1,  0,  0, 0,  0, // correlate(x2, x3)
    0,  0,  0,  0, -3, 13,  6,  9, 12, 17, 18,  3, 11,  6,  4,  0,  0, 0,  0, // correlate(x2, x4)
    0,  0,  0,  0,  9,  0, 14,  5,-12,  4, -7,-12,-10, -3,  0,  0,  0, 0,  0, // correlate(x2, x5)
    0,  0,  0,  0,  1, -3,  5, -2, -5, 12, -5, -2,  5, -3,  1,  0,  0, 0,  0, // correlate(x3, x3)
    0,  0,  0,  0,  1, -7, 12, -8, -2, -3, -2,  4,-13,  6, -4,  0,  0, 0,  0, // correlate(x3, x4)
    0,  0,  0,  0, -3,  8,-12, 10, -3,-17, 18, -3,  1,  3,  0,  0,  0, 0,  0, // correlate(x3, x5)
    0,  0,  0,  0, -4, 22,-17, 30,-27, 56,-27, 30,-17, 22, -4,  0,  0, 0,  0, // correlate(x4, x4)
    0,  0,  0,  0, 12,-14, 33,-35, 28,-24, 22,-22,-11,  3,  0,  0,  0, 0,  0, // correlate(x4, x5)
    0,  0,  0,  0,  0, -9, -6,  8,-24, 63,-24,  8, -6, -9,  0,  0,  0, 0,  0  // correlate(x5, x5)
    });
    std::array<double, 19*15> phaseXCsRef({
    0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,
   -6.151232e-02,-1.271464e-01,9.234959e-02,5.188365e-02,1.007947e-01,-1.384503e-01,-4.599154e-02,-6.692472e-02,3.462207e-01,3.382016e-01,4.211022e-01,-1.708914e-01,3.141261e-01,2.171429e-01,-5.416231e-01,1.130543e-01,1.694097e-02,1.927867e-01,-5.206377e-02,
   -6.110641e-02,6.078961e-02,-3.636048e-02,7.143334e-02,-5.332973e-02,5.446901e-02,-1.684043e-01,-1.622949e-01,-4.044898e-01,-1.501134e-01,-5.443236e-01,2.201561e-01,4.897521e-01,-3.524817e-01,9.411259e-03,-9.483153e-02,1.368126e-01,-4.441473e-02,2.932639e-02,
   1.513612e-01,-8.707414e-03,4.476860e-02,6.066890e-02,-2.312901e-01,-1.121518e-01,8.754047e-02,1.274589e-02,3.712141e-01,6.228320e-01,-1.835419e-01,1.252373e-01,-5.147736e-02,4.936928e-01,-1.558544e-01,-1.796304e-01,7.367992e-02,-1.181570e-01,-2.930616e-03,
   -4.657576e-02,-3.346773e-02,-9.391717e-03,-3.150363e-03,1.719709e-02,6.837526e-03,9.243102e-02,1.206705e-01,1.504667e-01,-2.162812e-01,1.558791e-01,-8.723456e-01,-2.205388e-01,1.788610e-01,-3.777590e-02,-1.540325e-01,-1.117871e-01,-4.039521e-03,-1.295668e-02,
   0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,
   -4.610714e-02,-1.039259e-01,-1.496421e-02,1.105277e-01,8.872632e-02,2.485514e-02,-2.575440e-01,-7.296426e-01,2.644207e-01,-4.018670e-01,-1.927628e-01,6.646649e-02,2.133414e-01,3.405788e-02,-1.344336e-01,-1.150401e-01,1.913617e-02,1.077269e-01,6.702860e-02,
   1.475951e-01,6.455310e-02,-4.838415e-03,-1.600762e-01,-4.171381e-01,3.854786e-01,2.620457e-02,1.320523e-01,2.039504e-01,3.606466e-01,3.711038e-01,-3.282225e-01,2.356607e-01,1.242011e-01,3.584132e-02,-6.641257e-02,-8.087082e-02,-2.460889e-01,2.163599e-01,
   1.322476e-01,-1.039455e-02,-1.003740e-03,-2.234629e-01,1.492313e-01,-1.002434e-01,4.766515e-01,1.419308e-01,-4.781173e-01,1.588244e-01,-2.570643e-01,-3.540900e-01,-3.792448e-01,6.728978e-02,3.736448e-02,1.053783e-02,-2.125248e-01,-8.238289e-02,-7.554904e-02,
   0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,
   1.597250e-01,-1.332154e-01,-1.194161e-01,-2.116534e-02,1.314730e-02,-2.572205e-02,5.814538e-01,-2.357973e-01,-3.709761e-01,-3.022459e-01,-9.454845e-02,1.743675e-02,-5.013380e-01,8.370718e-02,-1.379571e-01,-4.836907e-02,9.174166e-02,1.222882e-01,-7.874913e-02,
  -1.309352e-02,9.409002e-02,-5.381068e-02,7.185475e-02,-9.770780e-02,1.677586e-01,-1.729898e-01,1.732512e-01,-2.963524e-01,-5.123604e-01,5.618361e-01,2.734105e-01,2.607573e-01,9.586955e-02,2.115252e-01,4.106354e-02,1.036794e-01,-9.210729e-03,1.004292e-01,
   0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,
   2.220108e-03,-1.109708e-01,-6.928753e-02,1.892797e-01,5.044686e-02,-2.024351e-02,3.270075e-01,-3.664379e-01,1.114898e-01,-9.370640e-02,1.633439e-01,-4.985934e-01,-5.989414e-01,3.595927e-02,3.197732e-02,-5.634580e-02,1.188118e-01,-5.935997e-02,-1.566497e-01,
0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+0});
    // Populate the parameter class
    CorrelogramParameters parameters;
    parameters.setNumberOfSamples(nSamples);
    parameters.setNumberOfPaddedSamples(nPaddedSamples);
    parameters.setCorrelationPairs(nSignals, ldoAutoCorrelograms);
    EXPECT_TRUE(parameters.isValid());
    // Create a correlation engine
    Correlograms<double> dcorr(parameters);
    EXPECT_TRUE(dcorr.isInitialized());
    EXPECT_EQ(dcorr.getInputSignalLength(), 6);
    EXPECT_EQ(dcorr.getCorrelogramLength(), 19);  // 2*nPaddedSamples - 1
    EXPECT_EQ(dcorr.getNumberOfCorrelograms(), 15);
    // Set the input signals 
    for (int i=0; i<nSignals; ++i)
    {
        double *signalPtr = signals.data() + 10*i;
        dcorr.setInputSignal(i, nSamples, signalPtr);
    }
    // Compute correlograms
    dcorr.computeCrossCorrelograms();
    // Get the correlogram
    double error = 0;
    for (int ixc=0; ixc<dcorr.getNumberOfCorrelograms(); ++ixc)
    {
        auto xcPtr = dcorr.getRawCorrelogramPointer(ixc);
        ippsNormDiff_Inf_64f(xcPtr, xcsRef.data()+ixc*19, 19, &error);
        EXPECT_LE(error, 1.e-12);
        xcPtr = dcorr.getProcessedCorrelogramPointer(ixc);
        ippsNormDiff_Inf_64f(xcPtr, xcsRef.data()+ixc*19, 19, &error);
        EXPECT_LE(error, 1.e-12);
    }
    // Compute phase correlograms
    dcorr.computePhaseCorrelograms();
    for (int ixc=0; ixc<dcorr.getNumberOfCorrelograms(); ++ixc)
    {
        auto xcPtr = dcorr.getRawCorrelogramPointer(ixc);
        ippsNormDiff_Inf_64f(xcPtr, phaseXCsRef.data()+ixc*19, 19, &error);
        EXPECT_LE(error, 1.e-7);
        xcPtr = dcorr.getProcessedCorrelogramPointer(ixc);
        ippsNormDiff_Inf_64f(xcPtr, phaseXCsRef.data()+ixc*19, 19, &error);
        EXPECT_LE(error, 1.e-7);
    }
}

TEST(testCorrelograms, generalTest)
{
    int nSignals = 5;
    int nSamples = 1000;
    int nPaddedSamples = 1024;
    bool ldoAutoCorrelograms = false;
    // Populate the parameter class and initialize correlation engine
    CorrelogramParameters parameters;
    parameters.setNumberOfSamples(nSamples);
    parameters.setNumberOfPaddedSamples(nPaddedSamples);
    parameters.setCorrelationPairs(nSignals, ldoAutoCorrelograms);
    EXPECT_TRUE(parameters.isValid());
    Correlograms<double> dcorr(parameters);
    // Create an RNG
    std::random_device rndDevice;
    std::mt19937 mersenneEngine {rndDevice()};  // Generates random numbers 
    std::uniform_real_distribution<> dist(-10.0, 10.0);
    auto gen = [&dist, &mersenneEngine](){return dist(mersenneEngine);};
    // Set the signals to correlate
    auto nxcs = parameters.getNumberOfCorrelationPairs();
    std::vector<std::vector<double>> signals(nSignals);
    for (int i=0; i<nSignals; ++i)
    {
        std::vector<double> x(nPaddedSamples, 0);
        std::generate(x.begin(), x.begin()+nSamples, gen); 
        signals[i] = x;
        dcorr.setInputSignal(i, nSamples, x.data());
    }
    // Compute cross-correlations
    dcorr.computeCrossCorrelograms();
    // Compare
    int lenxc = dcorr.getCorrelogramLength();
    int ixc = 0; 
    for (int i=0; i<nSignals; ++i)
    {
        for (int j=i+1; j<nSignals; ++j)
        {
            auto xcPair = parameters.getCorrelationPair(ixc);
            auto xc1 = signals[xcPair.first];
            auto xc2 = signals[xcPair.second];
            auto xcRef = computeReferenceCorrelogram(nPaddedSamples,
                                                     xc1.data(), xc2.data(),
                                                     false);
            EXPECT_EQ(lenxc, static_cast<int> (xcRef.size()));
            auto xcPtr = dcorr.getRawCorrelogramPointer(ixc);
            double error;
            ippsNormDiff_Inf_64f(xcPtr, xcRef.data(), lenxc, &error); 
            EXPECT_LE(error, 1.e-10);
            ixc = ixc + 1;
        }
    }
    EXPECT_EQ(ixc, nxcs);
    // Compute phase-correlations
    dcorr.computePhaseCorrelograms();
    // Compare
    ixc = 0;
    for (int i=0; i<nSignals; ++i)
    {
        for (int j=i+1; j<nSignals; ++j)
        {
            auto xcPair = parameters.getCorrelationPair(ixc);
            auto xc1 = signals[xcPair.first];
            auto xc2 = signals[xcPair.second];
            auto xcRef = computeReferenceCorrelogram(nPaddedSamples,
                                                     xc1.data(), xc2.data(),
                                                     true);
            EXPECT_EQ(lenxc, static_cast<int> (xcRef.size()));
            auto xcPtr = dcorr.getRawCorrelogramPointer(ixc);
            double error;
            ippsNormDiff_Inf_64f(xcPtr, xcRef.data(), lenxc, &error);
            EXPECT_LE(error, 1.e-10);
            ixc = ixc + 1;
        }
    }
}

std::vector<double>
computeReferenceCorrelogram(const int npts,
                            const double signal1[],
                            const double signal2[],
                            const bool doPhase)
{
    // Initialize transform
    int lenxc = 2*npts - 1;
    RTSeis::Utilities::Transforms::DFTRealToComplex dft;
    dft.initialize(lenxc);
    // Transform
    int lenft = dft.getTransformLength();
    std::vector<std::complex<double>> dft1(lenft);
    std::vector<std::complex<double>> dft2(lenft);
    std::vector<std::complex<double>> xcft(lenft);
    std::complex<double> *ptr = dft1.data();
    dft.forwardTransform(npts, signal1, lenft, &ptr);
    ptr = dft2.data();
    dft.forwardTransform(npts, signal2, lenft, &ptr);
    // Correlate
    for (int i=0; i<lenft; ++i)
    {
        xcft[i] = dft1[i]*std::conj(dft2[i]);
    }
    if (doPhase)
    {
        for (int i=0; i<lenft; ++i)
        {
            xcft[i] = xcft[i]/std::max(1.e-14, std::abs(xcft[i]));
        }
    }
    // Inverse transform
    std::vector<double> xcWork(lenxc, 0);
    double *rptr = xcWork.data();
    dft.inverseTransform(lenft, xcft.data(), lenxc, &rptr);
    // Shuffle the end of the transform to the acausal part
    auto xc = RTSeis::Utilities::Transforms::DFTUtilities::fftShift(xcWork); 
    return xc;
}

}
