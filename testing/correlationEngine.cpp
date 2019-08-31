#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include <ipps.h>
#include "xcloc/correlationEngineParameters.hpp"
#include "xcloc/correlationEngine.hpp"
#include "xcloc/correlogramPostProcessorParameters.hpp"
#include "xcloc/correlogramPostProcessor.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace XCLoc;

TEST(testCorrelationEngine, parameters)
{
    CorrelationEngineParameters parameters;

    EXPECT_NO_THROW(parameters.setNumberOfSamples(100));
    EXPECT_EQ(parameters.getNumberOfSamples(), 100);
    EXPECT_EQ(parameters.getNumberOfPaddedSamples(), 100);

    EXPECT_NO_THROW(parameters.setNumberOfPaddedSamples(120));
    EXPECT_EQ(parameters.getNumberOfPaddedSamples(), 120);

    parameters.setMKLFloatingPointAccuracy(MKLFloatingPointAccuracy::LOW_ACCURACY);
    EXPECT_EQ(parameters.getMKLFloatingPointAccuracy(),
              MKLFloatingPointAccuracy::LOW_ACCURACY);

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
    CorrelationEngineParameters copyParms(parameters);
    EXPECT_EQ(copyParms.getNumberOfSamples(), 100);
    EXPECT_EQ(copyParms.getNumberOfPaddedSamples(), 120);
    EXPECT_EQ(copyParms.getMKLFloatingPointAccuracy(),
              MKLFloatingPointAccuracy::LOW_ACCURACY);
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

TEST(testCorrelationEngine, correlograms)
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
        0,-3,-4, 5,-2,  3,  0, 0,  0,  0  // signal 5
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
    // Populate the parameter class
    CorrelationEngineParameters parameters;
    parameters.setNumberOfSamples(nSamples);
    parameters.setNumberOfPaddedSamples(nPaddedSamples);
    parameters.setCorrelationPairs(nSignals, ldoAutoCorrelograms);
    EXPECT_TRUE(parameters.isValid());
    // Create a correlation engine
    CorrelationEngine<double> dcorr(parameters);
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
    for (int ixc=0; ixc<dcorr.getNumberOfCorrelograms(); ++ixc)
    {
        auto xcPtr = dcorr.getCorrelogramPointer(ixc);
        double error; 
        ippsNormDiff_Inf_64f(xcPtr, xcsRef.data()+ixc*19, 19, &error);
        EXPECT_LE(error, 1.e-12);
    }
}

TEST(testCorrelogramPostProcessor, parameters)
{
    CorrelogramPostProcessorParameters parameters;
    EXPECT_FALSE(parameters.isValid());
    
    parameters.setNoFiltering();
    EXPECT_TRUE(parameters.isValid());
    EXPECT_EQ(parameters.getFilteringType(),
              CorrelogramFilteringType::NO_FILTERING);
    parameters.clear();

    // fir envelope
    EXPECT_FALSE(parameters.isValid());
    EXPECT_NO_THROW(parameters.setFIREnvelopeFiltering(199));
    EXPECT_EQ(parameters.getFilteringType(),
              CorrelogramFilteringType::FIR_ENVELOPE_FILTERING);
    EXPECT_EQ(parameters.getFIREnvelopeFilterLength(), 199);
    // deal with odd number
    EXPECT_NO_THROW(parameters.setFIREnvelopeFiltering(200));
    EXPECT_EQ(parameters.getFIREnvelopeFilterLength(), 201);
    EXPECT_TRUE(parameters.isValid());

    // Test copy constructor
    CorrelogramPostProcessorParameters copyParms(parameters);
    EXPECT_EQ(copyParms.getFilteringType(),
              CorrelogramFilteringType::FIR_ENVELOPE_FILTERING);
    EXPECT_EQ(copyParms.getFIREnvelopeFilterLength(), 201);
}

TEST(testCorrelogramPostProcessor, processor)
{
    CorrelogramPostProcessorParameters parameters;
    parameters.setFIREnvelopeFiltering(301);

    CorrelogramPostProcessor<double> dproc;

} 

}
