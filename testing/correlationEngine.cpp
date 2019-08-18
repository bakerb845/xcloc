#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include "xcloc/correlationEngineParameters.hpp"
#include "xcloc/correlationEngine.hpp"
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
    dcorr.computeCrossCorrelograms();
}

}
