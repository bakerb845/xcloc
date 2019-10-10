#include <cstdio>
#include <cstdlib>
#include <string>
#include "xcloc/waveformIdentifier.hpp"
#include <gtest/gtest.h>
namespace
{
using namespace XCLoc;
TEST(WaveformIdentifier, waveid)
{
    WaveformIdentifier waveid("UU", "DUG", "P");
    WaveformIdentifier waveidCopy(waveid);
    EXPECT_TRUE(waveid == waveidCopy);
    EXPECT_EQ(waveid.getNetwork(), "UU");
    EXPECT_EQ(waveid.getStation(), "DUG");
    EXPECT_EQ(waveid.getPolarization(), "P");
    waveidCopy.setPolarization("S");
    EXPECT_TRUE(waveid != waveidCopy);
    EXPECT_TRUE( waveid > waveidCopy || waveidCopy > waveid);
}
}
