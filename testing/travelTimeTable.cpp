#include <cstdio>
#include <cstdlib>
#include <string>
#include "xcloc/travelTimeTableName.hpp"
#include <gtest/gtest.h>
namespace
{
using namespace XCLoc;

TEST(TravelTimeTables, name)
{
    /// Test c'tor -> which calls subsidiary functions
    TravelTimeTableName nameRef("FK", "STA1", "P");
    TravelTimeTableName name(nameRef);
    EXPECT_TRUE(name.isValid());
    EXPECT_EQ(name.getNetwork(), "FK");
    EXPECT_EQ(name.getStation(), "STA1");
    EXPECT_EQ(name.getPhase(), "P");
}


}
