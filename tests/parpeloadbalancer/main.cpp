#include "loadBalancerMasterTest.h"

#include <gtest/gtest.h>
#include <cstdlib>
#include <ctime>

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
