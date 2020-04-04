#include "amiciSimulationRunnerTest.h"
#include "multiConditionDataProviderTest.h"
#include "multiConditionProblemTest.h"
#include "simulationResultWriterTest.h"
#include "hierarchicalOptimizationTest.h"

#include <gtest/gtest.h>

#include <cstdlib>
#include <ctime>

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
