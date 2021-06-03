#include <parpecommon/parpeConfig.h>

#include <multiStartOptimizationTest.h>
#include <minibatchOptimizationTest.h>
#include <optimizationResultWriterTest.h>
#include <optimizationOptionsTest.h>
#include <optimizationProblemTest.h>

#ifdef PARPE_ENABLE_IPOPT
#include "localOptimizationIpoptTest.h"
#endif

#ifdef PARPE_ENABLE_FIDES
#include "localOptimizationFidesTest.h"
#endif

#ifdef PARPE_ENABLE_CERES
#include "localOptimizationCeresTest.h"
#endif

#ifdef PARPE_ENABLE_TOMS611
#include "localOptimizationToms611Test.cpp"
#endif

#ifdef PARPE_ENABLE_FSQP
#include "localOptimizationFsqpTest.cpp"
#endif

#include <gtest/gtest.h>

#include <cstdlib>
#include <ctime>


int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
