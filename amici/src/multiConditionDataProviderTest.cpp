#include "MultiConditionDataProvider.h"
#include "testingMisc.h"

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include <amici/model.h>
// #include "../tests/cpputest/testfunctions.h" // for Modell_Test

/**
 * @brief Mock MultiConditionDataProvider
 */
/*
class MultiConditionDataProviderTest : public parpe::MultiConditionDataProvider {
  public:
    MultiConditionDataProviderTest() {

        model = std::unique_ptr<amici::Model>(new amici::Model_Test(0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                                                    amici::AMICI_O2MODE_NONE,std::vector<realtype>(10),
                                                                    std::vector<realtype>(),std::vector<int>(),
                                                                    std::vector<realtype>(),std::vector<int>()));
    }

    int getNumberOfConditions() const override { return numConditions; }

//    int getNumConditionSpecificParametersPerSimulation() const override {
//        return numCondSpecParamPerSim;
//    }

    int numConditions = 1;
};
*/

// clang-format off
TEST_GROUP(multiConditionDataProvider){
    void setup() {
        parpe::initHDF5Mutex();
    }

    void teardown() {
    }
};
// clang-format on


/**
 * @brief Test mapping simulation<->optimization parameters with no
 * condition-specific parameters
 */
//TEST(multiConditionDataProvider, testDataProviderParameterMapping) {
//    MultiConditionDataProviderTest dp;
//    dp.numConditions = 1;
//    dp.numCondSpecParamPerSim = 0;

//    CHECK_EQUAL(10, dp.getNumOptimizationParameters());
//}

/**
 * @brief Test mapping simulation<->optimization parameters with
 * condition-specific parameters
 */
//TEST(multiConditionDataProvider,
//     testDataProviderParameterMappingConditionSpec) {
//    MultiConditionDataProviderTest dp;
//    dp.numConditions = 3;
//    dp.numCondSpecParamPerSim = 2;

//    CHECK_EQUAL(14, dp.getNumOptimizationParameters());
//}
