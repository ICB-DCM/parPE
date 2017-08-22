#include "MultiConditionDataProvider.h"
#include "testingMisc.h"

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include <amici_model.h>

/**
 * @brief Mock MultiConditionDataProvider
 */
class MultiConditionDataProviderTest : public MultiConditionDataProvider {
  public:
    MultiConditionDataProviderTest() {
        model = new Model(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                          23, 24, 25, AMICI_O2MODE_NONE);
    }

    int getNumberOfConditions() const { return numConditions; }

    int getNumConditionSpecificParametersPerSimulation() const {
        return numCondSpecParamPerSim;
    }

    int numConditions = 1;

    int numCondSpecParamPerSim = 0;

    ~MultiConditionDataProviderTest() { delete model; }
};

TEST_GROUP(multiConditionDataProvider){void setup(){initHDF5Mutex();
}

void teardown() { destroyHDF5Mutex(); }
}
;

/**
 * @brief Test mapping simulation<->optimization parameters with no
 * condition-specific parameters
 */
TEST(multiConditionDataProvider, testDataProviderParameterMapping) {
    MultiConditionDataProviderTest dp;
    dp.numConditions = 1;
    dp.numCondSpecParamPerSim = 0;

    CHECK_EQUAL(10, dp.getNumCommonParameters());

    CHECK_EQUAL(10, dp.getNumOptimizationParameters());

    CHECK_EQUAL(10,
                dp.getIndexOfFirstConditionSpecificOptimizationParameter(0));

    CHECK_EQUAL(10,
                dp.getIndexOfFirstConditionSpecificOptimizationParameter(1));
}

/**
 * @brief Test mapping simulation<->optimization parameters with
 * condition-specific parameters
 */
TEST(multiConditionDataProvider,
     testDataProviderParameterMappingConditionSpec) {
    MultiConditionDataProviderTest dp;
    dp.numConditions = 3;
    dp.numCondSpecParamPerSim = 2;

    CHECK_EQUAL(8, dp.getNumCommonParameters());

    CHECK_EQUAL(14, dp.getNumOptimizationParameters());

    CHECK_EQUAL(8, dp.getIndexOfFirstConditionSpecificOptimizationParameter(0));

    CHECK_EQUAL(10,
                dp.getIndexOfFirstConditionSpecificOptimizationParameter(1));
}
