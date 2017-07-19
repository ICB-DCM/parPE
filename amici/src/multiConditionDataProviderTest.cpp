#include "testingMisc.h"
#include "MultiConditionDataProvider.h"

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

/** @brief Dummy function so we don't need to link against an AMICI model.
 */
UserData getUserData() {
    return UserData();
}

/**
 * @brief Mock MultiConditionDataProvider
 */
class MultiConditionDataProviderTest : public MultiConditionDataProvider {
public:
    MultiConditionDataProviderTest()
    {
        *const_cast<int*>(&modelDims.np) = 10;
    }

    int getNumberOfConditions() const {
        return numConditions;
    }

    int getNumConditionSpecificParametersPerSimulation() const {
        return numCondSpecParamPerSim;
    }

    int numConditions = 1;

    int numCondSpecParamPerSim = 0;

    ~MultiConditionDataProviderTest() {}
};


TEST_GROUP(multiConditionDataProvider)
{
    void setup() {
        initHDF5Mutex();
    }

    void teardown() {
        destroyHDF5Mutex();
    }
};

/**
 * @brief Test mapping simulation<->optimization parameters with no condition-specific parameters
 */
TEST(multiConditionDataProvider, testDataProviderParameterMapping) {
    MultiConditionDataProviderTest dp;
    dp.numConditions = 1;
    dp.numCondSpecParamPerSim = 0;

    CHECK_EQUAL(10, dp.getNumCommonParameters());

    CHECK_EQUAL(10, dp.getNumOptimizationParameters());

    CHECK_EQUAL(10, dp.getIndexOfFirstConditionSpecificOptimizationParameter(0));

    CHECK_EQUAL(10, dp.getIndexOfFirstConditionSpecificOptimizationParameter(1));
}


/**
 * @brief Test mapping simulation<->optimization parameters with condition-specific parameters
 */
TEST(multiConditionDataProvider, testDataProviderParameterMappingConditionSpec) {
    MultiConditionDataProviderTest dp;
    dp.numConditions = 3;
    dp.numCondSpecParamPerSim = 2;

    CHECK_EQUAL(8, dp.getNumCommonParameters());

    CHECK_EQUAL(14, dp.getNumOptimizationParameters());

    CHECK_EQUAL(8, dp.getIndexOfFirstConditionSpecificOptimizationParameter(0));

    CHECK_EQUAL(10, dp.getIndexOfFirstConditionSpecificOptimizationParameter(1));
}
