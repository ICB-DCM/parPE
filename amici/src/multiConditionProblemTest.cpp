#include "multiConditionProblem.h"
#include "testingMisc.h"

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

// mock amici
int runAmiciSimulation(UserData *, ExpData const *, ReturnData *) { return 0; }

/**
 * @brief Mock MultiConditionProblem
 */
class MultiConditionProblemTest : public parpe::MultiConditionProblem {
  public:
    MultiConditionProblemTest() {}
    void addSimulationGradientToObjectiveFunctionGradient(
        const double *simulationGradient, double *objectiveFunctionGradient,
        int numCommon, int numConditionSpecificParams,
        int firstIndexOfCurrentConditionsSpecificOptimizationParameters) {
        MultiConditionProblem::
            addSimulationGradientToObjectiveFunctionGradientConditionSpecificParameters(
                simulationGradient, objectiveFunctionGradient, numCommon,
                numConditionSpecificParams,
                firstIndexOfCurrentConditionsSpecificOptimizationParameters);
    }
};

TEST_GROUP(multiConditionProblem){void setup(){parpe::initHDF5Mutex();
}

void teardown() { parpe::destroyHDF5Mutex(); }
}
;

/**
* @brief Check gradient aggregation
 */
TEST(multiConditionProblem, testAggregateGradientAllCommon) {
    MultiConditionProblemTest p;

    const double simulationGradient[] = {1.0, 2.0, 3.0, 4.0};
    const double objectiveFunctionGradientExpected[] = {-1, -1, -1, -1};
    double objectiveFunctionGradient[] = {-1, -1, -1, -1};

    int numCommon = 4;
    int numConditionSpecificParams = 0;
    int firstIndexOfCurrentConditionsSpecificOptimizationParameters = 0;

    // should not alter objectiveFunctionGradient
    p.addSimulationGradientToObjectiveFunctionGradient(
        simulationGradient, objectiveFunctionGradient, numCommon,
        numConditionSpecificParams,
        firstIndexOfCurrentConditionsSpecificOptimizationParameters);

    for (int i = 0; i < numCommon; ++i)
        CHECK_EQUAL(objectiveFunctionGradientExpected[i],
                    objectiveFunctionGradient[i]);
}

/**
 * @brief Check gradient aggregation
 */
TEST(multiConditionProblem, testAggregateGradientSpecific) {
    MultiConditionProblemTest p;

    const int numCommon = 2;
    const int numConditionSpecificParams = 2;
    const int numConditions = 2;
    const int firstIndexOfCurrentConditionsSpecificOptimizationParameters = 4;

    double objectiveFunctionGradient[] = {-1.0, -1.0, -1.0, -1.0, 2.0, 3.0};
    const double simulationGradient[] = {1.0, 2.0, 3.0, 4.0};
    const double objectiveFunctionGradientExpected[] = {-1.0, -1.0, -1.0,
                                                        -1.0, -1.0, -1.0};
    // should not alter objectiveFunctionGradient
    p.addSimulationGradientToObjectiveFunctionGradient(
        simulationGradient, objectiveFunctionGradient, numCommon,
        numConditionSpecificParams,
        firstIndexOfCurrentConditionsSpecificOptimizationParameters);

    for (int i = 0; i < numCommon + numConditions * numConditionSpecificParams;
         ++i)
        CHECK_EQUAL(objectiveFunctionGradientExpected[i],
                    objectiveFunctionGradient[i]);
}
