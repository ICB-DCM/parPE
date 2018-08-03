#include "multiConditionProblem.h"
#include "testingMisc.h"

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

// mock amici
namespace amici {
void runAmiciSimulation(Solver &solver, const ExpData *edata,
                       ReturnData *rdata, Model &model);
}

/**
 * @brief Mock MultiConditionProblem
 */
class AmiciSummedGradientFunctionTest : public parpe::AmiciSummedGradientFunction<int> {
  public:
    AmiciSummedGradientFunctionTest() = default;
//    void addSimulationGradientToObjectiveFunctionGradient(
//        const double *simulationGradient, double *objectiveFunctionGradient,
//        int numCommon, int numConditionSpecificParams,
//        int firstIndexOfCurrentConditionsSpecificOptimizationParameters) {

//        AmiciSummedGradientFunction::addSimulationGradientToObjectiveFunctionGradientConditionSpecificParameters(
//                simulationGradient, objectiveFunctionGradient, numCommon,
//                numConditionSpecificParams,
//                firstIndexOfCurrentConditionsSpecificOptimizationParameters);
//    }
};

// clang-format off
TEST_GROUP(multiConditionProblem){
    void setup() {
        parpe::initHDF5Mutex();
    }

    void teardown() {
    }
};
// clang-format on


/**
* @brief Check gradient aggregation
 */
//TEST(multiConditionProblem, testAggregateGradientAllCommon) {
//    AmiciSummedGradientFunctionTest p;

//    const double simulationGradient[] = {1.0, 2.0, 3.0, 4.0};
//    const double objectiveFunctionGradientExpected[] = {-1, -1, -1, -1};
//    double objectiveFunctionGradient[] = {-1, -1, -1, -1};

//    int numCommon = 4;
//    int numConditionSpecificParams = 0;
//    int firstIndexOfCurrentConditionsSpecificOptimizationParameters = 0;

//    // should not alter objectiveFunctionGradient
//    p.addSimulationGradientToObjectiveFunctionGradient(
//        simulationGradient, objectiveFunctionGradient, numCommon,
//        numConditionSpecificParams,
//        firstIndexOfCurrentConditionsSpecificOptimizationParameters);

//    for (int i = 0; i < numCommon; ++i)
//        CHECK_EQUAL(objectiveFunctionGradientExpected[i],
//                    objectiveFunctionGradient[i]);
//}

/**
 * @brief Check gradient aggregation
 */
//TEST(multiConditionProblem, testAggregateGradientSpecific) {
//    AmiciSummedGradientFunctionTest p;

//    const int numCommon = 2;
//    const int numConditionSpecificParams = 2;
//    const int numConditions = 2;
//    const int firstIndexOfCurrentConditionsSpecificOptimizationParameters = 4;

//    double objectiveFunctionGradient[] = {-1.0, -1.0, -1.0, -1.0, 2.0, 3.0};
//    const double simulationGradient[] = {1.0, 2.0, 3.0, 4.0};
//    const double objectiveFunctionGradientExpected[] = {-1.0, -1.0, -1.0,
//                                                        -1.0, -1.0, -1.0};
//    // should not alter objectiveFunctionGradient
//    p.addSimulationGradientToObjectiveFunctionGradient(
//        simulationGradient, objectiveFunctionGradient, numCommon,
//        numConditionSpecificParams,
//        firstIndexOfCurrentConditionsSpecificOptimizationParameters);

//    for (int i = 0; i < numCommon + numConditions * numConditionSpecificParams;
//         ++i)
//        CHECK_EQUAL(objectiveFunctionGradientExpected[i],
//                    objectiveFunctionGradient[i]);
//}

// TODO: check proper starting points are used
