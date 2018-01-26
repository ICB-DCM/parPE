#include <bits/stl_tree.h>

#include "MultiConditionDataProvider.h"
#include "multiConditionProblem.h"
#include "multiConditionProblemResultWriter.h"
#include "testingMisc.h"
#include "hdf5Misc.h"

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

// clang-format off
TEST_GROUP(multiConditionProblemResultWriter){
    void setup() {
        parpe::initHDF5Mutex();
    }

    void teardown() {
    }
};
// clang-format on

TEST(multiConditionProblemResultWriter, testResultWriter) {
    parpe::JobIdentifier id;
    parpe::MultiConditionProblemResultWriter w("deleteme.h5", true, id);

    w.setRootPath("/bla/");

    w.logLocalOptimizerIteration(1, NULL, 0, 2, NULL, 1);

    w.logLocalOptimizerObjectiveFunctionEvaluation(NULL, 0, 1, NULL, 1, 2, 3);

    parpe::logSimulation(w.getFileId(), "/test", NULL, 1, NULL, 1, 1, 2, NULL, NULL, 0, NULL, 1, 2, 0);

    w.saveLocalOptimizerResults(1, NULL, 0, 12, 0);
}

// IGNORE_TEST(multiConditionProblemResultWriter, testResultWriter) {

//    // TODO steadys tate example
//    UserData udata(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,AMI_SCALING_LN,
//    AMI_O2MODE_NONE);
//    MultiConditionDataProvider dataprovider("********.h5") ;
//    OptimizationProblem *problem = new MultiConditionProblem(&dataprovider);

//    double parameters[problem->numOptimizationParameters];
//    double objectiveFunctionValue = 3.0;
//    int numFunctionCalls = 12;
//    double timeElapsed = 123.4;

//    logLocalOptimizerObjectiveFunctionEvaluation(problem, parameters,
//    objectiveFunctionValue, NULL, numFunctionCalls, timeElapsed);

//    // TODO check

//    double objectiveFunctionGradient[problem->numOptimizationParameters];

//    logLocalOptimizerObjectiveFunctionEvaluation(problem, parameters,
//    objectiveFunctionValue, objectiveFunctionGradient,
//                                                                  numFunctionCalls,
//                                                                  timeElapsed);

//    // TODO check

//    flushResultWriter();

//    int numIterations = 13;
//    JobIdentifier path = {1};
//    logLocalOptimizerIteration(path, numIterations, parameters,
//                                    objectiveFunctionValue,
//                                    objectiveFunctionGradient,
//                                    timeElapsed,
//                                    problem->numOptimizationParameters,
//                                    1,//int alg_mod,
//                                    1,2,//double inf_pr, double inf_du,
//                                    1,//double mu,
//                                    1,//double d_norm,
//                                    1,//double regularization_size,
//                                    1,1,//double alpha_du, double alpha_pr,
//                                    1);//int ls_trials);

//    logSimulation(path, parameters, objectiveFunctionValue,
//    objectiveFunctionGradient,
//                       timeElapsed, problem->numOptimizationParameters, 0,
//                       NULL, NULL, NULL, 1, 10);

//    double totalTime = 1234;
//    saveTotalWalltime(totalTime);

//    saveLocalOptimizerResults(problem, path, objectiveFunctionValue,
//                                   parameters,
//                                   totalTime,
//                                   10);

//    delete problem;
//}
