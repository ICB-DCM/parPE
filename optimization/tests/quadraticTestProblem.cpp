#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "optimizationOptions.h"
#include "quadraticTestProblem.h"
#include <cmath>
#include <cstdio>

namespace parpe {

QuadraticTestProblem::QuadraticTestProblem() : OptimizationProblem(1) {
    optimizationOptions.maxOptimizerIterations = 12;
    optimizationOptions.optimizer = OPTIMIZER_IPOPT;

    parametersMin_[0] = -1e5;
    parametersMax_[0] = 1e5;
}

int QuadraticTestProblem::evaluateObjectiveFunction(const double *parameters,
                                                    double *objFunVal,
                                                    double *objFunGrad) {
    if (objFunGrad) {
        mock().actualCall("testObjGrad");

        objFunVal[0] = pow(parameters[0] + 1.0, 2) + 42.0;
        objFunGrad[0] = 2.0 * parameters[0] + 2.0;

        if(printDebug)
            printf("g: x: %f f(x): %f f'(x): %f\n", parameters[0], *objFunVal,
                    objFunGrad[0]);
    } else {
        mock().actualCall("testObj");

        *objFunVal = pow(parameters[0] + 1.0, 2) + 42.0;

        if(printDebug)
            printf("f: x: %f f(x): %f\n", parameters[0], *objFunVal);
    }
    return 0;
}

void QuadraticTestProblem::logOptimizerFinished(double optimalCost,
                                                const double *optimalParameters,
                                                double masterTime,
                                                int exitStatus) {
    mock().actualCall("logFinish").withIntParameter("exitStatus", exitStatus);

    this->optimalCost = optimalCost;
    this->optimalParameter = optimalParameters[0];

    if(printDebug)
        printf("finished: f(x*): %f x*: %f t: %fs exit: %d\n", optimalCost, optimalParameters[0],
                masterTime, exitStatus);
}


std::unique_ptr<OptimizationProblem> QuadraticOptimizationMultiStartProblem::getLocalProblemImpl(
        int multiStartIndex) {
    return std::unique_ptr<OptimizationProblem>(new QuadraticTestProblem());
}

} // namespace parpe
