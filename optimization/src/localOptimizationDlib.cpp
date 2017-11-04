#include <localOptimizationDlib.h>
#include "logging.h"
#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include <dlib/optimization.h>
#include <dlib/matrix.h>

namespace parpe {

typedef dlib::matrix<double,0,1> column_vector;
typedef dlib::matrix<double> general_matrix;

column_vector dlibColumnVectorFromDoubleArray(double const *src, int len) {
    column_vector colVec(len);
    std::copy(src, src + len, colVec.begin());
    return colVec;
}

int parpe::OptimizerDlibLineSearch::optimize(OptimizationProblem *problem)
{
    int status = 0;

    int numParams = problem->getNumOptimizationParameters();
    column_vector startingPoint = dlibColumnVectorFromDoubleArray(problem->getInitialParameters(), numParams);
    column_vector min = dlibColumnVectorFromDoubleArray(problem->getParametersMin(), numParams);
    column_vector max = dlibColumnVectorFromDoubleArray(problem->getParametersMax(), numParams);

    clock_t timeBegin = clock();
    clock_t timeIterBegin = clock();

    double finalFVal = dlib::find_min_box_constrained(
                dlib::lbfgs_search_strategy(10),
                dlib::objective_delta_stop_strategy(
                    1e-9, problem->getOptimizationOptions().maxOptimizerIterations),
                [&problem](const column_vector& x){
        // objective function
        static __thread int numFunctionCalls = 0;
        clock_t timeBegin = clock();

        double fval = NAN;
        problem->evaluateObjectiveFunction(x.begin(), &fval, nullptr);

        if(problem->getOptimizationOptions().printToStdout)
            std::cout<<std::endl<<numFunctionCalls++<<": "<<fval;

        clock_t timeEnd = clock();
        double wallTime = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;

        problem->logObjectiveFunctionEvaluation(x.begin(), fval, nullptr,
                                                numFunctionCalls, wallTime);

        return fval;

    },
    [&problem, &timeIterBegin](const column_vector& x){
        // objective function gradient
        static __thread int numFunctionCalls = 0;
        clock_t timeBegin = clock();

        double fVal = NAN;
        column_vector fGrad(problem->getNumOptimizationParameters());
        problem->evaluateObjectiveFunction(x.begin(), &fVal, fGrad.begin());

        if(problem->getOptimizationOptions().printToStdout)
            std::cout<<" g"<<numFunctionCalls++;

        clock_t timeEnd = clock();
        double wallTime = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;
//        double wallTimeIter = (double)(timeEnd - timeIterBegin) / CLOCKS_PER_SEC;
        timeIterBegin = clock();

        problem->logObjectiveFunctionEvaluation(x.begin(), fVal, fGrad.begin(), numFunctionCalls, wallTime);
        problem->intermediateFunction(
                    0, numFunctionCalls, fVal, 0, 0, 0, 0, 0, 0, 0, 0);
        return fGrad;
    },
    startingPoint, min, max);

    clock_t timeEnd = clock();
    double wallTime = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;

    problem->logOptimizerFinished(finalFVal, startingPoint.begin(), wallTime, 0);

    return status;
}

} // namespace parpe
