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

    dlib::find_min_box_constrained(
                dlib::lbfgs_search_strategy(10),
                dlib::objective_delta_stop_strategy(
                    1e-9, problem->getOptimizationOptions().maxOptimizerIterations),
                [&problem](const column_vector& x){
        // objective function
        static __thread int numFunctionCalls = 0;

        double fval = NAN;
        problem->evaluateObjectiveFunction(x.begin(), &fval, nullptr);

        if(problem->getOptimizationOptions().printToStdout)
            std::cout<<std::endl<<numFunctionCalls++<<": "<<fval;

        return fval;

    },
    [&problem](const column_vector& x){
        // objective function gradient
        static __thread int numFunctionCalls = 0;

        double unusedFVal = NAN;
        column_vector fGrad(problem->getNumOptimizationParameters());
        problem->evaluateObjectiveFunction(x.begin(), &unusedFVal, fGrad.begin());

        if(problem->getOptimizationOptions().printToStdout)
            std::cout<<" g"<<numFunctionCalls++;

        return fGrad;
    },
    startingPoint, min, max);

    return status;
}

} // namespace parpe
