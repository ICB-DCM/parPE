#include <localOptimizationDlib.h>
#include "logging.h"
#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include <dlib/optimization.h>
#include <dlib/matrix.h>

namespace parpe {

typedef dlib::matrix<double,0,1> column_vector;
typedef dlib::matrix<double> general_matrix;

int parpe::OptimizerDlibLineSearch::optimize(OptimizationProblem *problem)
{
    int status = 0;

    column_vector startingPoint(problem->getNumOptimizationParameters());
    auto tmpStartingPoint = std::unique_ptr<double const[]>(problem->getInitialParameters());
    std::copy(tmpStartingPoint.get(), tmpStartingPoint.get() + problem->getNumOptimizationParameters(),
              startingPoint.begin());

    column_vector min(problem->getNumOptimizationParameters());
    column_vector max(problem->getNumOptimizationParameters());
    std::fill(min.begin(), min.end(), -2.0);
    std::fill(max.begin(), max.end(), 2.0);
//    std::cout<<startingPoint<<"/"<<min<<"/"<<max<<std::endl;

    // TODO toColumnVector  for constraints
    dlib::find_min_box_constrained(dlib::lbfgs_search_strategy(10),
                                   dlib::objective_delta_stop_strategy(1e-9),
                                   [&problem](const column_vector& x){
        static int iter = 0;
        double result = NAN;
        problem->evaluateObjectiveFunction(x.begin(), &result, nullptr);
        std::cout<<std::endl<<iter++<<": "<<result;
        return result;

    },
    [&problem](const column_vector& x){
        double unusedResult = NAN;
        column_vector grad(problem->getNumOptimizationParameters());
        problem->evaluateObjectiveFunction(x.begin(), &unusedResult, grad.begin());
        std::cout<<"g";
        return grad;
    },
    startingPoint, min, max); // TODO: constraints from problem

    return status;
}

} // namespace parpe
