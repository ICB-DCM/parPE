#include <parpeoptimization/localOptimizationDlib.h>

#include <parpecommon/logging.h>
#include <parpeoptimization/optimizationOptions.h>
#include <parpeoptimization/optimizationProblem.h>

#include <dlib/optimization.h>
#include <dlib/matrix.h>

typedef dlib::matrix<double,0,1> column_vector;
typedef dlib::matrix<double> general_matrix;

namespace gsl {

template <
    typename T,
    long num_rows,
    long num_cols,
    typename mem_manager,
    typename layout
    >
span<const T> make_span(dlib::matrix<T,num_rows,num_cols, mem_manager,layout> const& matrix) {
    return span<const T>(matrix.begin(), matrix.size());
}

template <
    typename T,
    long num_rows,
    long num_cols,
    typename mem_manager,
    typename layout
    >
span<T> make_span(dlib::matrix<T,num_rows,num_cols, mem_manager,layout> &matrix) {
    return span<T>(matrix.begin(), matrix.size());
}

}

namespace parpe {


column_vector dlibColumnVectorFromDoubleArray(double const *src, int len) {
    column_vector colVec(len);
    std::copy(src, src + len, colVec.begin());
    return colVec;
}

std::tuple<int, double, std::vector<double> >
parpe::OptimizerDlibLineSearch::optimize(OptimizationProblem *problem)
{
    int numParams = problem->cost_fun_->numParameters();

    column_vector startingPoint(numParams);
    problem->fillInitialParameters(gsl::make_span(startingPoint));

    column_vector min(numParams);
    problem->fillParametersMin(gsl::make_span(startingPoint));

    column_vector max(numParams);
    problem->fillParametersMax(gsl::make_span(startingPoint));


    auto optimizationController = problem->getReporter();

    optimizationController->starting(gsl::make_span(startingPoint));

    double finalFVal = dlib::find_min_box_constrained(
        dlib::lbfgs_search_strategy(10),
        dlib::objective_delta_stop_strategy(
            1e-9, problem->getOptimizationOptions().maxOptimizerIterations),
        [&optimizationController](const column_vector& x){
            // objective function
            double fval = NAN;
            optimizationController->evaluate(
                gsl::make_span(x),
                fval, gsl::span<double>(nullptr, 0));

            return fval;

        },
        [&optimizationController](const column_vector& x){
            // objective function gradient

            double fVal = NAN;
            column_vector fGrad(optimizationController->numParameters());
            optimizationController->evaluate(
                gsl::make_span(x), fVal, gsl::make_span(fGrad));
            return fGrad;
        },
        startingPoint, min, max);

    optimizationController->finished(
        finalFVal, gsl::make_span(startingPoint), 0);

    return std::make_tuple(
        0, finalFVal,
        std::vector<double>(startingPoint.begin(), startingPoint.end()));
}

} // namespace parpe
