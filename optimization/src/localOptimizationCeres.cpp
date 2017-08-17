#include "localOptimizationCeres.h"
#include "optimizationOptions.h"
#include <ceres/ceres.h>

class CeresWrapper : public ceres::FirstOrderFunction {

  public:
    CeresWrapper(OptimizationProblem *problem) : problem(problem) {}

    /**
     * @brief Evaluate cost function
     * @param parameters
     * @param cost
     * @param gradient If not NULL, evaluate gradient
     * @return true on success, false otherwise
     */
    virtual bool Evaluate(const double *parameters, double *cost,
                          double *gradient) const {
        bool status =
            problem->evaluateObjectiveFunction(parameters, cost, gradient);

        return status == 0;
    }

    virtual int NumParameters() const {
        return problem->numOptimizationParameters;
    }

  private:
    OptimizationProblem *problem;
};

class MyIterationCallback : public ceres::IterationCallback {
  public:
    MyIterationCallback(OptimizationProblem *problem) : problem(problem) {}

    virtual ceres::CallbackReturnType
    operator()(const ceres::IterationSummary &summary) {
        switch (problem->intermediateFunction(
            0, summary.iteration, summary.cost, 0, 0, 0, 0, 0, 0, 0, 0)) {
        case 0:
            return ceres::SOLVER_CONTINUE;
        default:
            return ceres::SOLVER_ABORT;
        }
    }

  private:
    OptimizationProblem *problem = NULL;
};

/**
 * @brief getLocalOptimumCeres
 * @param problem
 * @return 1 on failure, 0 on success
 */
int getLocalOptimumCeres(OptimizationProblem *problem) {
    double *parameters = (double *)malloc(sizeof(*parameters) *
                                          problem->numOptimizationParameters);
    double *startingPoint = problem->getInitialParameters();
    if (startingPoint) {
        // copy, because will be update each iteration
        memcpy(parameters, startingPoint,
               sizeof(*parameters) * problem->numOptimizationParameters);
    } else {
        getRandomStartingpoint(problem->parametersMin, problem->parametersMax,
                               problem->numOptimizationParameters, parameters);
    }

    ceres::GradientProblemSolver::Options options;
    options.minimizer_progress_to_stdout =
        problem->optimizationOptions->printToStdout;
    options.max_num_iterations =
        problem->optimizationOptions->maxOptimizerIterations;
    //    options.gradient_tolerance = 1e-18;
    options.function_tolerance =
        problem->optimizationOptions->functionTolerance;

    MyIterationCallback callback(problem);
    options.callbacks.push_back(&callback);

    // Use quadratic approximation to not require gradient for line search
    // NOTE: WOLFE line_search_type will always require gradient
    // options.line_search_interpolation_type = ceres::QUADRATIC;

    ceres::GradientProblemSolver::Summary summary;

    ceres::GradientProblem ceresProblem(new CeresWrapper(problem));
    ceres::Solve(options, ceresProblem, parameters, &summary);

    //    std::cout<<summary.FullReport();

    problem->logOptimizerFinished(summary.final_cost, parameters,
                                  summary.total_time_in_seconds,
                                  summary.termination_type);

    free(parameters);

    return summary.termination_type == ceres::FAILURE ||
           summary.termination_type == ceres::USER_FAILURE;
}
