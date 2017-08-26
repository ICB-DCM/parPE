#include "localOptimizationCeres.h"
#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include <ceres/ceres.h>

class MyCeresFirstOrderFunction : public ceres::FirstOrderFunction {

  public:
    MyCeresFirstOrderFunction(OptimizationProblem *problem)
        : problem(problem) {}

    /**
     * @brief Evaluate cost function
     * @param parameters
     * @param cost
     * @param gradient If not NULL, evaluate gradient
     * @return true on success, false otherwise
     */
    virtual bool Evaluate(const double *parameters, double *cost,
                          double *gradient) const override {
        bool status =
            problem->evaluateObjectiveFunction(parameters, cost, gradient);

        return status == 0;
    }

    virtual int NumParameters() const override {
        return problem->getNumOptimizationParameters();
    }

  private:
    OptimizationProblem *problem;
};

class MyIterationCallback : public ceres::IterationCallback {
  public:
    MyIterationCallback(OptimizationProblem *problem) : problem(problem) {}

    virtual ceres::CallbackReturnType
    operator()(const ceres::IterationSummary &summary) override {
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

OptimizerCeres::OptimizerCeres() {}

int OptimizerCeres::optimize(OptimizationProblem *problem) {
    double *parameters = new double[problem->getNumOptimizationParameters()];
    double *startingPoint = problem->getInitialParameters();
    if (startingPoint) {
        // copy, because will be update each iteration
        memcpy(parameters, startingPoint,
               sizeof(*parameters) * problem->getNumOptimizationParameters());
    } else {
        getRandomStartingpoint(
            problem->getParametersMin(), problem->getParametersMax(),
            problem->getNumOptimizationParameters(), parameters);
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

    ceres::GradientProblem ceresProblem(new MyCeresFirstOrderFunction(problem));

    ceres::GradientProblemSolver::Summary summary;
    ceres::Solve(options, ceresProblem, parameters, &summary);

    //    std::cout<<summary.FullReport();

    problem->logOptimizerFinished(summary.final_cost, parameters,
                                  summary.total_time_in_seconds,
                                  summary.termination_type);

    delete[] parameters;

    return summary.termination_type == ceres::FAILURE ||
           summary.termination_type == ceres::USER_FAILURE;
}
