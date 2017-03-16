#include "localOptimizationCeres.hpp"

#include "ceres/ceres.h"

class CeresWrapper : public ceres::FirstOrderFunction {

public:
    CeresWrapper(OptimizationProblem *problem) : problem(problem) {}

    virtual bool Evaluate(const double* parameters,
                          double* cost,
                          double* gradient) const {
        bool status;

        if(gradient)
            status = problem->objectiveFunctionGradient(problem, parameters, cost, gradient);
        else
            status = problem->objectiveFunction(problem, parameters, cost);

        return status == 0;
    }

    virtual int NumParameters() const { return problem->numOptimizationParameters; }

private:
    OptimizationProblem *problem;
};

int getLocalOptimumCeres(OptimizationProblem *problem)
{
    double *parameters = (double *) malloc(sizeof(*parameters) * problem->numOptimizationParameters);
    // copy, because will be update each iteration
    memcpy(parameters, problem->initialParameters, sizeof(*parameters) * problem->numOptimizationParameters);

    ceres::GradientProblemSolver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = problem->maxOptimizerIterations;
//    options.gradient_tolerance = 1e-18;
//    options.function_tolerance = 1e-18;

    ceres::GradientProblemSolver::Summary summary;

    ceres::GradientProblem ceresProblem(new CeresWrapper(problem));
    ceres::Solve(options, ceresProblem, parameters, &summary);

//    std::cout<<summary.FullReport();

    problem->logOptimizerFinished(problem, summary.final_cost, parameters, summary.total_time_in_seconds, summary.termination_type);

    free(parameters);

    return summary.termination_type > 0;
}
