#include "localOptimizationCeres.h"

#include "ceres/ceres.h"

class CeresWrapper : public ceres::FirstOrderFunction {

public:
    CeresWrapper(OptimizationProblem *problem) : problem(problem) {}

    virtual bool Evaluate(const double* parameters,
                          double* cost,
                          double* gradient) const {
        bool status = problem->evaluateObjectiveFunction(parameters, cost, gradient);

        return status == 0;
    }

    virtual int NumParameters() const { return problem->numOptimizationParameters; }

private:
    OptimizationProblem *problem;
};

int getLocalOptimumCeres(OptimizationProblem *problem)
{
    double *parameters = (double *) malloc(sizeof(*parameters) * problem->numOptimizationParameters);
    if(problem->initialParameters) {
        // copy, because will be update each iteration
        memcpy(parameters, problem->initialParameters, sizeof(*parameters) * problem->numOptimizationParameters);
    } else {
        getRandomStartingpoint(problem->parametersMin, problem->parametersMax, problem->numOptimizationParameters, parameters);
    }

    ceres::GradientProblemSolver::Options options;
    options.minimizer_progress_to_stdout = problem->printToStdout;
    options.max_num_iterations = problem->maxOptimizerIterations;
//    options.gradient_tolerance = 1e-18;
//    options.function_tolerance = 1e-18;

    ceres::GradientProblemSolver::Summary summary;

    ceres::GradientProblem ceresProblem(new CeresWrapper(problem));
    ceres::Solve(options, ceresProblem, parameters, &summary);

//    std::cout<<summary.FullReport();

    problem->logOptimizerFinished(summary.final_cost, parameters, summary.total_time_in_seconds, summary.termination_type);

    free(parameters);

    return summary.termination_type > 0;
}
