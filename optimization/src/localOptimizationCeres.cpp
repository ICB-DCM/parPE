#include "localOptimizationCeres.h"
#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include <ceres/ceres.h>
#include <logging.h>

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

ceres::GradientProblemSolver::Options getCeresOptions(
                     OptimizationProblem *problem) {
    ceres::GradientProblemSolver::Options options;

    options.minimizer_progress_to_stdout =
        problem->optimizationOptions->printToStdout;
    options.max_num_iterations =
        problem->optimizationOptions->maxOptimizerIterations;
    //    options.gradient_tolerance = 1e-18;
    options.function_tolerance =
        problem->optimizationOptions->functionTolerance;

    // Use quadratic approximation to not require gradient for line search
    // NOTE: WOLFE line_search_type will always require gradient
    // options.line_search_interpolation_type = ceres::QUADRATIC;

    problem->optimizationOptions->for_each(setCeresOption, &options);

    return options;
}

int OptimizerCeres::optimize(OptimizationProblem *problem) {
    std::vector<double> parameters(problem->getNumOptimizationParameters());

    double *startingPoint = problem->getInitialParameters();
    if (startingPoint) {
        // copy, because will be update each iteration
        memcpy(parameters.data(), startingPoint,
               sizeof(double) * parameters.size());
    } else {
        problem->fillInitialParameters(parameters.data());
    }

    ceres::GradientProblem ceresProblem(new MyCeresFirstOrderFunction(problem));

    ceres::GradientProblemSolver::Options options = getCeresOptions(problem);
    MyIterationCallback callback(problem);
    options.callbacks.push_back(&callback);

    ceres::GradientProblemSolver::Summary summary;
    ceres::Solve(options, ceresProblem, parameters.data(), &summary);

    //    std::cout<<summary.FullReport();

    problem->logOptimizerFinished(summary.final_cost, parameters.data(),
                                  summary.total_time_in_seconds,
                                  summary.termination_type);

    return summary.termination_type == ceres::FAILURE ||
           summary.termination_type == ceres::USER_FAILURE;
}


void setCeresOption(const std::pair<const std::string, const std::string> &pair, void* arg) {
    // for iterating over OptimizationOptions

    auto options = static_cast<ceres::GradientProblemSolver::Options*>(arg);

    const std::string &key = pair.first;
    const std::string &val = pair.second;

    // TODO: set enums from string

    if(key == "line_search_direction_type") {
        options->line_search_direction_type =
                static_cast<ceres::LineSearchDirectionType>(std::stoi(val));
    } else if(key == "line_search_type") {
        options->line_search_type =
                static_cast<ceres::LineSearchType>(std::stoi(val));
    } else if(key == "nonlinear_conjugate_gradient_type") {
        options->nonlinear_conjugate_gradient_type =
                static_cast<ceres::NonlinearConjugateGradientType>(std::stoi(val));
    } else if(key == "max_lbfgs_rank") {
        options->max_lbfgs_rank = std::stoi(val);
    } else if(key == "use_approximate_eigenvalue_bfgs_scaling") {
        options->use_approximate_eigenvalue_bfgs_scaling = std::stoi(val);
    } else if(key == "line_search_interpolation_type") {
        options->line_search_interpolation_type =
                static_cast<ceres::LineSearchInterpolationType>(std::stoi(val));
    } else if(key == "min_line_search_step_size") {
        options->min_line_search_step_size = std::stod(val);
    } else if(key == "line_search_sufficient_function_decrease") {
        options->line_search_sufficient_function_decrease = std::stod(val);
    } else if(key == "max_line_search_step_contraction") {
        options->max_line_search_step_contraction = std::stod(val);
    } else if(key == "min_line_search_step_contraction") {
        options->min_line_search_step_contraction = std::stod(val);
    } else if(key == "max_num_line_search_step_size_iterations") {
        options->max_num_line_search_step_size_iterations = std::stoi(val);
    } else if(key == "max_num_line_search_direction_restarts") {
        options->max_num_line_search_direction_restarts = std::stoi(val);
    } else if(key == "line_search_sufficient_curvature_decrease") {
        options->line_search_sufficient_curvature_decrease = std::stod(val);
    } else if(key == "max_line_search_step_expansion") {
        options->max_line_search_step_expansion = std::stod(val);
    } else if(key == "max_num_iterations") {
        options->max_num_iterations = std::stoi(val);
    } else if(key == "max_solver_time_in_seconds") {
        options->max_solver_time_in_seconds = std::stod(val);
    } else if(key == "function_tolerance") {
        options->function_tolerance = std::stod(val);
    } else if(key == "gradient_tolerance") {
        options->gradient_tolerance = std::stod(val);
    } else if(key == "parameter_tolerance") {
        options->parameter_tolerance = std::stod(val);
    } else if(key == "logging_type") {
        options->logging_type = static_cast<ceres::LoggingType>(std::stoi(val));
    } else if(key == "minimizer_progress_to_stdout") {
        options->minimizer_progress_to_stdout = std::stoi(val);
    } else {
        logmessage(LOGLVL_WARNING, "Ignoring unknown optimization option %s.", key);
    }
}
