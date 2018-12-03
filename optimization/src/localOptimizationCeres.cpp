#include "localOptimizationCeres.h"
#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include <logging.h>
#include <misc.h>

#include <ceres/ceres.h>

// !! Don't use. Leads to race conditions. Also: unable to assign sinks to specific ceres instances.
#undef PARPE_CERES_MINIGLOG_REDIRECT
#ifdef PARPE_CERES_MINIGLOG_REDIRECT
#include <ceres/../../../../internal/ceres/miniglog/glog/logging.h>
#endif

namespace parpe {

void setCeresOption(const std::pair<const std::string, const std::string> &pair, ceres::GradientProblemSolver::Options* options);

#ifdef PARPE_CERES_MINIGLOG_REDIRECT
/**
 * @brief LogSinkAdapter redirectsceres miniglog output to logging.cpp.
 */
class LogSinkAdapter : public google::LogSink {
public:
    LogSinkAdapter() {
        id = counter++;
    }
    void send(google::LogSeverity severity,
                      const char* full_filename,
                      const char* base_filename,
                      int line,
                      const struct tm* tm_time,
                      const char* message,
                      size_t message_len) override {
        // Map log levels
        loglevel lvl = LOGLVL_INFO;
        switch (severity) {
        case google::INFO:
            lvl = LOGLVL_INFO;
            break;
        case google::WARNING:
            lvl = LOGLVL_WARNING;
            break;
        case google::ERROR:
            lvl = LOGLVL_ERROR;
            break;
        case google::FATAL:
            lvl = LOGLVL_CRITICAL;
            break;
        }

        parpe::logmessage(lvl, "ceres #%d: %s", id, message);

    }
    void WaitTillSent() override {}

private:
    /** count instantiations */
    static int counter;
    /** prefix for logging output to identifiy concurrent optimizer runs */
    int id = 0;
};

int LogSinkAdapter::counter = 0;
#endif

/**
 * @brief Adapter class for parpe::OptimizationProblem and ceres::FirstOrderFunction
 */
class MyCeresFirstOrderFunction : public ceres::FirstOrderFunction {

  public:
    MyCeresFirstOrderFunction(OptimizationProblem *problem, OptimizationReporter *reporter)
        : problem(problem), reporter(reporter) {
        numParameters = problem->costFun->numParameters();

        // bounds are not natively supported by CERES; a naive check is currently implemented which fails function evaluation if parameters are out of bounds
        parametersMin.resize(numParameters);
        problem->fillParametersMin(parametersMin);
        parametersMax.resize(numParameters);
        problem->fillParametersMax(parametersMax);
    }

    /**
     * @brief Evaluate cost function
     * @param parameters
     * @param cost
     * @param gradient If not NULL, evaluate gradient
     * @return true on success, false otherwise
     */
    bool Evaluate(const double *parameters, double *cost,
                          double *gradient) const override {

        // Naive bounds check: report failure if not within
        if(!withinBounds(numParameters, parameters,
                         parametersMin.data(), parametersMax.data()))
            return false;

        auto result = reporter->evaluate(gsl::make_span<double const>(parameters, numParameters),
                                         *cost,
                                         gsl::make_span<double>(gradient, gradient?numParameters:0));

        return result == functionEvaluationSuccess;
    }

    int NumParameters() const override {
        return numParameters;
    }

  private:
    OptimizationProblem *problem;
    int numParameters = 0;

    // non-owning
    OptimizationReporter* reporter;

    std::vector<double> parametersMin;
    std::vector<double> parametersMax;
};


/**
 * @brief Callback functor for to be called between ceres iterations
 */
class MyIterationCallback : public ceres::IterationCallback {
  public:
    // Non-owning
    explicit MyIterationCallback(OptimizationReporter *reporter) : reporter(reporter) {}

    ceres::CallbackReturnType
    operator()(const ceres::IterationSummary &summary) override {

        // TODO: print here

        int status = reporter->iterationFinished(gsl::span<double const>(), summary.cost, gsl::span<double const>());
        switch (status) {
        case 0:
            return ceres::SOLVER_CONTINUE;
        default:
            return ceres::SOLVER_ABORT;
        }
    }

  private:
    OptimizationReporter *reporter = nullptr;
};


ceres::GradientProblemSolver::Options getCeresOptions(
                     OptimizationProblem *problem) {
    ceres::GradientProblemSolver::Options options;

    // don't: use vlog which we can redirect more easily
    //options.minimizer_progress_to_stdout =
    //    problem->getOptimizationOptions().printToStdout;
    options.max_num_iterations =
        problem->getOptimizationOptions().maxOptimizerIterations;

    problem->getOptimizationOptions().for_each<ceres::GradientProblemSolver::Options*>(setCeresOption, &options);

    return options;
}


std::tuple<int, double, std::vector<double> > OptimizerCeres::optimize(OptimizationProblem *problem) {
#ifdef PARPE_CERES_MINIGLOG_REDIRECT
    // Redirect ceres output (actually it's copied; can't remove ceres sink)
    LogSinkAdapter log;
    google::AddLogSink(&log);
#endif

    std::vector<double> parameters(problem->costFun->numParameters());
    problem->fillInitialParameters(parameters);

    auto reporter = problem->getReporter();
    // GradientProblem takes ownership of
    ceres::GradientProblem ceresProblem(new MyCeresFirstOrderFunction(problem, reporter.get()));

    ceres::GradientProblemSolver::Options options = getCeresOptions(problem);
    // Can use reporter from unique_ptr here, since callback will be destroyed first
    MyIterationCallback callback(reporter.get());
    options.callbacks.push_back(&callback);

    ceres::GradientProblemSolver::Summary summary;

    reporter->starting(gsl::span<double>(parameters));

    ceres::Solve(options, ceresProblem, parameters.data(), &summary);
    reporter->finished(summary.final_cost,
                       gsl::span<double>(parameters),
                       summary.termination_type);

    //    std::cout<<summary.FullReport();

#ifdef PARPE_CERES_MINIGLOG_REDIRECT
    google::RemoveLogSink(&log); // before going out of scope
#endif

    return std::tuple<int, double, std::vector<double> >(summary.termination_type == ceres::FAILURE ||
           summary.termination_type == ceres::USER_FAILURE, summary.final_cost, parameters);
}

/**
 * @brief Set string options to ceres options object.
 *
 * Used for iterating over OptimizationOptions map.
 * @param pair key => value pair
 * @param options
 */
void setCeresOption(const std::pair<const std::string, const std::string> &pair,
                    ceres::GradientProblemSolver::Options* options) {
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
        // NOTE: this currently disabled until shippable CI has ceres 1.12 debian package available. this option is not available in the current 1.11 package and breaks CI.
//    } else if(key == "parameter_tolerance") {
//        options->parameter_tolerance = std::stod(val);
    } else if(key == "logging_type") {
        options->logging_type = static_cast<ceres::LoggingType>(std::stoi(val));
    } else if(key == "minimizer_progress_to_stdout") {
        options->minimizer_progress_to_stdout = std::stoi(val);
    } else {
        logmessage(LOGLVL_WARNING, "Ignoring unknown optimization option %s.", key.c_str());
        return;
    }

    logmessage(LOGLVL_DEBUG, "Set optimization option %s to %s.", key.c_str(), val.c_str());
}


} // namespace parpe
