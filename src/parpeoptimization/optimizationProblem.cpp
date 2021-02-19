#include <parpeoptimization/optimizationProblem.h>

#include <parpecommon/logging.h>
#include <parpecommon/misc.h>
#include <parpecommon/parpeException.h>
#include <parpeoptimization/optimizationOptions.h>
#include <parpeoptimization/optimizer.h>
#include <parpeoptimization/minibatchOptimization.h>
#include <parpeoptimization/optimizationResultWriter.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cassert>
#include <numeric>
#include <random>
#include <memory>
#include <utility>

namespace parpe {


int getLocalOptimum(OptimizationProblem *problem) {
    // TODO how to make this nicer? minibatchOptimizer should not inherit
    // from Optimizer since they have different interfaces, so we can not
    // use the same factory method
    auto options = problem->getOptimizationOptions();
    if (options.optimizer == optimizerName::OPTIMIZER_MINIBATCH_1) {
        auto minibatchProblem =
                dynamic_cast<MinibatchOptimizationProblem<int>*>(problem);
        if (!minibatchProblem)
            throw ParPEException("Minibatch optimizer selected but given "
                                 "optimization problem cannot be solved by "
                                 "minibatch optimizer");
        auto status = runMinibatchOptimization(minibatchProblem);
        return std::get < 0 > (status);
    }

    auto optimizer = std::unique_ptr < Optimizer > (
                problem->getOptimizationOptions().createOptimizer());
    if (!optimizer)
        throw ParPEException("Invalid optimizer selected. Did you compile "
                             "parPE with support for the selected optimizer?");
    auto status = optimizer->optimize(problem);
    return std::get < 0 > (status);
}


void *getLocalOptimumThreadWrapper(void *optimizationProblemVp) {
    auto problem = static_cast<OptimizationProblem *>(optimizationProblemVp);
    auto *result = new int;
    *result = getLocalOptimum(problem);
    return result;
}

void optimizationProblemGradientCheckMultiEps(OptimizationProblem *problem,
                                      int numParameterIndicesToCheck
                                      ) {
    // set eps
    std::vector<double> multi_eps {1e-1, 1e-3, 1e-5, 1e-7, 1e-9};

    // setting the number of parameters to the minimum of
    // numParamaterIndicesToCheck and dimension of the problem
    int numParameters = problem->cost_fun_->numParameters();
    numParameterIndicesToCheck =
        std::min(numParameterIndicesToCheck, numParameters);
    // choose random parameters to check
    std::vector<int> parameterIndices(numParameterIndicesToCheck);
    std::iota(parameterIndices.begin(), parameterIndices.end(), 0);
    std::random_device rd;
    std::mt19937 g(rd());

    optimizationProblemGradientCheckMultiEps(problem, parameterIndices, multi_eps);
}

void optimizationProblemGradientCheckMultiEps(OptimizationProblem *problem,
                                              gsl::span<const int> parameterIndices,
                                              gsl::span<double> multi_eps){
    //get a random theta
    std::vector<double> theta(problem->cost_fun_->numParameters());
    problem->fillInitialParameters(theta);
    double fc = 0; // f(theta)
    //evaluate the objective function at theta and get analytical gradient
    std::vector<double> gradient(theta.size());
    problem->cost_fun_->evaluate(theta, fc, gradient);

    std::vector<double> thetaTmp(theta);

    for(int curInd : parameterIndices) {
        double eps_best = 0.0;
        double regRelErr_best = 0.0;
        double absErr_best = 0.0;
        double fd_c_best = 0.0;

        // getting the analytical gradient of current index
        double curGrad = gradient[curInd];

        for(double epsilon : multi_eps) {
            double fb = 0.0; // f(theta - eps)
            double ff = 0.0; // f(theta + eps)

            thetaTmp[curInd] = theta[curInd] + epsilon;
            problem->cost_fun_->evaluate(gsl::span<double>(thetaTmp), ff,
                                         gsl::span<double>());

            thetaTmp[curInd] = theta[curInd] - epsilon;
            problem->cost_fun_->evaluate(gsl::span<double>(thetaTmp), fb,
                                         gsl::span<double>());
            // calculating the finite difference
            double fd_c = (ff - fb) / (2 * epsilon);
            //reverting thetaTmp back to original
            thetaTmp[curInd] = theta[curInd];

            double abs_err = std::fabs(curGrad - fd_c);
            double regRelError = std::fabs(abs_err / (fd_c + epsilon));

            // comparing results with current best epsilon
            // if better, replace. Also replace if no eps currently saved.
            if((regRelError < regRelErr_best) || (regRelErr_best == 0)){
                eps_best = epsilon;
                regRelErr_best = regRelError;
                fd_c_best = fd_c;
                absErr_best = abs_err;
            }
        }
        loglevel ll = LOGLVL_INFO;
        if (fabs(regRelErr_best) > 1e-3)
            ll = LOGLVL_WARNING;
        if (fabs(regRelErr_best) > 1e-2)
            ll = LOGLVL_ERROR;
        logmessage(ll, "%5d g: %12.6g  fd_c: %12.6g  |Δ/fd_c|: %.6e  |Δ|: %12.6g"
                       "  ϵ: %12.6g ",
                   curInd, curGrad, fd_c_best, regRelErr_best, absErr_best, eps_best);
    }
}

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      int numParameterIndicesToCheck,
                                      double epsilon) {
    int numParameters = problem->cost_fun_->numParameters();
    numParameterIndicesToCheck =
            std::min(numParameterIndicesToCheck, numParameters);
    // choose random parameters to check
    std::vector<int> parameterIndices(numParameterIndicesToCheck);
    std::iota(parameterIndices.begin(), parameterIndices.end(), 0);
    std::random_device rd;
    std::mt19937 g(rd());
    //std::shuffle(parameterIndices.begin(), parameterIndices.end(), g);

    optimizationProblemGradientCheck(problem, parameterIndices, epsilon);
}


void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      gsl::span<const int> parameterIndices,
                                      double epsilon) {
    double fc = 0; // f(theta)
    std::vector<double> theta(problem->cost_fun_->numParameters());
    problem->fillInitialParameters(theta);

    std::vector<double> gradient(theta.size());
    problem->cost_fun_->evaluate(theta, fc, gradient);

    std::vector<double> thetaTmp(theta);

    // printf("Index\tGradient\tfd_f\t\t(delta)\t\tfd_c\t\t(delta)\t\tfd_b\t\t(delta)\n");

    for (int curInd : parameterIndices) {
        double fb = 0, ff = 0; // f(theta + eps) , f(theta - eps)

        thetaTmp[curInd] = theta[curInd] + epsilon;
        problem->cost_fun_->evaluate(gsl::span<double>(thetaTmp), ff,
                                   gsl::span<double>());

        thetaTmp[curInd] = theta[curInd] - epsilon;
        problem->cost_fun_->evaluate(gsl::span<double>(thetaTmp), fb,
                                   gsl::span<double>());

        // double fd_f = (ff - fc) / epsilon;

        // double fd_b = (fc - fb) / epsilon;

        double fd_c = (ff - fb) / (2 * epsilon);

        thetaTmp[curInd] = theta[curInd];

        double curGrad = gradient[curInd];

        // check
        //        char status[] = " ";
        //        if(!((fd_c >= fd_f && fd_c <= fd_b)
        //             || (fd_c <= fd_f && fd_c >= fd_b)))
        //            status[0] = '!';
        //        if(!((curGrad >= fd_f && curGrad <= fd_b)
        //             || (curGrad <= fd_f && curGrad >= fd_b)))
        //            status[0] = '!';

        //        printf("%5d%s g: %12.6g fd_f: %12.6g (Δ%12.6g) fd_c: %12.6g (Δ%12.6g) fd_b: %12.6g (Δ%12.6g)",
        //               curInd, status, curGrad,
        //               fd_f, curGrad - fd_f,
        //               fd_c, curGrad - fd_c,
        //               fd_b, curGrad - fd_b);
        //        printf("fb: %12.6g fc: %12.6g ff: %12.6g\n", fb, fc, ff);

        double reg = 1e-5;
        double regRelError = (curGrad - fd_c) / (ff + reg);
        loglevel ll = LOGLVL_INFO;
        if (fabs(regRelError) > 1e-3)
            ll = LOGLVL_WARNING;
        if (fabs(regRelError) > 1e-2)
            ll = LOGLVL_ERROR;

        logmessage(ll, "%5d g: %12.6g  fd_c: %12.6g  Δ/ff: %.6e  f: %12.6g",
                   curInd, curGrad, fd_c, regRelError, ff);
    }
}

OptimizationProblem::OptimizationProblem(
        std::unique_ptr<GradientFunction> costFun,
        std::unique_ptr<Logger> logger)
    : cost_fun_(std::move(costFun)), logger_(std::move(logger)) {

}

const OptimizationOptions &OptimizationProblem::getOptimizationOptions() const {
    return optimization_options_;
}

void OptimizationProblem::setOptimizationOptions(const OptimizationOptions &options) {
    optimization_options_ = options;
}

std::unique_ptr<OptimizationReporter> OptimizationProblem::getReporter() const {
    return std::make_unique < OptimizationReporter > (
                cost_fun_.get(), std::make_unique < Logger > (*logger_));
}

void OptimizationProblem::fillInitialParameters(gsl::span<double> buffer) const {
    int numParameters = cost_fun_->numParameters();
    std::vector<double> parametersMin(numParameters);
    std::vector<double> parametersMax(numParameters);
    fillParametersMin(parametersMin);
    fillParametersMax(parametersMax);

    fillArrayRandomDoubleIndividualInterval(parametersMin, parametersMax,
                                            buffer);
}

OptimizationReporter::OptimizationReporter(GradientFunction *gradFun,
                                           std::unique_ptr<Logger> logger) :
        OptimizationReporter(gradFun, nullptr, std::move(logger)) {
    default_logger_prefix_ = this->logger_->getPrefix();
}

OptimizationReporter::OptimizationReporter(GradientFunction *gradFun,
                                           std::unique_ptr<OptimizationResultWriter> rw,
                                           std::unique_ptr<Logger> logger) :
        result_writer_(std::move(rw)), logger_(std::move(logger)) {
    setGradientFunction(gradFun);
    default_logger_prefix_ = this->logger_->getPrefix();
}

FunctionEvaluationStatus OptimizationReporter::evaluate(gsl::span<const double> parameters,
                                                        double &fval,
                                                        gsl::span<double> gradient,
                                                        Logger *logger,
                                                        double *cpuTime) const {
    double functionCpuSec = 0.0;
    if (cpuTime)
        *cpuTime = 0.0;

    if (beforeCostFunctionCall(parameters) != 0)
        return functionEvaluationFailure;

    if (gradient.data()) {
        if (!have_cached_gradient_
                || !std::equal(parameters.begin(), parameters.end(),
                               cached_parameters_.begin())) {
            // Have to compute anew
            cached_status_ = grad_fun_->evaluate(
                        parameters, cached_cost_, cached_gradient_,
                        logger ? logger : this->logger_.get(), &functionCpuSec);
            have_cached_cost_ = true;
            have_cached_gradient_ = true;
        }
        // recycle old result
        std::copy(cached_gradient_.begin(), cached_gradient_.end(), gradient.begin());
        fval = cached_cost_;
    } else {
        if (!have_cached_cost_
                || !std::equal(parameters.begin(), parameters.end(),
                               cached_parameters_.begin())) {
            // Have to compute anew
            cached_status_ = grad_fun_->evaluate(
                        parameters, cached_cost_, gsl::span<double>(),
                        logger ? logger : this->logger_.get(), &functionCpuSec);
            have_cached_cost_ = true;
            have_cached_gradient_ = false;
        }
        fval = cached_cost_;
    }

    // update cached parameters
    cached_parameters_.resize(num_parameters_);
    std::copy(parameters.begin(), parameters.end(), cached_parameters_.begin());

    cpu_time_iteration_sec_ += functionCpuSec;
    cpu_time_total_sec_ += functionCpuSec;
    if (cpuTime)
        *cpuTime = functionCpuSec;

    if (afterCostFunctionCall(
                parameters, cached_cost_,
                gradient.data() ? cached_gradient_ : gsl::span<double>()) != 0)
        return functionEvaluationFailure;

    return cached_status_;
}

int OptimizationReporter::numParameters() const {
    return grad_fun_->numParameters();
}

void OptimizationReporter::printObjectiveFunctionFailureMessage() const {
    if (logger_)
        logger_->logmessage(LOGLVL_ERROR, "Objective function evaluation failed!");
}

bool OptimizationReporter::starting(gsl::span<const double> initialParameters) const {
    // If this is called multiple times (happens via IpOpt), don't do anything
    if (started_)
        return false;

    wall_timer_.reset();

    if (result_writer_)
        result_writer_->starting(initialParameters);

    started_ = true;

    logger_->setPrefix(default_logger_prefix_ + "i" + std::to_string(num_iterations_));

    return false;
}

bool OptimizationReporter::iterationFinished(gsl::span<const double> parameters,
                                             double objectiveFunctionValue,
                                             gsl::span<double const> objectiveFunctionGradient) const {
    double wallTimeIter = wall_timer_.getRound();
    double wallTimeOptim = wall_timer_.getTotal();

    if (logger_)
        logger_->logmessage(LOGLVL_INFO,
                           "iter: %d cost: %g time_iter: wall: %gs cpu: %gs time_optim: wall: %gs cpu: %gs",
                           num_iterations_, objectiveFunctionValue, wallTimeIter,
                           cpu_time_iteration_sec_, wallTimeOptim,
                           cpu_time_total_sec_);

    if (result_writer_)
        result_writer_->logOptimizerIteration(
                num_iterations_, parameters.empty() ? cached_parameters_ : parameters, objectiveFunctionValue,
                // This might be misleading, the gradient could evaluated at other parameters if there was a line search inbetween
                objectiveFunctionGradient.empty() ? cached_gradient_ : objectiveFunctionGradient, wallTimeIter,
                cpu_time_iteration_sec_);
    ++num_iterations_;

    logger_->setPrefix(default_logger_prefix_ + "i" + std::to_string(num_iterations_));

    cpu_time_iteration_sec_ = 0.0;

    return false;
}

bool OptimizationReporter::beforeCostFunctionCall(gsl::span<const double> /*parameters*/) const {
    ++num_function_calls_;
    //timeCostEvaluationBegin = clock();

    return false;
}

bool OptimizationReporter::afterCostFunctionCall(gsl::span<const double> parameters,
                                                 double objectiveFunctionValue,
                                                 gsl::span<const double> objectiveFunctionGradient) const {
    double wallTime = wall_timer_.getTotal();

    if (!std::isfinite(objectiveFunctionValue))
        printObjectiveFunctionFailureMessage();

    if (result_writer_) {
        result_writer_->logObjectiveFunctionEvaluation(
                    parameters, objectiveFunctionValue,
                    objectiveFunctionGradient, num_iterations_, num_function_calls_,
                    wallTime);
    }
    return false;
}

void OptimizationReporter::finished(double optimalCost,
                                    gsl::span<const double> parameters,
                                    int exitStatus) const {
    double timeElapsed = wall_timer_.getTotal();

    if ((optimalCost <= cached_cost_ || cached_parameters_.empty()) && !parameters.empty()) {
        cached_cost_ = optimalCost;
        cached_parameters_.assign(parameters.begin(), parameters.end());
    } else if (cached_cost_ > optimalCost && parameters.empty()) {
        // the optimal value is not from the cached parameters and we did not get
        // the optimal parameters from the optimizer. since we don't know them, rather set to nan
        if (logger_)
            logger_->logmessage(LOGLVL_INFO, "cachedCost != optimalCost && parameters.empty()");
        cached_parameters_.assign(cached_parameters_.size(), NAN);
        cached_cost_ = optimalCost;
    } // else: our cached parameters were better. use those

    if (logger_)
        logger_->logmessage(LOGLVL_INFO, "Optimizer status %d, final llh: %e, "
                                        "time: wall: %f cpu: %f.", exitStatus,
                           cached_cost_, timeElapsed, cpu_time_total_sec_);

    if (result_writer_)
        result_writer_->saveOptimizerResults(cached_cost_, cached_parameters_,
                                           timeElapsed, cpu_time_total_sec_,
                                           exitStatus);
}

double OptimizationReporter::getFinalCost() const {
    return cached_cost_;
}

const std::vector<double> &OptimizationReporter::getFinalParameters() const {
    return cached_parameters_;
}

void OptimizationReporter::setGradientFunction(GradientFunction *gradFun) const {
    this->grad_fun_ = gradFun;
    num_parameters_ = gradFun->numParameters();
    cached_gradient_.resize(num_parameters_);
}

std::vector<std::string> OptimizationReporter::getParameterIds() const
{
    return grad_fun_->getParameterIds();
}

void OptimizationProblemImpl::fillParametersMin(gsl::span<double> buffer) const {
    std::copy(parametersMin.begin(), parametersMin.end(), buffer.begin());
}

void OptimizationProblemImpl::fillParametersMax(gsl::span<double> buffer) const {
    std::copy(parametersMax.begin(), parametersMax.end(), buffer.begin());
}

void OptimizationProblemImpl::setParametersMin(std::vector<double> parametersMin) {
    this->parametersMin = std::move(parametersMin);
}

void OptimizationProblemImpl::setParametersMax(std::vector<double> parametersMax) {
    this->parametersMax = std::move(parametersMax);
}

void OptimizationProblemImpl::setInitialParameters(std::vector<double> initial) {
    parametersStart = std::move(initial);
}

void OptimizationProblemImpl::fillInitialParameters(gsl::span<double> buffer) const {
    if (!parametersStart.empty()) {
        std::copy(parametersStart.begin(), parametersStart.end(), buffer.begin());
    } else {
        OptimizationProblem::fillInitialParameters(buffer);
    }
}

} // namespace parpe
