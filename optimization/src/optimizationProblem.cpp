#include "optimizationProblem.h"
#include "localOptimizationCeres.h"
#include "localOptimizationIpopt.h"
#include "logging.h"
#include "misc.h"
#include "optimizationOptions.h"
#include <optimizer.h>
#include <parpeException.h>
#include <minibatchOptimization.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cassert>
#include <numeric>
#include <random>
#include <memory>

namespace parpe {


/**
 * @brief getLocalOptimum
 * @param problem
 * @return int indicating status. 0: success, != 0: failure
 */

int getLocalOptimum(OptimizationProblem *problem) {
    // TODO how to make this nicer? minibatchOptimizer should not inherit
    // from Optimizer since they have different interfaces, so we can not
    // use the same factory method
    auto options = problem->getOptimizationOptions();
    if(options.optimizer == optimizerName::OPTIMIZER_MINIBATCH_1) {
        auto minibatchProblem = dynamic_cast<MinibatchOptimizationProblem<int>*>(problem);
        if(!minibatchProblem)
            throw ParPEException("Minibatch optimizer selected but given optimization problem cannot be solved by minibatch optimizer");
        auto status = runMinibatchOptimization(minibatchProblem);
        return std::get<0>(status);
    }

    auto optimizer = std::unique_ptr<Optimizer>(problem->getOptimizationOptions().createOptimizer());
    if(!optimizer)
        throw ParPEException("Invalid optimizer selected. Did you compile parPE with support for the selected optimizer?");
    auto status = optimizer->optimize(problem);
    return std::get<0>(status);
}

/**
 * @brief getLocalOptimumThreadWrapper wrapper for using getLocalOptimum with
 * pThreads.
 * @param problem
 * @return Pointer to int indicating status. 0: success, != 0: failure
 */

void *getLocalOptimumThreadWrapper(void *optimizationProblemVp) {
    auto problem = static_cast<OptimizationProblem *>(optimizationProblemVp);
    auto *result = new int;
    *result = getLocalOptimum(problem);
    return result;
}

//void runOptimizationsParallel(OptimizationProblem **problems,
//                              int numProblems) {
//    runInParallelAndWaitForFinish(getLocalOptimumThreadWrapper,
//                                  reinterpret_cast<void **>(problems), numProblems);
//}

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      int numParameterIndicesToCheck, double epsilon) {
    int numParameters = problem->costFun->numParameters();
    numParameterIndicesToCheck = std::min(numParameterIndicesToCheck, numParameters);
    // choose random parameters to check
    std::vector<int> parameterIndices(numParameters);
    std::iota(parameterIndices.begin(), parameterIndices.end(), 0);
    std::random_device rd;
    std::mt19937 g(rd());
    //std::shuffle(parameterIndices.begin(), parameterIndices.end(), g);

    optimizationProblemGradientCheck(problem, parameterIndices, epsilon);
}

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      gsl::span<const int> parameterIndices, double epsilon) {
    double fc = 0; // f(theta)
    std::vector<double> theta(problem->costFun->numParameters());
    problem->fillInitialParameters(theta);

    std::vector<double> gradient(theta.size());
    problem->costFun->evaluate(theta, fc, gradient);

    std::vector<double> thetaTmp(theta);

    // printf("Index\tGradient\tfd_f\t\t(delta)\t\tfd_c\t\t(delta)\t\tfd_b\t\t(delta)\n");

    for (int i = 0; i < parameterIndices.size(); ++i) {
        int curInd = parameterIndices[i];
        double fb = 0, ff = 0; // f(theta + eps) , f(theta - eps)

        thetaTmp[curInd] = theta[curInd] + epsilon;
        problem->costFun->evaluate(gsl::span<double>(thetaTmp), ff, gsl::span<double>());

        thetaTmp[curInd] = theta[curInd] - epsilon;
        problem->costFun->evaluate(gsl::span<double>(thetaTmp), fb, gsl::span<double>());

        double fd_f = (ff - fc) / epsilon;

        double fd_b = (fc - fb) / epsilon;

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
        double regRelError = (curGrad - fd_c) / (fd_c + reg);
        loglevel ll = LOGLVL_INFO;
        if(fabs(regRelError) > reg)
            ll = LOGLVL_WARNING;
        if(fabs(regRelError) > 0.1)
            ll = LOGLVL_ERROR;

        logmessage(ll, "%5d g: %12.6g  fd_c: %12.6g  Δ/fd_c: %.6e  f: %12.6g",
               curInd, curGrad, fd_c, regRelError, fc);

    }
}


OptimizationProblem::OptimizationProblem(std::unique_ptr<GradientFunction> costFun, std::unique_ptr<Logger> logger)
    : costFun(std::move(costFun)),
      logger(std::move(logger))
{

}


const OptimizationOptions &OptimizationProblem::getOptimizationOptions() const
{
    return optimizationOptions;
}

void OptimizationProblem::setOptimizationOptions(const OptimizationOptions &options)
{
    optimizationOptions = options;
}

std::unique_ptr<OptimizationReporter> OptimizationProblem::getReporter() const
{
    return std::make_unique<OptimizationReporter>(costFun.get(), std::make_unique<Logger>(*logger));
}


void OptimizationProblem::fillInitialParameters(gsl::span<double> buffer) const {
    int numParameters = costFun->numParameters();
    std::vector<double> parametersMin(numParameters);
    std::vector<double> parametersMax(numParameters);
    fillParametersMin(parametersMin);
    fillParametersMax(parametersMax);

    fillArrayRandomDoubleIndividualInterval(parametersMin.data(), parametersMax.data(),
                                            numParameters, buffer.data());
}


OptimizationReporter::OptimizationReporter(GradientFunction *gradFun, std::unique_ptr<Logger> logger)
    : OptimizationReporter(gradFun, nullptr, std::move(logger))
{
    defaultLoggerPrefix = this->logger->getPrefix();
}

OptimizationReporter::OptimizationReporter(GradientFunction *gradFun,
        std::unique_ptr<OptimizationResultWriter> rw, std::unique_ptr<Logger> logger)
    : resultWriter(std::move(rw)),
      logger(std::move(logger))
{
    setGradientFunction(gradFun);
    defaultLoggerPrefix = this->logger->getPrefix();
}

FunctionEvaluationStatus OptimizationReporter::evaluate(
        gsl::span<const double> parameters,
        double &fval,
        gsl::span<double> gradient,
        Logger *logger,
        double *cpuTime) const
{  
    double functionCpuSec = 0.0;
    if(cpuTime)
        *cpuTime = 0.0;

    if(beforeCostFunctionCall(parameters) != 0)
        return functionEvaluationFailure;

    if(gradient.data()) {
        if (!haveCachedGradient || !std::equal(parameters.begin(), parameters.end(),
                                               cachedParameters.begin())) {
            // Have to compute anew
            cachedStatus = gradFun->evaluate(parameters, cachedCost, cachedGradient, logger, &functionCpuSec);
            haveCachedCost = true;
            haveCachedGradient = true;
        }
        // recycle old result
        std::copy(cachedGradient.begin(), cachedGradient.end(), gradient.begin());
        fval = cachedCost;
    } else {
        if (!haveCachedCost || !std::equal(parameters.begin(), parameters.end(),
                                           cachedParameters.begin())) {
            // Have to compute anew
            cachedStatus = gradFun->evaluate(parameters, cachedCost, gsl::span<double>(), logger, &functionCpuSec);
            haveCachedCost = true;
            haveCachedGradient = false;
        }
        fval = cachedCost;
    }

    // update cached parameters
    cachedParameters.resize(numParameters_);
    std::copy(parameters.begin(), parameters.end(), cachedParameters.begin());

    cpuTimeIterationSec += functionCpuSec;
    cpuTimeTotalSec += functionCpuSec;
    if(cpuTime)
        *cpuTime = functionCpuSec;

    if(afterCostFunctionCall(parameters, cachedCost, gradient.data() ? cachedGradient : gsl::span<double>()) != 0)
        return functionEvaluationFailure;

    return cachedStatus;
}

int OptimizationReporter::numParameters() const {
    return gradFun->numParameters();
}

void OptimizationReporter::printObjectiveFunctionFailureMessage() const
{
    if(logger)
        logger->logmessage(LOGLVL_ERROR, "Objective function evaluation failed!");
}

bool OptimizationReporter::starting(gsl::span<const double> initialParameters) const
{
    // If this is called multiple times (happens via IpOpt), don't do anything
    if(started)
        return false;

    wallTimer.reset();
//    timeOptimizationBegin = clock();
//    timeIterationBegin = timeOptimizationBegin;
//    timeCostEvaluationBegin = timeOptimizationBegin;

    if(resultWriter)
        resultWriter->starting(initialParameters);

    started = true;

    logger->setPrefix(defaultLoggerPrefix + "i" + std::to_string(numIterations));

    return false;
}

bool OptimizationReporter::iterationFinished(gsl::span<const double> parameters,
                                             double objectiveFunctionValue,
                                             gsl::span<double const> objectiveFunctionGradient) const
{
    double wallTimeIter = wallTimer.getRound(); //(double)(clock() - timeIterationBegin) / CLOCKS_PER_SEC;
    double wallTimeOptim = wallTimer.getTotal(); //double)(clock() - timeOptimizationBegin) / CLOCKS_PER_SEC;

    if(logger)
        logger->logmessage(LOGLVL_INFO, "iter: %d cost: %g time_iter: wall: %gs cpu: %gs time_optim: wall: %gs cpu: %gs",
                           numIterations, objectiveFunctionValue,
                           wallTimeIter, cpuTimeIterationSec, wallTimeOptim, cpuTimeTotalSec);

    if(resultWriter)
        resultWriter->logOptimizerIteration(
                    numIterations,
                    parameters.empty() ? cachedParameters : parameters,
                    objectiveFunctionValue,
                    objectiveFunctionGradient.empty() ? cachedGradient : objectiveFunctionGradient, // This might be misleading, the gradient could evaluated at other parameters if there was a line search inbetween
                    wallTimeIter, cpuTimeIterationSec, logGradientEachIteration);
    ++numIterations;

    logger->setPrefix(defaultLoggerPrefix + "i" + std::to_string(numIterations));

    cpuTimeIterationSec = 0.0;

    return false;
}

bool OptimizationReporter::beforeCostFunctionCall(gsl::span<const double>  /*parameters*/) const
{
    ++numFunctionCalls;
    //timeCostEvaluationBegin = clock();

    return false;
}

bool OptimizationReporter::afterCostFunctionCall(gsl::span<const double> parameters,
                                                 double objectiveFunctionValue,
                                                 gsl::span<const double> objectiveFunctionGradient) const
{
    //clock_t timeCostEvaluationEnd = clock();

    double wallTime = wallTimer.getTotal();//(double)(timeCostEvaluationEnd - timeCostEvaluationBegin) / CLOCKS_PER_SEC;

    if(!std::isfinite(objectiveFunctionValue))
        printObjectiveFunctionFailureMessage();

    if(resultWriter) {
        resultWriter->logObjectiveFunctionEvaluation(parameters,
                                                     objectiveFunctionValue,
                                                     objectiveFunctionGradient,
                                                     numIterations,
                                                     numFunctionCalls,
                                                     wallTime,
													 logGradientEachFunctionEvaluation,
													 logParametersEachFunctionEvaluation);
    }
    return false;
}

void OptimizationReporter::finished(double optimalCost, gsl::span<const double> parameters, int exitStatus) const
{
    double timeElapsed = wallTimer.getTotal();

    if((optimalCost <= cachedCost || cachedParameters.empty()) && !parameters.empty()) {
        cachedCost = optimalCost;
        cachedParameters.assign(parameters.begin(), parameters.end());
    } else if(cachedCost > optimalCost && parameters.empty()) {
        // the optimal value is not from the cached parameters and we did not get
        // the optimal parameters from the optimizer. since we don't know them, rather set to nan
        if(logger) logger->logmessage(LOGLVL_INFO, "cachedCost != optimalCost && parameters.empty()");
        cachedParameters.assign(cachedParameters.size(), NAN);
        cachedCost = optimalCost;
    } // else: our cached parameters were better. use those

    if(logger)
        logger->logmessage(LOGLVL_INFO, "Optimizer status %d, final llh: %e, time: wall: %f cpu: %f.",
               exitStatus, cachedCost, timeElapsed, cpuTimeTotalSec);


    if(resultWriter)
        resultWriter->saveOptimizerResults(cachedCost, cachedParameters,
                                           timeElapsed, cpuTimeTotalSec, exitStatus);
}

void OptimizationReporter::setLoggingEachIteration(bool logGradient) const
{
	logGradientEachIteration = logGradient;
}

void OptimizationReporter::setLoggingEachFunctionEvaluation(bool logGradient, bool logParameters) const
{
    logGradientEachFunctionEvaluation = logGradient;
    logParametersEachFunctionEvaluation = logParameters;
}

double OptimizationReporter::getFinalCost() const
{
    return cachedCost;
}

const std::vector<double> &OptimizationReporter::getFinalParameters() const
{
    return cachedParameters;
}

void OptimizationReporter::setGradientFunction(GradientFunction *gradFun) const
{
    this->gradFun = gradFun;
    numParameters_ = gradFun->numParameters();
    cachedGradient.resize(numParameters_);
}

} // namespace parpe
