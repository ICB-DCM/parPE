#include "optimizationProblem.h"
#include "localOptimizationCeres.h"
#include "localOptimizationIpopt.h"
#include "logging.h"
#include "misc.h"
#include "optimizationOptions.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <optimizer.h>
#include <iostream>
#include <cassert>
#include <numeric>
#include <random>

namespace parpe {


/**
 * @brief getLocalOptimum
 * @param problem
 * @return int indicating status. 0: success, != 0: failure
 */

int getLocalOptimum(OptimizationProblem *problem) {
    Optimizer *optimizer = problem->getOptimizationOptions().createOptimizer();
    auto status = optimizer->optimize(problem);
    delete optimizer;
    return std::get<0>(status);
}

/**
 * @brief getLocalOptimumThreadWrapper wrapper for using getLocalOptimum with
 * pThreads.
 * @param problem
 * @return Pointer to int indicating status. 0: success, != 0: failure
 */

void *getLocalOptimumThreadWrapper(void *optimizationProblemVp) {
    OptimizationProblem *problem = (OptimizationProblem *)optimizationProblemVp;
    int *result = new int;
    *result = getLocalOptimum(problem);
    return result;
}

void runOptimizationsParallel(const OptimizationProblem **problems,
                              int numProblems) {
    runInParallelAndWaitForFinish(getLocalOptimumThreadWrapper,
                                  (void **)problems, numProblems);
}

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

    optimizationProblemGradientCheck(problem, parameterIndices.data(),
                                     numParameterIndicesToCheck, epsilon);
}

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      const int parameterIndices[],
                                      int numParameterIndices, double epsilon) {
    double fc = 0; // f(theta)
    std::vector<double> theta(problem->costFun->numParameters());
    problem->fillInitialParameters(theta.data());

    std::vector<double> gradient(theta.size());
    problem->costFun->evaluate(theta.data(), fc, gradient.data());

    std::vector<double> thetaTmp(theta);

    // printf("Index\tGradient\tfd_f\t\t(delta)\t\tfd_c\t\t(delta)\t\tfd_b\t\t(delta)\n");

    for (int i = 0; i < numParameterIndices; ++i) {
        int curInd = parameterIndices[i];
        double fb = 0, ff = 0; // f(theta + eps) , f(theta - eps)

        thetaTmp[curInd] = theta[curInd] + epsilon;
        problem->costFun->evaluate(thetaTmp.data(), ff, nullptr);

        thetaTmp[curInd] = theta[curInd] - epsilon;
        problem->costFun->evaluate(thetaTmp.data(), fb, nullptr);

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

        logmessage(ll, "%5d g: %12.6g  fd_c: %12.6g  Δ/fd_c: %.6e",
               curInd, curGrad, fd_c, regRelError);

    }
}


OptimizationProblem::OptimizationProblem(std::unique_ptr<GradientFunction> costFun)
    : costFun(std::move(costFun))
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
    return std::unique_ptr<OptimizationReporter>(new OptimizationReporter(costFun.get()));
}


void OptimizationProblem::fillInitialParameters(double *buffer) const {
    int numParameters = costFun->numParameters();
    std::vector<double> parametersMin(numParameters);
    std::vector<double> parametersMax(numParameters);
    fillParametersMin(parametersMin.data());
    fillParametersMax(parametersMax.data());

    fillArrayRandomDoubleIndividualInterval(parametersMin.data(), parametersMax.data(),
                                            numParameters, buffer);
}


OptimizationReporter::OptimizationReporter(GradientFunction *gradFun)
    : OptimizationReporter(gradFun, nullptr)
{

}

OptimizationReporter::OptimizationReporter(
        GradientFunction *gradFun,
        std::unique_ptr<OptimizationResultWriter> rw)
    : resultWriter(std::move(rw)),
      gradFun(gradFun)
{
    numParameters_ = gradFun->numParameters();
    cachedGradient.resize(numParameters_);
}

FunctionEvaluationStatus OptimizationReporter::evaluate(const double * const parameters, double &fval, double *gradient) const {
    if(beforeCostFunctionCall(numParameters_, parameters) != 0)
        return functionEvaluationFailure;

    if(gradient) {
        if (!haveCachedGradient || !std::equal(parameters, parameters + numParameters_,
                                               finalParameters.begin())) {
            // Have to compute anew
            cachedErrors = gradFun->evaluate(parameters, cachedCost, cachedGradient.data());
            std::copy(cachedGradient.begin(), cachedGradient.end(), gradient);
            haveCachedCost = true;
            haveCachedGradient = true;
        } else {
            // recycle old result
            std::copy(cachedGradient.begin(), cachedGradient.end(), gradient);
            fval = cachedCost;
        }
    } else {
        if (!haveCachedCost || !std::equal(parameters, parameters + numParameters_,
                                           finalParameters.begin())) {
            // Have to compute anew
            cachedErrors = gradFun->evaluate(parameters, cachedCost, nullptr);
            haveCachedCost = true;
            haveCachedGradient = false;
        }
        fval = cachedCost;
    }

    // update cached parameters
    finalParameters.resize(numParameters_);
    std::copy(parameters, parameters + numParameters_, finalParameters.begin());

    if(afterCostFunctionCall(numParameters_, parameters, cachedCost, gradient?cachedGradient.data():nullptr) != 0)
        return functionEvaluationFailure;

    return cachedErrors == 0 ? functionEvaluationSuccess : functionEvaluationFailure;
}

int OptimizationReporter::numParameters() const {
    return gradFun->numParameters();
}

bool OptimizationReporter::starting(int numParameters, const double * const initialParameters) const
{
    // If this is called multiple times (happens via IpOpt), don't do anything
    if(started)
        return false;

    wallTimer.reset();
//    timeOptimizationBegin = clock();
//    timeIterationBegin = timeOptimizationBegin;
//    timeCostEvaluationBegin = timeOptimizationBegin;

    if(resultWriter)
        resultWriter->starting(numParameters, initialParameters);

    started = true;

    return false;
}

bool OptimizationReporter::iterationFinished(const double * const parameters, double objectiveFunctionValue, const double * const objectiveFunctionGradient) const
{
    double wallTimeIter = wallTimer.getRound(); //(double)(clock() - timeIterationBegin) / CLOCKS_PER_SEC;
    double wallTimeOptim = wallTimer.getTotal(); //double)(clock() - timeOptimizationBegin) / CLOCKS_PER_SEC;

    logmessage(LOGLVL_INFO, "iter: %d cost: %g time_iter: %gs time_optim: %gs", numIterations, objectiveFunctionValue, wallTimeIter, wallTimeOptim);

    if(resultWriter)
        resultWriter->logLocalOptimizerIteration(numIterations, finalParameters.data(),
                                                 numParameters_, objectiveFunctionValue,
                                                 cachedGradient.data(), // This might be misleading, the gradient could evaluated at other parameters if there was a line search inbetween
                                                 wallTimeIter);
    ++numIterations;

    return false;
}

bool OptimizationReporter::beforeCostFunctionCall(int numParameters, const double * const parameters) const
{
    ++numFunctionCalls;
    //timeCostEvaluationBegin = clock();

    return false;
}

bool OptimizationReporter::afterCostFunctionCall(int numParameters, const double * const parameters, double objectiveFunctionValue, const double * const objectiveFunctionGradient) const
{
    //clock_t timeCostEvaluationEnd = clock();

    double wallTime = wallTimer.getTotal();//(double)(timeCostEvaluationEnd - timeCostEvaluationBegin) / CLOCKS_PER_SEC;

    if(resultWriter) {
        // TODO: problem.costfuj.getLastHiearchicalParameers lastFullParameterVector?
        resultWriter->logLocalOptimizerObjectiveFunctionEvaluation(parameters, numParameters, objectiveFunctionValue,
                                                                   objectiveFunctionGradient, numIterations, numFunctionCalls, wallTime);
    }
    return false;
}

void OptimizationReporter::finished(double optimalCost, const double *optimalParameters, int exitStatus) const
{
    double timeElapsed = wallTimer.getTotal();

    if(resultWriter)
        resultWriter->saveLocalOptimizerResults(optimalCost, optimalParameters, numParameters_, timeElapsed, exitStatus);
}

double OptimizationReporter::getFinalCost() const
{
    return cachedCost;
}

std::vector<double> OptimizationReporter::getFinalParameters() const
{
    return finalParameters;
}

} // namespace parpe
