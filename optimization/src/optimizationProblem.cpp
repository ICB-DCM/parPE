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

    printf("Index\tGradient\tfd_f\t\t(delta)\t\tfd_c\t\t(delta)\t\tfd_b\t\t("
           "delta)\n");

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
        char status[] = " ";
        if(!((fd_c >= fd_f && fd_c <= fd_b)
             || (fd_c <= fd_f && fd_c >= fd_b)))
            status[0] = '!';
        if(!((curGrad >= fd_f && curGrad <= fd_b)
             || (curGrad <= fd_f && curGrad >= fd_b)))
            status[0] = '!';


        printf("%d\t%s\tg: %f\tfd_f: %f\t(%f)\tfd_c: %f\t(%f)\tfd_b: %f\t(%f)",
               curInd, status, curGrad,
               fd_f, curGrad - fd_f,
               fd_c, curGrad - fd_c,
               fd_b, curGrad - fd_b);
        printf("\t\tfb: %f\tfc: %f\tff: %f\t\n", fb, fc, ff);
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
    return std::unique_ptr<OptimizationReporter>(new OptimizationReporter());
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


OptimizationReporter::OptimizationReporter()
{

}

OptimizationReporter::OptimizationReporter(std::unique_ptr<OptimizationResultWriter> rw) : resultWriter(std::move(rw))
{

}

bool OptimizationReporter::starting(int numParameters, const double * const initialParameters)
{
    // If this is called multiple times (happens via IpOpt), don't do anything
    if(started)
        return false;

    this->numParameters = numParameters;
    timeOptimizationBegin = clock();
    timeIterationBegin = timeOptimizationBegin;
    timeCostEvaluationBegin = timeOptimizationBegin;

    if(resultWriter)
        resultWriter->starting(numParameters, initialParameters);

    started = true;

    return false;
}

bool OptimizationReporter::iterationFinished(int numParameters, const double * const parameters, double objectiveFunctionValue, const double * const objectiveFunctionGradient)
{
    double wallTimeIter = (double)(clock() - timeIterationBegin) / CLOCKS_PER_SEC;
    double wallTimeOptim = (double)(clock() - timeOptimizationBegin) / CLOCKS_PER_SEC;

    logmessage(LOGLVL_INFO, "iter: %d cost: %g time_iter: %gs time_optim: %gs", numIterations, objectiveFunctionValue, wallTimeIter, wallTimeOptim);

    if(resultWriter)
        resultWriter->logLocalOptimizerIteration(numIterations, parameters, numParameters, objectiveFunctionValue, objectiveFunctionGradient,
                                                 wallTimeIter);
    ++numIterations;

    return false;
}

bool OptimizationReporter::beforeCostFunctionCall(int numParameters, const double * const parameters)
{
    ++numFunctionCalls;
    timeCostEvaluationBegin = clock();

    return false;
}

bool OptimizationReporter::afterCostFunctionCall(int numParameters, const double * const parameters, double objectiveFunctionValue, const double * const objectiveFunctionGradient)
{
    clock_t timeCostEvaluationEnd = clock();

    double wallTime = (double)(timeCostEvaluationEnd - timeCostEvaluationBegin) / CLOCKS_PER_SEC;

    if(resultWriter)
        resultWriter->logLocalOptimizerObjectiveFunctionEvaluation(parameters, numParameters, objectiveFunctionValue,
                                                                   objectiveFunctionGradient, numIterations, numFunctionCalls, wallTime);
    return false;
}

void OptimizationReporter::finished(double optimalCost, const double *optimalParameters, int exitStatus)
{
    clock_t timeEnd = clock();
    double timeElapsed = (double)(timeEnd - timeOptimizationBegin) / CLOCKS_PER_SEC;

    if(resultWriter)
        resultWriter->saveLocalOptimizerResults(optimalCost, optimalParameters,  numParameters, timeElapsed, exitStatus);
}

} // namespace parpe
