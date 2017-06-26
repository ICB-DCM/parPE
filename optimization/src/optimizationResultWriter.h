#ifndef OPTIMIZATIONRESULTWRITER_H
#define OPTIMIZATIONRESULTWRITER_H

#include "optimizationProblem.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string>

class OptimizationProblem;

class OptimizationResultWriter
{
public:
    OptimizationResultWriter();

    // TODO: remove problem here; only required for numParams
    OptimizationResultWriter(OptimizationProblem *problem, hid_t file_id);

    OptimizationResultWriter(OptimizationProblem *problem, const char *filename, bool overwrite);

    int initResultHDFFile(const char *filename, bool overwrite);

    void closeResultHDFFile();

    void logLocalOptimizerObjectiveFunctionEvaluation(const double *parameters,
                                                      double objectiveFunctionValue,
                                                      const double *objectiveFunctionGradient,
                                                      int numFunctionCalls,
                                                      double timeElapsedInSeconds);

    // TODO add likelihood for individual experiments
    void logLocalOptimizerIteration(int numIterations,
                                    double *theta,
                                    double objectiveFunctionValue,
                                    const double *gradient,
                                    double timeElapsedInSeconds,
                                    int nTheta,
                                    int alg_mod,
                                    double inf_pr, double inf_du,
                                    double mu,
                                    double d_norm,
                                    double regularization_size,
                                    double alpha_du, double alpha_pr,
                                    int ls_trials);

    void logSimulation(const double *theta,
                       double llh, const double *gradient,
                       double timeElapsedInSeconds,
                       int nTheta, int numStates,
                       double *states, double *stateSensi,
                       double *y, int jobId,
                       int iterationsUntilSteadystate);

    void flushResultWriter();

    void saveTotalWalltime(const double timeInSeconds);

    void saveLocalOptimizerResults(double finalNegLogLikelihood,
                                   const double *optimalParameters,
                                   double masterTime,
                                   int exitStatus);

    ~OptimizationResultWriter();

    OptimizationProblem *problem;
    hid_t file_id;
    std::string rootPath = "/";
};

#endif // OPTIMIZATIONRESULTWRITER_H
