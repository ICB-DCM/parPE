#ifndef LOGGER_H
#define LOGGER_H

#include "MultiConditionDataProvider.h"
#include "optimizationProblem.h"
#include "optimizationResultWriter.h"

namespace parpe {

// TODO get rid of JobIdentifier
// TODO merge with simulation resultwriter?
class MultiConditionProblemResultWriter : public OptimizationResultWriter {

public:
    MultiConditionProblemResultWriter();

    MultiConditionProblemResultWriter(hid_t file_id);

    MultiConditionProblemResultWriter(hid_t file_id, JobIdentifier id);

    MultiConditionProblemResultWriter(std::string filename, bool overwrite,
                                      JobIdentifier id);

    MultiConditionProblemResultWriter(MultiConditionProblemResultWriter const& other);

    void logLocalOptimizerObjectiveFunctionEvaluation(
            const double *parameters, int numParameters, double objectiveFunctionValue,
            const double *objectiveFunctionGradient, int numIterations, int numFunctionCalls,
            double timeElapsedInSeconds);

    std::string getIterationPath(int iterationIdx) override;

    std::string getSimulationPath(JobIdentifier id);

    std::string getOptimizationPath() override;

    JobIdentifier getJobId();

    void setJobId(JobIdentifier id);

    void logSimulation(JobIdentifier id, const double *theta, double llh,
                       const double *gradient, double timeElapsedInSeconds,
                       int nTheta, int numStates, double *states,
                       double *stateSensi, int numY, double *y, int jobId,
                       int iterationsUntilSteadystate, int status);

    void printObjectiveFunctionFailureMessage() const;

    void saveLocalOptimizerResults(double finalNegLogLikelihood,
                                           const double *optimalParameters,
                                           int numParameters, double masterTime,
                                           int exitStatus);


private:

    bool logLineSearch = false;
    JobIdentifier id;
    hid_t file_id = 0;
};

} // namespace parpe

#endif
