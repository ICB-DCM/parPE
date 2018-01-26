#ifndef LOGGER_H
#define LOGGER_H

#include "MultiConditionDataProvider.h"
#include "optimizationProblem.h"
#include "optimizationResultWriter.h"

namespace parpe {

// TODO get rid of JobIdentifier
// TODO merge with simulation resultwriter? the only thing this one is doing extra is tracking
// job identifier and printing 2 messages
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

    std::string getIterationPath(int iterationIdx) const override;

    std::string getOptimizationPath() const override;

    JobIdentifier getJobId();

    void setJobId(JobIdentifier id);

    void printObjectiveFunctionFailureMessage() const;

    void saveLocalOptimizerResults(double finalNegLogLikelihood,
                                           const double *optimalParameters,
                                           int numParameters, double masterTime,
                                           int exitStatus) const override;


private:

    JobIdentifier id;
    hid_t file_id = 0;
};

} // namespace parpe

#endif
