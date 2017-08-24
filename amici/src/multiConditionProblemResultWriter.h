#ifndef LOGGER_H
#define LOGGER_H

#include "MultiConditionDataProvider.h"
#include "optimizationProblem.h"
#include "optimizationResultWriter.h"
#include <stdbool.h>

class MultiConditionProblemResultWriter : public OptimizationResultWriter {

  public:
    MultiConditionProblemResultWriter();

    MultiConditionProblemResultWriter(hid_t file_id, JobIdentifier id);

    MultiConditionProblemResultWriter(std::string filename, bool overwrite,
                                      JobIdentifier id);

    std::string getIterationPath(int iterationIdx);

    std::string getSimulationPath(JobIdentifier id);

    std::string getOptimizationPath();

    JobIdentifier getJobId();

    void setJobId(JobIdentifier id);

    void logSimulation(JobIdentifier id, const double *theta, double llh,
                       const double *gradient, double timeElapsedInSeconds,
                       int nTheta, int numStates, double *states,
                       double *stateSensi, int numY, double *y, int jobId,
                       int iterationsUntilSteadystate, int status);

    bool logLineSearch = false;

  protected:
    JobIdentifier id;
};

#endif
