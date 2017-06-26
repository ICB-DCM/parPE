#ifndef LOGGER_H
#define LOGGER_H

#include <stdbool.h>
#include "optimizationProblem.h"
#include "MultiConditionDataProvider.h"
#include "optimizationResultWriter.h"


class MultiConditionProblemResultWriter : public OptimizationResultWriter {

public:
    MultiConditionProblemResultWriter();

    MultiConditionProblemResultWriter(hid_t file_id, JobIdentifier id);

    MultiConditionProblemResultWriter(const char *filename, bool overwrite, JobIdentifier id);

    void logSimulation(const double *theta,
                       double llh, const double *gradient,
                       double timeElapsedInSeconds,
                       int nTheta, int numStates,
                       double *states, double *stateSensi,
                       double *y, int jobId,
                       int iterationsUntilSteadystate);

    JobIdentifier path;
};

#endif
