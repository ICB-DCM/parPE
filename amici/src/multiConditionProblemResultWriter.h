#ifndef LOGGER_H
#define LOGGER_H

#include <stdbool.h>
#include "optimizationProblem.h"
#include "MultiConditionDataProvider.h"
#include "optimizationResultWriter.h"


class MultiConditionProblemResultWriter : public OptimizationResultWriter {

public:
    MultiConditionProblemResultWriter();

    MultiConditionProblemResultWriter(OptimizationProblem *problem, hid_t file_id, JobIdentifier id);

    MultiConditionProblemResultWriter(OptimizationProblem *problem, const char *filename, bool overwrite, JobIdentifier id);

    JobIdentifier path;
};

#endif
