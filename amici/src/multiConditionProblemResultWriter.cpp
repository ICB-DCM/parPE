#include "multiConditionProblemResultWriter.h"


MultiConditionProblemResultWriter::MultiConditionProblemResultWriter() : OptimizationResultWriter()
{

}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(OptimizationProblem *problem, hid_t file_id, JobIdentifier id) : OptimizationResultWriter(problem, file_id)
{
    this->path = path;
    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/iteration/%d/genotype/%d/experiment/%d/",
            id.idxLocalOptimization, id.idxLocalOptimizationIteration, 0, id.idxConditions);
    rootPath = fullGroupPath;
}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(OptimizationProblem *problem, const char *filename, bool overwrite, JobIdentifier id) : OptimizationResultWriter(problem, filename, overwrite)
{
    this->path = path;
    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/iteration/%d/genotype/%d/experiment/%d/",
            id.idxLocalOptimization, id.idxLocalOptimizationIteration, 0, id.idxConditions);
    rootPath = fullGroupPath;
}
